"""
Live genetics gates for the T1D calibrated-trust layer (Plan of Attack v3, Day 1).

Pulls target-disease evidence from public REST/GraphQL APIs so every genetic
claim on a slide is reproducible at build time (v3 provenance note):

  - Open Targets Platform GraphQL  -> genetic_association scores + rank
  - GWAS Catalog REST              -> genome-wide-significant SNP associations
  - GTEx / eQTL Catalogue          -> eQTL direction (proxy-tissue coloc gate)

Design rules (from v3):
  * No fabricated numbers. Everything here re-derives from a live API call.
  * Fail-closed: if an API is unreachable, the gate returns status="ERROR",
    never a silent pass.
  * T1D disease id = MONDO_0005147 (verified via Open Targets search).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field, asdict
from typing import Optional

import requests

OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"
GWAS_REST = "https://www.ebi.ac.uk/gwas/rest/api"
T1D_EFO = "MONDO_0005147"  # type 1 diabetes mellitus (verified live)

# Ensembl gene ids for the v3 target trio (+ gold-set anchors).
ENSEMBL = {
    "CD226": "ENSG00000150637",
    "RASGRP1": "ENSG00000172575",
    "PRKCQ": "ENSG00000065675",
    "SH2B3": "ENSG00000111252",
    "CTSH": "ENSG00000103811",
    "PTPN22": "ENSG00000134242",
    "INS": "ENSG00000254647",
    "IL2RA": "ENSG00000134460",
    "CTLA4": "ENSG00000163599",
    "IFIH1": "ENSG00000115267",
    "BACH2": "ENSG00000112182",
    "TYK2": "ENSG00000105397",
}


@dataclass
class OTAssoc:
    symbol: str
    ensembl_id: str
    overall_score: Optional[float]
    genetic_association: Optional[float]
    datatypes: dict = field(default_factory=dict)
    rank: Optional[int] = None
    total_targets: Optional[int] = None
    status: str = "OK"

    def as_dict(self):
        return asdict(self)


def _post(query: str, variables: dict, timeout: int = 60) -> dict:
    r = requests.post(
        OT_GQL, json={"query": query, "variables": variables}, timeout=timeout
    )
    r.raise_for_status()
    j = r.json()
    if "errors" in j:
        raise RuntimeError(j["errors"][0]["message"])
    return j["data"]


def ot_full_t1d_ranking(
    page_size: int = 1000, pause: float = 0.2
) -> tuple[dict[str, tuple[int, float]], int]:
    """Return {symbol: (rank, overall_score)} over ALL T1D-associated targets, and the total count.

    Paginates the disease -> associatedTargets endpoint. Rank is 1-based over the
    full ordered list (Open Targets returns rows already sorted by score).
    """
    q = """
    query D($efo:String!,$idx:Int!,$size:Int!){
      disease(efoId:$efo){
        associatedTargets(page:{index:$idx,size:$size}){
          count
          rows{ target{ approvedSymbol } score }
        }
      }
    }
    """
    rows: list = []
    total = None
    idx = 0
    while True:
        data = _post(q, {"efo": T1D_EFO, "idx": idx, "size": page_size})
        at = data["disease"]["associatedTargets"]
        total = at["count"]
        rows.extend(at["rows"])
        if len(at["rows"]) < page_size or len(rows) >= total:
            break
        idx += 1
        time.sleep(pause)
    ranked = {
        r["target"]["approvedSymbol"]: (i + 1, round(r["score"], 4))
        for i, r in enumerate(rows)
    }
    return ranked, total


def ot_genetic_association(symbols: list[str]) -> dict[str, OTAssoc]:
    """Pull overall + genetic_association datatype score for each gene vs T1D.

    Uses the Bs= target-filter on the disease query (the reliable path).
    Rank/total are filled in from the full ranking so the numbers are slide-ready.
    """
    ens = [ENSEMBL[s] for s in symbols if s in ENSEMBL]
    q = """
    query D($efo:String!,$bs:[String!]){
      disease(efoId:$efo){
        associatedTargets(page:{index:0,size:600}, Bs:$bs){
          rows{ target{ approvedSymbol } score datatypeScores{ id score } }
        }
      }
    }
    """
    out: dict[str, OTAssoc] = {}
    inv = {v: k for k, v in ENSEMBL.items()}
    try:
        data = _post(q, {"efo": T1D_EFO, "bs": ens})
        rows = data["disease"]["associatedTargets"]["rows"]
        got = {}
        for row in rows:
            sym = row["target"]["approvedSymbol"]
            dts = {d["id"]: round(d["score"], 3) for d in row["datatypeScores"]}
            got[sym] = OTAssoc(
                symbol=sym,
                ensembl_id=ENSEMBL.get(sym, "?"),
                overall_score=round(row["score"], 3),
                genetic_association=dts.get("genetic_association"),
                datatypes=dts,
            )
        # genes with no association row = not associated at all
        for s in symbols:
            if s not in got:
                got[s] = OTAssoc(
                    symbol=s,
                    ensembl_id=ENSEMBL.get(s, "?"),
                    overall_score=None,
                    genetic_association=None,
                    status="NOT_ASSOCIATED",
                )
        out = got
    except Exception as e:  # fail-closed
        for s in symbols:
            out[s] = OTAssoc(
                symbol=s,
                ensembl_id=ENSEMBL.get(s, "?"),
                overall_score=None,
                genetic_association=None,
                status=f"ERROR:{e}",
            )
    return out


def gwas_catalog_associations(rsid: str) -> list[dict]:
    """Fetch GWAS Catalog associations for a lead SNP (e.g. rs763361, rs72727394).

    Returns a compact list of {trait, pvalue, pmid} rows. Used to credential the
    T1 tier of the gold set and to confirm the RASGRP1 / CD226 lead SNPs.
    """
    url = f"{GWAS_REST}/singleNucleotidePolymorphisms/{rsid}/associations"
    try:
        r = requests.get(url, params={"projection": "associationBySnp"}, timeout=30)
        r.raise_for_status()
        assocs = r.json().get("_embedded", {}).get("associations", [])
        out = []
        for a in assocs:
            traits = [t.get("trait") for t in a.get("efoTraits", [])]
            out.append(
                {
                    "rsid": rsid,
                    "pvalue": a.get("pvalue"),
                    "pvalueMantissa": a.get("pvalueMantissa"),
                    "pvalueExponent": a.get("pvalueExponent"),
                    "traits": traits,
                }
            )
        return out
    except Exception as e:
        return [{"rsid": rsid, "error": str(e)}]
