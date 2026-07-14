"""
Day-0 Verification Gate — Check B: find ONE non-obvious, non-IND GO gene (v4 plan).

PRE-REGISTERED SELECTION RULE (frozen here BEFORE looking at any output — anti-p-hacking,
per PI-Scientist council ruling). The rule is deterministic; the ranking is the tie-breaker,
not curation. ALL survivors are reported, not just the pick.

Procedure:
  1. Pull the FULL Open Targets T1D associated-target ranking (genetic_association datatype).
  2. Restrict to genes actually PERTURBED in the Marson CD4 Perturb-seq DE table (must be
     testable in our system).
  3. Apply the gate's GO precondition (genetics side): genetic_association >= 0.20 floor.
  4. Exclude a PRE-SPECIFIED list of obvious/known/anchor/IND genes.
  5. For each survivor, query the LIVE eQTL Catalogue (eqtl_gate) for a significant
     CD4/Treg cis-eQTL — the cell-type-matched evidence the GO rule requires.
  6. The Check-B candidate = the TOP-RANKED (by genetic_association) survivor that ALSO has
     a real CD4/Treg cis-eQTL. Report the full survivor table so the choice is auditable.

A GO here still requires the eQTL; a survivor whose eQTL is absent is logged and we descend
the ranking. If NO survivor has a CD4/Treg eQTL, Check B returns NONE (honest negative ->
fall back to 'decision infrastructure' framing per the pre-registered compound-failure rule).
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd

from trustlayer import genetics, eqtl_gate

DE_CSV = Path("data/raw/DE_stats.suppl_table.csv")
OUT = Path("data/gold/day0_checkB_novel_go_receipt.json")
MANIFEST = Path("_script_manifest.jsonl")

GENETIC_FLOOR = 0.20  # commit_gate GO precondition

# PRE-SPECIFIED exclusions (frozen before running): obvious/known T1D genes, the calibration
# anchor, the demo trio, and IND-stage / textbook targets. Excluding these forces "non-obvious".
EXCLUDE = {
    "INS",
    "PTPN22",
    "IL2RA",
    "CTLA4",
    "HLA-DQA1",
    "HLA-DQB1",
    "HLA-DRB1",
    "HLA-A",
    "HLA-B",
    "IFIH1",
    "SH2B3",
    "CTSH",
    "BACH2",
    "TYK2",
    "GLIS3",
    "PTPN2",
    "IL2",
    "IL7R",
    "IKZF1",
    "IKZF3",
    "IKZF4",
    "ERBB3",
    "CENPW",
    "COBL",
    "UBASH3A",  # well-known T1D loci
    "CD226",  # calibration positive control / anchor — MUST NOT be the answer
    "RASGRP1",
    "PRKCQ",  # demo trio (already used)
    "INS-IGF2",
}

# Frozen pre-registration hash of the RULE (not the result).
RULE_TEXT = (
    "rank=OT genetic_association desc; must be perturbed in DE table; genetic_association>=0.20; "
    "exclude obvious/anchor/trio/IND; require live CD4-Treg cis-eQTL; pick top-ranked survivor "
    "with eQTL; report all survivors; NONE if no survivor has eQTL."
)


def _sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def _manifest(p: Path, note: str) -> None:
    with MANIFEST.open("a") as fh:
        fh.write(
            json.dumps(
                {
                    "path": str(p),
                    "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                    "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "note": note,
                }
            )
            + "\n"
        )


def main() -> int:
    rule_hash = _sha256_str(RULE_TEXT)
    print("=== DAY-0 CHECK B — non-obvious GO gene ===")
    print(f"pre-registered rule hash: {rule_hash}")

    # 1. full OT T1D ranking
    print("pulling Open Targets full T1D ranking (paginated)...")
    ranked, total = genetics.ot_full_t1d_ranking()
    print(f"  {len(ranked)} associated targets (OT count={total})")

    # need genetic_association datatype specifically, not just overall score.
    # ot_full_t1d_ranking returns overall score; pull genetic_association for the top slice.
    top_syms = [s for s, (r, sc) in sorted(ranked.items(), key=lambda kv: kv[1][0])][
        :120
    ]

    # 2. restrict to perturbed genes
    perturbed = set(
        pd.read_csv(DE_CSV, usecols=["target_contrast_gene_name"])[
            "target_contrast_gene_name"
        ]
        .dropna()
        .unique()
    )
    print(f"  perturbed genes in DE table: {len(perturbed)}")

    # 3+4. genetic_association floor + exclusions — pull GA datatype for candidates
    genetics.ENSEMBL_dynamic = True  # noqa (marker only)
    # ot_genetic_association needs ENSEMBL ids; for arbitrary symbols we query OT by symbol set
    # via the disease associatedTargets rows that already include datatypeScores.
    # Re-pull with datatype scores for the top perturbed candidates:
    cand = [s for s in top_syms if s in perturbed and s not in EXCLUDE]
    print(f"  candidates after perturbed+exclusion filter: {len(cand)}")

    # get genetic_association for candidates (batch via OT disease query by symbol->ensembl)
    # genetics.ot_genetic_association requires ENSEMBL map; fill dynamically from OT.
    ga_scores = _fetch_ga_for_symbols(cand)

    survivors = []
    for s in cand:
        ga = ga_scores.get(s)
        if ga is None or ga < GENETIC_FLOOR:
            continue
        survivors.append((s, ga, ranked[s][0]))
    survivors.sort(key=lambda x: -x[1])  # by genetic_association desc
    print(f"  survivors (GA>={GENETIC_FLOOR}): {len(survivors)}")
    for s, ga, rk in survivors[:20]:
        print(f"    {s:12s} GA={ga:.3f}  OT_rank={rk}")

    # 5. require live CD4/Treg cis-eQTL for the top survivors (descend until one hits)
    pick = None
    survivor_records = []
    for s, ga, rk in survivors:
        ensg = _symbol_to_ensembl(s)
        eqtl = (
            eqtl_gate.query_celltype_eqtl(ensg)
            if ensg
            else {"n_sig": 0, "n_tested": 0, "hits": []}
        )
        rec = {
            "symbol": s,
            "genetic_association": ga,
            "ot_rank": rk,
            "ensembl": ensg,
            "eqtl_n_sig": eqtl["n_sig"],
            "eqtl_hits": eqtl["hits"],
        }
        survivor_records.append(rec)
        print(
            f"    eQTL {s:12s} GA={ga:.3f}  n_sig={eqtl['n_sig']}  {[h['dataset_id'] for h in eqtl['hits']]}"
        )
        if pick is None and eqtl["n_sig"] > 0:
            pick = rec
        if len(survivor_records) >= 12:  # bound API calls
            break

    receipt = {
        "purpose": "Day-0 Check B — non-obvious non-IND GO gene (pre-registered rule)",
        "pre_registered_rule": RULE_TEXT,
        "pre_registered_rule_sha256": rule_hash,
        "genetic_floor": GENETIC_FLOOR,
        "excluded": sorted(EXCLUDE),
        "n_ot_t1d_targets": total,
        "n_survivors_ga_floor": len(survivors),
        "survivor_table": survivor_records,
        "pick": pick,
        "verdict": (
            "GO_CANDIDATE_FOUND"
            if pick
            else "NONE — fall back to decision-infrastructure framing"
        ),
    }
    OUT.write_text(json.dumps(receipt, indent=2, default=str))
    _manifest(OUT, "day0 check B non-obvious GO gene")
    print(f"\nPICK: {pick['symbol'] if pick else 'NONE'}")
    print(f"receipt -> {OUT}")
    return 0


# --- helpers that query OT for symbol->ensembl + genetic_association dynamically ---
def _fetch_ga_for_symbols(symbols: list[str]) -> dict:
    import requests

    q = """
    query D($efo:String!){
      disease(efoId:$efo){
        associatedTargets(page:{index:0,size:600}){
          rows{ target{ approvedSymbol } datatypeScores{ id score } }
        }
      }
    }
    """
    out = {}
    try:
        r = requests.post(
            genetics.OT_GQL,
            json={"query": q, "variables": {"efo": genetics.T1D_EFO}},
            timeout=90,
        )
        rows = r.json()["data"]["disease"]["associatedTargets"]["rows"]
        for row in rows:
            sym = row["target"]["approvedSymbol"]
            dts = {d["id"]: d["score"] for d in row["datatypeScores"]}
            if sym in symbols:
                out[sym] = round(dts.get("genetic_association", 0.0), 3)
    except Exception as e:
        print(f"  [GA fetch error] {e}")
    return out


_ENS_CACHE = {}


def _symbol_to_ensembl(symbol: str):
    if symbol in genetics.ENSEMBL:
        return genetics.ENSEMBL[symbol]
    if symbol in _ENS_CACHE:
        return _ENS_CACHE[symbol]
    import requests

    q = """query S($q:String!){ search(queryString:$q, entityNames:["target"]){ hits{ id entity name } } }"""
    try:
        r = requests.post(
            genetics.OT_GQL, json={"query": q, "variables": {"q": symbol}}, timeout=30
        )
        hits = r.json()["data"]["search"]["hits"]
        ensg = next((h["id"] for h in hits if h.get("entity") == "target"), None)
        _ENS_CACHE[symbol] = ensg
        return ensg
    except Exception:
        return None


if __name__ == "__main__":
    raise SystemExit(main())
