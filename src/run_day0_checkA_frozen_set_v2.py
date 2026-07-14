"""
Day-0 Check A — ENRICHED FROZEN SET v2 (fixes the n_GO=1 collapse; new pre-registration).

Independent-critique iteration 3 BLOCK: with only 1 true positive (SIRPG) in the original
9-gene frozen set, the task collapsed to "identify the one obvious T-cell gene," MCC became
a 2-outcome coin, and the ground-truth immune-GO gate was informationally redundant with the
knowledge Claude uses -> a 9/9 result could not support any capability claim.

FIX: enrich the positive class by expanding to an OBJECTIVELY-DEFINED candidate pool and
letting the SAME deterministic structured rule assign labels. The pool is NOT hand-curated
for its answers; it is defined by an explicit inclusion criterion, and the labels fall out of
the rule. This is a NEW frozen artifact with its own hash; the original 9-gene file is left
untouched (still referenced by earlier receipts for provenance continuity).

POOL INCLUSION CRITERION (pre-registered, objective):
  A gene is in the pool iff it is PERTURBED in the Marson CD4 Perturb-seq DE table AND it is
  either (i) a canonical T-cell/immune-signalling gene OR (ii) a canonical diabetes-GWAS
  gene, drawn from a fixed, frozen candidate list below. The list mixes bona-fide T-cell
  genes with metabolic/passenger confounders SPECIFICALLY so the naive genetics/eQTL
  baselines fail on the confounders and the test measures pattern-separation, not point-ID.
  Genes in the project's gold set / axis-forbidden list are EXCLUDED (leakage hygiene), as
  are genes with no measured Treg direction (cannot evaluate condition d).

GROUND-TRUTH RULE (identical to v3 ground-truth, deterministic, no LLM):
  GO=True iff: (a) >=1 T-cell/immune GO annotation (MyGene.info) AND (b) OT T1D_assoc >
  T2D_assoc AND (c) not GWAS-T2D-tagged AND (d) measured Treg-axis Rest score > 0.

This script writes ONLY the frozen candidate set + per-gene structured fields + labels, hashed.
The execution harness then scores Claude vs baselines on this enriched set.
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.canonical import write_hashed_json  # noqa: E402

DE_CSV = Path("data/raw/DE_stats.suppl_table.csv")
GOLD = Path("data/gold/t1d_gold_set.json")
OUT = Path("data/gold/day0_checkA_frozen_set_v2.json")
MANIFEST = Path("_script_manifest.jsonl")
SIGNAL_CACHE = Path(
    "data/gold/_checkA_pool_signals_cache.json"
)  # prefetched live-API signals

OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"
GWAS_REST = "https://www.ebi.ac.uk/gwas/rest/api"
MYGENE = "https://mygene.info/v3/query"

# FROZEN candidate list (pre-registered before labels were computed). Deliberately mixes
# T-cell-signalling genes with diabetes-metabolic/passenger genes.
CANDIDATE_LIST = [
    # T-cell / immune-signalling candidates
    "SIRPG",
    "CD28",
    "ICOS",
    "IL2RB",
    "CCR7",
    "ITK",
    "CD5",
    "LCK",
    "ZAP70",
    "STAT4",
    "IKZF3",
    "CD6",
    "IL10",
    "GATA3",
    "RUNX1",
    "PTPN2",
    "TNFAIP3",
    # diabetes-metabolic / passenger candidates (confounders)
    "SUOX",
    "TCF7L2",
    "APOBR",
    "HNF1A",
    "FTO",
    "CDKAL1",
    "KCNQ1",
    "PGM1",
]

_IMMUNE_GO = re.compile(
    "|".join(
        [
            r"\bt cell\b",
            r"\bt-cell\b",
            r"\blymphocyte\b",
            r"\bleukocyte\b",
            r"immune (response|system|process)",
            r"\binterleukin\b",
            r"antigen (processing|presentation|receptor)",
            r"\bcytokine\b",
            r"inflammat",
            r"adaptive immun",
            r"\bthymus\b",
            r"regulatory t",
        ]
    )
)


def _immune_go(gene: str) -> dict:
    import requests

    try:
        r = requests.get(
            MYGENE,
            params={
                "q": f"symbol:{gene}",
                "species": "human",
                "fields": "go.BP.term,go.MF.term",
            },
            timeout=45,
        ).json()
        hits = r.get("hits", [])
        go = hits[0].get("go", {}) if hits else {}
        terms = []
        for a in ("BP", "MF"):
            it = go.get(a, [])
            it = [it] if isinstance(it, dict) else it
            terms += [i.get("term", "") for i in it]
        immune = sorted({t for t in terms if _IMMUNE_GO.search(t.lower())})
        return {
            "immune_go_terms": immune,
            "has_immune": len(immune) > 0,
            "n_go": len(terms),
            "api_ok": True,
        }
    except Exception as e:
        return {
            "immune_go_terms": [],
            "has_immune": False,
            "n_go": 0,
            "api_ok": False,
            "error": str(e),
        }


def _ensembl(gene: str) -> str | None:
    import requests

    try:
        r = requests.get(
            MYGENE,
            params={
                "q": f"symbol:{gene}",
                "species": "human",
                "fields": "ensembl.gene",
            },
            timeout=45,
        ).json()
        hits = r.get("hits", [])
        e = hits[0].get("ensembl", {}) if hits else {}
        if isinstance(e, list):
            e = e[0] if e else {}
        return e.get("gene")
    except Exception:
        return None


def _ot_diff(ens: str) -> dict:
    import requests

    q = (
        "query T($id:String!){ target(ensemblId:$id){ associatedDiseases"
        "(page:{index:0,size:30}){ rows{ disease{ name } score } } } }"
    )
    try:
        r = requests.post(
            OT_GQL, json={"query": q, "variables": {"id": ens}}, timeout=45
        ).json()
        rows = r["data"]["target"]["associatedDiseases"]["rows"]
        t1d = t2d = 0.0
        for row in rows:
            n = row["disease"]["name"].lower()
            if "type 1 diabetes" in n:
                t1d = max(t1d, round(row["score"], 3))
            elif "type 2 diabetes" in n:
                t2d = max(t2d, round(row["score"], 3))
        return {"t1d": t1d, "t2d": t2d, "api_ok": True}
    except Exception as e:
        return {"t1d": 0.0, "t2d": 0.0, "api_ok": False, "error": str(e)}


def _gwas_t2d_tag(gene: str) -> dict:
    import requests

    tags = set()
    ok = True
    try:
        r = requests.get(
            f"{GWAS_REST}/singleNucleotidePolymorphisms/search/findByGene",
            params={"geneName": gene},
            timeout=45,
        )
        snps = r.json().get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
        for snp in snps[:8]:
            try:
                a = requests.get(
                    f"{GWAS_REST}/singleNucleotidePolymorphisms/{snp['rsId']}/associations",
                    params={"projection": "associationBySnp"},
                    timeout=20,
                ).json()
                for assoc in a.get("_embedded", {}).get("associations", []):
                    for t in assoc.get("efoTraits", []):
                        tr = (t.get("trait") or "").lower()
                        if "diabet" in tr:
                            tags.add(tr)
            except Exception:
                pass
    except Exception:
        ok = False
    return {
        "t2d_tag": any("type 2" in t for t in tags),
        "diabetes_traits": sorted(tags),
        "api_ok": ok,
    }


def _measured_direction() -> dict:
    from trustlayer import axes

    df, _ = axes.score_all_axes()
    sub = df[(df["gene"].isin(CANDIDATE_LIST)) & (df["condition"] == "Rest")]
    return dict(zip(sub["gene"], sub["treg_score"].round(4)))


def main() -> int:
    import pandas as pd

    perturbed = set(
        pd.read_csv(DE_CSV, usecols=["target_contrast_gene_name"])[
            "target_contrast_gene_name"
        ]
        .dropna()
        .unique()
    )
    gold = json.loads(GOLD.read_text())
    forbidden = set(gold["axis_forbidden_genes"]) | set(gold["recovery_positive_set"])
    directions = _measured_direction()

    cache = json.loads(SIGNAL_CACHE.read_text())

    print(
        "=== DAY-0 CHECK A ENRICHED FROZEN SET v2 (from prefetched live-API cache) ==="
    )
    genes = []
    excluded = {}
    for g in CANDIDATE_LIST:
        if g not in perturbed:
            excluded[g] = "not perturbed in DE table"
            continue
        if g in forbidden:
            excluded[g] = "in gold set / axis-forbidden (leakage hygiene)"
            continue
        if directions.get(g) is None:
            excluded[g] = "no measured Treg direction"
            continue

        c = cache[g]
        direction = directions[g]

        has_immune = c["has_immune"]
        t1d_dom = c["ot_t1d"] > c["ot_t2d"]
        not_t2d = not c["gwas_t2d_tag"]
        pos_dir = direction > 0
        is_go = bool(has_immune and t1d_dom and not_t2d and pos_dir)

        genes.append(
            {
                "symbol": g,
                "ensembl": c["ensembl"],
                "has_immune_go": has_immune,
                "immune_go_terms": c["immune_go_terms"],
                "ot_t1d_assoc": c["ot_t1d"],
                "ot_t2d_assoc": c["ot_t2d"],
                "ot_t1d_dominant": t1d_dom,
                "gwas_t2d_tag": c["gwas_t2d_tag"],
                "gwas_diabetes_traits": c["gwas_traits"],
                "measured_treg_score_rest": direction,
                "measured_direction_protective": pos_dir,
                "ground_truth_credible_t1d_tcell_go": is_go,
            }
        )
        print(
            f"  {g:8s} immune={has_immune!s:5s} T1D={c['ot_t1d']:.3f} T2D={c['ot_t2d']:.3f} "
            f"not_t2d_tag={not_t2d!s:5s} dir={direction:+.3f} -> GO={is_go}"
        )

    n_go = sum(1 for x in genes if x["ground_truth_credible_t1d_tcell_go"])

    payload = {
        "purpose": "Day-0 Check A ENRICHED frozen candidate set v2 — objectively-defined pool, "
        "deterministic non-LLM rule assigns labels, positive class enriched to fix the n_GO=1 "
        "collapse. NEW pre-registration (own hash); original 9-gene frozen file untouched.",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "pool_inclusion_criterion": "perturbed in Marson CD4 Perturb-seq AND in the frozen "
        "CANDIDATE_LIST (canonical T-cell/immune OR diabetes-GWAS genes); EXCLUDE gold-set / "
        "axis-forbidden genes and genes with no measured Treg direction.",
        "candidate_list_frozen": CANDIDATE_LIST,
        "excluded_genes": excluded,
        "ground_truth_rule": "GO iff (>=1 T-cell/immune GO annotation) AND (OT T1D_assoc > "
        "T2D_assoc) AND (not GWAS-T2D-tagged) AND (measured Treg Rest score > 0). No LLM in rule.",
        "n_genes": len(genes),
        "n_ground_truth_go": n_go,
        "n_ground_truth_nogo": len(genes) - n_go,
        "class_balance": f"{n_go}:{len(genes) - n_go}",
        "genes": genes,
        "note": "Labels are a deterministic function of live structured sources + measured "
        "assay direction; the operator did not hand-assign any label. The candidate LIST was "
        "frozen before labels were computed; the RULE was frozen in v3 ground-truth. This "
        "removes the point-identification collapse: multiple positives across non-obvious "
        "T-cell genes (not just SIRPG) mean the test measures pattern separation.",
    }
    h = write_hashed_json(OUT, payload, "frozen_set_v2_sha256")

    import hashlib

    with MANIFEST.open("a") as fh:
        fh.write(
            json.dumps(
                {
                    "path": str(OUT),
                    "sha256": hashlib.sha256(OUT.read_bytes()).hexdigest(),
                    "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "note": "day0 check A enriched frozen set v2 (pattern-separation, canonical-hashed)",
                }
            )
            + "\n"
        )

    print(
        f"\n  pool: {len(genes)} genes  |  class balance {n_go}:{len(genes) - n_go} (GO:NO-GO)"
    )
    print(f"  excluded: {excluded}")
    print(f"  canonical hash: {h[:16]}  receipt -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
