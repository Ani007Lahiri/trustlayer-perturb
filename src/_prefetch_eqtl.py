"""Fast parallel eQTL n_sig prefetch for the enriched pool -> data/gold/_checkA_eqtl_cache.json.
Uses the same CD4/Treg dataset list as trustlayer.eqtl_gate, short timeouts, parallel."""

from __future__ import annotations
import json, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import requests
from trustlayer import eqtl_gate

CAND = [
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
    "SUOX",
    "TCF7L2",
    "APOBR",
    "HNF1A",
    "FTO",
    "CDKAL1",
    "KCNQ1",
    "PGM1",
]
OUT = Path("data/gold/_checkA_eqtl_cache.json")
DS = eqtl_gate.CD4_TREG_DATASETS
API = "https://www.ebi.ac.uk/eqtl/api/v2"

# original-9 eqtl values already known (from frozen reasoning set)
ORIG = json.loads(Path("data/gold/day0_checkA_frozen_reasoning_set.json").read_text())
KNOWN = {g["symbol"]: g["eqtl_n_sig"] for g in ORIG["reasoning_dependent_genes"]}


def ensembl(g):
    try:
        r = requests.get(
            "https://mygene.info/v3/query",
            params={"q": f"symbol:{g}", "species": "human", "fields": "ensembl.gene"},
            timeout=20,
        ).json()
        h = r.get("hits", [])
        e = h[0].get("ensembl", {}) if h else {}
        e = e[0] if isinstance(e, list) and e else e
        return (e or {}).get("gene")
    except Exception:
        return None


def nsig_for(ensg):
    if not ensg:
        return 0
    n = 0
    for ds in DS:
        try:
            url = f"{API}/datasets/{ds}/associations?gene_id={ensg}&size=1&nlog10p={eqtl_gate.SIG_NLOG10P}"
            r = requests.get(url, timeout=12)
            if r.status_code == 200 and r.json():
                n += 1
        except Exception:
            pass
    return n


def one(g):
    if g in KNOWN:
        return g, KNOWN[g]
    return g, nsig_for(ensembl(g))


if __name__ == "__main__":
    print("prefetching eqtl n_sig (parallel)...")
    out = {}
    with ThreadPoolExecutor(max_workers=8) as ex:
        for g, n in ex.map(one, CAND):
            out[g] = n
            print(f"  {g:8s} eqtl_n_sig={n}", flush=True)
    OUT.write_text(json.dumps(out, indent=2))
    print(f"wrote {OUT}")
