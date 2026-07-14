"""
eQTL gate — LIVE cell-type-matched eQTL evidence (v4 fix, replaces hand-typed booleans).

The believe/veto money-shot keys on ONE boolean per gene: `celltype_matched_eqtl`
(does the gene have a significant cis-eQTL in CD4/Treg/T-cells?). In v3 that field was a
hand-typed dict literal — the deepest residual circularity the audit caught. This module
derives it from the eQTL Catalogue (github.com/eQTL-Catalogue) via a live query.

Design:
  * `query_celltype_eqtl(ensembl_id)` hits the eQTL Catalogue REST API for every
    CD4/Treg/T-cell gene-expression dataset and returns the number with a significant
    (p <= 1e-5) cis-eQTL.
  * `load_or_query()` prefers the committed receipt (`data/gold/eqtl_gate_receipt.json`,
    produced by this module) so the pipeline stays laptop-runnable and offline-reproducible,
    and falls back to a live query when `--refresh` is passed.

The point is not that the API is called at demo time; it is that the boolean is now a
recorded, reproducible query result with a citable dataset list — not an assertion.
"""

from __future__ import annotations

import json
import time
import urllib.request
from pathlib import Path
from typing import Optional

EQTL_API = "https://www.ebi.ac.uk/eqtl/api/v2"
RECEIPT = Path("data/gold/eqtl_gate_receipt.json")
SIG_NLOG10P = 5.0  # p <= 1e-5, standard cis-eQTL significance

# Curated CD4/Treg/T-cell datasets (eQTL Catalogue QTD accessions, ge quant).
CD4_TREG_DATASETS = [
    "QTD000031",  # BLUEPRINT T-cell
    "QTD000036",  # Bossini-Castillo_2019 Treg_naive
    "QTD000464",  # Schmiedel_2018 Treg_memory
    "QTD000469",  # Schmiedel_2018 Treg_naive
    "QTD000479",  # Schmiedel_2018 CD4_T-cell_naive
    "QTD000484",  # Schmiedel_2018 CD4_T-cell_anti-CD3-CD28
    "QTD000600",  # Perez_2022 T4
    "QTD000611",  # OneK1K CD4_CTL
    "QTD000612",  # OneK1K CD4_Naive
    "QTD000613",  # OneK1K CD4_TCM
    "QTD000614",  # OneK1K CD4_TEM
]

TRIO_ENSEMBL = {
    "CD226": "ENSG00000150637",
    "RASGRP1": "ENSG00000172575",
    "PRKCQ": "ENSG00000065675",
}


def _get(url: str, tries: int = 3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(1.5)
    return None


def query_celltype_eqtl(ensembl_id: str, datasets=None) -> dict:
    """Count CD4/Treg/T-cell datasets with a significant cis-eQTL for `ensembl_id`.

    Returns {n_tested, n_sig, hits:[{dataset_id, nlog10p}]}. An empty hit list means
    "not tested / not significant", never "no eQTL anywhere" (cis window is +/-1 Mb).
    """
    datasets = datasets or CD4_TREG_DATASETS
    hits = []
    for did in datasets:
        url = (
            f"{EQTL_API}/datasets/{did}/associations"
            f"?gene_id={ensembl_id}&nlog10p={SIG_NLOG10P}&size=50"
        )
        try:
            rows = _get(url) or []
        except Exception:
            rows = []
        if rows:
            best = max(rows, key=lambda x: x.get("nlog10p", 0))
            hits.append({"dataset_id": did, "nlog10p": round(best.get("nlog10p", 0), 1)})
    return {"n_tested": len(datasets), "n_sig": len(hits), "hits": hits}


def celltype_matched_eqtl(ensembl_id: str) -> bool:
    """The single load-bearing boolean: >=1 CD4/Treg dataset with a significant cis-eQTL."""
    return query_celltype_eqtl(ensembl_id)["n_sig"] > 0


def load_or_query(symbol: str, refresh: bool = False) -> dict:
    """Return the recorded eQTL ruling for `symbol`.

    Prefers the committed receipt (offline-reproducible); with refresh=True, re-queries
    the live API and returns fresh evidence. Either way the boolean is a query RESULT.
    """
    if not refresh and RECEIPT.exists():
        rec = json.loads(RECEIPT.read_text())
        t = rec["trio"][symbol]
        return {
            "symbol": symbol,
            "celltype_matched_eqtl": t["celltype_matched_eqtl"],
            "n_datasets_sig": t.get("n_datasets_sig"),
            "evidence": t.get("evidence"),
            "source": "committed eqtl_gate_receipt.json (live query, date recorded)",
        }
    ensg = TRIO_ENSEMBL[symbol]
    q = query_celltype_eqtl(ensg)
    return {
        "symbol": symbol,
        "celltype_matched_eqtl": q["n_sig"] > 0,
        "n_datasets_sig": q["n_sig"],
        "evidence": f"live query: {q['n_sig']}/{q['n_tested']} CD4/Treg datasets with sig cis-eQTL",
        "source": "eQTL Catalogue REST API (live)",
    }


if __name__ == "__main__":
    import sys

    refresh = "--refresh" in sys.argv
    for sym in TRIO_ENSEMBL:
        r = load_or_query(sym, refresh=refresh)
        print(f"{sym:8s} celltype_matched_eqtl={r['celltype_matched_eqtl']!s:5s}  "
              f"({r['n_datasets_sig']} sig)  [{r['source']}]")
