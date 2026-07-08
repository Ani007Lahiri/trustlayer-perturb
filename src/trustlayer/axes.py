"""
Day-4: program axes (Treg + orthogonal controls) for the specificity control.
Critique-corrected (Gate-1 Day-4): interferon dropped as a clean control (IFIH1 is a
gold gene AND an IFN-pathway gene); replaced with cholesterol/metabolic + ribosome-
biogenesis. Adds a magnitude-only baseline. All axes gold-disjoint + per-row on-target
exclusion. Frozen gene-set lists (pre-registered before compute).
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np

H5AD = Path("data/raw/GWCD4i.DE_stats.h5ad")
GOLD = Path("data/gold/t1d_gold_set.json")

# ---- FROZEN, pre-registered gene sets (curated; gold-disjoint enforced at runtime) ----
AXES = {
    "treg": {  # stable/suppressive (POS) minus ex-Treg/effector (NEG)
        "pos": [
            "IKZF2",
            "IL2RB",
            "IL10",
            "TNFRSF18",
            "TNFRSF4",
            "IKZF4",
            "LRRC32",
            "ENTPD1",
            "NT5E",
            "TGFB1",
            "LGALS1",
            "CCR8",
            "TIGIT",
            "FOXP3",
            "IL2RA",
            "CTLA4",
        ],
        "neg": [
            "IFNG",
            "TBX21",
            "IL17A",
            "RORC",
            "IL21",
            "CSF2",
            "IL23R",
            "RUNX3",
            "STAT4",
            "IL12RB2",
            "CCR6",
        ],
    },
    "cellcycle": {  # proliferation program (orthogonal control)
        "pos": [
            "MKI67",
            "PCNA",
            "TOP2A",
            "CCNB1",
            "CCNB2",
            "CDK1",
            "CCNA2",
            "AURKB",
            "BUB1",
            "CENPA",
            "CENPE",
            "FOXM1",
            "PLK1",
            "TYMS",
            "MCM2",
            "MCM6",
        ],
        "neg": ["CDKN1A", "CDKN1B", "CDKN2A", "RB1"],
    },
    "cholesterol": {  # cholesterol/metabolic (orthogonal replacement for interferon)
        "pos": [
            "HMGCR",
            "HMGCS1",
            "LDLR",
            "SREBF2",
            "SQLE",
            "DHCR7",
            "DHCR24",
            "MVD",
            "FDPS",
            "IDI1",
            "MSMO1",
            "INSIG1",
            "CYP51A1",
            "ACAT2",
        ],
        "neg": ["ABCA1", "ABCG1", "CYP7A1"],
    },
    "apoptosis": {  # apoptosis program (orthogonal control)
        "pos": [
            "BAX",
            "BAK1",
            "CASP3",
            "CASP7",
            "CASP9",
            "BBC3",
            "PMAIP1",
            "BID",
            "APAF1",
            "FAS",
            "FADD",
            "TNFRSF10B",
        ],
        "neg": ["BCL2", "BCL2L1", "MCL1", "BIRC5", "XIAP"],
    },
    "ribosome": {  # ribosome biogenesis (orthogonal control)
        "pos": [
            "RPL3",
            "RPL5",
            "RPL7",
            "RPS3",
            "RPS6",
            "RPS19",
            "NPM1",
            "FBL",
            "NCL",
            "BOP1",
            "GNL3",
            "RRS1",
            "WDR12",
            "UTP20",
        ],
        "neg": [],
    },
}
# interferon retained ONLY as a labeled shared-mechanism control (IFIH1 is a gold gene)
IFN_SHARED_MECHANISM = {
    "pos": [
        "ISG15",
        "MX1",
        "MX2",
        "OAS1",
        "OAS2",
        "OAS3",
        "IFIT1",
        "IFIT3",
        "IRF7",
        "STAT1",
        "STAT2",
        "RSAD2",
        "USP18",
        "IFI6",
        "IFI44",
        "DDX58",
    ],
    "neg": [],
}


def _decode_cat(h5obj) -> np.ndarray:
    if isinstance(h5obj, h5py.Group) and "categories" in h5obj and "codes" in h5obj:
        cats = np.array(
            [c.decode() if isinstance(c, bytes) else c for c in h5obj["categories"][:]]
        )
        return cats[h5obj["codes"][:]]
    return np.array([a.decode() if isinstance(a, bytes) else a for a in h5obj[:]])


def score_all_axes() -> tuple["object", dict]:
    """Compute every axis score + trans_effect_magnitude per perturbation-condition row.

    Returns (DataFrame[gene,condition,<axis>_score...,trans_effect_magnitude], manifest).
    All axes: signed contrast (POS mean - NEG mean), gold genes stripped, perturbed gene's
    own column excluded per row.
    """
    import pandas as pd

    gold = json.loads(GOLD.read_text())
    forbidden = set(gold["axis_forbidden_genes"])

    axis_defs = {
        k: {
            "pos": [g for g in v["pos"] if g not in forbidden],
            "neg": [g for g in v["neg"] if g not in forbidden],
        }
        for k, v in AXES.items()
    }
    axis_defs["ifn_shared"] = {
        "pos": [g for g in IFN_SHARED_MECHANISM["pos"] if g not in forbidden],
        "neg": [g for g in IFN_SHARED_MECHANISM["neg"] if g not in forbidden],
    }

    with h5py.File(H5AD, "r") as f:
        gene = _decode_cat(f["obs"]["target_contrast_gene_name"])
        cond = _decode_cat(f["obs"]["culture_condition"])
        var_names = _decode_cat(f["var"]["gene_name"])
        col = {g: i for i, g in enumerate(var_names)}

        # resolve columns per axis (only measured genes)
        cols = {}
        used = {}
        for name, d in axis_defs.items():
            pc = np.array([col[g] for g in d["pos"] if g in col])
            nc = np.array([col[g] for g in d["neg"] if g in col])
            cols[name] = (pc, nc)
            used[name] = {
                "pos_used": [g for g in d["pos"] if g in col],
                "neg_used": [g for g in d["neg"] if g in col],
            }

        z = f["layers"]["zscore"]
        n = z.shape[0]
        scores = {name: np.empty(n) for name in axis_defs}
        trans_mag = np.empty(n)
        CHUNK = 1000
        for s in range(0, n, CHUNK):
            e = min(s + CHUNK, n)
            block = np.nan_to_num(z[s:e, :], nan=0.0)
            full_sq = np.sum(block * block, axis=1)
            for r in range(s, e):
                prow = block[r - s]
                tcol = col.get(gene[r])
                ontarget_sq = prow[tcol] ** 2 if tcol is not None else 0.0
                trans_mag[r] = np.sqrt(max(full_sq[r - s] - ontarget_sq, 0.0))
                for name, (pc, nc) in cols.items():
                    pcx = pc[pc != tcol] if tcol is not None else pc
                    ncx = nc[nc != tcol] if tcol is not None else nc
                    pm = prow[pcx].mean() if pcx.size else 0.0
                    nm = prow[ncx].mean() if ncx.size else 0.0
                    scores[name][r] = pm - nm

    df = pd.DataFrame(
        {"gene": gene, "condition": cond, "trans_effect_magnitude": trans_mag}
    )
    for name in axis_defs:
        df[f"{name}_score"] = scores[name]

    # ---- assertions: no axis uses a gold gene ----
    for name, u in used.items():
        leaked = (set(u["pos_used"]) | set(u["neg_used"])) & forbidden
        assert not leaked, f"LEAK: axis {name} uses gold genes {leaked}"

    manifest = {
        "axes_used": used,
        "gold_disjoint": True,
        "per_row_ontarget_exclusion": True,
        "ifn_note": "ifn_shared retained ONLY as a shared-mechanism positive-leaning "
        "control (IFIH1 is a gold gene AND an IFN-pathway gene); NOT a "
        "clean orthogonal negative.",
    }
    return df, manifest
