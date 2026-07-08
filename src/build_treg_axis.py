"""
Day-2: build the Treg program axis (Plan v3 Fix 3 dependency), LEAKAGE-CORRECTED.

A UCell-style signed contrast score per perturbation over a curated Treg-program gene set,
computed on the trans zscore matrix. Score direction: a perturbation that pushes cells
TOWARD a stable-Treg program scores positive; toward ex-Treg/effector, negative.

Hard constraints (Gate-1 critique):
  1. The program gene set must NOT include any gene in axis_forbidden_genes (gold set).
  2. When scoring perturbation p, the perturbed gene p's own readout column is EXCLUDED
     (per-row on-target exclusion) so a perturbation targeting a program gene can't score
     high trivially via its own knockdown.

The gene set is canonical Treg biology (curated, cited in comments), deliberately DISJOINT
from the gold set so downstream recovery is name-level non-circular. The axis-swap control
(Day 4) tests biological-subspace disjointness against >=3 orthogonal non-Treg axes.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

H5AD = Path("data/raw/GWCD4i.DE_stats.h5ad")
GOLD = Path("data/gold/t1d_gold_set.json")
OUT = Path("data/interim/treg_axis_scores.parquet")
MANIFEST = Path("data/gold/treg_axis_manifest.json")

# Canonical stable-Treg / suppressive program (POSITIVE pole).
# Curated from core Treg biology (Sakaguchi lineage TFs, IL2-STAT5 axis, suppressive
# effectors). FOXP3 handled as gold-set T4 -> excluded if forbidden.
TREG_POS = [
    "IKZF2",
    "IL2RB",
    "IL10",
    "CTLA4",
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
]
# Ex-Treg / effector-Th1/Th17 destabilization program (NEGATIVE pole).
TREG_NEG = [
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
]


def _read_categorical(h5obj) -> np.ndarray:
    if isinstance(h5obj, h5py.Group) and "categories" in h5obj and "codes" in h5obj:
        cats = h5obj["categories"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
        codes = h5obj["codes"][:]
        return cats[codes]
    arr = h5obj[:]
    return np.array([a.decode() if isinstance(a, bytes) else a for a in arr])


def main() -> int:
    gold = json.loads(GOLD.read_text())
    forbidden = set(gold["axis_forbidden_genes"])

    # ---- HARD CONSTRAINT 1: strip forbidden (gold) genes from the program ----
    pos = [g for g in TREG_POS if g not in forbidden]
    neg = [g for g in TREG_NEG if g not in forbidden]
    dropped = sorted((set(TREG_POS) | set(TREG_NEG)) & forbidden)

    with h5py.File(H5AD, "r") as f:
        target_gene = _read_categorical(f["obs/target_contrast_gene_name"])
        condition = _read_categorical(f["obs/culture_condition"])
        var_names = _read_categorical(f["var/gene_name"])
        col = {g: i for i, g in enumerate(var_names)}
        pos_cols = np.array([col[g] for g in pos if g in col])
        neg_cols = np.array([col[g] for g in neg if g in col])
        pos_used = [g for g in pos if g in col]
        neg_used = [g for g in neg if g in col]

        z = f["layers/zscore"]
        n = z.shape[0]
        score = np.empty(n)
        CHUNK = 1000
        for s in range(0, n, CHUNK):
            e = min(s + CHUNK, n)
            block = np.nan_to_num(z[s:e, :], nan=0.0)
            for r in range(s, e):
                # HARD CONSTRAINT 2: exclude the perturbed gene's own column
                tcol = col.get(target_gene[r])
                prow = block[r - s]
                pc = pos_cols[pos_cols != tcol] if tcol is not None else pos_cols
                nc = neg_cols[neg_cols != tcol] if tcol is not None else neg_cols
                pos_mean = prow[pc].mean() if pc.size else 0.0
                neg_mean = prow[nc].mean() if nc.size else 0.0
                score[r] = pos_mean - neg_mean  # signed contrast

    df = pd.DataFrame(
        {"gene": target_gene, "condition": condition, "treg_axis_score": score}
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT)

    # ---- assertions (fail-closed) ----
    assert not (set(pos_used) & forbidden), "LEAK: forbidden gene in POS program"
    assert not (set(neg_used) & forbidden), "LEAK: forbidden gene in NEG program"

    manifest = {
        "axis": "Treg stable/suppressive (POS) minus ex-Treg/effector (NEG)",
        "pos_genes_used": pos_used,
        "neg_genes_used": neg_used,
        "dropped_because_in_gold_set": dropped,
        "per_row_ontarget_exclusion": True,
        "forbidden_disjoint": True,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2))

    print("=" * 62)
    print("TREG PROGRAM AXIS  (leakage-corrected)")
    print("=" * 62)
    print(f"  POS genes used ({len(pos_used)}): {pos_used}")
    print(f"  NEG genes used ({len(neg_used)}): {neg_used}")
    print(f"  dropped (in gold set): {dropped}")
    print(f"  per-row on-target exclusion: ON")
    print(
        f"\n  score dist: mean={score.mean():+.3f} std={score.std():.3f} "
        f"[{score.min():+.2f}, {score.max():+.2f}]"
    )
    # show trio
    for g in ["CD226", "RASGRP1", "PRKCQ"]:
        sub = df[df.gene == g]
        if len(sub):
            print(
                f"  {g:8s} treg_axis (by cond): "
                + ", ".join(
                    f"{r['condition']}={r['treg_axis_score']:+.3f}"
                    for _, r in sub.iterrows()
                )
            )
    print(f"\n  -> {OUT}\n  -> {MANIFEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
