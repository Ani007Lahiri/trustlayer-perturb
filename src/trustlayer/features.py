"""
Day-2: assemble the LEAKAGE-SAFE feature/label table for the base predictor.

Target (leakage-corrected per Gate-1 critique):
    trans_effect_magnitude(p,c) = || zscore vector EXCLUDING the perturbed gene ||_2
    primary target = log1p(trans_effect_magnitude)

Features (all independent of the DE readout that produces the target):
    - target_baseMean        (pre-perturbation expression; admissible)
    - n_cells_target         (assay depth)
    - condition one-hot      (Rest / Stim8hr / Stim48hr)
    - [optional] ESM-2 embedding of the target protein (added separately if it helps)

FORBIDDEN as features (outputs of the same DE computation as the target):
    ontarget_effect_size, n_up_genes, n_down_genes, n_total_de_genes, n_downstream,
    ontarget_significant, ontarget_effect_category, n_total_genes_category
A build-time assertion enforces this disjointness.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

H5AD = Path("data/raw/GWCD4i.DE_stats.h5ad")
SPLITS = Path("data/gold/frozen_splits.json")
OUT = Path("data/interim/day2_feature_table.parquet")

FORBIDDEN_FEATURES = {
    "ontarget_effect_size",
    "n_up_genes",
    "n_down_genes",
    "n_total_de_genes",
    "n_downstream",
    "ontarget_significant",
    "ontarget_effect_category",
    "n_total_genes_category",
}


def _read_categorical(h5obj) -> np.ndarray:
    """Decode an anndata-encoded categorical (categories + codes) or plain string dataset."""
    if isinstance(h5obj, h5py.Group) and "categories" in h5obj and "codes" in h5obj:
        cats = h5obj["categories"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
        codes = h5obj["codes"][:]
        return cats[codes]
    arr = h5obj[:]
    return np.array([a.decode() if isinstance(a, bytes) else a for a in arr])


def compute_trans_effect_magnitude() -> pd.DataFrame:
    """Stream the zscore layer, compute per-row trans-only L2 norm (excluding the
    perturbed gene's own column). Returns a DataFrame with obs keys + target."""
    with h5py.File(H5AD, "r") as f:
        target_gene = _read_categorical(f["obs/target_contrast_gene_name"])
        condition = _read_categorical(f["obs/culture_condition"])
        baseMean = f["obs/target_baseMean"][:]
        n_cells = f["obs/n_cells_target"][:]
        var_names = _read_categorical(f["var/gene_name"])
        gene_to_col = {g: i for i, g in enumerate(var_names)}

        z = f["layers/zscore"]  # (33983, 10282)
        n_obs = z.shape[0]
        full_sq = np.empty(n_obs)  # ||z||^2 over ALL genes
        ontarget_sq = np.zeros(n_obs)  # perturbed gene's own z^2 (to subtract)

        CHUNK = 1000
        for start in range(0, n_obs, CHUNK):
            end = min(start + CHUNK, n_obs)
            block = z[start:end, :]
            block = np.nan_to_num(block, nan=0.0)
            full_sq[start:end] = np.sum(block * block, axis=1)
            # subtract the on-target column for each row where present
            for r in range(start, end):
                col = gene_to_col.get(target_gene[r])
                if col is not None:
                    ontarget_sq[r] = block[r - start, col] ** 2

        trans_sq = np.clip(full_sq - ontarget_sq, 0, None)
        trans_mag = np.sqrt(trans_sq)

    df = pd.DataFrame(
        {
            "gene": target_gene,
            "condition": condition,
            "target_baseMean": baseMean,
            "n_cells_target": n_cells,
            "trans_effect_magnitude": trans_mag,
            "y": np.log1p(trans_mag),
        }
    )
    return df


def build_feature_table() -> pd.DataFrame:
    df = compute_trans_effect_magnitude()
    # condition one-hot
    for c in ["Rest", "Stim8hr", "Stim48hr"]:
        df[f"cond_{c}"] = (df["condition"] == c).astype(int)
    df["log_n_cells"] = np.log1p(df["n_cells_target"].astype(float))
    df["log_baseMean"] = np.log1p(df["target_baseMean"].clip(lower=0).astype(float))

    # attach split fold
    splits = json.loads(SPLITS.read_text())
    fold = {}
    for g in splits["train_genes"]:
        fold[g] = "train"
    for g in splits["calib_genes"]:
        fold[g] = "calib"
    for g in splits["test_genes"]:
        fold[g] = "test"
    df["fold"] = df["gene"].map(fold)

    # ---- leakage assertion: no forbidden feature columns present ----
    feature_cols = [
        "log_baseMean",
        "log_n_cells",
        "cond_Rest",
        "cond_Stim8hr",
        "cond_Stim48hr",
    ]
    assert not (set(feature_cols) & FORBIDDEN_FEATURES), (
        "LEAK: forbidden feature in matrix"
    )
    df.attrs["feature_cols"] = feature_cols
    return df


def feature_cols() -> list[str]:
    return ["log_baseMean", "log_n_cells", "cond_Rest", "cond_Stim8hr", "cond_Stim48hr"]


if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df = build_feature_table()
    df.to_parquet(OUT)
    print(f"feature table: {df.shape}  -> {OUT}")
    print(df[["gene", "condition", "trans_effect_magnitude", "y", "fold"]].head())
    print("\nfold counts (rows):", df["fold"].value_counts().to_dict())
    print(
        "target y: mean=%.3f std=%.3f min=%.3f max=%.3f"
        % (df["y"].mean(), df["y"].std(), df["y"].min(), df["y"].max())
    )
