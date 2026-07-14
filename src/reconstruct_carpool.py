"""Reconstruct CaRPool-seq (GSE213957) pseudobulk singles+doubles -> synergy table.

Follows data/gold/scc_carpool_replication_prereg.json (sha 9d48478e) exactly:
  - RNA (GEX) channel only; CP10K -> log1p per cell.
  - Cells assigned to perturbation groups from the QC'd metadata (GenePair/Guide.Class).
  - HVG (top-2000) selected on NT control cells ONLY (leakage guard).
  - effect e_g = pseudobulk(g) - pseudobulk(NT), on HVG.
  - For each testable double A+B (both singles present):
      true_mag = ||e_AB||, additive_mag = ||e_A + e_B||, sum_mag = ||e_A|| + ||e_B||.
  - min_cells_per_group = 25.

Writes:
  data/interim/carpool_synergy_table.parquet  (pair,a,b,true_mag,additive_mag,synergy_mag,synergy_frac,cos_true_add)
  data/interim/carpool_synergy_feats.parquet  (pair,cos_singles,abs_cos,min_mag,max_mag,sum_mag,mag_ratio,synergy_mag,synergy_frac)
  data/interim/carpool_synergy_index.json
"""

from __future__ import annotations
import gzip
import json
import os

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
DDIR = "data/external/carpool"
LANES = [1, 2, 3, 4]
N_HVG = 2000
MIN_CELLS = 25


def load_gex():
    """Load & concatenate the 4 GEX lanes into a single (cells x genes) CSR matrix,
    with lane-prefixed barcodes matching the metadata index (L{lane}_{bc})."""
    mats, barcodes = [], []
    genes = None
    for lane in LANES:
        base = f"{DDIR}/THP1-CaRPool-seq_and_HEK293FTstabRNA.GEXGDO{lane}"
        with gzip.open(f"{base}.matrix.mtx.gz", "rb") as f:
            m = mmread(f).tocsc()  # genes x cells
        with gzip.open(f"{base}.barcodes.tsv.gz", "rt") as f:
            bc = [f"L{lane}_{ln.strip()}" for ln in f]
        if genes is None:
            with gzip.open(f"{base}.features.tsv.gz", "rt") as f:
                genes = [ln.split("\t")[1] for ln in f]  # gene symbol
        mats.append(m.T.tocsr())  # cells x genes
        barcodes.extend(bc)
        print(f"  lane {lane}: {m.shape[1]} cells", flush=True)
    X = sparse.vstack(mats).tocsr()
    return X, np.array(barcodes), np.array(genes)


def main():
    print("loading GEX lanes...", flush=True)
    X, barcodes, genes = load_gex()
    print(f"  total {X.shape[0]} cells x {X.shape[1]} genes", flush=True)

    meta = pd.read_csv(f"{DDIR}/meta.tsv.gz", sep="\t", index_col=0)
    # intersect metadata (QC'd) cells with GEX barcodes
    bc_to_row = {b: i for i, b in enumerate(barcodes)}
    keep = [b in bc_to_row for b in meta.index]
    meta = meta[keep]
    rows = np.array([bc_to_row[b] for b in meta.index])
    Xm = X[rows]
    print(f"  matched {Xm.shape[0]} QC cells to GEX", flush=True)

    # CP10K -> log1p
    lib = np.asarray(Xm.sum(axis=1)).ravel()
    lib = np.where(lib > 0, lib, 1.0)
    Xcp = Xm.multiply(1e4 / lib[:, None]).tocsr()
    Xcp.data = np.log1p(Xcp.data)

    gclass = meta["Guide.Class"].astype(str).values
    gpair = meta["GenePair"].astype(str).values

    # --- control cells & HVG on controls only (leakage guard) ---
    ctrl_mask = gclass == "NT"
    print(f"  n control cells: {ctrl_mask.sum()}", flush=True)
    Xctrl = Xcp[ctrl_mask]
    ctrl_mean = np.asarray(Xctrl.mean(axis=0)).ravel()
    # variance on controls
    ctrl_sq = np.asarray(Xctrl.multiply(Xctrl).mean(axis=0)).ravel()
    ctrl_var = np.maximum(ctrl_sq - ctrl_mean**2, 0.0)
    hvg = np.argsort(ctrl_var)[::-1][:N_HVG]
    hvg = np.sort(hvg)
    print(f"  selected {len(hvg)} HVG on control cells only", flush=True)

    def group_effect(mask):
        if mask.sum() < MIN_CELLS:
            return None
        pb = np.asarray(Xcp[mask].mean(axis=0)).ravel()[hvg]
        return pb - ctrl_mean[hvg]

    # --- single-gene effects ---
    single_effect = {}
    sing_mask = gclass == "Single"
    for gp in np.unique(gpair[sing_mask]):
        parts = gp.split("_")
        g = [p for p in parts if p != "NT"]
        if not g:
            continue
        gene = g[0]
        e = group_effect(sing_mask & (gpair == gp))
        if e is not None:
            single_effect[gene] = e
    print(
        f"  built {len(single_effect)} single-gene effects (>= {MIN_CELLS} cells)",
        flush=True,
    )

    # --- doubles ---
    dual_mask = gclass == "Dual"
    rows_tab, rows_feat = [], []
    n_dropped_cells, n_dropped_single = 0, 0
    for gp in np.unique(gpair[dual_mask]):
        parts = gp.split("_")
        if len(parts) != 2:
            continue
        a, b = parts
        if a not in single_effect or b not in single_effect:
            n_dropped_single += 1
            continue
        e_ab = group_effect(dual_mask & (gpair == gp))
        if e_ab is None:
            n_dropped_cells += 1
            continue
        e_a, e_b = single_effect[a], single_effect[b]
        additive = e_a + e_b
        true_mag = float(np.linalg.norm(e_ab))
        additive_mag = float(np.linalg.norm(additive))
        synergy_vec = e_ab - additive
        synergy_mag = float(np.linalg.norm(synergy_vec))
        synergy_frac = synergy_mag / true_mag if true_mag > 0 else np.nan
        cos_true_add = float(
            np.dot(e_ab, additive)
            / (np.linalg.norm(e_ab) * np.linalg.norm(additive) + 1e-12)
        )
        mag_a, mag_b = float(np.linalg.norm(e_a)), float(np.linalg.norm(e_b))
        cos_singles = float(np.dot(e_a, e_b) / (mag_a * mag_b + 1e-12))
        pair = f"{a}+{b}"
        rows_tab.append(
            dict(
                pair=pair,
                a=a,
                b=b,
                true_mag=true_mag,
                additive_mag=additive_mag,
                synergy_mag=synergy_mag,
                synergy_frac=synergy_frac,
                cos_true_add=cos_true_add,
            )
        )
        rows_feat.append(
            dict(
                pair=pair,
                cos_singles=cos_singles,
                abs_cos=abs(cos_singles),
                min_mag=min(mag_a, mag_b),
                max_mag=max(mag_a, mag_b),
                sum_mag=mag_a + mag_b,
                mag_ratio=min(mag_a, mag_b) / (max(mag_a, mag_b) + 1e-12),
                synergy_mag=synergy_mag,
                synergy_frac=synergy_frac,
            )
        )

    tab = pd.DataFrame(rows_tab)
    feat = pd.DataFrame(rows_feat)
    os.makedirs("data/interim", exist_ok=True)
    tab.to_parquet("data/interim/carpool_synergy_table.parquet")
    feat.to_parquet("data/interim/carpool_synergy_feats.parquet")

    idx = {
        "n_ctrl": int(ctrl_mask.sum()),
        "n_single_genes": len(single_effect),
        "n_testable_doubles": len(tab),
        "n_dropped_missing_single": int(n_dropped_single),
        "n_dropped_low_cells": int(n_dropped_cells),
        "n_hvg": int(len(hvg)),
        "min_cells_per_group": MIN_CELLS,
        "single_genes": sorted(single_effect.keys()),
    }
    json.dump(idx, open("data/interim/carpool_synergy_index.json", "w"), indent=2)
    print(json.dumps(idx, indent=2))
    print(
        f"\nsynergy_frac range: {tab['synergy_frac'].min():.3f} - {tab['synergy_frac'].max():.3f}"
    )
    print(f"testable doubles: {len(tab)}")


if __name__ == "__main__":
    main()
