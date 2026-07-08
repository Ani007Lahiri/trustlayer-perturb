"""
Day-5: derive REAL per-gene trust scores for the trio from the conformal machinery
(replaces the Day-1 placeholder). Critique-corrected: provenance stamped per gene,
and the runner asserts trust is NON-BINDING for the demo trio.

trust(gene) = 1 - normalized(mean ensemble spread over that gene's test rows), in [0,1].
  - RASGRP1: donor-blocked spread (present in all 6 donor pairs) = guaranteed estimator.
  - CD226 / PRKCQ: pooled-proxy spread (absent/sparse in donor-pair data) = different scale.

Normalization: spread is min-max scaled against the full pooled-test spread distribution
so trust is on a 0..1 scale; the SAME normalizer is used for all pooled-proxy genes.
RASGRP1's donor-blocked spread is normalized against the donor-blocked spread distribution.
Both normalizers + the raw spread are recorded so the provenance is auditable.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

PREDS = Path("data/interim/base_predictor_preds.parquet")
FEATURE_TABLE = Path("data/interim/day2_feature_table.parquet")
SEED = 20260708


def _pooled_ensemble_spread_by_gene(genes: list[str], n_models: int = 8):
    """Fit a small diverse bootstrap ensemble on the pooled TRAIN fold, predict on all
    rows, return per-row spread + the spread distribution (for normalization)."""
    df = pd.read_parquet(FEATURE_TABLE)
    fc = ["log_baseMean", "log_n_cells", "cond_Rest", "cond_Stim8hr", "cond_Stim48hr"]
    tr = df[df.fold == "train"]
    Xtr, ytr = tr[fc].values, tr["y"].values
    n = len(ytr)
    preds = []
    for s in range(n_models):
        rng = np.random.default_rng(SEED + s)
        idx = rng.integers(0, n, n)
        m = HistGradientBoostingRegressor(
            max_iter=200,
            learning_rate=0.05,
            max_depth=4,
            l2_regularization=1.0,
            random_state=SEED + s,
        ).fit(Xtr[idx], ytr[idx])
        preds.append(m.predict(df[fc].values))
    spread_all = np.std(np.vstack(preds), axis=0)
    df = df.assign(spread=spread_all)
    lo, hi = float(spread_all.min()), float(spread_all.max())
    per_gene = {}
    for g in genes:
        sub = df[df.gene == g]
        if len(sub):
            per_gene[g] = float(sub["spread"].mean())
    return per_gene, (lo, hi)


def _donor_blocked_spread_for_gene(gene: str):
    """Mean donor-blocked ensemble spread for a gene present in all donor pairs.
    Returns (mean_spread, (lo,hi) normalizer over the pooled donor-blocked test spread)."""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from trustlayer.conformal import (
        H5MU,
        DONORS,
        pairs_excluding,
        pairs_containing,
        _stack,
        _ensemble_spread,
    )
    import h5py

    gene_spreads = []
    all_spreads = []
    with h5py.File(H5MU, "r") as f:
        for d in DONORS:
            calib = _stack(f, pairs_excluding(d))
            test = _stack(f, pairs_containing(d))
            fold_seed = SEED + hash(d) % 100000
            spreads = _ensemble_spread(calib, test, n_models=8, base_seed=fold_seed)
            all_spreads.append(spreads)
            mask = test["gene"] == gene
            if mask.any():
                gene_spreads.extend(spreads[mask].tolist())
    alls = np.concatenate(all_spreads)
    lo, hi = float(alls.min()), float(alls.max())
    return (float(np.mean(gene_spreads)) if gene_spreads else None), (lo, hi)


def _norm(x, lo, hi):
    if hi <= lo:
        return 0.5
    return float(np.clip((x - lo) / (hi - lo), 0, 1))


def derive_trust() -> dict:
    """Return per-gene trust for the trio with full provenance."""
    pooled, pool_norm = _pooled_ensemble_spread_by_gene(["CD226", "RASGRP1", "PRKCQ"])
    out = {}
    # CD226 + PRKCQ: pooled proxy
    for g in ["CD226", "PRKCQ"]:
        sp = pooled.get(g)
        trust = 1.0 - _norm(sp, *pool_norm) if sp is not None else None
        out[g] = {
            "trust": round(trust, 4) if trust is not None else None,
            "raw_spread": round(sp, 5) if sp is not None else None,
            "normalizer": [round(pool_norm[0], 5), round(pool_norm[1], 5)],
            "provenance": "pooled-proxy (absent/sparse in donor-pair data)",
        }
    # RASGRP1: donor-blocked (present in all 6 pairs)
    db_sp, db_norm = _donor_blocked_spread_for_gene("RASGRP1")
    trust = 1.0 - _norm(db_sp, *db_norm) if db_sp is not None else None
    out["RASGRP1"] = {
        "trust": round(trust, 4) if trust is not None else None,
        "raw_spread": round(db_sp, 5) if db_sp is not None else None,
        "normalizer": [round(db_norm[0], 5), round(db_norm[1], 5)],
        "provenance": "donor-blocked (present in all 6 donor pairs)",
    }
    return out


if __name__ == "__main__":
    import json

    t = derive_trust()
    print(json.dumps(t, indent=2))
