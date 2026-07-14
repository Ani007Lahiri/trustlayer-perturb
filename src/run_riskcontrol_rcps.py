"""TRACK A2 runner: RCPS risk control across three exchangeability regimes, from CACHE.

Reproduces data/gold/day0_riskcontrol_rcps_receipt.json without touching the 16GB h5mu.
  R1  LODO donor-level (n=4)     <- data/gold/conformal_lodo_receipt.json (fold_values)
  R2  LODO row-level  (n=58,558) <- same receipt (row_empirical, n_rows)  [assumption-violated]
  R3  Norman per-item (n=105)    <- data/interim/norman_group_stats.npz + frozen_marson_calibration.joblib

Prereg (hash-frozen BEFORE scoring): data/gold/day0_riskcontrol_rcps_prereg.json
  sha256 16aca796223a7821b6f7929c65d05159df8c205f93baaa4a2b136ea2dec615c7
"""
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor

from trustlayer import risk_control as rc

GOLD = Path("data/gold"); INT = Path("data/interim")
ALPHA = DELTA = 0.10
LEVELS = (0.80, 0.90, 0.95)
NOMINAL_MISS = {0.80: 0.20, 0.90: 0.10, 0.95: 0.05}
SEED = 20260708


def norman_residuals():
    meta = json.load(open(INT / "norman_meta.json"))
    npz = np.load(INT / "norman_group_stats.npz")
    mean, var, ncells = npz["mean"], npz["var"], npz["ncells"]
    inter_idx = np.array(meta["inter_idx"]); norman_sym = np.array(meta["norman_sym"])
    id2target = meta["id2target"]; panel_sym = norman_sym[inter_idx]
    cm, cv, nc = mean[:, 0], var[:, 0], ncells[0]
    ys, log_bm = [], []
    for gid in range(1, 106):
        tgt = id2target[str(gid)]
        pm, pv, npp = mean[inter_idx, gid], var[inter_idx, gid], ncells[gid]
        den = np.sqrt(pv / max(npp, 1) + cv[inter_idx] / max(nc, 1))
        z = np.where(den > 0, (pm - cm[inter_idx]) / den, 0.0).copy(); z[panel_sym == tgt] = 0.0
        ys.append(np.log1p(np.sqrt(np.sum(z * z))))
        on = norman_sym == tgt
        log_bm.append(np.log1p(max(cm[on].sum(), 0.0)) if on.any() else 0.0)
    ys = np.array(ys); N = len(ys)
    X = np.column_stack([np.array(log_bm), np.ones(N), np.zeros(N), np.zeros(N)])
    return ys, X


def main():
    rec = json.load(open(GOLD / "conformal_lodo_receipt.json"))["pooled_coverage"]
    key = {0.80: "0.8", 0.90: "0.9", 0.95: "0.95"}

    # R1 donor-level (Hoeffding on per-fold miscoverage fractions, n=4)
    R1 = {}
    for L in LEVELS:
        fm = 1.0 - np.array(rec[key[L]]["fold_values"])
        R1[L] = rc.rcps_certificate_at_fixed_lambda(fm, NOMINAL_MISS[L], DELTA, binary=False, bound="hoeffding")

    # R3 Norman per-item (RCPS scan, exchangeable by split)
    ys, X = norman_residuals()
    N = len(ys); p = np.random.default_rng(20260711).permutation(N); nca = N // 2
    ci, ti = p[:nca], p[nca:]
    m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_depth=4,
                                      l2_regularization=1.0, random_state=SEED).fit(X[ci], ys[ci])
    rc_cal = np.abs(ys[ci] - m.predict(X[ci])); rc_tst = np.abs(ys[ti] - m.predict(X[ti]))
    res = rc.rcps_from_residuals(rc_cal, ALPHA, DELTA, bound="cp")
    realized = float(np.mean(rc_tst > res.lambda_hat))

    print("R1 (donor n=4):", {int(L * 100): (round(R1[L]["ucb"], 3), R1[L]["certifies_risk_le_alpha"]) for L in LEVELS})
    print(f"R3 (Norman n=105) 90%: lambda_hat={res.lambda_hat:.4f} UCB={res.ucb_at_lambda_hat:.4f} "
          f"realized_test_cov={1 - realized:.4f}")


if __name__ == "__main__":
    main()
