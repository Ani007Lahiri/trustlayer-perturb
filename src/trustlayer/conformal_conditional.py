"""
A3 fix: principled CONDITIONAL-coverage conformal methods.

The disclosed break (project A3): under donor-blocked LODO split-conformal
(conformal.py), MARGINAL coverage is valid but CONDITIONAL coverage on the
largest TRUE-effect-magnitude stratum collapses (~0.60-0.68 at nominal 0.90).
Mechanism: the leakage-safe base predictor is near-constant (corr(pred,y)~0.07
in the modest regime), so |residual| ~ |y - const| and the single global
quantile under-covers BOTH effect-magnitude tails while over-covering the
middle. A globally-valid interval is silently mis-calibrated where effects
are biggest.

This module implements four distribution-free fixes, each conditioning the
interval width on something other than a single global quantile. All operate
on the COMMITTED residuals/predictions (same base predictor, same folds); only
the calibration step changes.

Methods
-------
baseline        : global split-conformal (single per-fold |resid| quantile).
                  Reproduces conformal.py.
mondrian_pred   : (a) DEPLOYABLE effect-magnitude Mondrian. Vovk's Mondrian /
                  class-conditional conformal (Vovk et al. 2003; Vovk 2012):
                  partition by a taxonomy known at test time -- here bins of the
                  PREDICTED magnitude -- and take a separate calibration quantile
                  per bin. Guarantees coverage conditional on the predicted-mag
                  bin. Only helps the TRUE-effect stratum insofar as pred tracks y.
mondrian_oracle : (a') ORACLE effect-magnitude Mondrian. Bins by the TRUE y.
                  Uses the label -> NOT deployable. Reported as the upper bound
                  ("ceiling") of bin-conditional conformal for this break.
normalized      : (b) normalized / localized conformal (Papadopoulos et al. 2008;
                  Lei et al. 2018, JASA). Nonconformity s = |y - pred| / sigma(x),
                  with sigma(x) a difficulty estimate (regressor fit on calibration
                  |resid|). Interval half-width = Q_s * sigma(x_test), so the band
                  widens where the model expects to be less certain. Distribution-
                  free marginal validity; adapts iff sigma(x) carries signal.
cqr             : (c) Conformalized Quantile Regression (Romano, Patterson, Candes,
                  NeurIPS 2019). Fit conditional-quantile regressors q_lo, q_hi at
                  alpha/2, 1-alpha/2; conformity E = max(q_lo(x)-y, y-q_hi(x));
                  Q = ceil((n+1)(1-alpha))/n empirical quantile of E; interval
                  [q_lo(x)-Q, q_hi(x)+Q]. Adapts the band shape to the features.

References
----------
Vovk, Lindsay, Nouretdinov, Gammerman (2003); Vovk (2012) "Conditional validity
  of inductive conformal predictors." (Mondrian CP)
Papadopoulos, Gammerman, Vovk (2008) "Normalized nonconformity measures."
Lei, G'Sell, Rinaldo, Tibshirani, Wasserman (2018) "Distribution-free predictive
  inference for regression", JASA. (locally-weighted / normalized conformal)
Romano, Patterson, Candes (2019) "Conformalized Quantile Regression", NeurIPS.
"""

from __future__ import annotations

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor


# ----------------------------------------------------------------- helpers
def conformal_quantile_level(n: int, level: float) -> float:
    """Finite-sample conformal rank as a fractional level for np.quantile(method='higher').

    Split-conformal covers with the ceil((n+1)*level)-th smallest score. We pass
    `level` to np.quantile(..., method='higher'), matching conformal.py's committed
    call exactly (so baseline reproduces the headline byte-for-byte)."""
    return level


def _bin_edges(x: np.ndarray, k: int) -> np.ndarray:
    """Interior quantile edges; deduplicated so near-constant x degrades gracefully."""
    e = np.quantile(x, np.linspace(0, 1, k + 1))
    return np.unique(e[1:-1])


def _digitize(x: np.ndarray, interior_edges: np.ndarray) -> np.ndarray:
    return np.digitize(x, interior_edges)


# ----------------------------------------------------------------- methods
def hw_baseline(calib_resid, calib_pred, test_pred, level, **_):
    """Global split-conformal: one quantile for every test row."""
    q = float(np.quantile(calib_resid, level, method="higher"))
    return np.full(len(test_pred), q)


def _mondrian_hw(calib_resid, calib_binvar, test_binvar, level, k):
    """Per-bin split-conformal quantile; test row uses its own bin's quantile.

    Bins are quantile bins of `*_binvar`. A bin too thin to support the level
    (or unseen at test) falls back to +inf (always cover) so we never silently
    under-cover a thin stratum."""
    edges = _bin_edges(calib_binvar, k)
    cb = _digitize(calib_binvar, edges)
    tb = _digitize(test_binvar, edges)
    q_by_bin = {}
    for b in np.unique(cb):
        r = calib_resid[cb == b]
        n = len(r)
        q_by_bin[b] = float(np.quantile(r, level, method="higher")) if (n + 1) * level <= n else np.inf
    return np.array([q_by_bin.get(b, np.inf) for b in tb])


def hw_mondrian_pred(calib_resid, calib_pred, test_pred, level, k=5, **_):
    """(a) DEPLOYABLE: bins by PREDICTED magnitude (known at test time)."""
    return _mondrian_hw(calib_resid, calib_pred, test_pred, level, k)


def hw_mondrian_oracle(calib_resid, calib_pred, test_pred, level, calib_y=None, test_y=None, k=5, **_):
    """(a') ORACLE: bins by TRUE y (uses the label; NON-deployable ceiling)."""
    return _mondrian_hw(calib_resid, calib_y, test_y, level, k)


def hw_normalized(calib_resid, calib_pred, test_pred, level, calib_X=None, test_X=None,
                  seed=20260711, **_):
    """(b) normalized/localized conformal. sigma(x) = difficulty regressor on |resid|.

    Half-width for test row j = Q_s * sigma(x_j), Q_s = quantile of s=|resid|/sigma
    on calibration. sigma floored at a small positive value for numerical safety."""
    diff = HistGradientBoostingRegressor(
        max_iter=300, learning_rate=0.05, max_depth=4, l2_regularization=1.0, random_state=seed
    ).fit(calib_X, calib_resid)
    floor = max(1e-3, 0.05 * float(np.mean(calib_resid)))
    sig_c = np.maximum(diff.predict(calib_X), floor)
    sig_t = np.maximum(diff.predict(test_X), floor)
    s = calib_resid / sig_c
    q_s = float(np.quantile(s, level, method="higher"))
    return q_s * sig_t


def cqr_intervals(calib_X, calib_y, test_X, level, seed=20260711):
    """(c) CQR (Romano et al. 2019). Returns (lo, hi) test interval bounds.

    q_lo, q_hi are pinball-loss GBM quantile regressors at alpha/2, 1-alpha/2;
    E = max(q_lo(x)-y, y-q_hi(x)); Q = (1-alpha) quantile of E (method='higher');
    interval = [q_lo(x)-Q, q_hi(x)+Q]. Distribution-free marginal validity."""
    alpha = 1.0 - level
    lo_q, hi_q = alpha / 2.0, 1.0 - alpha / 2.0
    common = dict(loss="quantile", max_iter=300, learning_rate=0.05, max_depth=4,
                  l2_regularization=1.0, random_state=seed)
    m_lo = HistGradientBoostingRegressor(quantile=lo_q, **common).fit(calib_X, calib_y)
    m_hi = HistGradientBoostingRegressor(quantile=hi_q, **common).fit(calib_X, calib_y)
    c_lo, c_hi = m_lo.predict(calib_X), m_hi.predict(calib_X)
    E = np.maximum(c_lo - calib_y, calib_y - c_hi)
    Q = float(np.quantile(E, level, method="higher"))
    t_lo, t_hi = m_lo.predict(test_X) - Q, m_hi.predict(test_X) + Q
    return t_lo, t_hi


# ----------------------------------------------------------------- scoring
def covered_and_halfwidth(method, fold, level, k=5, seed=20260711):
    """Return (covered_bool[n_test], halfwidth[n_test]) for a method on one fold.

    `fold` = dict with calib/test arrays: X,y,pred,resid (calib) and X,y,pred (test).
    Symmetric methods define covered = |y_test-pred_test| <= hw. CQR defines covered
    from its two-sided interval and reports an effective half-width = width/2."""
    c, t = fold["calib"], fold["test"]
    test_err = np.abs(t["y"] - t["pred"])
    if method == "cqr":
        lo, hi = cqr_intervals(c["X"], c["y"], t["X"], level, seed=seed)
        covered = (t["y"] >= lo) & (t["y"] <= hi)
        hw = (hi - lo) / 2.0
        return covered, hw
    fn = {
        "baseline": hw_baseline,
        "mondrian_pred": hw_mondrian_pred,
        "mondrian_oracle": hw_mondrian_oracle,
        "normalized": hw_normalized,
    }[method]
    hw = fn(c["resid"], c["pred"], t["pred"], level,
            calib_X=c["X"], test_X=t["X"], calib_y=c["y"], test_y=t["y"], k=k, seed=seed)
    covered = test_err <= hw
    return covered, hw


DEPLOYABLE = ["baseline", "mondrian_pred", "normalized", "cqr"]
ALL_METHODS = ["baseline", "mondrian_pred", "mondrian_oracle", "normalized", "cqr"]
