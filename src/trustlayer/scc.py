"""Synergy-Conditioned Conformal (SCC) — reusable core engine.

Rebuilt from the frozen pre-registrations (data/gold/scc_prereg.json,
data/gold/scc_powered_prereg.json) after the original ad-hoc runner was found to
have never been persisted to disk. This module is dataset-agnostic: it operates on
a synergy table (one row per DOUBLE perturbation) and reproduces the frozen Norman
numbers bit-for-bit (see tests/test_scc_engine.py).

Method (verbatim from prereg sha 7c97f9d6 / 77bcaea9):
  - For each double A+B: additive baseline effect = e_A + e_B (singles only);
    true effect = e_AB (the measured double). Residual = |true_mag - additive_mag|.
  - Leakage-safe synergy prior s_hat = sum of the two single-effect magnitudes
    (SINGLES ONLY, never the double), standardized (z-score).
  - LOPO (leave-one-pair-out): 131 folds. Each double held out once; conformal
    calibration on the other 130 (60/40 predictor/calib split within them).
    Interval = additive_pred +/- width. width constant (vanilla) or sigma(s_hat)*Q (SCC).
  - Coverage indicator per pair (1 = true_mag in interval).
  - Primary estimand: logistic cov ~ method + synergy + method:synergy; the
    interaction coeff (SCC vs vanilla slope of coverage-vs-synergy) > 0.
  - Inference: permutation of the synergy covariate (B=2000), one-sided p; and
    pair-level cluster bootstrap 95% CI on the interaction coeff.

The SCC sigma(s_hat) is an isotonic fit of |residual| on s_hat over the calibration
pairs (higher expected synergy -> wider interval). Vanilla uses a constant width =
the (1-alpha) empirical quantile of calibration residuals. Both are re-scaled so the
marginal calibration matches at the fold level (standard split-conformal quantile).
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from sklearn.isotonic import IsotonicRegression

# statsmodels is used for the logistic interaction test; imported lazily so the
# reconstruction/feature code has no hard dependency on it.


def standardize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    mu, sd = x.mean(), x.std()
    return (x - mu) / (sd if sd > 0 else 1.0)


@dataclass
class SCCConfig:
    alpha: float = 0.10  # nominal miscoverage -> 90% intervals
    calib_frac: float = 0.40  # 60/40 predictor/calib split within the 130 train pairs
    B_perm: int = 2000
    B_boot: int = 2000
    seed: int = 0
    high_split: float = (
        0.50  # pre-registered top-50% (median) split for the composite bar
    )
    nominal: float = field(init=False)

    def __post_init__(self):
        self.nominal = 1.0 - self.alpha


def _fold_widths(
    resid_cal: np.ndarray, shat_cal: np.ndarray, shat_test: float, alpha: float
):
    """Return (vanilla_width, scc_width) for a held-out test pair.

    vanilla: constant conformal width = (1-alpha) empirical quantile of |resid|.
    scc: isotonic sigma(s_hat) normalizes residuals; width = sigma(s_hat_test)*Q,
         where Q = (1-alpha) quantile of normalized scores |resid|/sigma(s_hat_cal).
    """
    n = len(resid_cal)
    # conformal quantile level with finite-sample correction
    q_level = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)

    # vanilla
    w_vanilla = np.quantile(resid_cal, q_level, method="higher")

    # scc: isotonic sigma(s_hat) fit on calibration residual magnitudes
    iso = IsotonicRegression(increasing=True, out_of_bounds="clip")
    # floor sigma to avoid divide-by-zero / degenerate normalization
    sig_cal = iso.fit(shat_cal, resid_cal).predict(shat_cal)
    sig_floor = max(1e-6, np.median(resid_cal) * 1e-3)
    sig_cal = np.maximum(sig_cal, sig_floor)
    norm_scores = resid_cal / sig_cal
    Q = np.quantile(norm_scores, q_level, method="higher")
    sig_test = max(sig_floor, float(iso.predict([shat_test])[0]))
    w_scc = sig_test * Q
    return float(w_vanilla), float(w_scc)


def run_lopo(table: pd.DataFrame, cfg: SCCConfig):
    """Leave-one-pair-out conformal coverage for vanilla vs SCC.

    `table` must have columns: true_mag, additive_mag, sum_mag (per double).
    Returns dict with per-pair shat, covV, covS arrays and marginal coverages.
    """
    rng = np.random.default_rng(cfg.seed)
    true_mag = table["true_mag"].to_numpy(float)
    add_mag = table["additive_mag"].to_numpy(float)
    resid_all = np.abs(true_mag - add_mag)
    shat = standardize(table["sum_mag"].to_numpy(float))
    n = len(table)

    covV = np.zeros(n, dtype=int)
    covS = np.zeros(n, dtype=int)

    for i in range(n):
        train_idx = np.array([j for j in range(n) if j != i])
        # 60/40 predictor/calib split within the 130 remaining pairs.
        # predictor set is not used to refit the additive baseline (it is a fixed
        # biological null e_A+e_B), so it only reserves rows out of calibration to
        # honor the prereg's split; calib = the 40% that set the conformal quantile.
        perm = rng.permutation(len(train_idx))
        n_cal = max(2, int(round(cfg.calib_frac * len(train_idx))))
        cal_idx = train_idx[perm[:n_cal]]

        resid_cal = resid_all[cal_idx]
        shat_cal = shat[cal_idx]
        w_v, w_s = _fold_widths(resid_cal, shat_cal, shat[i], cfg.alpha)

        # coverage: is the measured true_mag within [additive +/- width]?
        # residual = |true - additive|, so covered iff residual <= width.
        covV[i] = int(resid_all[i] <= w_v)
        covS[i] = int(resid_all[i] <= w_s)

    return {
        "shat": shat,
        "covV": covV,
        "covS": covS,
        "cov_vanilla": float(covV.mean()),
        "cov_scc": float(covS.mean()),
        "nominal": cfg.nominal,
        "n": n,
    }


def _logit_interaction(cov: np.ndarray, method: np.ndarray, syn: np.ndarray) -> float:
    """Interaction coeff of logistic cov ~ const + method + synergy + method:synergy.

    Uses a tiny L2 (ridge) regularization so near-separated small-n folds converge
    stably; this only shrinks coefficients toward 0 (conservative for a >0 test).
    """
    import warnings
    import statsmodels.api as sm

    X = np.column_stack([np.ones_like(method, float), method, syn, method * syn])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            res = sm.Logit(cov, X).fit_regularized(alpha=1e-4, L1_wt=0.0, disp=0)
        except Exception:
            res = sm.Logit(cov, X).fit(disp=0, method="bfgs", maxiter=200)
    return float(res.params[3])


def interaction_test(lopo: dict, cfg: SCCConfig):
    """Primary estimand: SCC-vs-vanilla coverage-vs-synergy interaction coeff,
    with synergy-permutation one-sided p and pair-level cluster-bootstrap 95% CI."""
    shat = lopo["shat"]
    covV = lopo["covV"]
    covS = lopo["covS"]
    n = len(shat)

    # stacked long form
    method = np.r_[np.zeros(n), np.ones(n)]
    syn = np.r_[shat, shat]
    cov = np.r_[covV, covS]
    obs = _logit_interaction(cov, method, syn)

    # permutation of the synergy covariate (same shuffle applied to both methods,
    # preserving the paired structure — shuffles only what synergy labels the pair).
    rng = np.random.default_rng(cfg.seed + 1)
    ge = 1  # +1 for observed (one-sided, >=)
    for _ in range(cfg.B_perm):
        pi = rng.permutation(n)
        syn_p = np.r_[shat[pi], shat[pi]]
        c = _logit_interaction(cov, method, syn_p)
        if c >= obs:
            ge += 1
    perm_p = ge / (cfg.B_perm + 1)

    # pair-level cluster bootstrap: resample PAIRS (both method rows move together)
    rng2 = np.random.default_rng(cfg.seed + 2)
    boots = []
    for _ in range(cfg.B_boot):
        samp = rng2.integers(0, n, size=n)
        method_b = np.r_[np.zeros(n), np.ones(n)]
        syn_b = np.r_[shat[samp], shat[samp]]
        cov_b = np.r_[covV[samp], covS[samp]]
        try:
            boots.append(_logit_interaction(cov_b, method_b, syn_b))
        except Exception:
            continue
    boots = np.array(boots)
    ci = [float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))]

    return {
        "interaction_coeff": obs,
        "perm_p_one_sided": perm_p,
        "boot_ci": ci,
        "B_perm": cfg.B_perm,
        "B_boot": cfg.B_boot,
    }


def composite_clearing(lopo: dict, cfg: SCCConfig, itest: dict):
    """Frozen composite within-tier clearing criterion (prereg 77bcaea9):
      c1: interaction coeff > 0 at one-sided perm p < 0.05
      c2: vanilla high-synergy (top-50%) coverage CI excludes nominal from below
      c3: SCC high-synergy coverage CI contains nominal
    Clears iff c1 AND c2 AND c3.
    """
    shat = lopo["shat"]
    covV = lopo["covV"]
    covS = lopo["covS"]
    thr = np.quantile(shat, 1 - cfg.high_split)
    hi = shat >= thr
    lo = ~hi

    def wilson_ci(k, nn, z=1.96):
        if nn == 0:
            return [np.nan, np.nan, np.nan]
        p = k / nn
        denom = 1 + z * z / nn
        centre = (p + z * z / (2 * nn)) / denom
        half = z * np.sqrt(p * (1 - p) / nn + z * z / (4 * nn * nn)) / denom
        return [float(p), float(centre - half), float(centre + half)]

    vH = wilson_ci(covV[hi].sum(), hi.sum())
    sH = wilson_ci(covS[hi].sum(), hi.sum())
    vL = wilson_ci(covV[lo].sum(), lo.sum())
    sL = wilson_ci(covS[lo].sum(), lo.sum())

    nominal = cfg.nominal
    c1 = (itest["interaction_coeff"] > 0) and (itest["perm_p_one_sided"] < 0.05)
    c2 = vH[2] < nominal  # vanilla high-syn upper CI below nominal
    c3 = sH[1] <= nominal <= sH[2]  # SCC high-syn CI contains nominal
    return {
        "vanilla_high": vH,
        "scc_high": sH,
        "vanilla_low": vL,
        "scc_low": sL,
        "c1": bool(c1),
        "c2": bool(c2),
        "c3": bool(c3),
        "clears_within_tier": bool(c1 and c2 and c3),
        "high_split": cfg.high_split,
        "nominal": nominal,
    }
