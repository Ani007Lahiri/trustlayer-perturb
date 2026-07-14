"""
TRACK A2: distribution-free RISK CONTROL layer (RCPS) on the committed conformal
machinery.

Motivation
----------
The committed donor-blocked LODO split-conformal (conformal.py) reports *empirical*
coverage (fold_mean +/- fold_std). Because LODO breaks exchangeability by design
(calibration = other-donor pairs, test = held-out donor), those coverage numbers are
an OBSERVATION, not a certificate. A finite calibration sample also carries sampling
error that a point estimate hides.

RCPS (Bates, Angelopoulos, Lei, Malik, Jordan; "Distribution-Free, Risk-Controlling
Prediction Sets", JACM 2021) turns a threshold parameter lambda (here the conformal
interval half-width) into a set-valued predictor with a finite-sample, distribution-free
high-probability guarantee:

    P( R(lambda_hat) <= alpha ) >= 1 - delta

where R(lambda) = E[ L(lambda) ] is the risk of the deployed interval. We take the
loss to be the miscoverage indicator L_i(lambda) = 1[ |y_i - pred_i| > lambda ], so
R(lambda) is the true miscoverage rate and the guarantee reads: with probability
>= 1-delta over the calibration draw, the deployed interval pred +/- lambda_hat
miscovers at most a fraction alpha of the (exchangeable) test population.

Why RCPS, not LTT
-----------------
R(lambda) is a SINGLE-parameter, monotone-decreasing risk over nested prediction sets
(bigger lambda -> more coverage -> less risk). RCPS is purpose-built for that case:
scan lambda downward and stop at the smallest lambda whose pointwise upper-confidence
bound (UCB) on the risk still satisfies UCB <= alpha. Learn-then-Test would treat lambda
as one hypothesis in a family and control FWER via multiple testing -- it certifies the
SAME guarantee here through a strictly-not-tighter route. RCPS is the assumption-matched
instrument for one monotone lambda.

The honest boundary
-------------------
RCPS controls CALIBRATION-SAMPLING error under exchangeability of the calibration units
with the test unit. It does NOT repair distribution shift. Under LODO the exchangeable
unit is the DONOR (n=4); RCPS applied at that faithful granularity is honest but its
finite-sample UCB is loose. Applied at the row level it is tight but rests on exactly the
exchangeability the project flags as violated. Both are reported; only the exchangeable
regime carries a claimed guarantee.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import beta as _beta
from scipy.stats import binom as _binom


# --------------------------------------------------------------- upper conf. bounds
def cp_upper(k: int, n: int, delta: float) -> float:
    """Clopper-Pearson (exact binomial) upper confidence limit on a Bernoulli mean.

    Returns the largest p such that P(Bin(n,p) <= k) >= delta, i.e. the one-sided
    (1-delta) upper limit on the miss-rate given k misses in n trials. This is the
    EXACT and tightest valid RCPS upper bound when the loss is binary {0,1} (the
    Bentkus/binomial special case of the RCPS concentration bound).
    """
    if k >= n:
        return 1.0
    return float(_beta.ppf(1.0 - delta, k + 1, n - k))


def hoeffding_upper(rhat: float, n: int, delta: float) -> float:
    """Hoeffding upper confidence bound, valid for ANY loss bounded in [0,1].

    Used when the per-unit loss is a FRACTION in [0,1] (e.g. a donor fold's miscoverage
    rate), not a single Bernoulli draw. Distribution-free, worst-case, hence loose at
    small n -- which is the honest point in the n=4 donor regime.
    """
    return float(min(1.0, rhat + np.sqrt(np.log(1.0 / delta) / (2.0 * n))))


def hoeffding_bentkus_upper(rhat: float, n: int, delta: float) -> float:
    """Hoeffding-Bentkus UCB from Bates et al. 2021 (their eq. for the p-value g).

    R+ = sup{ R in [Rhat,1] : g(Rhat, R) >= delta }, where
        g(Rhat,R) = min( exp(-n h1(Rhat,R)) ,  e * P(Bin(n,R) <= ceil(n Rhat)) )
    and h1(a,b) = a log(a/b) + (1-a) log((1-a)/(1-b)) is the Bernoulli KL.
    Tighter than plain Hoeffding, especially for small risks. Reported as a cross-check.
    """
    if rhat >= 1.0:
        return 1.0
    k = int(np.ceil(n * rhat))

    def h1(a: float, b: float) -> float:
        if b >= 1.0:
            return np.inf if a < 1.0 else 0.0
        if b <= 0.0:
            return np.inf if a > 0.0 else 0.0
        if a <= 0:
            return -np.log(1 - b)
        if a >= 1:
            return -np.log(b)
        return a * np.log(a / b) + (1 - a) * np.log((1 - a) / (1 - b))

    def g(R: float) -> float:
        if R <= rhat:
            return 1.0
        hoef = np.exp(-n * h1(rhat, R))
        bent = np.e * _binom.cdf(k, n, R)
        return min(hoef, bent)

    # largest R with g(R) >= delta, via bisection on [rhat, 1]
    lo, hi = rhat, 1.0
    if g(hi) >= delta:
        return 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) >= delta:
            lo = mid
        else:
            hi = mid
    return float(hi)


# --------------------------------------------------------------- RCPS core
@dataclass
class RCPSResult:
    alpha: float
    delta: float
    bound: str
    lambda_hat: float | None          # smallest certified half-width (None => none certifiable)
    certified: bool
    ucb_at_lambda_hat: float | None
    empirical_risk_at_lambda_hat: float | None
    lambda_grid: list = field(default_factory=list)
    ucb_curve: list = field(default_factory=list)
    emp_curve: list = field(default_factory=list)
    n_units: int = 0
    note: str = ""


def _ucb_binary(errors: np.ndarray, lam: float, delta: float, bound: str) -> tuple[float, float]:
    """Return (empirical_risk, UCB) for binary miscoverage loss at threshold lam."""
    miss = errors > lam
    k, n = int(miss.sum()), int(len(errors))
    rhat = k / n
    if bound == "cp":
        return rhat, cp_upper(k, n, delta)
    if bound == "hoeffding":
        return rhat, hoeffding_upper(rhat, n, delta)
    if bound == "hb":
        return rhat, hoeffding_bentkus_upper(rhat, n, delta)
    raise ValueError(bound)


def rcps_from_residuals(
    calib_errors: np.ndarray,
    alpha: float = 0.10,
    delta: float = 0.10,
    bound: str = "cp",
    lambda_grid: np.ndarray | None = None,
) -> RCPSResult:
    """RCPS on per-item calibration residuals (binary miscoverage loss).

    Scans lambda over a grid (default: sorted unique |residual| values, the only points
    at which the empirical risk changes). Because R(lambda) is monotone decreasing in
    lambda, lambda_hat = the SMALLEST lambda whose UCB stays <= alpha. Returns the full
    UCB curve for plotting and the certified half-width.
    """
    calib_errors = np.asarray(calib_errors, float)
    if lambda_grid is None:
        lambda_grid = np.unique(np.concatenate([[0.0], np.sort(calib_errors)]))
    emp, ucb = [], []
    for lam in lambda_grid:
        r, u = _ucb_binary(calib_errors, lam, delta, bound)
        emp.append(r)
        ucb.append(u)
    ucb = np.array(ucb)
    emp = np.array(emp)
    # smallest lambda with UCB <= alpha (grid ascending); monotone so this is the RCPS pick
    ok = np.where(ucb <= alpha)[0]
    if len(ok):
        j = ok[0]
        return RCPSResult(
            alpha, delta, bound, float(lambda_grid[j]), True,
            float(ucb[j]), float(emp[j]),
            [float(x) for x in lambda_grid], [float(x) for x in ucb],
            [float(x) for x in emp], len(calib_errors),
            "RCPS: smallest half-width whose finite-sample UCB on miscoverage <= alpha.",
        )
    return RCPSResult(
        alpha, delta, bound, None, False, None, None,
        [float(x) for x in lambda_grid], [float(x) for x in ucb],
        [float(x) for x in emp], len(calib_errors),
        "No half-width on the grid could be certified at (alpha,delta).",
    )


def rcps_certificate_at_fixed_lambda(
    per_unit_losses: np.ndarray,
    alpha: float,
    delta: float,
    binary: bool,
    bound: str = "cp",
) -> dict:
    """Certificate at a DEPLOYED lambda: given per-unit losses (binary misses, or
    per-fold miscoverage FRACTIONS in [0,1]), return the (1-delta) UCB on the risk and
    whether it certifies R <= alpha. This is the donor-fold (n=4) / pooled-row use case,
    where lambda is fixed to the committed frozen quantile and we ask what risk is
    certifiable from the available exchangeable units.
    """
    per_unit_losses = np.asarray(per_unit_losses, float)
    n = len(per_unit_losses)
    rhat = float(per_unit_losses.mean())
    if binary:
        k = int(round(per_unit_losses.sum()))
        ucb = cp_upper(k, n, delta) if bound == "cp" else (
            hoeffding_upper(rhat, n, delta) if bound == "hoeffding"
            else hoeffding_bentkus_upper(rhat, n, delta))
    else:
        # fractional [0,1] losses -> Hoeffding (or HB); CP is not valid for non-binary
        ucb = hoeffding_upper(rhat, n, delta) if bound != "hb" else hoeffding_bentkus_upper(rhat, n, delta)
    return {
        "n_units": n,
        "empirical_risk": round(rhat, 4),
        "ucb": round(float(ucb), 4),
        "alpha": alpha,
        "delta": delta,
        "bound": bound if binary else ("hoeffding" if bound != "hb" else "hb"),
        "certifies_risk_le_alpha": bool(ucb <= alpha),
        "vacuous": bool(ucb >= 0.99 or (ucb - rhat) > alpha),
    }
