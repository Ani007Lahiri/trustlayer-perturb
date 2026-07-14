"""TrustLayer: calibrated trust for perturbation-effect prediction.

The public API is deliberately small:
    tl = TrustLayer(alpha=0.10)
    tl.calibrate(residual_scores_cal)      # split-conformal calibration
    lo_hi = tl.interval(y_pred, scale)     # prediction interval at 1-alpha
    tl.covers(y_true, y_pred)              # bool: inside the conformal set
    tl.should_commit(trust, genetic_assoc, eqtl_dir)  # default-deny gate

Design note (honest): the conformal guarantee is EXCHANGEABILITY-conditional.
Under distribution shift (e.g. leave-one-donor-out) exchangeability is void and
coverage is NOT guaranteed; `shift_gap()` quantifies the empirical undercoverage
so the caller can ABSTAIN. See the pre-registered calibration-law receipt.
"""
import numpy as np

MIN_GENETIC_ASSOCIATION = 0.20
MIN_TRUST_SCORE = 0.50

class TrustLayer:
    def __init__(self, alpha=0.10):
        self.alpha = float(alpha)
        self._q = None

    def calibrate(self, scores_cal):
        """scores_cal: 1D array of nonconformity scores on a calibration set."""
        s = np.asarray(scores_cal, float)
        # finite-sample split-conformal quantile
        n = len(s); level = np.ceil((n + 1) * (1 - self.alpha)) / n
        self._q = float(np.quantile(s, min(level, 1.0), method="higher"))
        return self._q

    def covers(self, score_test):
        if self._q is None: raise RuntimeError("call calibrate() first")
        return np.asarray(score_test, float) <= self._q

    def empirical_coverage(self, scores_test):
        return float(self.covers(scores_test).mean())

    def shift_gap(self, scores_test, nominal=None):
        """Undercoverage gap = nominal - empirical. Positive = undercovering."""
        nominal = (1 - self.alpha) if nominal is None else nominal
        return float(nominal - self.empirical_coverage(scores_test))

    def should_commit(self, trust_score, genetic_association, eqtl_direction_ok):
        """Default-DENY commit gate. GO only if ALL floors are cleared."""
        return bool(
            trust_score >= MIN_TRUST_SCORE
            and genetic_association >= MIN_GENETIC_ASSOCIATION
            and eqtl_direction_ok
        )
