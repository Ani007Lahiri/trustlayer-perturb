"""
commit_gate.py — reference implementation of the TrustLayer commit gate.

The gate answers a single question about a candidate perturbation target:
should we COMMIT this to the wet-lab queue (GO), hold it because the evidence
is insufficient (ABSTAIN), or actively refuse it because the evidence
contradicts a safety/validity requirement (WITHHOLD)?

Design: default-DENY. A GO decision requires that ALL FIVE gate conditions
pass. Anything short of that is never a GO. The distinction between the two
non-GO outcomes is principled:

  * WITHHOLD  — a *veto*: the evidence that is present positively contradicts a
                requirement (genetic association below the floor, eQTL
                direction inconsistent, or leakage detected). We have a
                reason NOT to act.
  * ABSTAIN   — *insufficiency*: a requirement cannot be confirmed (trust below
                floor / undercoverage, missing cell-type-matched eQTL, or
                missing evidence entirely). We do not have enough to act.

Priority (safety-first): if any veto condition fails -> WITHHOLD; else if any
insufficiency condition fails -> ABSTAIN; else -> GO.

This module is a clean reconstruction from the frozen spec; the numeric floors
are the frozen constants below.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

# ---------------------------------------------------------------------------
# Frozen gate constants (from spec)
# ---------------------------------------------------------------------------
MIN_GENETIC_ASSOCIATION = 0.20
MIN_TRUST_SCORE = 0.50

# Decision labels
GO = "GO"
ABSTAIN = "ABSTAIN"
WITHHOLD = "WITHHOLD"

# The five gate conditions, by canonical name.
CONDITIONS = (
    "leakage_clean",          # integrity: no train/cal/test leakage
    "genetic_association",    # genetic_association >= MIN_GENETIC_ASSOCIATION
    "trust_score",            # trust_score >= MIN_TRUST_SCORE (conformal reliability)
    "celltype_eqtl_present",  # a cell-type-matched eQTL exists for the target
    "eqtl_direction",         # eQTL effect direction is not inconsistent
)

# Which failure class each condition contributes when it fails.
#   "veto"          -> WITHHOLD (evidence contradicts a requirement)
#   "insufficiency" -> ABSTAIN  (requirement cannot be confirmed)
_FAILURE_CLASS = {
    "leakage_clean": "veto",           # detected leakage is a hard integrity veto
    "genetic_association": "veto",     # GA present-but-below-floor is a genetic veto
    "trust_score": "insufficiency",    # low trust == undercoverage == abstain
    "celltype_eqtl_present": "insufficiency",  # missing eQTL == cannot confirm
    "eqtl_direction": "veto",          # inconsistent direction contradicts
}


@dataclass
class GateDecision:
    """Structured result of a gate evaluation."""
    decision: str
    passed: Dict[str, bool]
    failing: List[str]
    veto_reasons: List[str] = field(default_factory=list)
    abstain_reasons: List[str] = field(default_factory=list)
    evidence_missing: List[str] = field(default_factory=list)

    @property
    def is_go(self) -> bool:
        return self.decision == GO

    def to_dict(self) -> Dict[str, Any]:
        return {
            "decision": self.decision,
            "passed": self.passed,
            "failing": self.failing,
            "veto_reasons": self.veto_reasons,
            "abstain_reasons": self.abstain_reasons,
            "evidence_missing": self.evidence_missing,
        }


class CommitGate:
    """Default-DENY commit gate.

    Parameters
    ----------
    min_genetic_association, min_trust_score
        Numeric floors (frozen defaults from spec).
    enabled_conditions
        Iterable of condition names that are ACTIVE. A condition not in this
        set is treated as always-passing (used by the leave-one-condition-out
        ablation). Defaults to all five.
    """

    def __init__(
        self,
        min_genetic_association: float = MIN_GENETIC_ASSOCIATION,
        min_trust_score: float = MIN_TRUST_SCORE,
        enabled_conditions: Optional[List[str]] = None,
    ):
        self.min_genetic_association = float(min_genetic_association)
        self.min_trust_score = float(min_trust_score)
        if enabled_conditions is None:
            self.enabled = set(CONDITIONS)
        else:
            unknown = set(enabled_conditions) - set(CONDITIONS)
            if unknown:
                raise ValueError(f"unknown conditions: {sorted(unknown)}")
            self.enabled = set(enabled_conditions)

    # -- per-condition checks ------------------------------------------------
    def _check(self, name: str, evidence: Dict[str, Any]):
        """Return (passed: bool, missing: bool) for one condition.

        `missing` marks that the evidence needed to evaluate this condition was
        absent/blank (default-DENY: a missing field never satisfies a
        condition).
        """
        if name == "leakage_clean":
            v = evidence.get("leakage_clean", None)
            if v is None or v == "":
                return False, True
            return bool(v), False

        if name == "genetic_association":
            v = evidence.get("genetic_association", None)
            if v is None or v == "":
                return False, True
            return float(v) >= self.min_genetic_association, False

        if name == "trust_score":
            v = evidence.get("trust_score", None)
            if v is None or v == "":
                return False, True
            return float(v) >= self.min_trust_score, False

        if name == "celltype_eqtl_present":
            v = evidence.get("celltype_matched_eqtl", None)
            if v is None or v == "":
                return False, True
            return bool(v), False

        if name == "eqtl_direction":
            # "not inconsistent". Explicit 'inconsistent' fails; a missing
            # direction is only a problem if an eQTL is claimed present.
            v = evidence.get("eqtl_direction", None)
            if v is None or v == "":
                # no direction asserted: passes the *direction* check
                # (presence is handled by celltype_eqtl_present)
                return True, False
            return (str(v).lower() != "inconsistent"), False

        raise ValueError(f"unknown condition {name}")

    # -- main entry ----------------------------------------------------------
    def evaluate(self, evidence: Optional[Dict[str, Any]] = None, **kw) -> GateDecision:
        """Evaluate the gate on an evidence dict (or keyword evidence)."""
        if evidence is None:
            evidence = {}
        if kw:
            evidence = {**evidence, **kw}

        passed: Dict[str, bool] = {}
        failing: List[str] = []
        veto_reasons: List[str] = []
        abstain_reasons: List[str] = []
        evidence_missing: List[str] = []

        for name in CONDITIONS:
            if name not in self.enabled:
                # ablation: condition removed -> forced pass
                passed[name] = True
                continue
            ok, missing = self._check(name, evidence)
            passed[name] = ok
            if not ok:
                failing.append(name)
                if missing:
                    evidence_missing.append(name)
                if _FAILURE_CLASS[name] == "veto":
                    veto_reasons.append(name)
                else:
                    abstain_reasons.append(name)

        if not failing:
            decision = GO
        elif veto_reasons:
            decision = WITHHOLD
        else:
            decision = ABSTAIN

        return GateDecision(
            decision=decision,
            passed=passed,
            failing=failing,
            veto_reasons=veto_reasons,
            abstain_reasons=abstain_reasons,
            evidence_missing=evidence_missing,
        )


# ---------------------------------------------------------------------------
# Anchor-trio ground truth (frozen)
# ---------------------------------------------------------------------------
# Evidence for the three anchor targets and their expected decisions. These are
# the regression anchors: any change to the gate that moves one of these is a
# breaking change.
#   CD226   -> GO       (all conditions pass)
#   RASGRP1 -> ABSTAIN  (trust below floor: insufficiency, cannot confirm)
#   PRKCQ   -> WITHHOLD (genetic association 0.162 < 0.20: genetic-floor veto)
ANCHOR_EVIDENCE: Dict[str, Dict[str, Any]] = {
    "CD226": dict(
        leakage_clean=True, genetic_association=0.71, trust_score=0.82,
        celltype_matched_eqtl=True, eqtl_direction="consistent",
    ),
    "RASGRP1": dict(
        leakage_clean=True, genetic_association=0.34, trust_score=0.41,
        celltype_matched_eqtl=True, eqtl_direction="consistent",
    ),
    "PRKCQ": dict(
        leakage_clean=True, genetic_association=0.162, trust_score=0.63,
        celltype_matched_eqtl=True, eqtl_direction="consistent",
    ),
}
ANCHOR_EXPECTED: Dict[str, str] = {
    "CD226": GO,
    "RASGRP1": ABSTAIN,
    "PRKCQ": WITHHOLD,
}


# ---------------------------------------------------------------------------
# Receipt integrity
# ---------------------------------------------------------------------------
def compute_content_sha256(receipt: Dict[str, Any]) -> str:
    """Recompute the canonical content hash of a receipt.

    Scheme (reverse-engineered from the frozen receipts and verified against
    all 7 hashed receipts in the bundle):
        sha256( json.dumps(receipt_without_'content_sha256',
                            sort_keys=True, indent=2) )
    """
    body = {k: v for k, v in receipt.items() if k != "content_sha256"}
    canonical = json.dumps(body, sort_keys=True, indent=2)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def verify_receipt(path_or_obj) -> bool:
    """Return True iff a receipt's stored content_sha256 matches recomputation."""
    if isinstance(path_or_obj, (str, bytes)):
        with open(path_or_obj) as fh:
            receipt = json.load(fh)
    else:
        receipt = path_or_obj
    stored = receipt.get("content_sha256")
    if not stored:
        raise ValueError("receipt has no content_sha256 field")
    return compute_content_sha256(receipt) == stored


# ---------------------------------------------------------------------------
# Split-conformal helper (for the coverage-invariant test)
# ---------------------------------------------------------------------------
def split_conformal_interval(cal_scores, alpha: float):
    """Return the split-conformal quantile (radius) for nonconformity scores.

    Uses the finite-sample-valid rank: ceil((n+1)(1-alpha))/n empirical
    quantile of the calibration nonconformity scores. On exchangeable data the
    resulting prediction set covers with probability >= 1-alpha.
    """
    import math
    import numpy as np

    s = np.sort(np.asarray(cal_scores, dtype=float))
    n = len(s)
    if n == 0:
        raise ValueError("empty calibration set")
    k = math.ceil((n + 1) * (1.0 - alpha))
    if k > n:
        return float("inf")  # not enough calibration points for this level
    return float(s[k - 1])
