"""
test_trustlayer.py — pytest suite for the TrustLayer commit gate.

Covers:
  (i)   anchor trio reproduces CD226->GO, RASGRP1->ABSTAIN, PRKCQ->WITHHOLD
  (ii)  default-deny: missing/blank evidence never yields GO
  (iii) conformal coverage invariant on synthetic exchangeable data
  (iv)  receipt-hash verification against the frozen receipt bundle
  (v)   counterfactual: PRKCQ with GA raised above the floor flips toward GO
"""
import glob
import json
import os

import numpy as np
import pytest

from commit_gate import (
    ABSTAIN,
    ANCHOR_EVIDENCE,
    ANCHOR_EXPECTED,
    CONDITIONS,
    GO,
    MIN_GENETIC_ASSOCIATION,
    WITHHOLD,
    CommitGate,
    compute_content_sha256,
    split_conformal_interval,
    verify_receipt,
)

HERE = os.path.dirname(os.path.abspath(__file__))
GOLD = os.path.join(HERE, "data", "gold")


# ---------------------------------------------------------------------------
# (i) Anchor trio
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("target", ["CD226", "RASGRP1", "PRKCQ"])
def test_anchor_trio_decisions(target):
    gate = CommitGate()
    decision = gate.evaluate(ANCHOR_EVIDENCE[target])
    assert decision.decision == ANCHOR_EXPECTED[target], (
        f"{target}: expected {ANCHOR_EXPECTED[target]}, got {decision.decision} "
        f"(failing={decision.failing})"
    )


def test_prkcq_vetoed_specifically_on_genetic_floor():
    """PRKCQ must be WITHHELD *because of* the genetic-association floor,
    not for some incidental reason."""
    gate = CommitGate()
    d = gate.evaluate(ANCHOR_EVIDENCE["PRKCQ"])
    assert d.decision == WITHHOLD
    assert "genetic_association" in d.veto_reasons
    assert ANCHOR_EVIDENCE["PRKCQ"]["genetic_association"] < MIN_GENETIC_ASSOCIATION


def test_rasgrp1_abstains_not_withholds():
    """RASGRP1 is an insufficiency (ABSTAIN), not a veto (WITHHOLD)."""
    gate = CommitGate()
    d = gate.evaluate(ANCHOR_EVIDENCE["RASGRP1"])
    assert d.decision == ABSTAIN
    assert d.veto_reasons == []


# ---------------------------------------------------------------------------
# (ii) Default-deny
# ---------------------------------------------------------------------------
def test_empty_evidence_is_not_go():
    gate = CommitGate()
    assert gate.evaluate({}).decision != GO


def test_none_evidence_is_not_go():
    gate = CommitGate()
    assert gate.evaluate(None).decision != GO


# The four fields whose ABSENCE must void a GO (presence-required evidence).
# eqtl_direction is deliberately excluded: its check is "not *inconsistent*",
# so a missing direction is not a failure — eQTL *presence* is enforced
# separately by celltype_matched_eqtl. An explicitly inconsistent direction is
# covered by test_all_five_conditions_are_individually_necessary_for_go.
PRESENCE_REQUIRED = [
    "leakage_clean",
    "genetic_association",
    "trust_score",
    "celltype_matched_eqtl",
]


@pytest.mark.parametrize("drop", PRESENCE_REQUIRED)
def test_dropping_any_presence_required_field_breaks_go(drop):
    """Take the passing CD226 evidence and remove one presence-required field
    at a time; the result must never remain GO (default-deny)."""
    gate = CommitGate()
    ev = dict(ANCHOR_EVIDENCE["CD226"])
    ev.pop(drop)
    assert gate.evaluate(ev).decision != GO, f"GO survived dropping {drop}"


def test_missing_eqtl_direction_alone_still_go():
    """A missing direction assertion is 'not inconsistent' -> does not void GO,
    provided eQTL presence is still confirmed. Documents the intended asymmetry."""
    gate = CommitGate()
    ev = dict(ANCHOR_EVIDENCE["CD226"])
    ev.pop("eqtl_direction")
    assert gate.evaluate(ev).decision == GO


@pytest.mark.parametrize("blank", ["", None])
def test_blank_values_are_not_go(blank):
    gate = CommitGate()
    ev = dict(ANCHOR_EVIDENCE["CD226"])
    ev["genetic_association"] = blank
    assert gate.evaluate(ev).decision != GO


def test_missing_field_is_flagged_as_missing_evidence():
    gate = CommitGate()
    ev = dict(ANCHOR_EVIDENCE["CD226"])
    ev.pop("trust_score")
    d = gate.evaluate(ev)
    assert "trust_score" in d.evidence_missing


# ---------------------------------------------------------------------------
# (iii) Conformal coverage invariant on synthetic exchangeable data
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("nominal", [0.80, 0.90, 0.95])
def test_split_conformal_coverage_on_exchangeable_data(nominal):
    """On exchangeable (i.i.d.) data, split-conformal prediction sets must
    cover the truth at >= the nominal level, up to Monte-Carlo tolerance.

    Setup: y = f(x) + noise. Point predictor = f_hat (a fit); nonconformity
    score = |y - f_hat(x)|. Calibrate the radius on a calibration split, then
    measure empirical coverage on a fresh test split. Averaged over many
    trials the empirical coverage concentrates on the nominal level.
    """
    rng = np.random.default_rng(0)
    alpha = 1.0 - nominal
    n_cal, n_test, n_trials = 500, 500, 200
    covs = []
    for _ in range(n_trials):
        # exchangeable residuals (same distribution in cal and test)
        cal_scores = np.abs(rng.standard_normal(n_cal))
        radius = split_conformal_interval(cal_scores, alpha)
        test_scores = np.abs(rng.standard_normal(n_test))
        covs.append(np.mean(test_scores <= radius))
    emp = float(np.mean(covs))
    # Split-conformal is a lower-bound guarantee; empirical mean should sit at
    # or just above nominal, and never fall materially below it.
    assert emp >= nominal - 0.02, f"under-coverage: {emp:.3f} < {nominal}"
    assert emp <= nominal + 0.05, f"gross over-coverage: {emp:.3f} >> {nominal}"


def test_conformal_radius_monotone_in_confidence():
    """Higher confidence -> wider (never narrower) interval."""
    rng = np.random.default_rng(1)
    scores = np.abs(rng.standard_normal(2000))
    r80 = split_conformal_interval(scores, 0.20)
    r90 = split_conformal_interval(scores, 0.10)
    r95 = split_conformal_interval(scores, 0.05)
    assert r80 <= r90 <= r95


# ---------------------------------------------------------------------------
# (iv) Receipt-hash verification (real integrity test)
# ---------------------------------------------------------------------------
def _hashed_receipts():
    out = []
    for p in sorted(glob.glob(os.path.join(GOLD, "*receipt*.json"))):
        with open(p) as fh:
            r = json.load(fh)
        if "content_sha256" in r:
            out.append(p)
    return out


def test_receipt_bundle_has_hashed_receipts():
    receipts = _hashed_receipts()
    assert len(receipts) >= 2, f"expected >=2 hashed receipts, found {len(receipts)}"


@pytest.mark.parametrize("path", _hashed_receipts()[:3])
def test_receipt_content_sha256_matches(path):
    """Recompute content_sha256 and assert it matches the stored field."""
    with open(path) as fh:
        receipt = json.load(fh)
    stored = receipt["content_sha256"]
    recomputed = compute_content_sha256(receipt)
    assert recomputed == stored, (
        f"{os.path.basename(path)}: hash mismatch\n"
        f"  stored     = {stored}\n  recomputed = {recomputed}"
    )
    assert verify_receipt(path) is True


def test_tampered_receipt_fails_verification():
    """A one-field mutation must break the hash (guards against a no-op check)."""
    receipts = _hashed_receipts()
    with open(receipts[0]) as fh:
        receipt = json.load(fh)
    tampered = dict(receipt)
    tampered["__injected__"] = "tamper"
    assert verify_receipt(tampered) is False


# ---------------------------------------------------------------------------
# (v) Counterfactual: gate is monotone, not hard-coded to gene names
# ---------------------------------------------------------------------------
def test_prkcq_flips_toward_go_when_ga_raised():
    """Raise PRKCQ's genetic association above the floor; the veto must clear.
    Because PRKCQ's other evidence already passes, it should flip to GO."""
    gate = CommitGate()
    before = gate.evaluate(ANCHOR_EVIDENCE["PRKCQ"])
    assert before.decision == WITHHOLD

    ev = dict(ANCHOR_EVIDENCE["PRKCQ"])
    ev["genetic_association"] = 0.25  # now >= 0.20 floor
    after = gate.evaluate(ev)
    assert after.decision == GO, (
        f"raising GA above floor did not flip PRKCQ to GO "
        f"(got {after.decision}, failing={after.failing})"
    )
    # and it must be strictly the genetic condition that changed
    assert "genetic_association" not in after.failing


def test_gate_is_monotone_in_genetic_association():
    """Sweeping GA from below to above the floor flips exactly at the floor,
    and never flips back — monotone, no hidden gene-name lookups."""
    gate = CommitGate()
    base = dict(ANCHOR_EVIDENCE["CD226"])  # everything else passes
    decisions = []
    for ga in [0.0, 0.10, 0.19, 0.20, 0.30, 0.90]:
        ev = dict(base)
        ev["genetic_association"] = ga
        decisions.append((ga, gate.evaluate(ev).decision))
    # below floor -> WITHHOLD (veto); at/above floor -> GO
    for ga, dec in decisions:
        if ga < MIN_GENETIC_ASSOCIATION:
            assert dec == WITHHOLD, (ga, dec)
        else:
            assert dec == GO, (ga, dec)


def test_all_five_conditions_are_individually_necessary_for_go():
    """From the passing CD226 case, breaking any single active condition
    removes the GO — confirms all five are load-bearing."""
    gate = CommitGate()
    breakers = {
        "leakage_clean": dict(leakage_clean=False),
        "genetic_association": dict(genetic_association=0.05),
        "trust_score": dict(trust_score=0.10),
        "celltype_eqtl_present": dict(celltype_matched_eqtl=False),
        "eqtl_direction": dict(eqtl_direction="inconsistent"),
    }
    assert set(breakers) == set(CONDITIONS)
    for cond, override in breakers.items():
        ev = dict(ANCHOR_EVIDENCE["CD226"])
        ev.update(override)
        assert gate.evaluate(ev).decision != GO, f"{cond} did not block GO"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
