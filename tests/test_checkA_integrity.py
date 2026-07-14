"""
Tests for Day-0 Check A integrity: canonical hashing round-trip + receipt self-consistency.

Guards against the hash-canonicalization bug the independent critique caught (a receipt that
hashes a sorted string but writes an unsorted file will not re-hash to its stored value).
All tests are offline (no API/model calls).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trustlayer.canonical import (  # noqa: E402
    compute_hash,
    write_hashed_json,
    verify_hashed_json,
)


def test_canonical_hash_roundtrip(tmp_path):
    """write_hashed_json -> verify_hashed_json must round-trip, incl. non-sorted input keys."""
    payload = {"z_last": 1, "a_first": {"nested": [3, 2, 1]}, "middle": "x"}
    p = tmp_path / "r.json"
    h = write_hashed_json(p, payload, "sha")
    ok, stored, recomputed = verify_hashed_json(p, "sha")
    assert ok, f"hash mismatch: {stored} != {recomputed}"
    assert stored == h


def test_hash_independent_of_input_key_order(tmp_path):
    """Same content in different key order must produce the same canonical hash."""
    a = {"b": 1, "a": 2, "c": {"y": 1, "x": 2}}
    b = {"c": {"x": 2, "y": 1}, "a": 2, "b": 1}
    assert compute_hash(a, "sha") == compute_hash(b, "sha")


def test_hash_detects_tampering(tmp_path):
    """Editing any content field after write must break verification."""
    p = tmp_path / "r.json"
    write_hashed_json(p, {"val": 10, "note": "orig"}, "sha")
    data = json.loads(p.read_text())
    data["val"] = 11  # tamper
    p.write_text(json.dumps(data, indent=2))
    ok, _, _ = verify_hashed_json(p, "sha")
    assert not ok, "tampering was not detected"


@pytest.mark.parametrize(
    "fname,key",
    [
        ("day0_checkA_frozen_set_v2.json", "frozen_set_v2_sha256"),
        ("day0_checkA_execution_receipt_v2.json", "execution_receipt_v2_sha256"),
        ("day0_checkA_sensitivity_receipt.json", "sensitivity_receipt_sha256"),
    ],
)
def test_committed_receipts_verify(fname, key):
    """The committed Check A v4 receipts must re-hash to their stored values."""
    path = ROOT / "data" / "gold" / fname
    if not path.exists():
        pytest.skip(f"{fname} not present")
    ok, stored, recomputed = verify_hashed_json(path, key)
    assert ok, f"{fname}: stored {stored[:16]} != recomputed {recomputed[:16]}"


def test_frozen_set_labels_are_deterministic_function():
    """Every GO label in the frozen set must equal its own rule: immune AND t1d_dom AND
    not-t2d-tagged AND positive-direction. Catches any hand-editing of labels."""
    path = ROOT / "data" / "gold" / "day0_checkA_frozen_set_v2.json"
    if not path.exists():
        pytest.skip("frozen set v2 not present")
    fs = json.loads(path.read_text())
    for g in fs["genes"]:
        expected = bool(
            g["has_immune_go"]
            and (g["ot_t1d_assoc"] > g["ot_t2d_assoc"])
            and (not g["gwas_t2d_tag"])
            and (g["measured_treg_score_rest"] > 0)
        )
        assert g["ground_truth_credible_t1d_tcell_go"] == expected, (
            f"{g['symbol']}: stored label != rule recomputation"
        )


def test_execution_outcome_matches_scores():
    """The receipt's outcome string must be consistent with its own MCC numbers + win rule."""
    path = ROOT / "data" / "gold" / "day0_checkA_execution_receipt_v2.json"
    if not path.exists():
        pytest.skip("execution receipt v2 not present")
    rc = json.loads(path.read_text())
    claude_mcc = rc["claude_score"]["mcc"]
    best = rc["best_baseline_mcc"]
    gap = round(claude_mcc - best, 4)
    if gap >= 0.15 and claude_mcc >= 0.50:
        expected = "WIN"
    elif claude_mcc < best:
        expected = "LOSS"
    else:
        expected = "TIE_NO_WIN"
    assert rc["outcome"] == expected, f"outcome {rc['outcome']} != derived {expected}"
