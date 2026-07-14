"""
Day-5 end-to-end tests: the believe/veto split, trust non-binding, and counterfactuals.
Verifies the critique-required properties:
  - trust is NEVER the binding constraint for the demo trio (genetics/eQTL carry it)
  - PRKCQ has the HIGHEST trust yet is still WITHHELD (trust doesn't drive decisions)
  - counterfactuals: RASGRP1+eQTL -> GO; PRKCQ+high-GA -> still WITHHELD (trust non-binding)
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from trustlayer.commit_gate import CommitGate, CriticVerdict, MIN_TRUST_SCORE  # noqa: E402

RECEIPT = Path("data/gold/pipeline_day5_receipt.json")
import pytest  # noqa: E402

pytestmark = pytest.mark.skipif(
    not RECEIPT.exists(), reason="run src/run_pipeline_day5.py first"
)


def _receipt():
    return json.loads(RECEIPT.read_text())


def test_believe_veto_split_matches_v3():
    r = _receipt()
    dec = {s["gene"]: s["decision"] for s in r["believe_veto_split"]}
    assert dec["CD226"] == "GO"
    assert dec["RASGRP1"] == "WITHHELD"
    assert dec["PRKCQ"] == "WITHHELD"


def test_trust_non_binding_for_trio():
    """Every trio gene clears the trust floor by margin -> trust is not the binding gate."""
    r = _receipt()
    assert r["trust_non_binding_for_trio"] is True
    for s in r["believe_veto_split"]:
        assert s["trust"] >= r["trust_margin_registered"]


def test_prkcq_highest_trust_but_withheld():
    """The strongest honesty check: PRKCQ has top trust yet is blocked by genetics."""
    r = _receipt()
    trust = {s["gene"]: s["trust"] for s in r["believe_veto_split"]}
    assert trust["PRKCQ"] == max(trust.values())
    dec = {s["gene"]: s["decision"] for s in r["believe_veto_split"]}
    assert dec["PRKCQ"] == "WITHHELD"


def test_rasgrp1_binding_is_eqtl_not_trust():
    r = _receipt()
    rasgrp1 = next(s for s in r["believe_veto_split"] if s["gene"] == "RASGRP1")
    assert any("eQTL" in b for b in rasgrp1["binding_constraints"])
    assert not any("trust_score" in b for b in rasgrp1["binding_constraints"])


def test_counterfactual_rasgrp1_eqtl_flips_to_go():
    r = _receipt()
    assert r["counterfactuals"]["RASGRP1_add_celltype_eqtl"] == "GO"


def test_counterfactual_prkcq_high_ga_flips_to_go():
    """v4 (live eQTL correction): PRKCQ genuinely HAS a CD4/Treg eQTL (3/11 datasets), so its
    ONLY binding constraint is the genetic floor. Raise GA above the floor -> flips to GO.
    This proves trust (0.883, the trio's HIGHEST) was never binding: PRKCQ is withheld purely
    on genetics. (Prior test asserted WITHHELD under the hand-typed eQTL=False, which the live
    query falsified.)"""
    r = _receipt()
    assert r["counterfactuals"]["PRKCQ_raise_GA_to_0.50"] == "GO"


def test_gate_is_deterministic():
    v = CriticVerdict(
        gene="X",
        genetic_association=0.8,
        genome_wide_sig_snp=True,
        celltype_matched_eqtl=True,
        eqtl_direction_consistent=True,
        proxy_tissue_only=False,
        trust_score=0.9,
        leakage_audit_passed=True,
    )
    g = CommitGate(out_dir="data/gold/nominations")
    assert (
        g.evaluate(v, "ANCHOR").content_hash() == g.evaluate(v, "ANCHOR").content_hash()
    )
