"""
Tests for the deterministic commit gate (Plan v3, Fix 4).

These encode the v3 target-hierarchy ruling as executable assertions:
  - CD226   (anchor)  -> GO       : genetics-secure, cell-type eQTL direction OK
  - RASGRP1 (bet)     -> WITHHELD : no cell-type-matched eQTL (proxy only) -> veto fires
  - PRKCQ   (control) -> WITHHELD : genetic_association 0.162 below floor
  - any leakage failure -> WITHHELD, always
Plus determinism (same input -> same hash).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from trustlayer.commit_gate import CommitGate, CriticVerdict  # noqa: E402


# ---- live-derived genetics values (from data/gold/genetics_gate_receipt.json) ----
def cd226_verdict(**over):
    v = dict(
        gene="CD226",
        genetic_association=0.834,
        genome_wide_sig_snp=True,
        celltype_matched_eqtl=True,
        eqtl_direction_consistent=True,
        proxy_tissue_only=False,
        trust_score=0.82,
        leakage_audit_passed=True,
        notes="rs763361 GoF; KD mimics protection",
    )
    v.update(over)
    return CriticVerdict(**v)


def rasgrp1_verdict(**over):
    # v3: risk allele RAISES RASGRP1 in LCL proxy (+0.57), but NO CD4/Treg eQTL exists.
    v = dict(
        gene="RASGRP1",
        genetic_association=0.506,
        genome_wide_sig_snp=True,
        celltype_matched_eqtl=False,
        eqtl_direction_consistent=None,
        proxy_tissue_only=True,
        trust_score=0.71,
        leakage_audit_passed=True,
        notes="rs72727394 p=4e-10; proxy-tissue only, cell-type-unconfirmed",
    )
    v.update(over)
    return CriticVerdict(**v)


def prkcq_verdict(**over):
    v = dict(
        gene="PRKCQ",
        genetic_association=0.162,
        genome_wide_sig_snp=False,
        celltype_matched_eqtl=False,
        eqtl_direction_consistent=None,
        proxy_tissue_only=True,
        trust_score=0.40,
        leakage_audit_passed=True,
        notes="sub-GWS; thyroid/allergy locus",
    )
    v.update(over)
    return CriticVerdict(**v)


def test_cd226_goes(tmp_path):
    gate = CommitGate(out_dir=tmp_path)
    nom = gate.evaluate(cd226_verdict(), role="ANCHOR")
    assert nom.decision == "GO", nom.reasons
    p = gate.commit(nom)
    assert "GO" in p.name and p.exists()


def test_rasgrp1_withheld_no_celltype_eqtl(tmp_path):
    """The money shot: the system refuses its OWN novel bet on real null data."""
    gate = CommitGate(out_dir=tmp_path)
    nom = gate.evaluate(rasgrp1_verdict(), role="BET")
    assert nom.decision == "WITHHELD"
    assert any("no cell-type-matched" in r for r in nom.reasons)
    p = gate.commit(nom)
    assert "BLOCKED" in p.name and p.exists()


def test_prkcq_withheld_below_genetic_floor(tmp_path):
    gate = CommitGate(out_dir=tmp_path)
    nom = gate.evaluate(prkcq_verdict(), role="CONTROL")
    assert nom.decision == "WITHHELD"
    assert any("genetic_association" in r and "floor" in r for r in nom.reasons)


def test_leakage_always_blocks(tmp_path):
    """Even a perfectly-anchored gene is blocked if the leakage audit fails."""
    gate = CommitGate(out_dir=tmp_path)
    nom = gate.evaluate(cd226_verdict(leakage_audit_passed=False), role="ANCHOR")
    assert nom.decision == "WITHHELD"
    assert any("leakage" in r for r in nom.reasons)


def test_low_trust_blocks_even_if_anchored(tmp_path):
    gate = CommitGate(out_dir=tmp_path)
    nom = gate.evaluate(cd226_verdict(trust_score=0.10), role="ANCHOR")
    assert nom.decision == "WITHHELD"
    assert any("trust_score" in r for r in nom.reasons)


def test_deterministic_hash():
    gate = CommitGate(out_dir="data/gold/nominations")
    a = gate.evaluate(cd226_verdict(), role="ANCHOR")
    b = gate.evaluate(cd226_verdict(), role="ANCHOR")
    assert a.content_hash() == b.content_hash()


def test_default_deny_on_missing_everything(tmp_path):
    gate = CommitGate(out_dir=tmp_path)
    v = CriticVerdict(
        gene="X",
        genetic_association=None,
        genome_wide_sig_snp=False,
        celltype_matched_eqtl=False,
        eqtl_direction_consistent=None,
        proxy_tissue_only=False,
        trust_score=None,
        leakage_audit_passed=False,
    )
    nom = gate.evaluate(v, role="BET")
    assert nom.decision == "WITHHELD"
    assert len([r for r in nom.reasons if r.startswith("BLOCK")]) >= 3
