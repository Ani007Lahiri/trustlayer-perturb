"""Tests for the frozen leakage-safe splits (Plan v3 anti-leakage rule)."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from trustlayer import splits  # noqa: E402

DATA = Path("data/raw/DE_stats.suppl_table.csv")
pytestmark = pytest.mark.skipif(not DATA.exists(), reason="DE_stats table not pulled")


def test_no_fold_overlap():
    m = splits.build()
    tr, ca, te = set(m.train_genes), set(m.calib_genes), set(m.test_genes)
    assert not (tr & ca) and not (tr & te) and not (ca & te)


def test_gold_forced_into_test_never_calibration():
    m = splits.build()
    assert set(m.gold_in_test).issubset(set(m.test_genes))
    assert not (set(m.gold_in_test) & set(m.calib_genes))


def test_deterministic_across_runs():
    a = splits.build()
    b = splits.build()
    assert a.content_hash == b.content_hash
    assert a.train_genes == b.train_genes


def test_covers_all_perturbations():
    m = splits.build()
    total = len(set(m.train_genes) | set(m.calib_genes) | set(m.test_genes))
    assert total == m.n_perturbations
