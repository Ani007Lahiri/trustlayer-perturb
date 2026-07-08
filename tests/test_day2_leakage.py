"""Day-2 leakage tests: feature/target independence + axis disjointness."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

DATA = Path("data/raw/GWCD4i.DE_stats.h5ad")
pytestmark = pytest.mark.skipif(not DATA.exists(), reason="DE_stats h5ad not pulled")

FORBIDDEN = {
    "ontarget_effect_size",
    "n_up_genes",
    "n_down_genes",
    "n_total_de_genes",
    "n_downstream",
    "ontarget_significant",
    "ontarget_effect_category",
    "n_total_genes_category",
}


def test_feature_cols_disjoint_from_forbidden():
    from trustlayer.features import feature_cols

    assert not (set(feature_cols()) & FORBIDDEN)


def test_axis_genes_disjoint_from_gold():
    manifest = json.loads(Path("data/gold/treg_axis_manifest.json").read_text())
    gold = json.loads(Path("data/gold/t1d_gold_set.json").read_text())
    forbidden = set(gold["axis_forbidden_genes"])
    used = set(manifest["pos_genes_used"]) | set(manifest["neg_genes_used"])
    assert not (used & forbidden), f"axis leaked gold genes: {used & forbidden}"


def test_base_predictor_beats_baselines():
    m = json.loads(Path("data/gold/base_predictor_metrics.json").read_text())
    # signal exists: model R2 > 0 and > shuffle on held-out test
    assert m["test"]["model"]["r2"] > 0
    assert m["test"]["model"]["r2"] > m["test"]["shuffle_baseline"]["r2"]
    # mean baseline ~ 0 (sanity)
    assert abs(m["test"]["mean_baseline"]["r2"]) < 0.02


def test_gold_genes_only_in_test_fold():
    """The base predictor must never have trained on a gold gene."""
    import pandas as pd

    df = pd.read_parquet("data/interim/day2_feature_table.parquet")
    gold = set(
        json.loads(Path("data/gold/t1d_gold_set.json").read_text())[
            "recovery_positive_set"
        ]
    )
    folds = df[df.gene.isin(gold)]["fold"].unique().tolist()
    assert folds == ["test"], f"gold genes leaked into non-test folds: {folds}"
