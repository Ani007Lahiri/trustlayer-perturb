"""
Day-2: fit the deliberately-simple base predictor + baselines (Plan v3 §1).

The base predictor is NOT the product — it exists so the Day-3 conformal layer has
residuals to calibrate. We only need it to beat mean/shuffle baselines (signal exists,
headroom remains). Leakage-safe features only (see features.py).

Outputs:
  data/interim/base_predictor_preds.parquet  (per-row predictions on all folds)
  data/gold/base_predictor_metrics.json      (train/calib/test R2, Spearman, baselines)
  models/base_predictor.joblib
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import r2_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.features import feature_cols  # noqa: E402

TABLE = Path("data/interim/day2_feature_table.parquet")
PREDS = Path("data/interim/base_predictor_preds.parquet")
METRICS = Path("data/gold/base_predictor_metrics.json")
MODEL = Path("models/base_predictor.joblib")
SEED = 20260708


def evaluate(y_true, y_pred) -> dict:
    return {
        "r2": float(r2_score(y_true, y_pred)),
        "spearman": float(spearmanr(y_true, y_pred).statistic),
        "rmse": float(np.sqrt(np.mean((y_true - y_pred) ** 2))),
        "n": int(len(y_true)),
    }


def main() -> int:
    df = pd.read_parquet(TABLE)
    fc = feature_cols()

    tr = df[df.fold == "train"]
    ca = df[df.fold == "calib"]
    te = df[df.fold == "test"]

    Xtr, ytr = tr[fc].values, tr["y"].values
    model = HistGradientBoostingRegressor(
        max_iter=300,
        learning_rate=0.05,
        max_depth=4,
        l2_regularization=1.0,
        random_state=SEED,
    )
    model.fit(Xtr, ytr)

    # ---- baselines ----
    global_mean = float(ytr.mean())
    rng = np.random.default_rng(SEED)
    # shuffle baseline: permute features (destroy feature->label link), refit
    Xtr_shuf = Xtr.copy()
    for j in range(Xtr_shuf.shape[1]):
        Xtr_shuf[:, j] = rng.permutation(Xtr_shuf[:, j])
    shuf_model = HistGradientBoostingRegressor(
        max_iter=300,
        learning_rate=0.05,
        max_depth=4,
        l2_regularization=1.0,
        random_state=SEED,
    ).fit(Xtr_shuf, ytr)

    metrics = {"features": fc, "seed": SEED, "global_mean_target": global_mean}
    for name, part in [("train", tr), ("calib", ca), ("test", te)]:
        X, y = part[fc].values, part["y"].values
        pred = model.predict(X)
        mean_pred = np.full_like(y, global_mean)
        shuf_pred = shuf_model.predict(X)
        metrics[name] = {
            "model": evaluate(y, pred),
            "mean_baseline": evaluate(y, mean_pred),
            "shuffle_baseline": evaluate(y, shuf_pred),
        }

    # residuals on calib (what the conformal layer will use) + full preds
    df["pred"] = model.predict(df[fc].values)
    df["residual"] = df["y"] - df["pred"]
    df[
        ["gene", "condition", "fold", "y", "pred", "residual", "trans_effect_magnitude"]
    ].to_parquet(PREDS)

    METRICS.parent.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, indent=2))
    MODEL.parent.mkdir(parents=True, exist_ok=True)
    try:
        import joblib

        joblib.dump(model, MODEL)
    except Exception:
        pass

    # ---- summary ----
    print("=" * 64)
    print("BASE PREDICTOR  (leakage-safe covariates; target=log1p trans-effect)")
    print("=" * 64)
    for part in ["train", "calib", "test"]:
        m = metrics[part]
        print(f"\n  [{part}]  n={m['model']['n']}")
        print(
            f"    model    : R2={m['model']['r2']:+.3f}  Spearman={m['model']['spearman']:+.3f}  RMSE={m['model']['rmse']:.3f}"
        )
        print(f"    mean     : R2={m['mean_baseline']['r2']:+.3f}  (should be ~0)")
        print(f"    shuffle  : R2={m['shuffle_baseline']['r2']:+.3f}  (should be <= 0)")
    beats = (
        metrics["test"]["model"]["r2"] > metrics["test"]["shuffle_baseline"]["r2"]
        and metrics["test"]["model"]["r2"] > 0
    )
    print(f"\n  Base predictor beats baselines on test (signal exists): {beats}")
    print(f"  -> {METRICS}\n  -> {PREDS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
