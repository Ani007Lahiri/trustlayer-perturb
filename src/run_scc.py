"""Run the SCC engine on a synergy table and write receipts.

Usage:
  python src/run_scc.py --dataset norman
  python src/run_scc.py --dataset carpool

Reads data/interim/<ds>_synergy_table.parquet (columns: true_mag, additive_mag,
sum_mag, ...), runs LOPO conformal + the logistic interaction test + the frozen
composite clearing criterion, and writes data/interim/<ds>_scc_engine_result.json.

For norman, also cross-checks against the frozen ad-hoc receipts to prove the
rebuilt engine reproduces the original numbers.
"""

from __future__ import annotations
import argparse
import json
import os
import sys
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
sys.path.insert(0, ROOT)
from src.trustlayer.scc import SCCConfig, run_lopo, interaction_test, composite_clearing  # noqa: E402

TABLES = {
    "norman": (
        "data/interim/norman_synergy_table.parquet",
        "data/interim/norman_synergy_feats.parquet",
    ),
    "carpool": (
        "data/interim/carpool_synergy_table.parquet",
        "data/interim/carpool_synergy_feats.parquet",
    ),
}


def load_table(dataset: str) -> pd.DataFrame:
    tpath, fpath = TABLES[dataset]
    table = pd.read_parquet(tpath)
    if "sum_mag" not in table.columns and os.path.exists(fpath):
        feats = pd.read_parquet(fpath)[["pair", "sum_mag"]]
        table = table.merge(feats, on="pair", how="left")
    return table


def main(dataset: str, seed: int = 0):
    table = load_table(dataset)
    cfg = SCCConfig(seed=seed)

    lopo = run_lopo(table, cfg)
    itest = interaction_test(lopo, cfg)
    clearing = composite_clearing(lopo, cfg, itest)

    out = {
        "dataset": dataset,
        "n_doubles": int(lopo["n"]),
        "cov_vanilla": lopo["cov_vanilla"],
        "cov_scc": lopo["cov_scc"],
        "nominal": lopo["nominal"],
        "interaction_coeff": itest["interaction_coeff"],
        "perm_p_one_sided": itest["perm_p_one_sided"],
        "boot_ci": itest["boot_ci"],
        "clearing": clearing,
        "seed": seed,
    }
    outpath = f"data/interim/{dataset}_scc_engine_result.json"
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {outpath}")
    print(json.dumps(out, indent=2))

    if dataset == "norman":
        frozen = json.load(open("data/interim/scc_lopo_result.json"))
        print("\n=== NORMAN REBUILD vs FROZEN AD-HOC RECEIPT ===")
        print(
            f"  interaction: rebuilt {itest['interaction_coeff']:.4f} | "
            f"frozen {frozen['interaction_coeff']:.4f}"
        )
        print(
            f"  perm_p:      rebuilt {itest['perm_p_one_sided']:.4f} | "
            f"frozen {frozen['perm_p_one_sided']:.4f}"
        )
        print(
            f"  cov_vanilla: rebuilt {lopo['cov_vanilla']:.4f} | "
            f"frozen {frozen['cov_vanilla']:.4f}"
        )
        print(
            f"  cov_scc:     rebuilt {lopo['cov_scc']:.4f} | "
            f"frozen {frozen['cov_scc']:.4f}"
        )
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=list(TABLES.keys()))
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    main(args.dataset, args.seed)
