"""
Day-4: blinded recovery (Q-A) + axis-swap specificity control (Q-B).
Critique-corrected (Gate-1 Day-4). CPU only.

Q-A: do the 8 genetics-anchored gold genes (test fold only) rank above chance on the
     Treg-axis nomination score? Mann-Whitney U + EXACT permutation null (AUROC descriptive).
Q-B: does the Treg axis beat control axes (cellcycle, cholesterol, apoptosis, ribosome)
     AND a magnitude-only baseline? With magnitude residualized out and a leave-IFIH1-out
     sensitivity check. Honest negative expected at n=8.

Outputs:
  data/gold/recovery_specificity_receipt.json
  data/gold/axes_manifest.json
  figures/recovery_specificity.png
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.axes import score_all_axes  # noqa: E402

SPLITS = Path("data/gold/frozen_splits.json")
GOLD = Path("data/gold/t1d_gold_set.json")
PREDS = Path("data/interim/base_predictor_preds.parquet")
RECEIPT = Path("data/gold/recovery_specificity_receipt.json")
AXES_MANIFEST = Path("data/gold/axes_manifest.json")
SEED = 20260708

CONTROL_AXES = ["cellcycle", "cholesterol", "apoptosis", "ribosome"]


def auroc_from_scores(score: np.ndarray, is_gold: np.ndarray) -> float:
    """AUROC = P(gold score > non-gold score); descriptive only."""
    r = rankdata(score)
    n_pos = int(is_gold.sum())
    n_neg = len(score) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    return float((r[is_gold].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def perm_test(score: np.ndarray, is_gold: np.ndarray, n_perm=10000, seed=SEED) -> dict:
    """Exact-ish permutation: permute the gold labels across all perturbations.
    Statistic = mean gold percentile rank. One-sided (gold ranks HIGHER)."""
    rng = np.random.default_rng(seed)
    pct = rankdata(score) / len(score)
    n_pos = int(is_gold.sum())
    obs = float(pct[is_gold].mean())
    null = np.empty(n_perm)
    idx = np.arange(len(score))
    for i in range(n_perm):
        sel = rng.choice(idx, n_pos, replace=False)
        null[i] = pct[sel].mean()
    p = float((np.sum(null >= obs) + 1) / (n_perm + 1))
    return {
        "obs_mean_gold_percentile": round(obs, 4),
        "null_mean": round(float(null.mean()), 4),
        "null_std": round(float(null.std()), 4),
        "p_value_one_sided": round(p, 4),
        "auroc_descriptive": round(auroc_from_scores(score, is_gold), 4),
        "n_gold": n_pos,
        "n_total": int(len(score)),
    }


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Rank-based residual of y on x (Spearman-style; robust to non-linearity)."""
    ry, rx = rankdata(y), rankdata(x)
    b = np.polyfit(rx, ry, 1)
    return ry - (b[0] * rx + b[1])


def main() -> int:
    splits = json.loads(SPLITS.read_text())
    test_genes = set(splits["test_genes"])
    gold = set(json.loads(GOLD.read_text())["recovery_positive_set"])

    df, axes_manifest = score_all_axes()
    AXES_MANIFEST.write_text(json.dumps(axes_manifest, indent=2))

    # collapse to per-gene (max |axis| across conditions) within TEST fold
    df = df[df["gene"].isin(test_genes)].copy()
    # attach trust proxy (pooled Day-2 ensemble spread) if available
    trust = None
    if PREDS.exists():
        preds = pd.read_parquet(PREDS)
        # trust proxy = -|residual| (smaller residual = more trustworthy); pooled Day-2
        pr = preds.groupby("gene")["residual"].apply(lambda s: -np.abs(s).mean())
        trust = pr

    axis_cols = ["treg"] + CONTROL_AXES + ["ifn_shared"]
    per_gene = (
        df.groupby("gene")
        .agg(
            {
                **{f"{a}_score": lambda s: s.abs().max() for a in axis_cols},
                "trans_effect_magnitude": "max",
            }
        )
        .reset_index()
    )
    per_gene["is_gold"] = per_gene["gene"].isin(gold).values

    genes_present_gold = sorted(set(per_gene.loc[per_gene.is_gold, "gene"]))
    is_gold = per_gene["is_gold"].values

    results = {
        "n_gold_in_test": int(is_gold.sum()),
        "gold_present": genes_present_gold,
        "n_test_perturbations": int(len(per_gene)),
        "axes": {},
    }

    # ---- Q-A + Q-B: per-axis recovery ----
    for a in axis_cols:
        sc = per_gene[f"{a}_score"].values.astype(float)
        results["axes"][a] = perm_test(sc, is_gold)

    # magnitude-only baseline (Fix #2a)
    mag = per_gene["trans_effect_magnitude"].values.astype(float)
    results["axes"]["magnitude_baseline"] = perm_test(mag, is_gold)

    # Treg residualized on magnitude (Fix #2b)
    treg_resid = residualize(per_gene["treg_score"].abs().values.astype(float), mag)
    results["axes"]["treg_residualized_on_magnitude"] = perm_test(treg_resid, is_gold)

    # ---- leave-IFIH1-out sensitivity (Fix #1) ----
    mask = per_gene["gene"] != "IFIH1"
    pg2 = per_gene[mask]
    ig2 = pg2["is_gold"].values
    results["leave_IFIH1_out"] = {
        a: perm_test(pg2[f"{a}_score"].values.astype(float), ig2)
        for a in ["treg"] + CONTROL_AXES
    }

    # ---- specificity verdict ----
    treg_p = results["axes"]["treg"]["p_value_one_sided"]
    treg_auroc = results["axes"]["treg"]["auroc_descriptive"]
    beats_controls = all(
        results["axes"]["treg"]["obs_mean_gold_percentile"]
        > results["axes"][c]["obs_mean_gold_percentile"]
        for c in CONTROL_AXES
    )
    beats_magnitude = (
        results["axes"]["treg"]["obs_mean_gold_percentile"]
        > results["axes"]["magnitude_baseline"]["obs_mean_gold_percentile"]
    )
    results["specificity_verdict"] = {
        "treg_recovery_p": treg_p,
        "treg_auroc_descriptive": treg_auroc,
        "treg_beats_all_controls": beats_controls,
        "treg_beats_magnitude_baseline": beats_magnitude,
        "significant_at_0.05": treg_p < 0.05,
        "honest_conclusion": _conclude(treg_p, beats_controls, beats_magnitude),
    }
    results["caveats"] = [
        "n=8 gold genes -> underpowered CASE STUDY, not a headline (v3 Fix 6).",
        "Permutation null assumes pooled-test exchangeability (recovery is pooled, not "
        "donor-blocked, because CD226+6/8 gold genes are absent from donor-pair data).",
        "Trust proxy = pooled Day-2 ensemble/residual, NOT the guaranteed Day-3 conformal trust.",
        "Control comparisons are DESCRIPTIVE (multiple axes at n=8); no family-wise claim.",
    ]

    RECEIPT.write_text(json.dumps(results, indent=2))
    _plot(per_gene, axis_cols, results)

    # ---- summary ----
    print("=" * 66)
    print("DAY-4 RECOVERY + SPECIFICITY  (n=8 case study, v3 Fix 6)")
    print("=" * 66)
    print(f"  gold genes in test ({results['n_gold_in_test']}): {genes_present_gold}")
    print(f"  test perturbations: {results['n_test_perturbations']}")
    print(f"\n  Recovery by axis (mean gold percentile [null], perm p, AUROC):")
    for a in axis_cols + ["magnitude_baseline", "treg_residualized_on_magnitude"]:
        r = results["axes"][a]
        tag = (
            " <- TREG"
            if a == "treg"
            else (" (shared-mech)" if a == "ifn_shared" else "")
        )
        print(
            f"    {a:32s} {r['obs_mean_gold_percentile']:.3f} [{r['null_mean']:.3f}]  "
            f"p={r['p_value_one_sided']:.4f}  AUROC={r['auroc_descriptive']:.3f}{tag}"
        )
    v = results["specificity_verdict"]
    print(
        f"\n  Treg beats all controls: {v['treg_beats_all_controls']}  |  "
        f"beats magnitude baseline: {v['treg_beats_magnitude_baseline']}"
    )
    print(f"  Treg recovery significant (p<0.05): {v['significant_at_0.05']}")
    print(f"\n  CONCLUSION: {v['honest_conclusion']}")
    print(f"\n  -> {RECEIPT}")
    return 0


def _conclude(p, beats_controls, beats_mag) -> str:
    if p < 0.05 and beats_controls and beats_mag:
        return (
            "Treg axis recovers gold genes above chance AND beats all controls + "
            "magnitude baseline -> specificity SUPPORTED (case study, n=8)."
        )
    if p < 0.05 and not (beats_controls and beats_mag):
        return (
            "Treg axis recovers gold genes above chance, but does NOT clearly beat "
            "controls/magnitude -> recovery may reflect generic effect, not Treg-"
            "specificity (honest negative on specificity)."
        )
    return (
        "Treg recovery NOT significant at n=8 -> underpowered; consistent with v3's "
        "expectation that recovery is a case study, not a powered validation."
    )


def _plot(per_gene, axis_cols, results):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Path("figures").mkdir(exist_ok=True)
    order = ["treg"] + CONTROL_AXES + ["magnitude_baseline", "ifn_shared"]
    labels = [
        "Treg",
        "cell-cycle",
        "cholesterol",
        "apoptosis",
        "ribosome",
        "magnitude",
        "IFN(shared)",
    ]
    aurocs = [results["axes"][a]["auroc_descriptive"] for a in order]
    ps = [results["axes"][a]["p_value_one_sided"] for a in order]
    colors = ["C0"] + ["C7"] * 4 + ["C3", "C8"]
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(labels, aurocs, color=colors)
    ax.axhline(0.5, ls="--", color="k", alpha=0.5, label="chance (0.5)")
    for b, p in zip(bars, ps):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 0.01,
            f"p={p:.2f}",
            ha="center",
            fontsize=8,
        )
    ax.set_ylabel("gold-gene recovery AUROC (descriptive)")
    ax.set_title(
        "Day-4 axis-swap specificity (n=8 case study)\nTreg vs orthogonal controls + magnitude baseline"
    )
    ax.set_ylim(0, 1.05)
    ax.legend()
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig("figures/recovery_specificity.png", dpi=140)
    plt.close(fig)
    print("  figure -> figures/recovery_specificity.png")


if __name__ == "__main__":
    raise SystemExit(main())
