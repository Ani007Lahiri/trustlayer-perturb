"""
Figures for the next-run items 2 & 3 (all numbers read from committed receipts).
Left  : reliability ceiling — model Spearman vs donor-block ceiling (matched quantity).
Right : rule-ablation flip matrix on the demo trio.
Writes figures/reliability_ceiling.png and figures/rule_ablation.png.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CEIL = Path("data/gold/reliability_ceiling_receipt.json")
ABL = Path("data/gold/rule_ablation_receipt.json")
FIGDIR = Path("figures")
FIGDIR.mkdir(exist_ok=True)


def fig_ceiling() -> None:
    r = json.loads(CEIL.read_text())
    p = r["PRIMARY_matched_quantity_comparison"]
    model = p["model_test_spearman"]
    ceil = p["ceiling_r_mean"]
    lo, hi = p["ceiling_r_ci95"]

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    # ceiling band
    ax.axhspan(lo, hi, color="#c9d7e8", alpha=0.7, label="donor-block ceiling 95% CI")
    ax.axhline(ceil, color="#2b5f9e", lw=2, label=f"ceiling Spearman = {ceil:.3f}")
    ax.bar(
        [0], [model], width=0.45, color="#e08214", label=f"model Spearman = {model:.3f}"
    )
    pct = p["model_pct_of_rank_ceiling_mean"]
    ax.annotate(
        f"{pct:.0f}% of ceiling\n(range {p['model_pct_of_rank_ceiling_range_over_ci'][0]:.0f}"
        f"–{p['model_pct_of_rank_ceiling_range_over_ci'][1]:.0f}% over CI)",
        xy=(0, model),
        xytext=(0.35, ceil + 0.03),
        fontsize=9,
        ha="left",
        arrowprops=dict(arrowstyle="->", color="#555"),
    )
    ax.set_xlim(-0.5, 1.4)
    ax.set_ylim(0, max(hi, model) + 0.12)
    ax.set_xticks([0])
    ax.set_xticklabels(["covariate\nmodel"])
    ax.set_ylabel("Spearman correlation\n(matched: log1p trans-only ‖zscore‖)")
    ax.set_title(
        "The model is near the assay's own reproducibility ceiling\n"
        "modest accuracy = noisy assay, not a weak model",
        fontsize=10,
    )
    ax.legend(fontsize=8, loc="lower right", framealpha=0.9)
    fig.tight_layout()
    out = FIGDIR / "reliability_ceiling.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"-> {out}")


def fig_ablation() -> None:
    r = json.loads(ABL.read_text())
    baseline = r["baseline_verdicts"]
    single = r["single_rule_ablation"]
    genes = list(baseline.keys())
    rules = list(single.keys())

    # matrix: rows=rules, cols=genes; 1 if that gene FLIPS when the rule is removed
    M = np.zeros((len(rules), len(genes)))
    for i, rule in enumerate(rules):
        for j, g in enumerate(genes):
            if g in single[rule]["changed"]:
                M[i, j] = 1

    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    ax.imshow(M, cmap="Oranges", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([f"{g}\n({baseline[g]})" for g in genes], fontsize=9)
    ax.set_yticks(range(len(rules)))
    ylabels = [f"{rule}\n[{single[rule]['label'].split(' ')[0]}]" for rule in rules]
    ax.set_yticklabels(ylabels, fontsize=9)
    for i, rule in enumerate(rules):
        for j, g in enumerate(genes):
            if M[i, j]:
                ax.text(
                    j,
                    i,
                    f"{baseline[g]}→GO",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="black",
                )
    ax.set_title(
        "Rule ablation: which rule flips which verdict\n"
        "genetics + eQTL are pivotal; trust & leakage are inactive",
        fontsize=10,
    )
    fig.tight_layout()
    out = FIGDIR / "rule_ablation.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"-> {out}")


if __name__ == "__main__":
    fig_ceiling()
    fig_ablation()
