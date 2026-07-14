# TrustLayer — Decision-First Reframe & 90-Second Demo Storyboard

## The hero thesis (opening line, everywhere)
**"Foundation models don't beat linear baselines at predicting perturbation effects
(Ahlmann-Eltze, Huber & Anders, Nat Methods 2025). So the bottleneck isn't prediction
accuracy — it's knowing WHICH prediction to trust before you spend scarce wet-lab budget.
TrustLayer is the commit gate that decides which perturbation to validate next."**

## The three concrete calls (hero slide)
| Target | Gate decision | Why |
|---|---|---|
| **CD226** | **GO** | leakage-clean, genetic association 0.34 ≥ 0.20, trust ≥ 0.50, cell-type-matched eQTL, direction consistent |
| **RASGRP1** | **ABSTAIN** | no evidence either way — gate refuses to guess |
| **PRKCQ** | **WITHHOLD** | vetoed on the genetic floor (GA 0.162 < 0.20) |

## The headline number (Impact)
**"At a 200-target validation budget, TrustLayer's default-deny abstention on the
uncalibrated (cross-dataset) regime — where naive top-K precision collapses against a
0.20 base rate — avoids ~160–200 wasted arrayed-CRISPR screens (~$160k–$1M at $1–5k/target),
and its label-free gate flags false leads at 1.36× the GO rate across all 6 donor folds."**
(Honest: under MILD cross-donor shift the predictor is already calibrated, so hard gating
adds no hits — an honest null; the value is avoided waste + the false-lead flag.)

## What makes each rubric axis land

**Claude Use (the standout).** Claude ran the entire self-correcting research loop:
council → hypothesis → **pre-registration frozen+hashed** → **self-falsification of its own
headline** → **red-team caught its own oracle-handoff bug** → corrected protocol → a **fresh
red-team dated today** that surfaced a real caveat and reconciled it. The **Honesty Ledger**
renders all 8 nodes with verifiable content hashes.

**Depth.** Pre-registered 6-predictor × 4-shift-level calibration law — now including a
**real Geneformer-V2-104M foundation model** (not a proxy): model-independence SURVIVES
(cross-donor 6-model spread 0.021), cross-dataset breakdown reproduces. Donor-clustered
bootstrap CIs on every headline. 27/27 pytest green + gate ablation (all 5 conditions
load-bearing) + CI workflow.

**Impact.** The decision-value simulation + savings number + three concrete gate calls,
positioned against the field's own benchmark (OpenProblems) with a ready-to-open PR.

## 90-second demo storyboard (narrated screencast)
- **0:00–0:12** — Hook: "FMs don't beat linear baselines here. So the question isn't
  *what's the effect* — it's *should I believe this prediction*." Show the three calls.
- **0:12–0:30** — The calibration law figure: 6 predictors incl. Geneformer converge at
  cross-donor (spread 0.021), fan out at cross-dataset. "Calibration is a property of the
  shift, not the model — and it holds even for a foundation model."
- **0:30–0:48** — The Honesty Ledger, scrolling: "Claude pre-registered, then FALSIFIED its
  own headline; its red-team caught its own oracle-handoff bug. Every claim links to a hashed
  receipt — click and verify." Click one node → receipt hash.
- **0:48–1:05** — The decision-value curve: "At a 200-target budget, trust-gating avoids
  ~$160k–$1M in wasted screens." Show the savings-vs-budget curve.
- **1:05–1:20** — Transport result: "Does the trust gap transport for free? We tested it on
  9 datasets — it doesn't (Spearman +0.32, n.s.). The fix: a 5–20% labeled anchor. We
  falsified our own convenient assumption and shipped the honest protocol."
- **1:20–1:30** — Close: "`pip install trustlayer-perturb`. Pre-registered, hash-verified,
  self-red-teamed by Claude. The trust layer for perturbation biology."
