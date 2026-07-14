# OpenProblems PR — ready to open (you own the GitHub account)

## Where
Repo: https://github.com/openproblems-bio/task_perturbation_prediction
(or the current perturbation-prediction task repo; confirm the active one)

## Branch + PR title
Branch:  `feat/conformal-coverage-metric`
Title:   `Add conformal coverage metric + optional interval-output extension (design issue + component)`

## PR body (paste as-is)
---
### Motivation
All current metrics on this task score **point-prediction accuracy**. None score whether a
method's stated **uncertainty is calibrated** — yet for perturbation prediction, where even
strong models do not reliably beat linear baselines (Ahlmann-Eltze, Huber & Anders,
Nat Methods 2025, doi:10.1038/s41592-025-02772-6), the actionable question is often
"*which predictions should I trust?*" not "*how accurate is the point estimate?*".

### What this PR proposes
A **conformal coverage** metric (coverage@0.90 + undercoverage gap). Because coverage is
undefined on a point estimate, it needs a small, backward-compatible extension to the
method-output contract: an **optional** `interval_lower` / `interval_upper` layer (or a
quantile-sample layer). Point-only methods are simply not scored on it.

### Contents
- `src/metrics/conformal_coverage/` — Viash component + script (in this PR).
- A reference split-conformal interval-emitting control method to exercise the metric.

### Evidence this metric surfaces something real
In a pre-registered study (5 architectures × a 4-level shift ladder, public perturb-seq),
conformal coverage@0.90 is model-independent (across-model spread 0.03) up to cross-donor
shift, but the spread blows up to 0.44 under cross-dataset shift — a failure mode invisible
to accuracy metrics. Full receipts + code: <trustlayer-perturb repo URL>.

### Ask
Opening as a **design issue first**: would maintainers accept the optional interval-output
extension before the component is merged? Happy to iterate on the API shape.
---

## Files to include in the branch
Copy from the release bundle: `release/openproblems_pr/coverage/config.vsh.yaml`,
`release/openproblems_pr/coverage/script.py`, and reference the design doc
`release/openproblems_pr/PR_DESIGN_ISSUE.md`.
