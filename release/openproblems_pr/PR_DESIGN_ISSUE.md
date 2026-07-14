# Design issue + PR: add a calibration/coverage metric to task_perturbation_prediction

## Motivation
All 5 current metrics on `task_perturbation_prediction` score POINT-prediction accuracy.
None score whether a method's stated uncertainty is *calibrated* — yet for perturbation
prediction, where even strong models do not reliably beat linear baselines
(Ahlmann-Eltze, Huber & Anders, Nat Methods 2025, doi:10.1038/s41592-025-02772-6), the
actionable question is often "when should I trust a prediction?" not "how accurate is it?".

## Proposal
Add a **conformal coverage** metric (coverage@0.90 + undercoverage gap). Because coverage
is undefined on a point estimate, this requires a small, backward-compatible extension to
the method-output contract: an OPTIONAL `interval_lower` / `interval_upper` layer (or a
quantile-sample layer). Point-only methods are simply not scored on it.

## What this PR contains
- `src/metrics/conformal_coverage/` (Viash component + script, unit-testable).
- A reference interval-emitting control method (split-conformal wrapper) to exercise it.

## Evidence this metric surfaces something real
In a pre-registered study (5 architectures × 4-level shift ladder, public perturb-seq),
conformal coverage@0.90 is model-independent (spread 0.03) up to cross-donor shift but the
across-model spread blows up to 0.44 under cross-dataset shift — a failure mode completely
invisible to accuracy metrics. Receipt + code: <trustlayer-perturb repo>.

## Status / ask
Opening as a design issue first: would maintainers accept the optional interval-output
extension? Happy to iterate on the API shape before the component is merged.
