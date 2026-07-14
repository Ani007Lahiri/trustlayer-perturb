# Preprint abstract (bioRxiv-ready)

**Title:** Calibrated trust for perturbation prediction: conformal coverage is
model-independent under mild shift and fails predictably under severe shift

**Abstract.** Deep models for genetic-perturbation-effect prediction do not reliably
outperform simple linear baselines on unseen and combinatorial perturbations
(Ahlmann-Eltze, Huber & Anders, 2025). This shifts the practical question from *how
accurate is a prediction* to *when should it be believed*. We study calibration rather
than accuracy: we wrap five architecturally distinct predictors (identity, mean, ridge,
gradient-boosted trees, deep ensemble) in split-conformal prediction and evaluate coverage
across a pre-registered four-level distribution-shift ladder (in-distribution, held-out
perturbation, cross-donor, cross-dataset) on public perturb-seq data. We find that
conformal coverage at the 90% level is remarkably **model-independent** (across-model
spread 0.03) and undercoverage is **mild** (≤0.056) up to cross-donor shift, but that
model-independence **breaks** under the most severe cross-dataset shift (spread 0.44), with
some predictors over-covering and others collapsing. A frozen pre-registration overturned
our own earlier, method-dependent claim of severe undercoverage — a caution about
protocol sensitivity in calibration benchmarking. We package the layer, its default-deny
commit gate, and the full diagnostic as open source, and propose a coverage metric for the
OpenProblems perturbation-prediction benchmark. Calibration is a property of the shift
regime more than of the model; the size of the coverage gap, and the point at which models
diverge, is the actionable trust signal.
