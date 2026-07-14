# TrustLayer — calibrated trust for perturbation prediction

A small, honest layer that answers *"should I believe this perturbation-effect
prediction?"* — not *"what is the effect?"*. It wraps any point predictor with
split-conformal intervals, a default-deny commit gate, and pre-registered
shift-robustness diagnostics.

## Install
```
pip install trustlayer-perturb   # (scaffold; publish from release/trustlayer)
```

## Quickstart
```python
from trustlayer import TrustLayer
tl = TrustLayer(alpha=0.10)          # target 90% coverage
tl.calibrate(scores_cal)             # nonconformity scores on a calibration split
cov = tl.empirical_coverage(scores_test)
gap = tl.shift_gap(scores_test)      # undercoverage under shift; abstain if large
go  = tl.should_commit(trust=0.62, genetic_association=0.34, eqtl_direction_ok=True)
```

## What the evidence says (pre-registered, sha256 `3799ac42…`)
Across **5 architecturally-distinct predictors × a 4-level shift ladder** on public
perturb-seq data, split-conformal coverage@0.90 is **model-independent (spread 0.03)**
up to cross-donor shift and undercoverage is **mild (≤0.056)**; model-independence
**breaks only under the most severe cross-dataset shift** (spread 0.44). A frozen
pre-registration *overturned* an earlier, method-dependent "severe undercoverage"
headline — the diagnostic caught our own overclaim. Full receipt + red-team report ship
with the repo.

## Honest scope
- The conformal guarantee is exchangeability-conditional; under shift it is a
  *diagnostic*, not a guarantee. `shift_gap()` is how you know to abstain.
- This is a methods/diagnostics contribution, not a wet-lab-validated discovery.

## License
MIT. See `CITATION.cff` for citation metadata (Zenodo DOI on first tagged release).
