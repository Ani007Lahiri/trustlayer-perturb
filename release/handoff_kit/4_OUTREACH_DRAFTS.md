# Outreach drafts — YOU send these from your own account. I do not send them, and any
# reply must actually arrive from a real person; do not treat these as adopter interest until
# a genuine response comes back.

## A. OpenProblems maintainers (after the PR is open)
Subject: Coverage/calibration metric for task_perturbation_prediction — design issue + component

Hi [maintainer name],

I opened a PR/design issue proposing an optional conformal-coverage metric for
task_perturbation_prediction: [PR URL]. The motivation is that all current metrics score
point accuracy, but since strong models don't reliably beat linear baselines here
(Ahlmann-Eltze/Huber 2025), the actionable question is often *which predictions to trust*.
In a pre-registered study I found conformal coverage@0.90 is model-independent up to
cross-donor shift but breaks under cross-dataset shift (spread 0.03 -> 0.44) — invisible to
accuracy metrics.

One concrete question: would an optional interval-output extension to the method contract be
acceptable in principle, so coverage can be scored without breaking point-only methods?

Happy to iterate on the API shape. Code + receipts: [repo/DOI].
Thanks, Anirudh

## B. A perturb-seq / T-cell-engineering lab (pick one whose work you know)
Subject: A "commit gate" for deciding which perturbation to validate next — worth a look?

Hi [PI name],

I built a small open tool that wraps any perturbation-effect predictor with calibrated
uncertainty and a default-deny "commit gate": it decides which predictions are trustworthy
enough to act on, and abstains under distribution shift unless you supply a few cheap
labeled anchor perturbations to recalibrate. In a retrospective simulation, trust-gated
selection of validation targets recovered [X — insert the decision-value headline once the
simulation lands] more true strong perturbations at a fixed budget than naive top-K by
predicted effect.

Given your group's perturb-seq work, one question: would a trust gate like this change which
targets you'd prioritize for arrayed validation? If it's useful I'd love feedback on the gate
criteria. `pip install trustlayer-perturb`; details + figures: [repo/DOI].

Thanks for your time, Anirudh

---
INTEGRITY NOTE (for the submission): only quote a named person's interest if they actually
reply expressing it. Do not represent these drafts as endorsements.
