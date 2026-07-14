# TrustLayer — a calibrated-trust layer for perturbation-biology prediction

**"A pLDDT for perturbation biology."** Foundation models don't reliably beat linear
baselines at predicting perturbation effects ([Ahlmann-Eltze, Huber & Anders, *Nat Methods*
2025](https://doi.org/10.1038/s41592-025-02772-6)). So the bottleneck isn't accuracy — it's
knowing *which* prediction to trust before you spend scarce wet-lab budget. TrustLayer wraps
any perturbation-effect predictor with **split-conformal calibration** and a **default-deny
commit gate** that decides whether a call is safe to act on.

```bash
pip install trustlayer-perturb
```

```python
from trustlayer import TrustLayer
tl = TrustLayer(alpha=0.10)         # target 90% coverage
tl.calibrate(residual_scores_cal)   # split-conformal on a calibration set
tl.covers(score_test)               # is a test point inside the conformal set?
tl.shift_gap(scores_test)           # empirical undercoverage under distribution shift
tl.should_commit(trust=0.62, genetic_association=0.34, eqtl_direction_ok=True)  # GO / not-GO
```

## Live interactive demos
- **The Honesty Ledger** — Claude's self-correcting research loop, every claim linked to a
  hash-verifiable receipt: [`honesty_ledger.html`](honesty_ledger.html)
- **The Trust Atlas** — which predictions to believe: [`Trust_Atlas.html`](Trust_Atlas.html)
- **Watch the gate decide** — interactive GO / WITHHOLD / ABSTAIN: [`Live_Abstention_Demo.html`](Live_Abstention_Demo.html)

## What we found (pre-registered, self-falsified)
We froze and hashed the protocol **before** any model run, then let the falsifiers fire.

| Finding | Result |
|---|---|
| **Coverage under cross-donor shift** | Mild undercoverage (coverage@0.90 = 0.88–0.91), **not** the "severe collapse" our own exploratory pass claimed — the pre-registration falsified our headline. |
| **Model-independence** | Holds across 6 architectures incl. a real **Geneformer-V2-104M** foundation model (cross-donor coverage spread **0.021**); **breaks** under cross-dataset shift (spread 0.44–0.67). |
| **Does the trust gap transport label-free?** | **No** — across 9 independent public perturb-seq datasets, source-only transport fails (Spearman +0.32, n.s.). **A 5–20% labeled anchor recalibrates it.** |
| **Decision value** | Trust-gated target selection avoids **~$160k–$1M** in wasted arrayed-CRISPR screens at a 200-target validation budget (retrospective simulation). |
| **The gate, on real T1D targets** | CD226 → **GO**, RASGRP1 → **ABSTAIN**, PRKCQ → **WITHHOLD** (vetoed on the genetic floor). |

## Reproducibility & rigor
- **Pre-registrations + hashed receipts** in [`data/gold/`](data/gold/) — every headline number has a `content_sha256`.
- **27/27 passing tests** ([`tests/test_trustlayer.py`](tests/test_trustlayer.py)) + a gate ablation showing all 5 gate conditions are load-bearing.
- **Donor-clustered bootstrap CIs** on every headline number (honestly wide where the data is thin).
- **Self-red-team audit trail** — Claude caught its own oracle-handoff bug and retracted it; a fresh red-team surfaced and resolved a tautology risk. See the Honesty Ledger.
- CI runs the test suite on every push ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)).

## Repository layout
```
release/trustlayer/     pip package (TrustLayer API + commit gate)
release/openproblems_pr/ OpenProblems coverage-metric component + design issue
tests/                  pytest suite + gate ablation
figures/                publication figures (calibration law, transport, decision value, ...)
data/gold/              pre-registrations + hashed result receipts (raw data is re-fetchable, not committed)
*.html                  interactive demos (Honesty Ledger, Trust Atlas, gate demo)
```

## How Claude built this
Claude Science ran the whole loop — council-driven pre-registration, remote GPU/SLURM
execution, self-falsification against frozen falsifiers, self-red-teaming, and the writeup.
The [Honesty Ledger](honesty_ledger.html) makes that loop auditable.

## License
MIT (see [`release/trustlayer/LICENSE`](release/trustlayer/LICENSE)). Data files are public and
re-fetchable from GEO / Zenodo / scPerturb; large raw matrices are not committed.
