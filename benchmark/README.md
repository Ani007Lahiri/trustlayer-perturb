# Calibration & Selective-Prediction Benchmark for Perturbation Prediction (public v2)

**The axis no existing perturbation benchmark reports.** OP3, PerturBench, PertEval-scFM, VCBench, and the Arc Virtual Cell Challenge all score *point accuracy*. None reports, as its primary axis, whether a predictor's per-prediction **confidence is calibrated** and whether **selective abstention** improves the risk taken on. This benchmark does — and it is **model-agnostic**: the same conformal trust gate wraps any predictor with `fit`/`predict`.

**v2 adds** a third public dataset (Datlinger 2017 CROP-seq Jurkat T-cell, GSE92872) and a fifth, structurally-different predictor — a **GEARS-style relational predictor**.

Pre-registered: `benchmark_prereg.json` (sha256 `3a41ee682ea903c1…`), metric defs frozen before running. Machinery is standard and cited — split conformal (Vovk; Angelopoulos & Bates 2021), locally-adaptive nonconformity (Papadopoulos 2008; Lei 2018). **Calibration/selective-prediction contribution, not a new accuracy SOTA.**

## The 5 predictors
- `mean_baseline`, `linear` (Ridge), `hist_gbr`, `random_forest` — generic feature regressors.
- `relational_gears_style` — **GEARS-STYLE relational**: predicts a held-out perturbation's effect as a co-expression-similarity-weighted average of *related* training perturbations (leave-one-perturbation-out); similarity from target-gene co-expression in **control (unperturbed) cells** — leakage-free. This captures GEARS' core inductive bias ("generalize from related perturbations"). **It is NOT GEARS itself** — no graph neural network, no GO graph, no cell-level training. Full GEARS needs a GPU + torch/PyG; this is the laptop-feasible analog of its mechanism, evaluated on the dataset with full control-cell data (Datlinger).

## The 3 datasets
Same target throughout: **y = log1p(‖z_trans‖)** (L2 norm of the Welch-z trans-response vector, on-target gene excluded). Calibration is assessed **within** each dataset — never by comparing raw y across datasets (different gene universes).

### Gladstone CD4+ T-cell Perturb-seq (T1D) — well-powered (committed folds)
| predictor | coverage 80/90/95 | cal-err | AURC | null | beats | R² | ρ |
|---|---|---|---|---|---|---|---|
| mean_baseline | 0.795/0.897/0.947 | 0.003 | 0.0947 | 0.1110 | ✓ | -0.000 | 0.000 |
| linear | 0.798/0.892/0.946 | 0.004 | 0.0896 | 0.1064 | ✓ | 0.082 | 0.353 |
| hist_gbr | 0.792/0.891/0.948 | 0.005 | 0.0885 | 0.1054 | ✓ | 0.099 | 0.357 |
| random_forest | 0.791/0.892/0.948 | 0.007 | 0.0896 | 0.1069 | ✓ | 0.074 | 0.351 |

### Norman 2019 K562 (public, GSE133344) — small-n (n=105, 200 splits)
| predictor | coverage 80/90/95 | cal-err | AURC | null | beats | R² | ρ |
|---|---|---|---|---|---|---|---|
| mean_baseline | 0.799/0.903/0.961 | 0.068 | 0.2667 | 0.2659 | ✗ | -0.041 | nan |
| linear | 0.794/0.896/0.955 | 0.067 | 0.2693 | 0.2705 | ✓ | -0.079 | -0.045 |
| hist_gbr | 0.772/0.884/0.953 | 0.071 | 0.2597 | 0.2674 | ✓ | -0.064 | 0.086 |
| random_forest | 0.394/0.500/0.596 | 0.264 | 0.3156 | 0.3265 | ✓ | -0.597 | -0.044 |

### Datlinger 2017 CROP-seq Jurkat T-cell (public, GSE92872) — small-n (n=64, 200 splits)
| predictor | coverage 80/90/95 | cal-err | AURC | null | beats | R² | ρ |
|---|---|---|---|---|---|---|---|
| mean_baseline | 0.805/0.909/0.970 | 0.077 | 0.1702 | 0.1716 | ✓ | -0.079 | 0.000 |
| linear | 0.782/0.887/0.953 | 0.080 | 0.1427 | 0.1418 | ✗ | 0.247 | 0.302 |
| hist_gbr | 0.805/0.909/0.970 | 0.077 | 0.1702 | 0.1714 | ✓ | -0.079 | 0.000 |
| random_forest | 0.533/0.654/0.773 | 0.200 | 0.1565 | 0.1556 | ✗ | 0.085 | 0.422 |
| relational_gears_style | 0.826/0.919/0.972 | 0.085 | 0.1307 | 0.1268 | ✗ | 0.385 | 0.607 |

## What the numbers say (honestly)
- **Model-agnostic calibration holds across contexts:** on the well-powered Gladstone data all predictors get valid coverage (±0.03); on both small-n public sets 3/4 generic predictors + the relational predictor stay valid within the small-n tolerance.
- **The benchmark has teeth:** RandomForest's calibration **breaks on both small-n datasets** (coverage 0.39 / 0.53 at nominal 0.80) — the same failure mode caught twice. A benchmark that never fails anything isn't measuring anything.
- **The GEARS-style relational predictor is the strongest real learner on the T-cell data** (R²=0.385, Spearman=0.607) — beating linear (0.25) and the tree models. Borrowing from *related* perturbations genuinely helps predict held-out T-cell perturbation effects, and the trust gate **calibrates it** (coverage 0.83/0.92/0.97). That is the question v2 was built to answer: the gate wraps a structurally different predictor, not just feature regressors.
- **Honest small-n caveat:** at n=64–105 the selective-prediction AURC-vs-null margins are within noise (several "✗"), so the selective axis is only convincingly demonstrated on the well-powered Gladstone data.
- **A confound was removed pre-run:** Datlinger's `log_n_cells` correlates −0.95 with y purely because Welch z scales with 1/√n_cells — a target-construction artifact, not signal. Dropping it took linear R² from a spurious 0.92 to an honest 0.25.

## Reproduce / extend
```bash
python benchmark/calibration_benchmark.py   # writes benchmark_receipt.json
```
```python
from benchmark.calibration_benchmark import load_gladstone, run_predictor_fixed_split
class MyModel:
    def fit(self, X, y): ...; return self
    def predict(self, X): return yhat
d,_ = load_gladstone()
res = run_predictor_fixed_split(*d['train'], *d['calib'], *d['test'], MyModel())
```

## Honest scope
Three datasets, five predictors, one relational (non-GNN) analog of GEARS — a **seed** benchmark, not yet the field standard. Earning "the first calibration benchmark for perturbation biology" outright still needs more public datasets, a real GNN predictor (GEARS/scGPT on GPU), and an open submittable leaderboard. Receipt sha256 `a0a16f317df03aa6…`.
