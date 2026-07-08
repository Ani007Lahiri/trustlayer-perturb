# Methodology — Day 2: Base Perturbation-Effect Predictor + Treg Program Axis

*Component methodology for the v3 build. The overall research question, hypothesis, and
success criteria are defined in `revised_plans/PLAN_OF_ATTACK_v3.html` and
`research-state.md`. This document specifies the Day-2 predictor + axis only, so the
Day-3 conformal layer (the actual product) has a defensible base to calibrate.*

## Research Question & Hypothesis (Day-2 scope)
- **Question:** Can a deliberately simple, borrowed base predictor produce per-perturbation
  effect predictions good enough that a calibrated trust layer on top is meaningful?
- **Hypothesis:** A pseudobulk-delta feature representation + ESM-2 target embedding +
  DEG priors predicts a scalar perturbation-effect summary (per perturbation-condition)
  with signal above a shuffle/mean baseline — WITHOUT competing on raw accuracy.
- **Non-goal (v3 §1):** We do NOT try to beat GEARS/scGPT/linear baselines. We only need
  a base predictor whose *residuals* the conformal layer can calibrate. Accuracy is a
  saturated game (Ahlmann-Eltze 2025); we compete on decision under uncertainty.
- **Success criteria (Day 2):**
  1. Predictor produces per-(gene,condition) predictions on the frozen splits with NO
     gold-set gene in train/calibration (leakage assertion passes).
  2. Prediction target correlates with a held-out mean/shuffle baseline at level that
     leaves headroom for a trust layer (i.e., not perfect, not zero — the interesting regime).
  3. Treg program axis is computable and is DISJOINT from every gold-set gene name.

## Data Sources
- `data/raw/GWCD4i.DE_stats.h5ad` (16.79 GB, verified): 33,983 (perturbation×condition)
  × 10,282 measured genes. Layers: `log_fc`, `zscore` (DESeq2). Dense, no NaN.
- `data/gold/frozen_splits.json`: hashed train/calib/test partition (gold forced to test).
- `data/gold/t1d_gold_set.json`: 8-gene blinded positive set + axis_forbidden_genes.
- ESM-2 embeddings: precomputed feature prior (v3 demotes sequence→function to a
  precomputed feature, NOT a trained chain). Source: fair-esm2 / HF facebook/esm2_t33.

## Prediction Task (precise) — LEAKAGE-CORRECTED (critique Gate-1, cycle 1)
For each perturbation p in condition c, the base predictor maps target-gene features to a
**summary of the perturbation's TRANS transcriptional effect** (downstream regulation only):
  - `trans_effect_magnitude(p,c)` = L2 norm of the zscore vector over measured genes
    **EXCLUDING the perturbed gene's own component** (i.e. trans genes only).
The primary regression target is `log1p(trans_effect_magnitude)`.

**Why the on-target gene is excluded (BLOCKING fix):** the pre-COMPUTE critique found that
including the perturbed gene in the norm creates target-feature circularity — the on-target
zscore is typically the single largest component of the norm, and any on-target-derived
feature (`ontarget_effect_size`) would then partially determine the target by construction
(`target = sqrt(ontarget² + Σ_trans zscore²)`). Restricting to trans genes removes this
mechanical leak AND is scientifically better: trans-effects (downstream gene regulation)
are exactly what a Treg lab triages on, not the trivial fact that a guide knocks down its
own target.

Secondary/robustness target: `n_downstream(p,c)` from obs (already a trans-only count) —
used only as a cross-check, never as a feature.

Rationale for a SCALAR target (simplicity preference): predicting the full 10,282-dim
delta is the crowded accuracy game we explicitly avoid. A calibrated scalar trans-effect
score is exactly what the trust/decision layer needs and is honestly evaluable at n_eff ~ tens.

## Features (borrowed substrate, v3 §1) — LEAKAGE-CORRECTED
Per target gene. **CRITICAL:** no feature may be derived from the same DE computation as
the target. Specifically `ontarget_effect_size`, `n_up/down_genes`, `n_total_de_genes`,
and `n_downstream` are FORBIDDEN as features (they are outputs of the DE that produces the
target). They may be used ONLY as QC filters, never as model inputs.
1. **ESM-2 embedding** (mean-pooled, 1280-d for esm2_t33) of the target protein — the
   sequence→function prior, precomputed and frozen. This is the primary feature: it is a
   property of the protein's sequence, fully independent of the DE readout.
2. **Perturbation delivery covariates** (independent of the DE readout): `n_cells_target`
   (assay depth for this perturbation), `target_baseMean` (baseline expression BEFORE
   perturbation), condition one-hot, `n_guides` / `single_guide_estimate` (assay design).
   baseMean is a pre-perturbation quantity, not a DE output, so it is admissible.
3. NO pseudobulk-delta feature of the target gene's own post-perturbation profile
   (would leak). Cross-gene network priors (Krogan PPI degree, etc.) may be added later as
   they are perturbation-independent.

**Leakage audit (build-time assertion):** assert the feature matrix column set is disjoint
from {ontarget_effect_size, n_up_genes, n_down_genes, n_total_de_genes, n_downstream,
ontarget_significant, ontarget_effect_category} before fitting.

## Treg Program Axis (v3 Fix 3 dependency) — LEAKAGE-CORRECTED
- UCell-style signed contrast score per perturbation over a curated Treg program gene set,
  computed on the measured-gene delta (zscore) matrix.
- **HARD CONSTRAINT 1:** the program gene set must not include ANY gene in
  `axis_forbidden_genes`. Enforced by assertion at build time. This makes recovery of the
  gold set non-circular (name-level disjoint) — the axis-swap control (Day 4) then tests
  biological-subspace disjointness.
- **HARD CONSTRAINT 2 (per-row on-target exclusion):** when scoring perturbation p, the
  perturbed gene p itself is excluded from the readout columns, mirroring the trans-only
  target. This prevents a perturbation that happens to target a program gene from scoring
  high trivially via its own on-target knockdown.

## Analysis Pipeline
### Step 1: Precompute ESM-2 embeddings for all perturbed target genes
- Method: map gene symbol → UniProt canonical sequence → ESM-2 mean embedding.
- Tools: gget/UniProt for sequences; fair-esm2 (esm2_t33_650M) for embeddings.
- Compute: embeddings for ~11.5k proteins. GPU-accelerated on Modal if local is slow;
  otherwise CPU batch (small model) — decide after timing a batch of 50.
- Output: `data/interim/esm2_target_embeddings.parquet` (gene → 1280-d vector).

### Step 2: Assemble the feature/label table on frozen splits
- Build (gene, condition) rows with features (Step 1 + obs-derived) and targets
  (effect_magnitude, n_downstream) read from the h5ad layers/obs.
- Assert: no gold gene in train/calib; gold genes only in test.
- Output: `data/interim/day2_feature_table.parquet`.

### Step 3: Fit the base predictor
- Model: gradient-boosted trees (HistGradientBoostingRegressor) — simple, CPU-fast,
  strong tabular baseline. (Deep ensemble for conformal comes Day 3.)
- Train on train fold; report performance on calib + test folds.
- Baselines: (a) predict global mean, (b) shuffle-target permutation.

### Step 4: Build + validate the Treg axis
- Compute UCell signed contrast per (gene,condition); assert forbidden-gene disjointness.
- Output: `data/interim/treg_axis_scores.parquet`.

## Controls & Validation
- **Negative control:** shuffle target labels → predictor R² should collapse to ~0.
- **Baseline:** global-mean predictor → the base predictor must beat it (else no signal).
- **Leakage assertion:** programmatic check that gold genes are absent from train/calib
  and from the axis gene set (fail-closed).

## Statistical Plan
- Primary metric: Spearman ρ and R² of predicted vs actual effect_magnitude on the
  held-out test fold (gold genes), with bootstrap 95% CI (donor structure caveated:
  true n_eff is tens, per data_facts_receipt.json).
- Report the mean/shuffle baselines as competing curves (Day 3 makes this the headline).

## Compute Requirements
- Step 1 (ESM-2): if CPU batch of 50 proteins > ~5 min, move to Modal (T4/A10, ~$0.5-2,
  <30 min) — will present cost estimate and get approval BEFORE any GPU spend.
- Steps 2-4: laptop CPU (feature table + HGBR + UCell are all light).

## Limitations & Assumptions
- Scalar effect target is a deliberate simplification; full-profile prediction is out of scope.
- DE_stats is aggregated across 4 donors → donor-blocked LODO (Fix 6) needs the
  by_donors.h5mu (16.87 GB), deferred to Day 3.
- ESM-2 embeddings are a frozen prior; we do not fine-tune (v3: precomputed feature).
- UCell axis gene set is curated from canonical Treg biology, not learned — and must be
  gold-set-disjoint, which slightly weakens it but is required for non-circularity.
