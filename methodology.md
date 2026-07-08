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

---

# Methodology — Day 3: Calibrated Trust Layer (the product)

*The Day-2 predictor exists to be calibrated. Day 3 builds the actual contribution:
a conformal trust score with a coverage guarantee, headlined via DONOR-blocked
leave-one-donor-out (LODO), per v3 Fix 6.*

## Research Question & Hypothesis (Day-3 scope)
- **Question:** Does a conformal prediction interval on the trans-effect prediction achieve
  its target coverage under an HONEST donor-blocked split, and does selective prediction
  (abstaining on low-confidence calls) monotonically reduce risk?
- **Hypothesis:** MAPIE split-conformal intervals calibrated on held-out donors achieve
  ~90% empirical coverage at nominal 90%, and the selective-risk (AURC) curve shows that
  filtering to high-trust predictions lowers error — i.e., the trust score is USABLE for
  triage.
- **Success criteria (critique-aligned):**
  1. The Clopper-Pearson interval on pooled donor-blocked coverage CONTAINS the nominal
     level (80/90/95), with the CP interval WIDTH reported as a first-class result. At
     n_eff ~ tens the CP interval is expected to be WIDE (~12 test points pooled) — so
     "contains nominal" is reported alongside width, not as a discrimination claim.
  2. Selective-risk curve below the no-selection baseline on held-out-donor test pairs,
     with the shuffle-residual control curve FLAT (guards circularity).
  3. Coverage assessed under DONOR-blocked LODO (harder, honest regime) AND
     perturbation-blocked (both reported; donor-block is the headline).

## Data Sources (added Day 3)
- `data/raw/GWCD4i.DE_stats.by_donors.h5mu` (16.87 GB): per-donor-pair DE modalities
  (one AnnData per disjoint donor pair). Enables donor-blocked calibration: calibrate on
  donor-pairs excluding a held-out donor, test on the held-out donor's pairs.
- Reuse: frozen_splits.json (perturbation blocking), base_predictor_preds.parquet (residuals).

## Analysis Pipeline
### Step 1: Deep ensemble (3-5 seeds) of the base predictor
- Refit the HGBR base predictor at 3-5 random seeds; ensemble mean = point prediction,
  ensemble spread = a heuristic uncertainty (feeds the conformal nonconformity score).
- Rationale (simplicity): a small ensemble is the lightest way to get a variance signal
  without deep nets. If spread adds nothing over absolute-residual CQR, drop it (log as
  negative result).

### Step 2: MAPIE split-conformal (SplitConformalRegressor, MAPIE 1.4.1)
- Nonconformity: absolute residual (and CQR variant if quantile base is added).
- Calibrate on the calibration fold; evaluate coverage + interval width on test.
- Target levels: 0.80, 0.90, 0.95.

### Step 3: DONOR-blocked LODO conformal (the headline, Fix 6) — CRITIQUE-CORRECTED
**Exchangeability fix (Gate-1 BLOCKING #1):** donor-pairs SHARE donors, so leaving out
only pairs *containing* d is necessary but NOT sufficient — calibration pairs would still
share the other donors with test pairs, breaking exchangeability. Correct construction:
  - Held-out donor d: **calibration set = pairs whose BOTH members != d** (drawn only from
    the other 3 donors); **test set = pairs containing d.**
  - Assert donor-disjointness: no member of any calibration pair equals d (fail-closed).
  - With 4 donors (6 disjoint pairs), holding out d leaves C(3,2)=3 clean calibration pairs
    and 3 test pairs.
- Pre-register effective-N per fold (independent perturbations), NOT raw 33,983.

### Step 4: Selective prediction / coverage curves — CRITIQUE-CORRECTED
**Circularity fix (Gate-1 BLOCKING #3):** three-way separation. Conformal quantiles/widths
are fit ONLY on the calibration set (Step 3). The trust score AND the selective-risk/AURC
curve are computed ONLY on the held-out-donor TEST pairs — never on calibration data; no
width/threshold is tuned on test pairs.
- Trust score = f(interval width, ensemble spread), evaluated on test pairs only.
- Sort test predictions by trust; plot risk vs coverage; compute AURC.
- Competing curve: mean/shuffle baseline (must be worse).

## Controls & Validation
- **Coverage calibration plot:** empirical vs nominal coverage (diagonal = perfect).
- **Negative control:** shuffle residuals -> selective-risk curve should be FLAT (no gain).
- **Donor-block vs perturbation-block:** report both; donor-block is the honest headline.

## Statistical Plan — CRITIQUE-CORRECTED (Gate-1 BLOCKING #2 + Q4)
- **Power honesty:** at n_eff ~ tens with 4 LODO folds, coverage estimates CANNOT
  discriminate 80 vs 90 vs 95%. So the claim is DEMOTED from "within-tolerance calibration
  at 3 levels" to **"the Clopper-Pearson interval contains the nominal level"**, and the CP
  interval WIDTH is reported as a first-class result (not a footnote). Coverage is POOLED
  across folds rather than claimed per-level-per-fold.
- **Selective risk:** AURC on held-out-donor test pairs. Report Clopper-Pearson/CP as the
  primary uncertainty statement. A donor-block bootstrap is noted but has near-zero
  resolution at 4 clusters (~35 distinct resamples) — reported as a limitation, not a
  headline CI. Pair-level bootstrap understates variance (pairs non-independent) — avoided
  as primary.
- Honesty: every interval reported with the effective-N caveat; recovery (Day 4) stays a
  CI-shown case study.

## Compute Requirements
- All CPU / laptop (small ensemble + MAPIE + curves). No GPU. No spend.

## Limitations & Assumptions
- Conformal validity assumes exchangeability WITHIN the calibration/test block; donor
  blocking is the honest way to approximate this given only 4 donors.
- The ensemble uncertainty is heuristic, not Bayesian; the conformal GUARANTEE is what
  carries the coverage claim, not the ensemble spread.

---

# Methodology — Day 4: Blinded Recovery + Axis-Swap Specificity Control

*The coverage headline (Day 3) proves calibration. Day 4 tests two orthogonal claims:
(A) does the trust layer, applied blind, RECOVER known T1D genes it never saw in
calibration? and (B) is the Treg axis SPECIFIC to Treg biology, or a generic
"big-effect" detector? Per v3 Fix 6, recovery is an underpowered (n=8) CI-shown case
study, NOT a headline.*

## Research Questions & Hypotheses
- **Q-A (recovery):** Ranking held-out perturbations by (Treg-axis effect x trust), do the
  8 genetics-anchored T1D genes (never in train/calibration) rank higher than chance?
- **Q-B (specificity):** Is the Treg-axis signal distinguishable from generic transcriptional
  magnitude? A gene that just "does a lot" should NOT automatically score high on the Treg
  axis unless it genuinely moves Treg biology.
- **H-A:** The 8 gold genes rank above the median of all perturbations on the combined
  Treg-axis-x-trust score, but with WIDE CIs (n=8) -> case study, not proof.
- **H-B:** The Treg axis correlates with gold-gene recovery MORE than >=3 orthogonal
  non-Treg control axes (e.g. cell-cycle, interferon, apoptosis) do. If a random/off-target
  axis recovers gold genes equally well, the Treg axis is not specific -> honest negative.

## Data Sources
- Pooled DE_stats (data/interim/day2_feature_table.parquet + treg_axis_scores.parquet):
  recovery runs on POOLED data because CD226 + 6/8 gold genes are absent from donor-pair
  data (by_donors_facts.json).
- frozen_splits.json: 8 gold genes are in the TEST fold only (never calibration).
- base_predictor_preds.parquet: trust proxy from Day 2 (upgraded conceptually by the
  Day-3 conformal machinery; here we use the ensemble-spread trust on pooled test).

## Analysis Pipeline
### Step 1: Blinded recovery ranking (Q-A) — CRITIQUE-CORRECTED
- For every perturbation in the TEST fold, compute a nomination score =
  |Treg-axis score| combined with trust (trust = POOLED Day-2 ensemble spread, stated
  plainly — NOT the guaranteed Day-3 conformal trust, which is donor-blocked and excludes
  6/8 gold genes).
- Rank test-fold perturbations; locate the 8 gold genes' ranks.
- **Primary metric (small-N appropriate):** Mann-Whitney U / rank-sum test of gold vs
  non-gold nomination percentiles, with an EXACT permutation null (permute the 8 gold
  labels across all test perturbations, >=10000 perms). AUROC reported as a descriptive
  monotone summary of U, NOT as the inferential statistic.
- Report the permutation p-value + CI. Explicitly n=8 underpowered case study.

### Step 2: Axis-swap specificity control (Q-B)  [the real scientific control] — CRITIQUE-CORRECTED
- Build >=3 orthogonal control axes on the SAME zscore matrix, each gold-disjoint with
  per-row on-target exclusion (same construction as the Treg axis):
    * cell-cycle / proliferation axis
    * **cholesterol/metabolic axis** (orthogonal replacement for interferon)
    * apoptosis axis
    * ribosome-biogenesis axis (4th, for margin)
- **BLOCKING FIX #1 — interferon is NOT a clean control:** IFIH1 (a gold gene) is a
  canonical IFN-pathway gene, so an IFN axis shares mechanism with the positives. IFN is
  therefore DROPPED as a negative control. If reported at all, it is labeled a
  "shared-mechanism positive-leaning control", and every axis also gets a
  leave-IFIH1-out sensitivity check.
- **BLOCKING FIX #2 — magnitude confound CONTROLLED (not just detected):**
    (a) MAGNITUDE-ONLY baseline axis = rank by trans_effect_magnitude alone; the Treg axis
        must beat this baseline, else recovery is a "big-effect" artifact.
    (b) Partial out magnitude: rank on residual of Treg-axis regressed on
        trans_effect_magnitude; re-run the rank-sum test on residualized scores.
- Specificity holds ONLY IF the Treg axis beats each control AND the magnitude-only baseline
  by a margin that survives the n=8 permutation test. Otherwise report an honest negative.

## Controls & Validation
- **Random-axis control:** a random signed gene set (matched size) -> expected null.
- **Magnitude-only baseline:** explicit competing axis (Fix #2a) the Treg axis must beat.
- **Leave-IFIH1-out sensitivity:** re-run all axes dropping IFIH1 (Fix #1).
- **Leakage:** gold genes never in calibration; all axes gold-disjoint; per-row on-target
  exclusion (all enforced + asserted).

## Statistical Plan — CRITIQUE-CORRECTED
- Primary: Mann-Whitney U rank-sum (gold vs non-gold) with EXACT permutation p-value
  (>=10000 label permutations). AUROC = descriptive monotone summary only.
- Specificity: Treg-vs-control and Treg-vs-magnitude-baseline rank-sum differences.
  Multiple axes -> comparisons are DESCRIPTIVE (pre-registered as such; at n=8 nothing is
  expected to reach corrected significance). Do NOT claim specificity if the permutation
  CI overlaps the null.
- Explicitly label n=8 as underpowered; report effect + CI width prominently.

## Compute Requirements
- CPU / laptop only (gene-set contrasts + ranking + bootstrap). No GPU, no spend.

## Limitations & Assumptions
- n=8 recovery is underpowered; this is a CASE STUDY illustrating the workflow, not a
  powered validation (v3 Fix 6). The headline remains Day-3 coverage.
- Control axes are curated, not exhaustive; a negative specificity result at n=8 does not
  disprove Treg relevance, only that this dataset/N cannot resolve it.
