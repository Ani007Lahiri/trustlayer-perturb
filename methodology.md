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

---

# Methodology — Day 5: Critic -> Structured Verdict -> commit_gate (end-to-end veto)

*Days 0-4 built the pieces: genetics gates, gold set, base predictor, Treg axis, the
Day-3 conformal trust, and the deterministic commit_gate (Day 1, 7 tests). Day 5 wires
them into ONE end-to-end pipeline where a critic produces a structured CriticVerdict per
target, and the deterministic gate is the SOLE writer of a nomination. The demo's
believe/veto split (CD226 GO / RASGRP1 WITHHELD / PRKCQ WITHHELD) now runs on REAL trust
numbers, not placeholders.*

## Research Question & Hypothesis
- **Q:** When the commit_gate is fed a critic verdict built from LIVE evidence
  (genetics receipt + Day-3 conformal trust + Day-2/4 leakage audit), does it reproduce
  the v3 believe/veto split deterministically -- and does the RASGRP1 veto fire for
  data-grounded reasons, not hand-tuned thresholds?
- **H:** CD226 -> GO (genetics-secure, cell-type eQTL direction-consistent, trust>=floor,
  leakage-clean). RASGRP1 -> WITHHELD (no cell-type eQTL + low cross-donor reproducibility).
  PRKCQ -> WITHHELD (genetic_association 0.162 below the 0.20 floor). All three
  deterministic (same inputs -> same content hash).

## Data Sources (all already produced, no new compute)
- data/gold/genetics_gate_receipt.json (Open Targets + GWAS, live).
- data/gold/conformal_lodo_receipt.json (Day-3 trust machinery).
- data/gold/data_facts_receipt.json (RASGRP1 cross-donor r=0.072, n_downstream=992).
- data/gold/frozen_splits.json (leakage audit inputs).

## Critic -> Verdict construction (the only new logic)
For each target, the critic assembles a CriticVerdict from EXISTING receipts (no LLM call
inside the gate; the gate stays deterministic per v3 Fix 4):
  - genetic_association  <- genetics_gate_receipt (live Open Targets)
  - genome_wide_sig_snp  <- GWAS Catalog lead-SNP evidence
  - celltype_matched_eqtl / eqtl_direction_consistent / proxy_tissue_only <- the v3 ruling
    (CD226 GoF cell-type-consistent; RASGRP1 proxy-only; PRKCQ unanchorable)
  - trust_score          <- a REAL number derived from the Day-3 conformal machinery
    (not a placeholder). Derivation stated explicitly and computed, not hand-typed.
  - leakage_audit_passed <- assert gene not in program axis AND splits frozen (hash match)

## Trust-score derivation (replaces the Day-1 placeholder) — CRITIQUE-CORRECTED
The Day-3 conformal layer yields, per perturbation, an ensemble spread. For the trio,
compute per-gene trust = 1 - normalized(mean ensemble spread over that gene's test rows),
clipped to [0,1]. Provenance differs by gene and is RECORDED per verdict:
  - RASGRP1: present in all 6 donor pairs -> donor-blocked spread (guaranteed machinery).
  - CD226 / PRKCQ: absent/sparse in donor-pair data -> POOLED proxy spread (NOT the
    guaranteed estimator; a different scale).

**BLOCKING FIX #1 + #2 (trust is a co-equal gate + incommensurable estimators):**
Because the 0.50 trust floor is applied uniformly but the estimators differ, and because
`commit_gate.evaluate()` makes trust a NECESSARY (AND) gate (trust<0.50 blocks regardless
of genetics), we do NOT claim "trust is secondary." Instead:
  (a) The runner computes and REPORTS each gene's real trust value up front.
  (b) It ASSERTS that for the demo trio, trust is NEVER the binding constraint: every gene
      either clears the 0.50 floor by a margin (so genetics/eQTL carry the GO/WITHHELD) OR
      is already blocked by genetics/eQTL. If this assertion fails, the DEMO CLAIM changes
      -- the threshold does NOT.
  (c) Each emitted decision lists its BINDING constraint(s), so a trust-driven flip can
      never masquerade as a genetics story.
  (d) Trust provenance (donor-blocked vs pooled-proxy) is stamped on every verdict; the
      uniform floor is applied only after confirming trust is non-binding for the trio.

## Analysis Pipeline
1. Load all receipts; build CriticVerdict per target from live evidence.
2. Run commit_gate.evaluate() + commit() -> GO or BLOCKED artifact per gene.
3. Assert the split matches v3 (CD226 GO, RASGRP1/PRKCQ WITHHELD) and is deterministic.
4. Emit an end-to-end lineage receipt: every gate decision -> the exact receipt fields
   and hashes it depended on (auditable provenance chain).

## Controls & Validation
- **Determinism:** run the gate twice; content hashes must match (already tested Day 1).
- **Counterfactual (shown, not hand-waved):** flip RASGRP1's celltype_matched_eqtl to True
  -> gate flips to GO. Demonstrates the veto is driven by the eQTL fact, not a hard-coded
  gene name. (A test, not a claim.)
- **Threshold provenance:** the 0.20 genetic floor and 0.50 trust floor are pre-registered
  in commit_gate.py; changing them is a code+test change (already the case).

## Statistical Plan
- No new statistics; Day 5 is integration + determinism. The trust numbers carry the
  Day-3 CIs already reported. The believe/veto split is a deterministic function, not an
  estimate.

## Compute Requirements
- CPU only; reads existing receipts. No GPU, no spend.

## Limitations & Assumptions — CRITIQUE-CORRECTED (honest framing of what the gate proves)
- **The commit_gate is a deterministic PROPAGATOR, not a discovery engine.** It proves that
  a fixed body of evidence yields a reproducible, auditable, fail-closed nomination with a
  content hash -- NOT that CD226 is a true target or RASGRP1 a true veto. This is an
  engineering/governance guarantee (determinism, default-deny, sole-writer, threshold
  provenance in code), not a biological validation.
- **Two of the three decisive fields are human-curated inputs**, so the believe/veto split
  is a WORKED EXAMPLE of the governance layer, not independent evidence the ruling is right:
    * eQTL cell-type-match / direction booleans encode the v3 human ruling (cited: rs763361
      GoF cell-type-consistent for CD226; RASGRP1 GTEx LCL-proxy-only; PRKCQ unanchorable).
      A fully-automated eQTL-direction call is future work.
    * trust estimator is provenance-mixed (donor-blocked for RASGRP1, pooled proxy for
      CD226/PRKCQ) and is confirmed non-binding for the trio before the floor is applied.
  Only the PRKCQ genetic-floor veto (GA=0.162 < 0.20) is driven by LIVE data independent of
  the desired outcome.
- The counterfactual test (flip RASGRP1 eQTL -> GO) proves the gate keys on the eQTL FACT,
  not the gene string. It does NOT prove the boolean was assigned correctly -- that rests on
  the cited human ruling.

---

# NEXT-RUN METHODOLOGY (2026-07-09) — Decision-System Validation (Items 1–3)

*Appended for the "From Calibrated Model to Decision System" run. Governing memo:
`Impact_and_GPU_Strategy.md`. All three items run LOCALLY on CPU ($0); no cluster/GPU.
Every output → `data/gold/` receipt + `_script_manifest.jsonl` sha256, matching the
Day 0–5 provenance discipline. Frozen seed 20260708. Frozen split hash 45ca2893cbe7e282.*

*REVISED after pre-COMPUTE critique (VERDICT: BLOCKING, 5 issues B1–B5 + 3 obs O1–O3).
The revisions below implement every fix. Re-critique required before COMPUTE.*

## Research question & hypothesis (REVISED — B1)
- **Q (item 1, reframed):** Does the model+trust layer add validated precision *BEYOND
  genetics alone*? (The original "gate beats naive" framing was circular: the truth channel
  shares its generative process with the gate's genetic-association floor — B1.)
- **H1 (lift-beyond-genetics):** *within the genetics-eligible stratum* (genes already
  passing the genetic-association floor), ranking by model+trust adds validated precision
  over ranking by genetics/effect-size alone, beyond a paired permutation null.
- **H0:** conditional on genetics, model+trust adds no precision (paired lift CI contains 0).
- **H2 (ceiling):** model rank performance is a large fraction of the between-donor
  reproducibility ceiling → "honest model vs noisy assay" (descriptive, donor-limited CI).
- **Success criteria (pre-registered):** H1 accepted ONLY if the PAIRED lift (model+trust −
  genetics-alone), at matched coverage, has a permutation-null CI excluding 0 at ≥1 coverage
  level, AND the effect is reported *stratified* by genetic-floor pass/fail (B1). A positive
  that lives only in the floor-pass stratum is reported as such, not as a general claim.

## Data-sizing pre-check (MANDATORY, run & receipt BEFORE analysis — B2)
Before any lift computation, compute and write to a receipt:
- count of TEST-fold genes with a truth label (split: gold-forced vs random-landed);
- count in the genetics-eligible stratum.
**Gate:** if TEST-fold truth positives < 15, item 1 is REPORTED AS A CASE STUDY with
explicit n honesty (it is NOT an "all-genes lift"); the memo's "all-genes" wording is
reconciled to "test-fold decision replay" (B2). No silent n=8 re-run.

## Data sources (on disk, no download)
- `data/raw/DE_stats.suppl_table.csv` — 33,983 rows; per-row `crossdonor_correlation_mean`,
  `ontarget_effect_size`, `n_downstream`, `target_baseMean`, `n_cells_target`.
- `data/raw/GWCD4i.DE_stats.by_donors.h5mu` — 6 donor-pairs × 2591 reproducible-core genes.
- `data/gold/frozen_splits.json` (hash 45ca2893cbe7e282) — leakage-safe folds.

## Item 1 — Lift BEYOND genetics (REVISED B1/B2/B3)
- **Scope:** TEST-fold genes only (leakage control). Report the data-sizing pre-check first.
- **Comparison (paired, B3):** two rankings over the SAME gene set — (A) model+trust
  (predicted effect × conformal trust), (B) genetics/effect-size alone. Truth = held-out
  causal-gene labels.
- **Truth channel + caveat (B1):** ClinVar pathogenic + GWAS hits, WITH the residual-leakage
  caveat stated explicitly (shared generative process with OT genetic_association). Primary
  reporting **partials out the genetics floor** (restrict to floor-pass stratum) so any lift
  is attributable to model+trust, not genetics. Also report a genetics-INDEPENDENT sensitivity
  using ClinVar entries NOT genome-wide-significant, if any exist at usable N.
- **Statistic (B3):** the PAIRED lift = precision_A@K − precision_B@K at matched coverage K.
  Null = permute SCORES and RE-SELECT top-K inside each of ≥2000 permutations
  (selection-aware); report the paired difference with its own permutation CI. Two separate
  CIs are NOT used to claim a difference.
- **Curve:** paired lift vs coverage K/N, with the permutation null band.

## Item 2 — Reliability ceiling (BETWEEN-DONOR PROXY; REVISED B4/O1/O2/O3)
- **Data correction:** no cell-level counts → no true split-half. Use
  `crossdonor_correlation_mean` + by_donors 2591-core as a between-donor proxy.
- **Framing (O1 — pick ONE):** report it as a **lower bound on the technical-reliability
  ceiling** IF the target is a single donor's delta; OR the *appropriate* ceiling if the
  prediction target is the donor-AVERAGED delta. State that between-donor r includes
  biological donor variation, so it is NOT interchangeable with within-perturbation
  technical split-half.
- **Axes (O2):** state the attenuation assumption (reliability² bounds R² only under
  additive-noise attenuation) and define the rank-ceiling derivation explicitly; do not
  equate Pearson cross-donor r with a Spearman ceiling without stating the mapping.
- **Min-support (O3):** PRE-REGISTER a FIXED threshold (n_cells_target ≥ median of the
  non-null-crossdonor subset, value frozen in the receipt) — not "chosen from the data" —
  and report a sensitivity curve across thresholds.
- **CI (B4):** donor-pair BLOCK bootstrap (resample the 6 pairs), NOT a perturbation
  bootstrap. State the CI is donor-limited (4 donors) and untightenable by perturbation
  count. The perturbation-bootstrap CI is never reported as primary.

## Item 3 — Rule ablation (REVISED B5)
- Ablate each `commit_gate` rule (leakage, genetic-assoc floor, trust floor, eQTL), singly
  AND pairwise.
- **Report FLIP COUNT FIRST.** A rule that flips zero test-set verdicts is labeled
  **"inactive — undetermined"**, NOT "not pivotal" (B5). (Expected: trust floor is
  non-binding for the trio → likely inactive on the small decision set.)
- **Pivotal** = a rule whose Δ(item-1 paired lift) has a CI (from the B3 paired null)
  **excluding 0**. Pre-register that the genetic floor will appear pivotal *trivially*
  because item-1 truth is genetics-confounded (B1) — so its "pivotalness" is discounted.

## Controls & validation
- **Positive control:** model+trust concentrates truth at low coverage within the stratum.
- **Negative control:** permute gene↔score → paired lift collapses into null band.
- **Leakage control:** scoring on TEST-fold genes only; caveat on truth/genetics overlap.

## Statistical plan
- Primary: PAIRED permutation null (≥2000, score-permute + re-select) on the lift; 95% CI.
- Curve over coverage; if a headline K is chosen, Bonferroni across reported K.
- Effect measure: paired precision lift (model+trust − genetics-alone) at matched coverage.

## Compute
- Platform: LOCAL CPU. Runtime: minutes. Cost: $0. (No cluster/GPU for items 1–3; scVI/G1
  is downstream, self-submitted by the user on their own account.)

## Limitations & assumptions (REVISED)
- **Item-1 truth is genetics-confounded** (B1): the DEFENSIBLE claim is "model+trust adds
  lift *conditional on* genetics," not "gate beats genetics." A fully genetics-independent
  truth (functional CRISPR screen) is future work / a separate data pull.
- **Item-1 power (B2):** driven by TEST-fold truth-positive count; if small, it is a case
  study, reported with n honesty — NOT a general lift claim.
- **Item-2 ceiling** is a between-donor proxy conflating biological + technical variance,
  with a donor-limited CI (4 donors); true split-half needs the 1.7 TB raw cells (deferred).
- n=4 donors → all donor-resolved CIs wide; directional only.
