# SCOPE: Chaining Beats Leaping — A Layered Test of the Genotype→Phenotype Map

## The one-sentence pitch
Most genotype→phenotype prediction "leaps" from DNA straight to a trait. We test whether
**chaining through intermediate molecular layers** — sequence → regulatory activity →
expression program → protein-network context — predicts a cellular phenotype better than
leaping, and whether the **omnigenic** core/peripheral split explains *which* genes chain
cleanly. Three Gladstone datasets each supply one layer, so the project is natively
integrative.

## Why this is a good hackathon bet
- **Uses the provided data as designed but connects it** — most teams will attack one
  dataset in isolation. Connecting layers is the differentiated, judge-friendly story.
- **Tests a *named theory*** (Pritchard's omnigenic model — and Pritchard is literally a
  co-author on the provided Perturb-seq dataset). That gives the empirical finding a spine.
- **Produces a discrete, reproducible finding + figures** — exactly the Research-track ask.
- **Fits your profile** — real ML (fine-tune a genomic sequence model), real GPU, real
  single-cell analysis, but scoped so a v1 result exists by mid-week.

## Research question (precise)
For a set of genes/regulatory elements measured in T cells:
1. **Q1 (chaining vs leaping):** Does a *layered* predictor — sequence → predicted
   regulatory activity → predicted expression change — explain more variance in the
   Perturb-seq transcriptomic phenotype than a *direct* sequence→phenotype model of equal
   capacity?
2. **Q2 (omnigenic structure):** Do "core" genes (large, direct trans-effects in
   Perturb-seq; hub position in the Krogan network) chain more cleanly than "peripheral"
   genes, as the omnigenic model predicts?

## Primary hypothesis (H1)
Layer-chained prediction explains significantly more phenotype variance than an
equal-capacity direct model, **and** the chaining advantage is concentrated in peripheral
genes (whose effects are indirect/network-mediated), while core genes are predictable even
by the direct model.

## Null hypothesis (H0)
The layered predictor explains no more variance than the direct model (Δ variance-explained
≤ noise floor from bootstrapping), and there is no core/peripheral difference in chaining
advantage.

## Pre-specified success criteria (define BEFORE any analysis)
- **Minimum success:** A trained DNA→regulatory-activity model with held-out performance
  reported honestly (Spearman/Pearson vs measured MPRA), + a Perturb-seq phenotype
  representation, + at least one clean figure comparing chained vs direct prediction with
  bootstrap CIs. Even a *null* result (chaining doesn't help) is a legitimate, publishable
  finding for this hackathon.
- **Target success:** Statistically significant Δ variance-explained (chaining > leaping,
  95% bootstrap CI excludes 0) on held-out genes.
- **Stretch success:** Q2 confirmed — chaining advantage significantly larger for
  peripheral than core genes; identification of specific network-mediated predictions.

## What each dataset contributes (layer map)
| Layer | Dataset | Mapping learned/measured |
|-------|---------|--------------------------|
| Sequence → regulatory | Pollard MPRA | Fine-tune seq model to predict element activity; in-silico mutagenesis for single-base effects |
| Perturbation → expression | Marson/Pritchard T-cell Perturb-seq | Gene KO → transcriptomic phenotype vector; define core vs peripheral empirically |
| Protein context | Krogan PPI network | Hub/periphery structure; network-mediated indirect effects; missing-member prediction |

## Scope boundaries (what we will NOT do)
- Not training a foundation model from scratch — fine-tune / probe pretrained genomic LMs.
- Not touching human whole-organism phenotype — "phenotype" = cellular transcriptomic state.
- Not claiming causal mechanism beyond what perturbation + held-out prediction support.

## Risk register & de-risking
- **Risk: sequence→MPRA model too weak to be a useful layer.** → Fallback: use measured
  MPRA activity (skip the learned layer) so the chaining test still runs on real values.
- **Risk: dataset access/licensing at hackathon.** → Pre-identify public equivalents
  (Replogle genome-wide Perturb-seq; ENCODE/lentiMPRA; BioGRID/STRING for PPI) so nothing
  blocks on a single gated download.
- **Risk: integration/ID-mapping friction (gene IDs across 3 datasets).** → Budget explicit
  time; use a single canonical gene ID space (Ensembl) early.

## Week plan (7 days)
- **Day 1 (kickoff):** Confirm dataset access. Finish literature review + reasoning. Lock
  methodology + pre-COMPUTE critique. Decide exact models/datasets.
- **Day 2:** Data ingestion + ID harmonization across all three layers. Build the
  Perturb-seq phenotype matrix; define core/peripheral labels.
- **Day 3:** Train/fine-tune sequence→regulatory model on MPRA; report held-out perf.
  In-silico saturation mutagenesis sanity check.
- **Day 4:** Build chained predictor + equal-capacity direct baseline. First variance-
  explained comparison with bootstrap CIs.
- **Day 5:** Q2 stratification (core vs peripheral); network-mediated effect analysis with
  Krogan PPI. Post-COMPUTE critique.
- **Day 6:** Figures (publication quality), robustness checks, negative-result honesty pass.
- **Day 7:** Writeup + reproducible repo + submission. Optional demo.

## Deliverables for submission
1. `finding`: chained-vs-direct variance-explained result (+ core/peripheral) with CIs.
2. `figures/`: 3–5 publication-quality figures.
3. Reproducible repo (data-fetch scripts, training, analysis, seeds).
4. Short writeup framing the result against the omnigenic model.
