# Angle 2 — Public datasets to broaden calibration & donor diversity (not more features)

**Frame:** the trust layer is a conformal/CQR gate on top of a HistGBR effect-size predictor, run under leave-one-donor-out (LODO). This session already proved an information limit: features predict the base predictor's |residual| at R²≈0.01 (vs 0.975 for the true effect). **Adding feature columns is a proven dead end.** New data helps *only* if it is (A) more of the SAME task (primary CD4⁺ T-cell perturbation with an effect-size DE readout) for calibration/coverage/recalibration, (B) donor-level expression variation to inform covariate-shift weighting under LODO, or (C) an independent T1D-relevance prior used as a gate signal (never as a residual-predicting feature). Every entry is tagged GOOD / DEAD-END accordingly. The external null (Zhu, AUROC 0.465) is ground truth: frozen transfer already failed → any gain must be shown after **recalibration** and validated on held-out donors, not asserted.

---

## TIER 1 — Direct calibration / recalibration assets (same task: primary CD4⁺ T-cell perturbation → effect-size DE)

These are the only sources that add *calibration points and donors for the exact task*. This is the "right kind of data."

| Rank | Dataset | What it is | Verdict | Laptop? |
|---|---|---|---|---|
| 1 | **Zhu et al. 2025** (bioRxiv 2025.12.23.696273) — *already on disk* | Genome-scale CRISPRi Perturb-seq, ~22M primary human CD4⁺ T cells, ~12,748 genes, **rest + stim** | **GOOD — recalibration.** Repurpose from blind-validation-only to a same-assay recalibration/second-context calibration set. It is the same cell type and task as Marson; the NULL means frozen calibration doesn't transfer → recalibrate CQR on Zhu donors and its rest/stim split, then re-measure LODO coverage on Zhu. Highest expected value because it is same-task AND local. | Yes for pseudobulk/DE stats (16.79GB h5ad on disk, donor-blocked); full 22M-cell matrix is cluster-only |
| 2 | **Schmidt et al. 2022**, *Science* 375, abj4008 (GEO; "CRISPRa Perturb-seq") | Genome-wide CRISPRa/i screens + CRISPRa Perturb-seq scRNA-seq of hits in **primary human CD4⁺/CD8⁺ T cells** (IL-2 / IFN-γ programs) | **GOOD (with recalibration) — perturbation-modality diversity + donors.** Adds gain-of-function perturbations (opposite direction to Marson KO/CRISPRi). Because modality differs, treat as a concept-shift candidate like Norman → recalibrate, don't drop in frozen. Tests whether the gate generalizes across CRISPRa vs CRISPRi. | Yes — the scRNA-seq arm is modest |
| 3 | **scPerturb** (Peidli et al. 2024, *Nat Methods* 21:531) — scperturb.org | 44 harmonized single-cell perturbation datasets, uniform QC, **E-distance effect-magnitude metric** | **PARTIAL GOOD — robustness/recalibration corpus.** Ready pool to stress-test covariate-shift-weighted conformal across contexts and to harvest the few primary-T-cell sets (Shifrut, Datlinger CROP-seq, Papalexi ECCITE-seq). Most members are K562/THP-1/RPE1 → cross-assay concept shift (same failure mode as Norman); useful as recalibration targets, **not** drop-in calibration or donor-diversity. E-distance gives a principled magnitude label. | Yes per-dataset (harmonized h5ad) |

---

## TIER 2 — Donor-diversity references for covariate-shift weighting (CD4 expression variation across many donors; NOT perturbation)

Not calibration points, but they parameterize donor shift for LODO — a genuine calibration lever, not volume.

| Rank | Dataset | What it is | Verdict | Laptop? |
|---|---|---|---|---|
| 4 | **Soskic et al. 2022**, *Nat Genet* 54:817 | 655,349 CD4⁺ T cells, **119 donors**, resting + 3 activation timepoints; 6,407 eGenes, 2,265 dynamic; 127 immune-disease-linked | **GOOD — donor-shift reference + genetic gate input.** Donor variation matched to CD4 activation states (directly comparable to Marson rest/stim) → improves covariate-shift weights and conditional coverage across donors. Dynamic eQTLs also feed Tier 3. | Yes — eQTL/summary tables; full matrix cluster-only |
| 5 | **OneK1K** (Yazar et al. 2022, *Science* 376:eabf3041) — onek1k.org | 1.27M PBMCs, **982 donors**, sc-eQTL in 14 immune cell types incl. CD4; 26,597 cis-eQTLs; 305 autoimmune loci resolved | **GOOD — donor-shift reference (population scale) + genetic gate input.** Largest donor panel for CD4 expression variation → strongest empirical prior for LODO donor-shift weighting. Its autoimmune-colocalizing eQTLs double as Tier-3 priors. | Yes — eQTL summary stats; 1.27M-cell matrix cluster-only |

---

## TIER 3 — Genetic gate inputs (independent T1D-relevance priors; use as a gate signal, NEVER a residual feature)

**Strict caveat:** feeding any of these into HistGBR as another column to predict effect size / |residual| = DEAD END (info-limit). The only defensible use is an *independent* prior that reweights or vetoes commits, validated by whether it improves gate **precision/AUROC on the Zhu null**, not by residual R².

| Rank | Dataset | What it is | Verdict | Laptop? |
|---|---|---|---|---|
| 6 | **Chiou et al. 2021**, *Nature* 594:398 | T1D GWAS (520,580 samples) + snATAC (131,554 nuclei); **fine-mapped credible sets**, T-cell-enriched cCREs, cCRE→gene links | **GOOD — the T1D-specific gate prior.** Which perturbed genes carry fine-mapped T1D causal support / T-cell regulatory evidence. Disease-matched, orthogonal to the base predictor. | Yes — credible-set/supp tables are tiny |
| 7 | **Open Targets Platform / Genetics** (Ochoa et al.; NAR 2025, rel 24.12) | Integrated, scored target–disease associations; genetics + colocalization; API + bulk download | **GOOD — easiest genetic prior to wire in today.** Per-gene T1D association score as a plausibility prior on perturbed genes. Lowest integration effort of all Tier 3. | Yes — API/download |
| 8 | **eQTL Catalogue** (Kerimov et al. 2021, *Nat Genet* 53:1290; 2023 update) — EBI | 31 uniformly reprocessed studies incl. **BLUEPRINT CD4⁺ T cells**; 1.7M fine-mapped associations; REST API | **GOOD — harmonized, CD4-specific colocalization inputs.** Best-engineered source; immune-cell eQTLs diverge strongly from whole blood, so use purified-CD4 QTLs (not GTEx whole blood). | Yes — API/summary stats |
| 9 | **DICE** (Schmiedel et al. 2018, *Cell* 175:1701; 2022 *Sci Immunol* sc-eQTL) — dice-database.org | Immune eQTL across 13–15 cell types incl. CD4 subsets (~91 donors); 2022 sc-eQTL in activated CD4 (89 donors, 19 subsets) | **GOOD — activation/cell-type-specific eQTL priors.** 41% of eGenes single-cell-type specific → cell-type-resolved priors. | Yes — eQTL tables |
| 10 | **ImmuNexUT** (Ota et al. 2021, *Cell* 184:3006) | 28 immune cell types, **337 patients + 79 controls**, 9,852 samples; eQTLs for 13,395 genes incl. **patient-only eQTLs** | **GOOD — disease-state genetic prior.** Unique autoimmune-patient eQTLs; broadest immune eQTL panel. | Yes — summary stats |

---

## DEAD ENDS (explicitly flagged)

- **Any Tier-3 resource used as a HistGBR feature column** to predict effect size or |residual| → dead end (R²≈0.01 proven this session). Genetic data is a *gate prior*, not a predictor input.
- **Norman K562 and Replogle 2022 (K562/RPE1) genome-scale Perturb-seq as drop-in calibration** → dead end for donor diversity: wrong cell type, and frozen transfer already failed on Norman (concept shift). Only valid as recalibration/robustness stress tests.
- **GTEx whole-blood eQTL as the immune prior** → weak: purified-immune-cell eQTLs share little with whole blood (shown in eQTL Catalogue). Use CD4-specific QTLs instead.
- **"More cells → better model" from large atlases (OneK1K/Zhu full matrix)** → volume, not calibration, unless the cells are used specifically for donor-shift weighting (Tier 2) or same-task recalibration (Tier 1).

---

## Validation protocol for any claimed gain
1. Never claim improvement from frozen transfer — recalibrate CQR on the new source's donors first (Norman/Zhu both proved frozen fails).
2. Report LODO **marginal + conditional (per-donor) coverage** on held-out donors of the new source, not Marson only.
3. For Tier-3 gate priors: A/B the gate **with vs without** the prior and report precision/AUROC on the **Zhu external null** — the prior earns its place only if it lifts separation there.
