# Wet-Lab Validation Protocol — Prospective Test of the Commit Gate
### CRISPRi knockdown of endogenous immune genes in primary human CD4+ T cells
**STATUS: PROPOSAL. This is a pre-registered experimental DESIGN a collaborator could execute. It contains NO results. All timeline and budget figures are ESTIMATES, explicitly labeled.**

---

## 0. Why this experiment exists
The trust layer makes falsifiable, per-gene commitments: a **GO** call asserts the model's predicted perturbation effect is trustworthy; a **WITHHOLD/ABSTAIN** call asserts it is not sufficiently supported. Every merit-tier gain so far has come from validation against *external* truth (the project's own rule: "the only thing that moves the tier is validation against truth you didn't build"). The honest computational ceiling is ~76–78; the only path into the 80s–90s is a **prospective wet-lab test** whose labels did not exist when the gate was built. This protocol specifies that test.

## 1. Objective & pre-registered hypothesis
**Primary hypothesis (H1, pre-registered, falsifiable):** Genes the gate calls **GO** exhibit larger and more reproducible on-target downstream transcriptional effects upon CRISPRi knockdown than genes it calls **WITHHOLD/ABSTAIN**.

**Primary endpoint:** separation (AUROC) between GO and WITHHOLD genes on a pre-specified downstream-effect-magnitude statistic (§4), measured blind to the gate label.
**Success criterion (a priori):** AUROC ≥ 0.70 with a 95% CI lower bound > 0.50.
**Failure / null criterion (a priori):** CI includes 0.50 → the gate does **not** prospectively separate — a **publishable negative** that bounds the method honestly (§8).
**Null is valuable:** this project reports nulls on the front page; a well-powered negative here is a real result, not a failure of the study.

## 2. Gene panel (pre-registered, powered)
**25 genes — 12 GO / 13 WITHHOLD**, balanced by design, drawn from the committed 98-gene eQTL/OpenTargets overlap set **plus the PRKCQ anchor** — a committed WITHHOLD on the genetic floor (GA=0.162<0.20) that sits *outside* the 98-gene atlas-overlap set and is added explicitly so the anchor trio is complete. Anchor trio: **CD226 (GO), RASGRP1 (WITHHOLD; role: ABSTAIN — no cis-eQTL in 120 T-cell datasets), PRKCQ (WITHHOLD; genetic floor)**.

**Selection rule (frozen):** the 3 anchors + the 11 highest-T1D-association GO genes + the 11 highest-effect WITHHOLD genes, so the contrast is not confounded with raw effect size. Gate labels are frozen from `routeA_zhu_frozen_scoring_table.json` (v704dd0ee) before any wet-lab work; PRKCQ's WITHHOLD call comes from the committed genetic-floor gate, not the atlas table.

| Gene | Gate call | T1D assoc | eQTL sig (datasets) | Gate's prediction |
|------|-----------|-----------|---------------------|-------------------|
| CD226 | **GO** | 0.533 | 4 | large/reproducible on-target effect |
| RASGRP1 | **WITHHOLD** | 0.32 | 0 | no committed prediction (withhold/abstain) |
| PRKCQ *(anchor; genetic floor, outside atlas-overlap set)* | **WITHHOLD** | 0.162 | 3 | no committed prediction (withhold/abstain) |
| INSR | **GO** | 0.629 | 2 | large/reproducible on-target effect |
| CD3E | **GO** | 0.564 | 1 | large/reproducible on-target effect |
| SIRPG | **GO** | 0.532 | 2 | large/reproducible on-target effect |
| TCF7L2 | **GO** | 0.527 | 2 | large/reproducible on-target effect |
| SUOX | **GO** | 0.516 | 6 | large/reproducible on-target effect |
| UBASH3A | **GO** | 0.484 | 6 | large/reproducible on-target effect |
| APOBR | **GO** | 0.472 | 3 | large/reproducible on-target effect |
| GLIS3 | **GO** | 0.463 | 3 | large/reproducible on-target effect |
| ITPR3 | **GO** | 0.446 | 1 | large/reproducible on-target effect |
| ZFP36L1 | **GO** | 0.42 | 2 | large/reproducible on-target effect |
| CDKAL1 | **GO** | 0.415 | 1 | large/reproducible on-target effect |
| NDUFAF3 | **WITHHOLD** | 0.397 | 0 | no committed prediction (withhold/abstain) |
| GATA3 | **WITHHOLD** | 0.325 | 0 | no committed prediction (withhold/abstain) |
| IL2RB | **WITHHOLD** | 0.378 | 0 | no committed prediction (withhold/abstain) |
| NDUFA8 | **WITHHOLD** | 0.394 | 0 | no committed prediction (withhold/abstain) |
| NDUFS8 | **WITHHOLD** | 0.399 | 0 | no committed prediction (withhold/abstain) |
| IL7R | **WITHHOLD** | 0.458 | 0 | no committed prediction (withhold/abstain) |
| NDUFA9 | **WITHHOLD** | 0.394 | 0 | no committed prediction (withhold/abstain) |
| PTPN2 | **WITHHOLD** | 0.479 | 0 | no committed prediction (withhold/abstain) |
| KEAP1 | **WITHHOLD** | 0.326 | 0 | no committed prediction (withhold/abstain) |
| NDUFC1 | **WITHHOLD** | 0.394 | 0 | no committed prediction (withhold/abstain) |
| NDUFB9 | **WITHHOLD** | 0.394 | 0 | no committed prediction (withhold/abstain) |

*(RASGRP1 is the teaching case: a large model-predicted effect, but its gate role is ABSTAIN — no cis-eQTL exists in 120 T-cell datasets — so the gate's prediction is "do not trust the effect," which this experiment tests directly. In the binary GO/WITHHOLD commit gate it is recorded as WITHHOLD.)*

**Power:** with 12 GO vs 13 WITHHOLD genes, a two-group AUROC test detects a true AUROC of 0.75 vs the 0.50 null at ~80% power, α=0.05 (two-sided), under moderate within-group variance (Hanley–McNeil approximation). If pilot variance is larger, the panel scales to 40 genes without redesign; the selection rule is pre-specified to extend by rank.


## 3. Experimental design
- **Perturbation:** CRISPRi via dCas9-KRAB, lentiviral delivery into primary human CD4+ T cells (method: Schmidt et al. 2022, *Science* 375:abj4008; Perturb-seq framework: Dixit 2016; Replogle 2022).
- **Guides:** 2–3 independent sgRNAs per gene targeting the TSS; **non-targeting** sgRNA controls and a **safe-harbor** (AAVS1) control; guide-level replication lets gene effects be separated from guide artifacts.
- **Donors:** ≥ 4 independent healthy donors (matches the project's 4-donor structure and its known **donor-generalization** concern; leave-one-donor-out is the analysis unit so donor is a modeled random effect, not a confound).
- **Conditions:** resting + TCR-stimulated (mirrors the atlas Rest/Stim8/Stim48 contexts — and lets the protocol test the project's own cross-context recalibration finding in the wet lab).
- **Randomization/blinding:** gate labels (GO/WITHHOLD) are **masked** during library prep, sequencing, and primary effect-size scoring; unblinded only at the pre-registered analysis step.

## 4. Readout
- **Assay:** single-cell RNA-seq (10x Genomics 3′) for genome-wide downstream effects; a targeted panel is the lower-cost fallback (§7).
- **QC gates:** ≥ 500 cells/guide passing standard scRNA QC (mito%, UMI, doublet removal); CRISPRi knockdown efficiency ≥ 50% at the mRNA level (§6).
- **Effect statistic (maps to the gate's prediction):** per gene, the magnitude of the downstream differential-expression response (number of significant downstream DE genes and/or aggregate signed program shift), computed with the **on-target gene zeroed** (the trivial-artifact guard carried directly from the computational sign-trust work). This is the quantity the gate commits on.

## 5. Analysis plan (pre-registered)
1. **Primary:** GO vs WITHHOLD AUROC on the §4 statistic + 95% CI (cluster bootstrap over donors). Evaluated once, blind-broken.
2. **Recalibration built in:** apply the project's per-assay recalibration (the frozen→recalibrated step that took the computational AUROC 0.47→0.71) to the wet-lab measurements; report both frozen and recalibrated separation, since the whole thesis is that calibration must be re-fit on the deployment assay.
3. **Coverage:** do the gate's conformal trust intervals cover the observed wet-lab effects at the nominal 80/90/95%? This is the direct prospective test of the calibration claim.
4. **Multiplicity:** Benjamini–Hochberg across genes; primary endpoint is the single AUROC, so no correction needed for it.
5. **Stopping rule:** fixed sample size (no interim peeking); pre-registered on OSF before donor 1.

## 6. Controls & rigor
- On-target knockdown confirmed per gene (qPCR or the CRISPRi guide's own on-target UMI drop) before a gene enters the effect analysis.
- Non-targeting + safe-harbor controls define the null distribution of "effect magnitude."
- Batch/donor confounds handled by donor as a random effect + batch covariates; guides nested within gene.
- Trivial-artifact guards inherited from the computational work: on-target column excluded from the downstream statistic; gene identity never used as a feature.

## 7. Timeline & budget — ESTIMATES ONLY
*(order-of-magnitude planning figures, not quotes; actuals depend on institution, core rates, and donor sourcing)*
- **Duration:** ~4–6 months (guide cloning + lenti prep ~4–6 wk; primary T-cell CRISPRi + expansion ~3–4 wk/donor batch, batched; scRNA-seq + sequencing ~3–4 wk; analysis ~4–6 wk).
- **Sequencing:** ~24 genes × ~3 guides × 4 donors × 2 conditions ≈ pooled Perturb-seq design; ~2–4 × 10x lanes → **~$25k–60k** sequencing (est.).
- **Reagents/lenti/culture:** **~$20k–40k** (est.).
- **Personnel:** ~0.5–1 FTE research associate for the run (est.).
- **Rough all-in:** **~$75k–150k** and one to two quarters, assuming an existing primary-T-cell CRISPRi-capable lab. *All figures are estimates and must be confirmed with the executing core.*

## 8. Expected outcomes & interpretation
- **Positive (AUROC ≥ 0.70, CI > 0.50):** the gate prospectively separates trustworthy from untrustworthy predictions against labels that did not exist when it was built — the tier-moving result (toward the 80s–90s on the internal scale) and the core claim of a publishable methods paper.
- **Negative (CI includes 0.50):** the gate does not prospectively separate on this readout. This **bounds the method honestly**, is directly publishable (a calibration/selective-prediction benchmark with a negative headline is exactly the kind of result the field now rewards), and tells the next iteration where the ceiling is.
- **Partial (GO effects real but recalibration required):** confirms the central finding — calibration must be re-fit per assay — now in the wet lab.

## 9. Feasibility & risks
- **Primary-T-cell CRISPRi efficiency** is the main technical risk; mitigated by per-gene knockdown QC and 2–3 guides/gene (established protocol, Schmidt 2022).
- **Donor variability** is expected and is the analysis unit (LODO), not a nuisance to hide.
- **Cost/scale**: the panel is deliberately sized to be powered yet affordable; the targeted-panel fallback roughly halves sequencing cost if budget-constrained.
- **What would make it fail:** knockdown too weak to produce measurable downstream effects (→ readout floor), or effect sizes so donor-variable that no method separates (→ the honest negative, still informative).

---
*Grounded in the committed gate calls (`routeA_zhu_frozen_scoring_table.json` v704dd0ee; anchor eQTL/genetic evidence from the committed eqtl_gate + OpenTargets receipts). Methods references: Schmidt et al. 2022 Science 375:abj4008 (CRISPRi in primary human T cells); Dixit et al. 2016 Cell (Perturb-seq); Replogle et al. 2022 Cell (genome-scale Perturb-seq). This document is a proposal; no experiments have been performed.*
