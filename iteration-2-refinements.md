# Iteration 2 — 20-agent deep dive: refinements & validation

> **⚠️ HISTORICAL AUDIT TRAIL — SUPERSEDED BY v2.1.** This file preserves the pre-v2.1 reasoning (CD226 = "direction-secure" anchor; RASGRP1 = hypermorphic novel bet). Both framings were later revised: CD226 is now the **hedged calibration flagship** (its KD-on-suppression direction is *contested* — Ma 2023 *Cell Reports* [PMID 37864795] vs Brusko 2023 *Diabetes*), and **PRKCQ is the primary novel bet** with RASGRP1 gated on a Day-1 eQTL check. For the current state see `revised_plans/PLAN_OF_ATTACK_v2.pdf` (§3a/§3b) and the reconciled `README.md` / `SCOPE.md` / `t1d-nomination-memo.md`. Kept as-is for the audit record.

*Synthesis of 20 parallel research agents across five clusters (trust/calibration methods · competitive+alternative angles · biology/nomination hardening · demo+Claude-use · impact/Gladstone). Verdict, then the specific refinements.*

---

## Verdict: HOLD the core, REFINE substantially (this is iteration, not a pivot)

The core bet — a **calibrated trust layer over CD4 T-cell perturbation predictions, anchored to Type 1 diabetes, driven by a self-critiquing agent council** — survived every stress test and came out validated as *genuinely orthogonal* to the field and *fresh as a headline*. But five refinements materially raise the ceiling, and one biology finding forces a restructure. The single most important upgrade is a **rebrand + a two-target nomination structure**:

> **"A pLDDT for perturbation biology"** — a calibrated trust layer that tells you which T-cell perturbation predictions to believe, proves it by recovering held-out T1D genes, and nominates the next experiment: a genetics-anchored **anchor (CD226)** that shows the pipeline finds real biology, and a genuinely novel **bet (RASGRP1)** it surfaces on its own.

Why this framing wins: AlphaFold didn't just predict structures — its **pLDDT confidence score is what made it usable** (people act only on high-confidence regions). Confidence-as-the-product is a proven, legible pattern in structure (pLDDT/PAE) and variant effect (AlphaMissense's explicit "uncertain" class) but **nobody leads with it in perturbation biology** — the mechanics exist as methods (PRESCRIBE, NeurIPS 2025) and the virtual-cell field openly calls calibrated uncertainty "the missing piece," yet no one has made it the headline identity. You'd be first. ([pLDDT](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/) · [AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492) · [PRESCRIBE](https://arxiv.org/abs/2510.07964))

---

## 1. The method is now rigorous, buildable, and defensible (Cluster A)

The "knows when it's wrong" identity now has a coherent, cheap, expert-proof stack — each layer independently recommended ADOPT:

| Layer | Method | Why | Source |
|---|---|---|---|
| Base predictor | Pseudobulk + **delta/residual** prediction + **ESM-2 protein embeddings** + DEG-frequency/mean-expression priors | The exact Arc-VCC-winning substrate; cheap to reproduce, frees your novelty budget for the trust layer | [Arc VCC wrap-up](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up) |
| Uncertainty | **Deep ensemble (3-5 seeds)** | Best calibration + shift-robustness per unit effort; disagreement doubles as a **mean-collapse detector** | [Ovadia 2019](https://arxiv.org/abs/1906.02530) |
| Guarantee | **Conformal**: split-CQR + **Mondrian** (per-slice) + **weighted conformal** for unseen-gene covariate shift | Distribution-free, finite-sample coverage; weighted CP makes "handles unseen genes" a *provable* property | [Weighted CP](https://arxiv.org/abs/1904.06019) · [CQR](https://arxiv.org/abs/1905.03222) |
| The demo metric | **Selective prediction**: risk-coverage curve + **AURC** vs a random-abstention null | THE judge-legible artifact: "error falls as we keep only high-trust predictions," formalized | [SelectiveNet](https://arxiv.org/abs/1901.09192) |
| OOD trust signal | **kNN-distance in the model's gene-embedding space** | Simplest defensible "distance from training"; PRESCRIBE shows such distance tracks per-gene accuracy | [PRESCRIBE](https://arxiv.org/html/2510.07964v1) |
| Proof of calibration | Reliability diagram + **proper score (Brier/NLL/CRPS)** + coverage-vs-sharpness + **per-slice (seen vs unseen) table** + constant baseline + bootstrap CIs | Pre-empts the two guaranteed attacks: "ECE is gameable" and "aggregate hides subgroup miscalibration" | [Kuleshov 2018](https://proceedings.mlr.press/v80/kuleshov18a/kuleshov18a.pdf) |

**REJECT evidential deep learning** (documented overconfidence — wrong tool when honesty is the point). **Evaluate on a held-out-perturbation (OOD) split, never IID** — otherwise your uncertainty looks great by rewarding confident, mean-collapsed, trivially-wrong predictions. Differentiate explicitly from PRESCRIBE/GEARS: "they filter by confidence; we formalize it as selective prediction with AURC in T1D, and prove per-slice calibration."

---

## 2. The biology restructure — the most important change (Cluster C)

**CD226's "novel target" claim is dead — but it becomes MORE useful as an anchor, not less.** Three findings force this:

- The exact mechanism (block CD226 → raise FOXP3 → treat autoimmunity) is **published, patented (Tsukuba/TNAX), and IND-cleared as Riverview RVW101** (ulcerative-colitis-first, July 2025). ([TNAX](https://tnaxbio.com/2025/07/30/news/))
- A **2023 *Cell Reports* paper directly contradicts** "KD stabilizes Tregs": Treg-specific CD226 deletion *increased* Treg numbers but *weakened* suppression via mTOR/Myc. ([Cell Reports 2023](https://www.cell.com/cell-reports/fulltext/S2211-1247(23)01318-9))
- CD226 is **essential for NK/CD8 tumor and antiviral surveillance** → systemic knockdown is risky for a pediatric indication. ([review](https://www.nature.com/articles/s41423-020-00633-0))

**So restructure the nomination into a two-target story that is actually stronger than one gene:**

| Role | Gene | What it does for the pitch |
|---|---|---|
| **The anchor** (proof the pipeline finds real biology) | **CD226** | Genetics-anchored (coding GoF variant rs763361), direction-secure, and **now independently de-risked by a competitor's IND** — cite that as *validation of the thesis*, and differentiate on modality (Treg-selective KD) + T1D indication. NOT a novel discovery. |
| **The novel bet** (the fresh, falsifiable discovery) | **RASGRP1** | *More novel* than CD226, independent T1D locus, **hypermorphic risk allele** (risk ↑RASGRP1/ERK, so KD mimics protection — same logic as CD226), and a **double mechanism**: peripheral loss ↑Treg expansion/function AND ↓RAS-ERK effector drive. Largely unprosecuted as a Treg CRISPRi target. ([GWAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC3272492/) · [risk↑expr](https://pmc.ncbi.nlm.nih.gov/articles/PMC6536009/)) |

This "recover a known anchor **and** make a novel falsifiable bet" structure is exactly what the demo judge said separates #1 from finalist. **Backup if RASGRP1's therapeutic window worries you: PRKCQ** (distinct T1D locus, druggable, well-tolerated loss). **RASGRP1 caveats to state:** complete loss = immunodeficiency/lymphoma → use *tunable/partial, Treg-restricted* CRISPRi; confirm the T1D-risk-allele eQTL direction in primary CD4/Treg first.

**Two more biology fixes:**
- **Program axis:** score perturbations with a **UCell signed contrast** — beneficial (stable-Treg module + Tr1 module) minus pathogenic (ex-Treg/Th1 *gated on low FOXP3* + Th17 + autoreactive BHLHE40/CSF2/TNF) — validated against the Genome Medicine 2024 **TMZ score**, TSDR methylation, and CITE-seq. Use UCell/AUCell (composition-robust), regress out cell-cycle/composition/donor/sex/HLA. ([axis atlas](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01300-z) · [UCell](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8271111/))
- **Genotype-by-perturbation → hypothesis only.** With n=4 donors, genotype is inseparable from donor identity (precedent needs 89-259 donors). Keep it as a **pre-registered, precedent-backed vignette** using within-donor allele-specific expression + external GTEx data — never a confirmatory p-value. ([n needed](https://www.nature.com/articles/s41588-025-02344-6))

---

## 3. Angle upgrades that raise the ceiling (Cluster B)

Keep the trust headline, but nest three things under it — each independently makes the story bigger without dilution:

- **Active learning = the payoff, not a rival headline.** "Which perturbation next" is largely owned by IterPert/GeneDisco. Nest it: *calibration is the contribution; picking the next experiment is the proof it's decision-grade.* Do **not** headline a retrospective budget curve (a known-weak form); frame it as a sanity check with random + prior baselines. ([IterPert](https://www.biorxiv.org/content/10.1101/2023.12.12.571389v1))
- **Causal regulator discovery = supporting mechanism, not sole headline.** There's a near-exact **CD4+ T-cell precedent** (Cell Genomics 2024, 84-TF causal GRN linked to immune GWAS). Scope it small/local, make it falsifiable by leave-one-perturbation-out interventional prediction, and **gate edges by the trust score** — that turns causal fragility into your differentiator. ([precedent](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00290-8))
- **Cross-disease transfer panel = cheap novelty + impact.** Keep T1D deep, then reuse the identical perturbation→program map and swap public GWAS for **RA / celiac / MS / lupus**. Reframe: *"a T1D-calibrated engine whose trust score predicts WHERE it transfers across CD4-driven autoimmunity."* Don't overclaim "one engine, all diseases" — disease-specific context is the signal the calibration detects. ([shared architecture](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10291128/))

---

## 4. Demo + Claude-use: concrete and buildable (Cluster D)

- **The money shot is triple-validated** — the last *three* Anthropic hackathon winners all built around adversarial self-verification. Lead with it (teaser at 0:00, full reshuffle by ~1:10), put **one blinded number on screen** ("recovered 8/10 known T1D genes it never saw"), cut architecture to <20s, reserve the final 2-3 hrs for the video. ([winners](https://claude.com/blog/meet-the-winners-of-our-claude-opus-4-8-build-day-hackathon))
- **Build the agent council on Claude Code + Agent SDK** (agents as `.claude/agents/*.md`, biology MCP connectors from the `anthropics/life-sciences` marketplace: Open Targets, ChEMBL, PubMed, 10x). Don't rebuild what Claude Science ships — **extend its reviewer into a vetoing, evidence-modality critic bench**. The Elo tournament is a solved pattern (open Swarms code) = an afternoon; spend the budget on the **veto gate + verification**. The critic's teeth = **tool-scoping**: give it `Read/PubMed` only, **no Write**. Show the `skeptic.md` frontmatter on a slide. Run the pipeline with/without the veto and show the **hallucination-rate delta**. ([subagents](https://docs.claude.com/en/docs/agent-sdk/subagents) · [Virtual Lab](https://github.com/zou-group/virtual-lab))
- **Reproducibility that's genuinely true:** scope "Reproduce" to a **deterministic eval pass, not GPU training** (bitwise GPU repro is a coin-flip liability). Content-address + canonically serialize the leaderboard (sorted rows, fixed floats, no timestamps) → SHA-256 → recompute on click → "VERIFIED — HASH MATCH." Build the provenance manifest first, UI second. ([DVC](https://dvc.org/doc/user-guide/project-structure/internal-files) · [PyTorch determinism](https://pytorch.org/docs/stable/notes/randomness.html))

---

## 5. Impact + Gladstone framing (Cluster E)

- **Build on their own data.** The hackathon ships **Marson's (Gladstone) T-cell Perturb-seq**; Marson explicitly names T1D and runs Foxp3/Treg CRISPR screens; Pollard (Gladstone data science) rewards **reusable open tools + hypothesis-prioritization**. Frame it as "we built on your data, and we ship a reusable tool, not a one-off." Stay in the autoimmunity/Treg lane; don't pitch beta-cell replacement or a risk dashboard; don't overclaim validation ("Lead with Integrity" is a stated Gladstone value). ([Marson](https://gladstone.org/people/alex-marson))
- **The impact anchor with hard numbers:** teplizumab delays onset a median of **~2 years — the *entire* ceiling of disease modification** — at ~$194k/course; T1D costs ~$500k lifetime and ~9 life-years per patient; 30-60% of kids present in DKA. The shared failure point of *every* current therapy is **Tregs destabilizing in the inflamed islet** (IL-2-deprivation apoptosis, TSDR hypermethylation, FOXP3 loss → ex-Tregs). Position the target as the missing **stability node** that makes transient teplizumab/IL-2/adoptive-Treg benefit *durable*. Go/no-go bar: beat the 24-month ceiling. ([teplizumab TN-10](https://www.nejm.org/doi/full/10.1056/NEJMoa1902226) · [Treg failure](https://pmc.ncbi.nlm.nih.gov/articles/PMC8275894/))
- **Name a real validation path.** CD226/RASGRP1 **RNP knockdown in primary human CD4/Treg → flow FOXP3/TIGIT + suppression assay** is routine, orderable, ~4-8 weeks, low-$thousands, and reproduces a published phenotype. Real labs: **Marson (UCSF/Gladstone), Levings (UBC)**; commercial: **Synthego/EditCo**. Present Tier 1 (knockdown+flow) as "the validation," CAR-Treg/in-vivo (Sonoma precedent) as a future milestone. Say "RNP knockout," not "CRISPRi screen," and don't claim any lab has committed. ([Marson CRISPR](https://www.science.org/doi/10.1126/science.abj4008))

---

## What to explicitly NOT do (traps the agents flagged)

- Don't train/fine-tune a bigger foundation model — scale didn't win Arc; priors + decision layer did.
- Don't game PDS by magnitude inflation (Arc is patching it) or chase MAE (structurally unbeatable vs the mean).
- Don't call CD226 a "novel target," or claim GxP as a result (n=4), or headline causal discovery or a retrospective active-learning curve.
- Don't attempt bitwise-reproducible GPU *training* on camera.
- Don't drift out of the autoimmunity/Treg lane or overclaim wet-lab validation.

---

## Updated one-liner & rubric map

> **"A pLDDT for perturbation biology: a calibrated trust layer that tells you which T-cell predictions to believe — proven by recovering held-out T1D genes, and cashed out as the next experiment to run. It recovers a genetics-anchored known target (CD226) and makes one novel, falsifiable bet (RASGRP1), with a critic that kills its own top pick on camera."**

- **Demo (30):** self-veto reshuffle to blinded recovered genes + a hash-verified reproduce.
- **Claude Use (25):** evidence-modality agent council + no-Write vetoing critic + with/without-veto ablation.
- **Impact (25):** build on Marson's data; the missing Treg-stability node vs teplizumab's 2-year ceiling; named validation path.
- **Depth (20):** conformal + selective-prediction + per-slice calibration on an OOD split; two mechanistically-defensible, genetics-anchored nominations.

---

## Key sources
Trust/method: [Weighted CP](https://arxiv.org/abs/1904.06019) · [CQR](https://arxiv.org/abs/1905.03222) · [SelectiveNet](https://arxiv.org/abs/1901.09192) · [PRESCRIBE](https://arxiv.org/abs/2510.07964) · [Kuleshov calibration](https://proceedings.mlr.press/v80/kuleshov18a/kuleshov18a.pdf) · [Ovadia ensembles](https://arxiv.org/abs/1906.02530) · [pLDDT](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/) · [AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492)
Competitive/angles: [Arc VCC](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up) · [IterPert](https://www.biorxiv.org/content/10.1101/2023.12.12.571389v1) · [CD4 causal GRN](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00290-8) · [autoimmune shared architecture](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10291128/)
Biology: [CD226 T1D](https://diabetesjournals.org/diabetes/article/72/11/1629/153546/) · [CD226 IND competitor](https://tnaxbio.com/2025/07/30/news/) · [CD226 contradiction](https://www.cell.com/cell-reports/fulltext/S2211-1247(23)01318-9) · [RASGRP1 GWAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC3272492/) · [program-axis atlas](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01300-z) · [GxP donor scale](https://www.nature.com/articles/s41588-025-02344-6)
Demo/Claude: [hackathon winners](https://claude.com/blog/meet-the-winners-of-our-claude-opus-4-8-build-day-hackathon) · [Agent SDK subagents](https://docs.claude.com/en/docs/agent-sdk/subagents) · [Virtual Lab](https://github.com/zou-group/virtual-lab) · [Claude Science](https://www.anthropic.com/news/claude-science-ai-workbench)
Impact/Gladstone: [Marson/Gladstone](https://gladstone.org/people/alex-marson) · [teplizumab TN-10](https://www.nejm.org/doi/full/10.1056/NEJMoa1902226) · [Treg failure in islet](https://pmc.ncbi.nlm.nih.gov/articles/PMC8275894/) · [Marson CRISPR in human T cells](https://www.science.org/doi/10.1126/science.abj4008)
