# README — Master index & current state

**Project:** an AI co-scientist for the *Built with Claude: Life Sciences* hackathon (Gladstone × Anthropic × Cerebral Valley, Research track, July 7–13 2026).
**Status:** planning complete through Iteration 2. This file is the single source of truth; the other files are detail/appendices and audit trail. Last consolidated after a 20-agent deep-dive.

> **Start here, then open `winning-plan.md` (execution) and `t1d-nomination-memo.md` (targets). `iteration-2-refinements.md` explains the latest reasoning.**

---

## The idea in one line

> **"A pLDDT for perturbation biology"** — a **calibrated trust layer** over CD4⁺ T-cell Perturb-seq that tells you *which* perturbation predictions to believe, proves it by recovering held-out Type 1 diabetes genes, and cashes out as the **next experiment to run** — recovering a genetics-anchored known target (**CD226**, the anchor) and making one genuinely novel, falsifiable bet (**RASGRP1**, the bet), with a **self-critiquing agent council** whose vetoing critic kills its own top pick on camera.

Four-beat arc: **predict → trust (calibrated "believe this / don't") → recommend next experiment → self-critique (live veto).**

Why "pLDDT for perturbation biology": AlphaFold's per-residue confidence (pLDDT) is what made it *usable* — people act only on high-confidence regions. Confidence-as-the-product is proven in structure (pLDDT/PAE) and variant effect (AlphaMissense's explicit "uncertain" class), but **nobody leads with it in perturbation biology** — so we're first to make it the headline.

---

## Locked decisions (do not relitigate without cause)

1. **Single focused story**, not a portfolio (the council was unanimous).
2. **Anchor disease = Type 1 diabetes** (proven T-cell-therapy axis via teplizumab; external genetics-derived core-gene validation set exists; clean CD4/Treg biology; avoids rivals' IBD/atopic-dermatitis).
3. **Headline identity = calibrated trust ("knows when it's wrong")**, framed as "a pLDDT for perturbation biology." Prediction accuracy and retrospective target-recovery are dead games — we compete on *decision under uncertainty*.
4. **Two-target nomination (v2.1-revised):** CD226 = genetics-anchored **anchor**, but reframed as the **hedged/moderate-trust calibration flagship** — its KD-on-suppression direction is *contested* in the literature (Ma 2023 *Cell Reports* vs Brusko 2023 *Diabetes*), so surfacing it with honest, hedged trust *is* the point. **PRKCQ is now the primary novel bet** (direction verified); **RASGRP1 only survives a Day-1 eQTL-direction go/no-go** (its hypermorphic-risk-allele premise is unconfirmed). See `revised_plans/` §3a.
5. **Build scope for the session that produced these docs = plan + deck**; the code is built during the event.
6. **Sequence→function chain is demoted to precomputed features** — not a trained model (vaporware in a week).

Two decisions still open (flagged as "fog"): recruit a biology partner vs. solo; which exact T1D GWAS locus becomes the demo's headline drill-down.

---

## The plan on one page

**Data (know it precisely).** [Marson/Pritchard genome-scale CD4⁺ T-cell Perturb-seq](https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1): every expressed gene, **CRISPRi = knockdown** (not knockout), ~22M cells, **4 donors**, **rest + stimulated**, GWAS-linked; public as GEO **GSE314342** / CZI. Plus **Krogan PPI** and **Pollard MPRA** used as *priors/features*.

**Method stack (rigorous, cheap, buildable):**
- Base predictor: **pseudobulk + delta/residual** prediction + **ESM-2 protein embeddings** + DEG-frequency/mean-expression priors (the Arc-VCC-winning substrate — borrow it cheaply).
- Uncertainty: **deep ensemble (3–5 seeds)** (best calibration + shift-robustness; disagreement flags mean-collapse).
- Guarantee: **conformal** — split-CQR + **Mondrian** (per-slice) + **weighted conformal** for unseen-gene covariate shift.
- Demo metric: **selective prediction** — risk-coverage curve + **AURC** vs a random-abstention null.
- OOD trust signal: **kNN-distance in the model's gene-embedding space**.
- Proof of calibration: reliability diagram + **proper score (Brier/NLL/CRPS)** + coverage-vs-sharpness + **per-slice (seen vs unseen) table** + constant baseline + bootstrap CIs. Evaluate on a **held-out-perturbation (OOD) split, never IID**. (Reject evidential deep learning — overconfident.)

**Program axis (what "beneficial" means).** UCell **signed contrast**: beneficial (stable-Treg module + Tr1 module) − pathogenic (ex-Treg/Th1 *gated on low FOXP3* + Th17 + autoreactive BHLHE40/CSF2/TNF); validate vs the Genome-Medicine-2024 **TMZ score**, TSDR methylation, CITE-seq. Use UCell/AUCell (composition-robust); regress out cell-cycle/composition/donor/sex/HLA.

**Nominations.** CD226 (anchor: coding *candidate-causal* risk variant rs763361 — risk-allele direction secure, but the **KD-on-Treg-suppression direction is CONTESTED**: Ma et al., *Cell Reports* 2023;42(10):113306 [PMID 37864795] finds Treg-conditional deletion *impairs* suppression and worsens GvHD, opposite to Brusko/Thirawatananond *Diabetes* 2023's NOD protection. **Do NOT call CD226 "direction-secure" — make it the hedged/moderate-trust calibration flagship** that cites the split on camera; a competitor is IND-cleared on the anti-DNAM-1 mechanism [TNAX asset]). RASGRP1 (bet: independent T1D locus, but the **hypermorphic-risk-allele direction is UNCONFIRMED** — not a significant RASGRP1 eQTL in any CD4/Treg resource, and germline loss *causes* autoimmunity/lymphoma; **treat as a Day-1 eQTL go/no-go, else swap to PRKCQ** as primary bet — PKC-θ blockade *enhances* Treg suppression, direction verified). **Landmines — do NOT nominate:** RNF20 (Marson's own hit, Nature 2020), RBPJ-NCOR (**Sakaguchi's** hit, Nature 2025 — a landmine, but NOT Marson's), SIRPG (wrong direction for KD), BACH2/IKZF4 (positive Treg regulators → KD destabilizes). **Genotype × perturbation = labeled hypothesis only** (n=4 donors can't support it). *(See `revised_plans/PLAN_OF_ATTACK_v2.pdf` §3a/§3b for the full v2.1 reconciliation.)*

**Angle upgrades (nested under the trust headline, don't dilute):** active learning = the *payoff* ("which perturbation next" proves the trust score is decision-grade); causal regulator discovery = *supporting* mechanism (CD4 T-cell precedent, gate edges by trust); cross-disease **transfer panel** (RA/celiac/MS) by swapping public GWAS on the same map = cheap novelty + impact.

**Agent council (creative Claude use).** Claude Code + Agent SDK, agents as `.claude/agents/*.md`, biology **MCP connectors** (Open Targets, ChEMBL, PubMed, 10x from `anthropics/life-sciences`). Evidence-**modality** agents (genetics / expression / network / literature) → **Elo tournament** (open Swarms code) → **vetoing critic in a fresh context with NO Write access** (tool-scoping = its teeth). Show the with/without-veto **hallucination-rate delta**.

**Demo (heaviest score).** Money shot = the critic **vetoes the pipeline's own #1**, live, in a visibly separate context window; leaderboard **reshuffles** to blinded, held-out **recovered T1D genes** (FOXP3/CTLA4/STAT1 from the genetics-derived set) with **one real recovery number** on screen; then a **hash-verified "Reproduce"** (deterministic eval, not GPU training). Lead with it (teaser at 0:00); reserve the final day for the video.

**Impact / Gladstone framing.** Build on **Marson's own data** (he's at Gladstone, names T1D, runs Foxp3/Treg CRISPR screens); ship a **reusable tool**, not a one-off (Pollard's value). The pitch: the missing **Treg-stability node** — Tregs collapse in the inflamed islet (IL-2-deprivation apoptosis, TSDR hypermethylation → ex-Tregs), the shared failure point of every current therapy — vs teplizumab's **~2-year** disease-modification ceiling ($194k/course; T1D ≈ $500k lifetime, ~9 life-years lost). Named **Tier-1 validation**: RNP knockdown in primary human CD4/Treg → flow FOXP3/TIGIT + suppression assay (~4–8 wks, low-$; Marson/Levings labs; Synthego/EditCo). Stay in the autoimmunity lane; don't overclaim validation ("Lead with Integrity" is a Gladstone value).

---

## Rubric map (weights unconfirmed publicly; robust either way)

| Axis (est. weight) | Our play |
|---|---|
| **Demo (30)** | Live self-veto reshuffle to blinded recovered genes + hash-verified reproduce |
| **Claude Use (25)** | Evidence-modality agent council + no-Write vetoing critic + with/without-veto ablation |
| **Impact (25)** | Build on Marson's data; the missing Treg-stability node vs teplizumab's 2-yr ceiling; named validation path |
| **Depth (20)** | Conformal + selective-prediction + per-slice calibration on an OOD split; two genetics-anchored nominations |

---

## What to explicitly NOT do (traps)

Don't train a bigger foundation model (scale didn't win Arc; priors + decision layer did) · don't game PDS by magnitude inflation or chase MAE · don't call CD226 "novel" · **don't call CD226 "direction-secure"** (the KD-on-suppression sign is contested — Ma 2023 vs Brusko 2023) · **don't present RASGRP1 as a confident bet** (direction unconfirmed — gate it or use PRKCQ) · don't claim genotype×perturbation as a result (n=4) · don't headline causal discovery or a retrospective active-learning curve · don't attempt bitwise-reproducible GPU *training* on camera · don't drift out of the autoimmunity/Treg lane or overclaim wet-lab validation.

## Factual corrections carried forward

CRISPRi = **knockdown**, not knockout · the Arc VCC overall winner was a ~100B-param model, real lesson = **metric-gaming**, not "scale wins" · **BioPlex ≠ Krogan** (Krogan = AP-MS host–pathogen interactomes) · closest prior art to differentiate from = **Ota/Pritchard 2025** (omnigenic × Perturb-seq, same lab/data — we do prospective calibrated nomination, not attribution) · **[v2.1] the CD226 Treg-deletion contradiction IS real** — Ma et al., *Cell Reports* 2023;42(10):113306, PMID 37864795 (verified via NCBI); an earlier audit wrongly called it nonexistent. Its title contains "mTOR/Myc." Use it, don't delete it · **[v2.1] RBPJ-NCOR is Sakaguchi's (Nature 2025), NOT Marson's** — only RNF20 is Marson's (Nature 2020) · **[v2.1] RASGRP1's hypermorphic risk allele is unconfirmed** in CD4/Treg eQTL data — the "risk↑RASGRP1" signal is a *different SNP in SLE*; germline loss causes autoimmunity/lymphoma.

---

## File index

| File | What it is | Status |
|---|---|---|
| **README.md** | This master index + current-state snapshot | **CURRENT — start here** |
| **SCOPE.md** | The refined project scope (what we're building) | **CURRENT** |
| **research-state.md** | Decisions log + state tracker | **CURRENT** |
| **winning-plan.md** | Execution plan: wayfinder tickets, MLE data-contract/eval-gates, TDD seams, demo, risks | **CURRENT** (Iteration-1 base + Iteration-2 banner) |
| **t1d-nomination-memo.md** | Target nominations (CD226 anchor + RASGRP1 bet), validation design, gget/PubMed command pack | **CURRENT** (Iteration-2 reconciled) |
| **iteration-2-refinements.md** | The 20-agent deep-dive that produced the latest refinements + rationale | **CURRENT — rationale/audit** |
| **improvement-research-report.md** | Initial competitive/landscape research (why the original idea was a trap) | **BACKGROUND — audit** |
| **pitch-deck.html** | The slide deck (open in a browser) | **DELIVERABLE** (rebrand pending) |

*Historical note:* `SCOPE.md` and `research-state.md` originally described a "Generalization Frontier / chained-vs-direct" idea that the adversarial council retired; they have been rewritten to the current direction. The chained-model concept survives only as *precomputed features*.

---

## Using this from Claude Science

Point Claude Science at this folder and start from `README.md`, then `winning-plan.md`. The plan already lists the exact method stack, data contract, eval gates, and agent architecture to build. When you begin building: freeze the held-out (perturbation- and donor-level) splits first, stand up the eval harness + baseline gauntlet before any modeling, and build the vetoing critic agent from hour one.
