# SCOPE — A pLDDT for Perturbation Biology (calibrated trust for T-cell Perturb-seq, anchored to Type 1 diabetes)

*This is the current, refined scope (post adversarial council + grilling + a 20-agent deep dive). It replaces the original "Generalization Frontier / chained-vs-direct" scope, which was retired — see `README.md` and `iteration-2-refinements.md` for why. The old chained-model idea survives only as precomputed features.*

## The one-sentence pitch
Perturbation-effect prediction is a saturated game (deep models don't beat linear baselines) and retrospective target-recovery is what everyone else does — so we don't compete on accuracy. We build a **calibrated trust layer over CD4⁺ T-cell Perturb-seq — "a pLDDT for perturbation biology"** — that tells you *which* predictions to believe, proves it by recovering held-out Type 1 diabetes genes, and cashes out as the **single next experiment to run**: a genetics-anchored known target (**CD226**, the anchor) plus one genuinely novel, falsifiable bet (**RASGRP1**, the bet). A **self-critiquing agent council** whose vetoing critic kills its own top pick — live, on camera — carries the demo.

## The four-beat arc
**predict → trust (calibrated "believe this / don't") → recommend the next experiment → self-critique (a fresh-context critic with veto power).**

## Why this wins (and why the original framing didn't)
- **Prediction accuracy is dead** ([Ahlmann-Eltze, *Nat Methods* 2025](https://www.nature.com/articles/s41592-025-02772-6); [Arc VCC 2025](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up)). So is retrospective ranking (what the rival teams do). The empty, winnable lane is **decision under uncertainty**.
- **Confidence-as-the-product is proven and legible** (AlphaFold pLDDT/PAE; AlphaMissense's "uncertain" class) but **nobody leads with it in perturbation biology** — the mechanics exist as methods (PRESCRIBE) while the virtual-cell field openly calls calibrated uncertainty "the missing piece." We're first to make it the headline.
- **Orthogonal to the state of the art:** none of the four Arc VCC 2025 winners built for calibration/trust/decision-making — so we borrow their accuracy substrate cheaply and spend the whole novelty budget on the trust + decision layer.
- **Recover a known anchor AND make a novel bet** is exactly the structure that separates a winner from a finalist.

## The datasets (and how each is used)
| Layer | Dataset | Role |
|---|---|---|
| Perturbation → expression (core) | [Marson/Pritchard genome-scale CD4⁺ T-cell Perturb-seq](https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1) — every gene, **CRISPRi knockdown**, ~22M cells, 4 donors, **rest + stim**, GWAS-linked (GEO GSE314342 / CZI) | The substrate; source of held-out perturbations and the T1D program axis |
| Protein context (prior) | **Krogan PPI** | Mechanism annotation + network-distance features; edges for causal support |
| Sequence → regulatory (prior) | **Pollard MPRA** | Precomputed variant/regulatory **features** only — NOT a trained sequence model |

## Research questions (refined)
- **RQ1 (trust):** Can a calibrated confidence (conformal + selective prediction) reliably separate trustworthy from untrustworthy CD4 perturbation-effect predictions, with *per-slice* calibration that holds on a held-out-perturbation (OOD) split?
- **RQ2 (validation):** Does the high-trust set recover the genetics-derived T1D "core genes" (FOXP3/CTLA4/STAT1) that were held out — a non-circular external check?
- **RQ3 (decision):** Does the trust score, used as an acquisition function, pick better next experiments than random/uncertainty baselines (active learning as the payoff)?
- **RQ4 (transfer):** Does the trust score predict *where* the engine transfers across CD4-driven autoimmune diseases (RA / celiac / MS)?
- **RQ5 (nomination):** Which knockdowns move CD4 cells toward the beneficial (stable-Treg) pole, with a mechanism and a trust score — yielding CD226 (anchor) and RASGRP1 (bet)?

## Method stack (rigorous, cheap, buildable)
- **Base predictor:** pseudobulk + **delta/residual** prediction + **ESM-2 protein embeddings** + DEG-frequency/mean-expression priors (borrowed from the Arc VCC winners).
- **Uncertainty:** **deep ensemble (3–5 seeds)** — best calibration/shift-robustness; ensemble disagreement flags mean-collapse. (Reject evidential DL — overconfident.)
- **Guarantee:** **conformal** — split-CQR + **Mondrian** (per-slice) + **weighted conformal** for unseen-gene covariate shift.
- **Demo metric:** **selective prediction** — risk-coverage curve + **AURC** vs a random-abstention null.
- **OOD trust signal:** **kNN-distance in the model's gene-embedding space**.
- **Calibration proof:** reliability diagram + **proper score (Brier/NLL/CRPS)** + coverage-vs-sharpness + **per-slice (seen vs unseen) table** + constant baseline + bootstrap CIs. **Evaluate on a held-out-perturbation (OOD) split, never IID.**

## The T1D "program axis" (what "beneficial" means)
A **UCell signed contrast**: beneficial (stable-Treg module: FOXP3/IL2RA/CTLA4/IKZF2/IKZF4/TNFRSF18/TIGIT/LRRC32 + Tr1 module: IL10/MAF/PRDM1/LAG3/ITGA2/…) **minus** pathogenic (ex-Treg/Th1: IFNG/TBX21/CXCR3 *gated on low FOXP3* + Th17: RORC/IL17A/IL23R + autoreactive: BHLHE40/CSF2/TNF). Validate against the [Genome-Medicine-2024 TMZ score](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01300-z), TSDR methylation, and CITE-seq. Use UCell/AUCell (composition-robust); regress out cell-cycle, composition, donor, sex, HLA.

## Nominations (see `t1d-nomination-memo.md` for the full case)
- **CD226 — the anchor, reframed as the hedged calibration flagship (v2.1).** Coding *candidate-causal* T1D risk variant (rs763361); the **risk-allele direction is secure, but the KD-on-Treg-suppression direction is CONTESTED**: [Ma et al., *Cell Reports* 2023;42(10):113306, PMID 37864795](https://pubmed.ncbi.nlm.nih.gov/37864795/) finds Treg-conditional CD226 deletion *impairs* suppression and worsens GvHD — the opposite of Brusko/Thirawatananond *Diabetes* 2023's NOD protection (same Foxp3-Cre design). **Do NOT present CD226 as "direction-secure."** Instead, have the pipeline surface it and assign it **moderate/hedged trust**, citing the Ma-vs-Thirawatananond split on camera — this makes CD226 the demo's proof that the trust layer is honest about contested biology. A competitor is IND-cleared on the anti-DNAM-1 mechanism (TNAX asset — describe the asset, not the "RVW101" code). Keep the human-Treg KD suppression experiment, framed as *resolving a real, cited split* (stabilization vs blunted activation at matched activation; TIGIT-dependence).
- **PRKCQ — the primary novel bet (v2.1); RASGRP1 = gated fallback.** RASGRP1's premise (hypermorphic risk allele → KD mimics protection) is **UNCONFIRMED** — the T1D lead SNP is not a significant RASGRP1 eQTL in any CD4/Treg resource (the "risk↑RASGRP1" signal is a *different SNP in SLE*), and germline loss *causes* autoimmunity/lymphoma. Treat RASGRP1's eQTL-direction check as a **hard Day-1 go/no-go**; if not cleanly hypermorphic in CD4/Treg, **PRKCQ is the primary bet** (PKC-θ blockade *enhances* Treg suppression — direction verified, druggable, loss tolerated). If RASGRP1 is kept, label it *"speculative, direction-unconfirmed"* and use tunable, Treg-restricted CRISPRi.
- **Do NOT nominate:** RNF20 (Marson's own hit, Nature 2020), RBPJ-NCOR (**Sakaguchi's** hit, Nature 2025 — a landmine, but NOT Marson's), SIRPG (wrong KD direction), BACH2/IKZF4 (positive Treg regulators). **Genotype × perturbation = labeled hypothesis only** (n=4 donors).

## Angle upgrades (nested under the trust headline — do not dilute)
- **Active learning = the payoff** — the trust score used to pick the next perturbation proves it's decision-grade (sanity-check curve with random + prior baselines; don't headline it).
- **Causal regulator discovery = supporting mechanism** — small/local, falsified by leave-one-perturbation-out interventional prediction, edges gated by trust (CD4 T-cell precedent exists).
- **Cross-disease transfer panel** — reuse the same map, swap public GWAS for RA/celiac/MS: "a T1D-calibrated engine that predicts where it transfers across autoimmunity."

## Creative Claude use — the agent council
Claude Code + Agent SDK; agents as `.claude/agents/*.md`; biology **MCP connectors** (Open Targets, ChEMBL, PubMed, 10x). Evidence-**modality** generators (genetics / expression / network / literature) → **Elo tournament** → **vetoing critic in a fresh context with NO Write access** (tool-scoping is the veto's teeth). Extend Claude Science's built-in reviewer into a *vetoing, multi-modality* critic bench; demo the with/without-veto **hallucination-rate delta**.

## The demo (heaviest score)
Money shot: the critic **vetoes the pipeline's own #1 candidate**, live, in a visibly separate context window; the leaderboard **reshuffles** to blinded, held-out **recovered T1D genes** with **one real recovery number** on screen; then a **hash-verified "Reproduce"** (deterministic eval — NOT GPU training). Lead with a 0:00 teaser; reserve the final day for the video.

## Impact & Gladstone framing
Build on **Marson's own dataset** (Gladstone; he names T1D and runs Foxp3/Treg CRISPR screens); ship a **reusable, auditable tool** (Pollard's value). The asset is the missing **Treg-stability node** — Tregs collapse in the inflamed islet (IL-2-deprivation apoptosis, TSDR hypermethylation → ex-Tregs), the shared failure point of every current therapy — measured against teplizumab's **~2-year** disease-modification ceiling. Named **Tier-1 validation**: RNP knockdown in primary human CD4/Treg → flow FOXP3/TIGIT + suppression assay.

## Scope boundaries (what we will NOT do)
- No foundation model trained/fine-tuned from scratch; no trained sequence→function chain (Pollard MPRA = features only).
- No headline claim on raw prediction accuracy, PDS/MAE, or a retrospective active-learning curve.
- No genotype×perturbation *result* (n=4 donors → hypothesis only).
- No beta-cell replacement, CGM, or risk-dashboard framing; stay in the autoimmunity/Treg lane.
- No wet-lab validation *claimed* — only a named, orderable validation *path*.

## Risk register
- **CD226 direction contested (Ma 2023 *Cell Reports* vs Brusko 2023 *Diabetes*)** → do NOT claim "direction-secure"; reframe CD226 as the **hedged/moderate-trust calibration flagship** that cites the split on camera (turns the landmine into the demo's honesty proof); monitor the TNAX competitor asset.
- **RASGRP1 direction unconfirmed** → Day-1 eQTL-direction go/no-go on the T1D lead SNP in CD4/Treg; swap to **PRKCQ** (direction verified) if it fails.
- **RASGRP1 toxicity window** → tunable, Treg-restricted CRISPRi; confirm eQTL direction first.
- **Calibration debunked on subgroups** → lead with proper scores + per-slice/OOD table + bootstrap CIs, not ECE.
- **Reproducibility theater** → scope "Reproduce" to a deterministic eval with content-addressed hashing, not GPU training.
- **Multi-agent token cost (~15×)** → cap Elo matches, Sonnet for judges/critics, Opus for orchestrator.

## Revised week plan
- **Day 1:** Lock the reframe; freeze perturbation- and donor-level held-out splits; stand up the agent council **including the critic from hour one**; write the recovery-of-known-genes eval first.
- **Day 2:** Harmonize data to Ensembl; build the base predictor (pseudobulk+delta+ESM-2+priors) and the UCell program axis.
- **Day 3:** Deep ensemble → conformal (CQR/Mondrian/weighted) → selective-prediction curve + AURC; per-slice/OOD calibration table.
- **Day 4:** Recover held-out T1D core genes (blinded); surface CD226 (anchor) + RASGRP1 (bet) with trust scores and mechanisms; red-team/critic veto pass.
- **Day 5:** Active-learning payoff curve; cross-disease transfer panel; causal-support edges gated by trust.
- **Day 6:** Finalize nominations + validation designs; reproducibility hash-match; figures; honesty pass.
- **Day 7:** Record the demo video (whole day); ship reproducible repo + writeup.

## Deliverables
1. A **calibrated trust layer** ("a pLDDT for perturbation predictions") with per-slice/OOD calibration proof.
2. **Blinded recovery** of genetics-derived T1D core genes as external validation.
3. **Two nominations** — CD226 (anchor) + RASGRP1 (bet) — each with mechanism, trust score, and a one-week validation design.
4. A **self-critiquing agent council** with a demonstrated with/without-veto hallucination-rate delta.
5. Reproducible, hash-verified repo + a demo-first 3-minute video.
