# Research State — "A pLDDT for perturbation biology" (calibrated trust for T-cell Perturb-seq, T1D)

*Live decision log + state tracker. Current as of **v3 BUILD, Day 0-1**. Supersedes Iteration 2 (which pre-dates the v3 PRKCQ reversal). See `revised_plans/PLAN_OF_ATTACK_v3.html` for the governing plan.*

## Current stage
BUILD — Day 2 COMPLETE. Base predictor + Treg axis built on real data, both critique gates PASS, 15/15 tests. Day-2 methodology written, PASSED pre-COMPUTE critique (1 fix cycle) and post-COMPUTE critique. Day 0-1 foundation COMPLETE.

## Day 2 results
- **Base predictor** (leakage-safe, 5 covariate features, HGBR): test R²=+0.096, Spearman=+0.356; mean baseline R²≈0; shuffle R²<0. Modest R² is the IDEAL regime — signal exists but leaves headroom for the trust layer to matter (v3 §1). Target = log1p(||trans zscore, on-target excluded||).
- **Treg program axis** (Fix 3 dependency): 13 POS + 10 NEG genes, gold-set-disjoint (CTLA4/FOXP3/IL2RA dropped), per-row on-target exclusion. RASGRP1 KD pushes toward Treg program (+0.85 Rest) — consistent with KD-mimics-protection thesis.
- **Leakage fix (critique-driven):** original ||zscore|| target contained the on-target gene → circularity with ontarget_effect_size feature. Fixed to trans-only target + forbidden-feature assertion. Critique confirmed "circularity genuinely severed."
- Provenance: _script_manifest.jsonl (11 output entries, sha256).

## Day 3 results (HEADLINE — v3 Fix 6)
- **Donor-blocked LODO conformal coverage** (calibrate on 3 donors, test on held-out 4th; 2591 cross-donor-reproducible genes): 80%→**0.793±0.010**, 90%→**0.895±0.005**, 95%→**0.949±0.005**. Every fold-std interval CONTAINS nominal. This is the "pLDDT for perturbation biology" proof — coverage holds under the hardest honest regime. Figure: figures/coverage_calibration.png.
- **Selective risk**: model AURC 0.11510 vs 2000-permutation null 0.11826±0.00047, **p=0.0005** — trust ordering is modest but statistically significant.
- **Honesty**: row-level CP is overconfident (non-independent rows); report fold-std across 4 donors as the interval. Recovery stays a pooled-data case study (CD226 + 6/8 gold genes absent from donor-pair data — the 2591-gene reproducible core is what donor-blocking can evaluate).
- Data: by_donors.h5mu (16.87GB, 6 pairs / 4 donors) pulled + verified.

## Day 4 (in progress)
Blinded recovery (pooled, n=8 case study) + axis-swap specificity control. Methodology PASSED critique after 1 fix cycle.
- Pre-COMPUTE Gate-1 (Day-4) cycle 1: **BLOCKINGx2** — (1) interferon control confounded (IFIH1 is a gold gene AND an IFN-pathway gene → shared mechanism, not orthogonal); (2) magnitude confound only detected, not controlled.
- cycle 2 (after fix): **PASS** — (1) IFN dropped, replaced with cholesterol/metabolic + ribosome-biogenesis controls, + leave-IFIH1-out sensitivity; (2) magnitude-only baseline axis the Treg axis must beat + residualize-on-magnitude. Primary metric → Mann-Whitney U + exact permutation null (AUROC descriptive only). Expect an honest negative at n=8 as the modal outcome.

## Day 4 results (HONEST NEGATIVE — correctly scoped)
- **Blinded recovery (n=8 case study):** the Treg axis does NOT recover the 8 genetics-anchored T1D genes above chance (AUROC=0.439, perm p=0.72). No control axis, IFN-shared, or magnitude baseline is significant either.
- **Specificity:** Treg beats no control and does not beat the magnitude baseline. Residualizing on magnitude and leaving IFIH1 out don't change the null.
- **Interpretation (critique-verified honest):** this is *absence of evidence* at n=8 (underpowered), NOT *evidence of absence* (method disproven). Recovery was always scoped as a case study (v3 Fix 6); the **Day-3 donor-blocked coverage remains the headline**. The project honestly reports its secondary validation as null rather than p-hacking a positive — which *is* the intellectual-honesty thesis in action.
- Both critique gates PASS. Figure: figures/recovery_specificity.png.
- Pre-COMPUTE Gate-1 (Day-4) cycle 1: BLOCKINGx2 (IFN control confounded by IFIH1; magnitude confound uncontrolled) → cycle 2 PASS (IFN dropped + cholesterol/ribosome controls + magnitude baseline + residualization; Mann-Whitney + exact perm null). Gate-2 PASS.

## Next (Day 5)
Wire critic → commit_gate; the believe/veto split gets a real trust_score from the conformal layer. The RASGRP1 veto is strengthened by two independent real-data facts now: (1) no cell-type-matched eQTL, (2) low cross-donor reproducibility (r=0.072, n_downstream=992).

## Critique history
- Pre-COMPUTE Gate-1 (Day-2 methodology), cycle 1: **BLOCKING** — target-feature circularity: effect_magnitude=||zscore|| contained the on-target gene, and ontarget_effect_size was a feature (target = sqrt(ontarget² + Σtrans²)). On-target is typically the largest norm component → non-negligible leak.
- Pre-COMPUTE Gate-1, cycle 2 (after fix): **PASS** — "core circularity genuinely severed." Fixes: (1) target = trans-only ||zscore|| (perturbed gene excluded); (2) all DE-derived quantities forbidden as features (build-time assertion); (3) Treg axis excludes perturbed gene per-row. Non-blocking note: target_baseMean predicting trans-effect is legitimate signal, not leakage.
- Post-COMPUTE Gate-2 (Day-2 outputs): **PASS** — every number traceable to script output; trans-only norm correct; IL17A->10 NEG genes honest (unmeasured gene dropped). Recommended _script_manifest.jsonl (added).
- Pre-COMPUTE Gate-1 (Day-3 conformal), cycle 1: **BLOCKING x3** — (1) donor-pair LODO didn't enforce full donor disjointness (calibration pairs shared donors with test -> exchangeability breach); (2) coverage at 80/90/95% unpowered at n_eff~tens/4 folds (over-claim); (3) selective-prediction circularity (same residuals calibrate + score + evaluate).
- Pre-COMPUTE Gate-1 (Day-3), cycle 2 (after fix): **PASS** — all 3 resolved. Fixes: (1) calibration = pairs with both members != held-out donor, fail-closed assertion; (2) demoted to "CP interval contains nominal" + report width, pool across folds; (3) three-way separation (quantiles on calib, trust+AURC on held-out-donor test only). Implementation notes to honor: recompute residuals PER-FOLD (not Day-2 global fit), n_eff counts perturbations within test pairs, CP interval will be wide (~12 pts).
- Post-COMPUTE Gate-2 (Day-3 results): coverage headline VERIFIED (donor-disjoint LODO, per-fold target, 3-way separation all correct; calibration sound). **1 BLOCKING** on selective-prediction: degenerate ensemble spread (rng re-seeded identically each fold) + AURC gain rested on a single shuffle. FIXED: diverse bootstrap (distinct seed per member+fold) + 2000-permutation null → p=0.0005. Stale comment removed. Claim softened to "modest but significant." Environment, live genetics gates, tiered cited gold set, deterministic commit_gate() veto, and frozen leakage-safe splits are all built and passing on REAL data (11/11 tests). The 16.79 GB Perturb-seq DE_stats delta file is downloaded, byte-verified against S3, and confirmed to open as 33,983×10,282 with log_fc+zscore layers. Ready for Day 2 (base predictor + UCell axis).

## Data (verified this session)
- `data/raw/GWCD4i.DE_stats.h5ad` — 16.79 GB, from CZI Virtual Cells Platform public S3 (MIT license). Byte-exact match to S3 content-length; opens 33983×10282, layers log_fc/zscore/p_value/adj_p_value/baseMean/lfcSE; all trio + gold genes present. Provenance: `data/gold/data_provenance_receipt.json`.
- `data/raw/*.suppl_table.csv` — DE_stats obs table, guide KD efficiency, sgRNA library, sample metadata, autoimmune enrichment (from GitHub mirror).
- Source found via: CZI VCP page https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq → S3 bucket `s3://genome-scale-tcell-perturb-seq/marson2025_data/`.
- NOTE for Day 3: donor-blocked LODO conformal (v3 Fix 6) needs `GWCD4i.DE_stats.by_donors.h5mu` (16.87 GB, same bucket) — not yet pulled; splits currently mark donor_blocking=UNAVAILABLE.

## v3 corrections applied this session (vs Iteration 2)
- **Target hierarchy reversed** (v3 Fix 1): CD226 = ANCHOR (genetic_association 0.834), RASGRP1 = novel BET (0.506), PRKCQ = DEMOTED cross-trait control (0.162). Verified live via Open Targets — hierarchy PASSES.
- **Gold set re-derived** (v3 Fix 2): retired MEOX1/CD1E (no T1D association) and LGALS3BP/CD247 (outside top-250); replaced with 8 GWS-confirmed T1D genes. INS dropped (not perturbed in CD4).
- **Veto rebuilt as deterministic `commit_gate()`** (v3 Fix 4): plain Python, default-deny, 7/7 tests pass. Replaces the falsifiable "no-Write tool scoping" claim.
- **Headline moved to calibration** (v3 Fix 6): real data confirms only 4,775/33,983 pairs have a cross-donor estimate, median r=0.41, 26% r<0.2 → effective-N is tens, recovery (n=8) is a CI-shown case study.

## Build artifacts produced (data/gold/)
- `genetics_gate_receipt.json` — live Open Targets + GWAS Catalog credentials for the trio (hierarchy PASS)
- `t1d_gold_set.json` / `.tsv` — tiered (T1-T4) per-gene cited gold set; n=8 blinded recovery positive set
- `data_facts_receipt.json` — effective-N + target coverage from the real DE_stats table
- `believe_veto_split.json` + `nominations/` — CD226 GO, RASGRP1 WITHHELD, PRKCQ WITHHELD (deterministic, lineage-backed)

## Research question (refined)
Can a **calibrated trust layer** over CD4⁺ T-cell perturbation-effect predictions tell you *which* predictions to believe — proven by recovering held-out Type 1 diabetes core genes and by picking the next experiment — and used to nominate a genetics-anchored target (CD226) plus one genuinely novel bet (RASGRP1)? Headline identity: "a pLDDT for perturbation biology." Prediction accuracy and retrospective ranking are the baselines we deliberately do NOT compete on.

## Context (hackathon)
- Event: Built with Claude: Life Sciences (Gladstone × Anthropic × Cerebral Valley), Research track, **built on Claude Science**. July 7–13 2026.
- Team ≤2. Prizes: 1st $30k / 2nd $10k / 3rd $5k + **$10k Gladstone Award** ("science that can overcome disease").
- Rubric (from user notes, publicly unconfirmed): Impact 25 / Claude Use 25 / Depth 20 / Demo 30 — plan is robust to reasonable variations.
- Provided Gladstone datasets: Marson/Pritchard CD4⁺ T-cell Perturb-seq (core), Pollard MPRA (feature prior), Krogan PPI (feature/mechanism prior).

## Core thesis
Perturbation prediction is saturated and retrospective ranking is crowded; the open lane is **decision under uncertainty**. Ship a calibrated trust score (conformal + selective prediction) that's provably calibrated per-slice on an OOD split, prove it by blinded recovery of the genetics-derived T1D core genes, cash it out as the next experiment, and surface CD226 (anchor) + RASGRP1 (bet). A self-critiquing agent council makes it a winning demo. Build on Marson's own data; ship a reusable tool; frame impact as the missing Treg-stability node vs teplizumab's 2-year ceiling.

## Competitive intel
- **Rivals (Discord):** Hyun (Cambridge) — IBD/autoimmune signature reversal; say181 (UCSF) — Th2/atopic-dermatitis "next-degrader." Both do retrospective Perturb-seq → ranked drug targets + GWAS check. We differentiate by being calibrated + prospective + decision-grade, and by anchoring to **T1D** (avoids their diseases).
- **State of the art (Arc VCC 2025):** all four winners used hybrid DL + classical stats + protein embeddings; none built for calibration/trust/decision → our angle is orthogonal. Borrow their substrate (pseudobulk+delta+ESM-2+DEG priors) as the cheap base model.
- **Closest prior art to differentiate from:** Ota/Pritchard 2025 (omnigenic × Perturb-seq, same lab/data — they do attribution; we do prospective calibrated nomination); PRESCRIBE (confidence for perturbation prediction — a method, not a headline product).
- **CD226 competition:** a competitor is IND-cleared on the CD226-block→raise-FOXP3 mechanism (Riverview RVW101) — cite as thesis-validation; differentiate on Treg-selective modality + T1D indication.

## Key decisions (locked)
1. Single focused story (not a portfolio). — *council, unanimous*
2. Disease anchor = **Type 1 diabetes**. — *grilling*
3. Headline = **calibrated trust ("knows when it's wrong")**, branded "a pLDDT for perturbation biology." — *council + iteration 2*
4. Nomination = **two-target**: CD226 (anchor) + RASGRP1 (bet); backup PRKCQ. — *iteration 2*
5. Session scope = **plan + deck**; build during the event. — *grilling*
6. Sequence→function chain **demoted to precomputed features** (not trained). — *council + ML red-team*
7. Method stack = base(pseudobulk+delta+ESM-2+priors) → deep ensemble → conformal(CQR/Mondrian/weighted) → selective-prediction(AURC) → kNN-OOD → per-slice/OOD calibration proof. — *iteration 2, Cluster A*
8. Angle upgrades nested under trust: active learning (payoff), causal discovery (support), cross-disease transfer (RA/celiac/MS). — *iteration 2, Cluster B*
9. Program axis = UCell signed contrast (stable-Treg + Tr1) − (ex-Treg/Th1 gated on low FOXP3 + Th17 + autoreactive). — *iteration 2, Cluster C*
10. Genotype × perturbation = **hypothesis only** (n=4 donors). — *iteration 2, Cluster C*
11. Agent council = evidence-modality agents + Elo tournament + no-Write vetoing critic; MCP connectors. — *iteration 2, Cluster D*
12. Reproducibility = deterministic-eval hash-match (not GPU training). — *iteration 2, Cluster D*
13. CRITICAL anti-leakage rule (carried forward): freeze **perturbation- and donor-level** held-out splits before any modeling.

## Experiment log
| Attempt | Method | Result | Status |
|---------|--------|--------|--------|
| (none yet — build begins on Claude Science) | | | |

## Critique history
- **Adversarial council (4 seats):** killed the "in-silico engineer / recover known targets" headline (wrong cell type, circular metric, pre-empted by Ota/Pritchard); converged on trust + next-experiment reframe. DONE.
- **Skeptical-immunologist red-team:** ranked CD226 top of the original shortlist; caught landmines (RNF20/RBPJ = Marson's own; SIRPG wrong direction). DONE.
- **20-agent deep dive (Iteration 2):** validated the trust angle as fresh + orthogonal; forced the CD226→anchor + RASGRP1→bet restructure; added the method stack, cross-disease panel, and Gladstone/impact framing. DONE.
- **Pre-COMPUTE critique of the build:** pending (run the critic on the eval design for leakage before modeling).

## What worked
- Reframing from accuracy/attribution → calibrated decision under uncertainty (orthogonal to all rivals + Arc winners).
- The "recover a known anchor AND make a novel bet" nomination structure.

## What didn't work (retired)
- Original "chained-vs-direct on held-out genes" headline (walks into the baseline buzzsaw).
- "In-silico T-cell engineer / recover RASA2" headline (wrong cell type; circular; pre-empted).
- CD226 as a "novel target" claim (competitor already at IND; 2023 contradiction).

## Open questions
- Recruit a biology partner (Marson/Levings-adjacent) vs. solo with parallel agents?
- Which exact T1D GWAS locus becomes the demo's headline drill-down?
- Confirm RASGRP1 T1D-risk-allele eQTL direction in primary CD4/Treg before committing it as the bet.
- Final blinded gold set for RQ2 (genetics-derived T1D core genes) + publication-time split details.

## Artifacts (current folder)
- `README.md`: DONE — master index + current-state snapshot (start here).
- `SCOPE.md`: DONE — refined scope (rewritten from the retired Generalization Frontier).
- `research-state.md`: DONE — this file.
- `winning-plan.md`: DONE — execution plan (wayfinder tickets, MLE eval-gates, TDD seams, demo); Iteration-2 banner added.
- `t1d-nomination-memo.md`: DONE — CD226 (anchor) + RASGRP1 (bet), validation design, gget/PubMed command pack.
- `iteration-2-refinements.md`: DONE — 20-agent deep-dive synthesis (rationale/audit).
- `improvement-research-report.md`: DONE — initial competitive/landscape research (background/audit).
- `pitch-deck.html`: DONE — deck (rebrand to "pLDDT for perturbation biology" + CD226-anchor/RASGRP1-bet slide pending).
- `figures/`: empty (build phase).
- `experiments.tsv`: missing (build phase).
