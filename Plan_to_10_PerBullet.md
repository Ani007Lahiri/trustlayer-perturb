# Plan to 10/10 — Per-Bullet Execution Map (Honest)

**Date:** 2026-07-12 · **Method:** one research agent per rubric-capping bullet (15 agents) → honest per-bullet plan → executed the movers.
**Reconciled internal merit: 78/100.** This document does NOT claim a path to 10/10 on every axis, because that path does not honestly exist in a compute-only hackathon. It states, per bullet, exactly what was done and what is structurally impossible.

## The unanimous finding

All 15 agents returned the same shape of answer: **no single bullet reaches a clean 10 in-scope.** The bullets split into three groups:

### Group 1 — STRUCTURALLY CAPPED (cannot reach 10 in a compute-only hackathon; honest ceiling stated)
| Bullet | Why 10 is impossible | Honest ceiling | Best in-scope substitute (done/available) |
|---|---|---|---|
| IMP1 Wet-lab | Cannot transfect/assay from a laptop | ~7 | Retrospective corroboration vs Norman's own fitness GI + Horlbeck 2018 GI map + Open Targets (available; not fabricated wet-lab) |
| IMP3 Different modality | Truly different modality (drug-combo viability) is off-GEO/ENA, needs network grant + new pipeline | capped | Document the one-modality-family boundary + pre-register NCI-ALMANAC as declared future work |
| IMP5 Adoption | External users acting over time cannot be manufactured pre-deadline | "adoptable", not "adopted" | Release package + CALL_FOR_SUBMISSIONS + Zenodo DOI staging (claim first-of-kind, not adopted) |
| CLA1 Novel capability | Core mechanism (recompute stat, reject on mismatch) IS statcheck (2015) | ~7–8 | Reframe as novel *application* + composition; 2×2 ablation shows reasoning-only auditors can't discriminate, receipt-gate does |
| CLA3 Reusable proven | New skill, no external track record | ~9.3–9.5 | 2nd external dogfood (GEARS claim) done; portable adapter + README shipped |
| DEP1 One-command repro | Two headline families depend on donor-blocked 16GB data not clonable | clean two-tier repro | run_all.sh offline tier verified today (SCC + benchmark reproduce from committed intermediates) |

### Group 2 — ALREADY EXECUTED by the agents (artifacts built; fold into deck/writeup)
| Bullet | Artifact built | What it does |
|---|---|---|
| IMP2 Concrete payoff | `decision_ledger.png` + receipt | Turns the meta-level trust layer into a concrete GO/WITHHOLD slot-allocation decision (CD226 GREENLIGHT vs RASGRP1 WITHHOLD) |
| IMP6 Effect size | `calibration_is_the_point.png` | Reframes modest R² as the *point*: a well-calibrated modest model beats an overconfident strong one (accuracy ⊥ coverage) |
| CLA2 Process→product | `auditor_catch_replay.html` | ~15s self-verifying replay of Claude's honesty-auditor catching Claude's OWN overclaim and forcing the retraction |
| DEP3 Citation stubs | `references.bib` fixed | Both [VERIFY] stubs resolved to real citations (Correia NeurIPS 2024 arXiv:2405.02140; Systema Nat Biotech 10.1038/s41587-025-02777-8) |
| DEP5 Narrative spine | `honesty_ratchet.png` + `README_SPINE.md` | Six-lane proposed→published provenance diagram; makes the through-line the first thing a judge sees |

### Group 3 — NEW SCIENCE EXECUTED THIS SESSION (the only bullet that could move the substantive tier)
- **IMP4 Absolute-coverage bar.** Pre-registered (sha `9b8887f3`) and tested 4 conformal methods (vanilla, SCC scalar prior, Mondrian synergy-tercile, conditional-quantile) against the absolute high-synergy coverage bar on both datasets. **Result (receipt sha `e852818b`): no conditional method DOMINATES the scalar SCC prior.** Mondrian helps on CaRPool but degrades marginal coverage on Norman; conditional-quantile ≈ vanilla. SCC satisfies the pre-registered bar on both datasets, but its POINT coverage stays below nominal (0.82 / 0.86) — the CI includes 0.90 largely because n≈130–140 is small. **Honest reading: "statistically consistent with nominal," not "demonstrably achieves nominal."** This ADDS rigor (3 alternatives pre-registered and beaten) and CONFIRMS the absolute bar is genuinely hard at this sample size — it is not a framing artifact. Does not change the 78 tier.

## What this means for the rubric (honest)

- **Depth & Execution → genuine ~10 is reachable**, but only via engineering the agents scoped (commit repo, verify two-tier repro, scar cleanup, stubs fixed ✓) — those touch git and await your approval.
- **Impact → ~9 ceiling.** IMP2/IMP6 framing done; IMP4 adds rigor; the last point needs wet-lab or adoption (impossible in-scope). The honest substitutes (retrospective corroboration, decision ledger) are built or available.
- **Claude Use → ~9–9.5 ceiling.** CLA2 replay + CLA3 external dogfood done; the gap to 10 is a novel *capability*, which the CLA1 agent showed honestly we do not have (it's statcheck's mechanism, composed).

**The honest bottom line: this pass raised the floor on three axes toward their true ceilings (~9 / ~9.5 / ~10-with-repo) and added one pre-registered rigor result (IMP4) — it did NOT manufacture a 10/10-on-everything, because that would require fabricating wet-lab, adoption, or capability the project does not have.** The remaining real points live in the Demo (30%) and in committing the repo.
