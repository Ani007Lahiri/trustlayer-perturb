# Positioning & Impact — how to raise the rating's QUALITY and land harder
*Synthesis of three web-grounded landscape studies: (1) past winning bio-AI hackathon projects, (2) the science×AI company landscape, (3) honest rating-quality levers. All sources named at bottom. Merit unchanged at 76.8 — this is positioning, not a new result.*

## The one honest through-line
Across all three studies the verdict is the same: **this project's substance is already at the winner line for its arena — the gap is packaging and positioning, not science.** The closest analog event (Arc's Virtual Cell Challenge) built its whole credibility on blind held-out data and on *openly reporting that models don't beat naive baselines*. This project's honest NULL (0.465) + pre-registered recovery (0.709) + self-harsh 76.8 is exactly that culture. In this arena, honesty is a **scarce, judge-visible signal** — most projects overclaim and get caught in Q&A. Do not sand it down.

---
## 1. What actually wins here (hackathon-winners study)
Scored against 10 winning elements distilled from 9 real winning projects (Arc VCC 2025, Evolved 2024's GenPlasmid/EvoCapsid, Recursion RxRx, Bio×AI 2025, Anthropic Build Day, AIxBio 2026, iGEM SUPER):
- **HAVE (4/10):** memorable quantified result (0.47→0.71); rigor/honest held-out eval (the strongest asset); uses the host's Gladstone data + Claude Science; deployability/continuity (a gate is infrastructure, not a one-off).
- **PARTIAL, cheaply addable (6/10):** live driveable demo · named beneficiary · reusable artifact · wet-lab path · novelty framing · balanced pitch.
- **MISSING (0/10):** nothing is absent. The gap is presentation.

**What separates a research-track winner from strong-but-unplaced:** a defensible result on data you didn't pick (✓ have — LODO); honesty about limits read as strength (✓ have); one crisp repeatable takeaway (✓ have the sentence); a concrete "who acts on it Monday" (⚠ must dramatize); scope honesty + a real next step (✓ the wet-lab protocol + detector).

**Top-3 placement moves (ranked by return/effort):**
1. **Ship a live, driveable gate demo** — *done this session:* `Live_Abstention_Demo.html` (type/pick a gene → GO/VERIFY/ABSTAIN, the money-shot toggle flips RASGRP1 from false-confident-green to correctly-grey). This was the single highest-leverage move and it's built.
2. **Package the "wraps any predictor" claim into a real left-behind artifact** — a pip-installable `trust-gate` wrapper + one-command LODO eval reproducing 0.465→0.709, with a README showing a stranger wrapping their own model in ~10 lines. GenPlasmid placed on exactly this (a dataset+benchmark left behind).
3. **Dramatize the beneficiary + state novelty honestly** — "A T1D lab has 200 candidate CD4+ knockouts, budget for ~40: the gate says COMMIT to 12, VERIFY 30, ABSTAIN the rest — here's the one plate that tests it." Pair with a one-slide honesty statement: *framing* novelty (calibrated default-deny gate for perturbation biology) + *borrowed* machinery (split-conformal) + self-scored 76.8. (Captured in the PITCH_PLAYBOOK.)

---
## 2. Where it sits in the field (science×AI landscape study)
- **The predictors it would wrap:** Arc Institute STATE + Virtual Cell Challenge, CZI Virtual Cells Platform, Recursion (merged w/ Exscientia), Noetik OCTO, Insitro, Cellarity, Vevo/Tahoe-100M. Enormous capital on *bigger predictors and bigger data*.
- **Who owns a confidence layer for perturbation prediction? Essentially no one.** Conformal prediction is mature and model-agnostic, and it's used for ADMET / drug-target tasks — but **not** packaged as a deploy/abstain gate for perturbation transcriptomics. The field explicitly names the gap: current perturbation benchmarks "treat predictions as point estimates"; calibration is at most an emerging leaderboard column, never a per-gene decision layer.
- **White-space verdict: real but narrow, honestly bounded.** Structural biology got pLDDT (per-residue 0–100, four adopted bands) and it was as load-bearing for adoption as accuracy itself — AlphaFold won the 2024 Nobel. Perturbation biology has **no equivalent**. The defensible claim is *integration/standardization + honest calibration*, **not** "invented uncertainty for biology."
- **Quantified pain (defensible numbers for the pitch):** ~50%+ of preclinical research is irreproducible (~$28B/yr US; Freedman 2015; Amgen replicated 6/53 landmark studies); ~90% clinical trial failure; the **FDA's Jan 2025 draft framework** for validating AI models in drug submissions means regulators are now asking exactly the "can we trust this output" question a gate answers.
- **Sharpest impact framing:** *"AlphaFold didn't win biology's trust because it was accurate — it won because every prediction shipped with a pLDDT telling you when to believe it. Perturbation biology is having its AlphaFold moment on accuracy but has no pLDDT. We built the missing layer — a model-agnostic GO/VERIFY/ABSTAIN gate that turns any virtual-cell predictor into a triage instrument, honest enough to report its own null and recover from it. The bottleneck has moved from making predictions to knowing which ones deserve a pipette."*

---
## 3. Honest levers on the rating (rating-quality study)
Of 11 available moves, **7 are impact/legibility-only, and only 3 can nudge merit without wet-lab** (a perturbation-native nonconformity score; a genuine shift-regime coverage guarantee; a *validated* — not proposed — concept-shift detector). Only the first two can plausibly reach the 80 novelty wall; the 90 wall needs the wet-lab protocol run. **Do not mistake the many easy impact wins for merit movement.**

**The single highest-leverage honest move — the minimal public benchmark.** Extend the *identical* protocol from one context to ≥2–3 public scPerturb/OP3 datasets × ≥2–3 wrapped predictors (a linear baseline is mandatory, per Ahlmann-Eltze; plus e.g. GEARS and one scFM), reporting **risk–coverage + calibration curves as the primary axis** — the thing no existing perturbation benchmark (OP3, PerturBench, PertEval-scFM, VCBench, VCC) reports — with standard field metrics (E-distance, DES/PDS) alongside for comparability, released as an open eval others can submit to.

Why this one: it (a) **earns the "field infrastructure" reframe** that's currently unearned, (b) hits the exact white space the literature flags, (c) inherits the proven adoption formula (standardized task + held-out data + named metric + open leaderboard — how scPerturb, OP3, and Arc VCC became infrastructure), and (d) turns any *future* genuine-merit move (novel nonconformity score; wet-lab validation) into a submission to your own benchmark rather than an orphan result.

**Honest classification:** impact-first, merit-adjacent (**~+1–3 pts, to the top of the 65–79 band, ~78–79**), **novelty-enabling — not a shortcut through the 80 wall.** A benchmark assembled from known primitives is not, by itself, method novelty. Said plainly so the number stays honest.

**The benchmark-reframe verdict:** "the first calibration/selective-prediction benchmark for perturbation biology" is **earned in principle, not yet in fact.** Claiming it *today* (one context, few predictors, single recovery result) would be the exact overclaim this project's ethos forbids. The honest present-tense voice is *"a single-context trust layer that could seed the first calibration benchmark"* — claim it outright only once the multi-dataset, multi-model public eval exists.

---
## Bottom line — the prioritized honest path
1. **Now (built/cheap, pure placement):** the live abstention demo (done), the dramatized beneficiary + honest-novelty slide (in the playbook), keep the null + 76.8 central. These move *placement*, not merit.
2. **Half-day (impact + a real +1–2 merit):** package the pip-installable `trust-gate` wrapper + one-command LODO repro.
3. **The one strategic move (impact-first, to ~78–79, novelty-enabling):** the minimal public multi-dataset/multi-model calibration benchmark — adopt the "toward the first…" framing now, claim it outright only when built.
4. **The honest ceiling (gated, out of scope):** a perturbation-native nonconformity method (the 80 wall) and the wet-lab validation of the gate (the 90 wall). Named, not faked.

*Sources (web-grounded, named): Arc Virtual Cell Challenge 2025 (Cell S0092-8674(25)00675-0; arcinstitute.org wrap-up) · Evolved 2024 / Lux Capital (GenPlasmid/OpenPlasmid, EvoCapsid, Brainstorm) · Recursion RxRx / MorphoLogic · Bio×AI Hackathon 2025 (bio.xyz) · Anthropic Claude Build Day (Tekton, Sim Francisco) · AIxBio 2026 (Apart Research) · Devpost/TAIKAI judging criteria · Arc STATE, CZI Virtual Cells, Noetik OCTO, Recursion, Vevo/Tahoe-100M · Ahlmann-Eltze, Huber & Anders 2025 Nat Methods 22:1657 (linear-baseline null) · scPerturb (Peidli 2024 Nat Methods 21:531) · Open Problems/OP3 (Luecken 2025 Nat Biotech) · conformal lineage (Angelopoulos & Bates; Tibshirani 2019) · Freedman 2015 PLoS Biol (~$28B/yr); Begley & Ellis 2012 Nature (6/53) · FDA draft AI guidance Jan 2025 · AlphaFold pLDDT / 2024 Nobel. Some benchmark-critique items are 2026 preprints (VCBench arXiv 2604.27646; virtual-cell-usefulness bioRxiv 2026) — flagged as most-recent gap-evidence.*
