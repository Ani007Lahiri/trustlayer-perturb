# Pitch Playbook — judge psychology, Q&A prep, and the memorable devices
*Research/Lab track (Gladstone × Anthropic "Built with Claude: Life Sciences"). Companion to the deck, the workflow graphics, and the Live Abstention Demo.*

## The one repeated line (say it 3×: open, at the recovery beat, close)
> **"It knows what it doesn't know."**

Hero number, shown as motion not text: **0.47 → 0.71** (flip the toggle live).

---
## Memorable devices (pick 1–2, use them consistently)
- **The tile** — one visual primitive (green = commit, amber = verify, grey hatched = abstain), the same object on every slide, in both widgets, in the video, on the handout. Borrowed straight from AlphaFold's pLDDT color language.
- **The circuit breaker** — "the gate is a circuit breaker for wrong biology: when confidence runs where coverage can't support it, it trips *on purpose*." A breaker tripping is universally understood as protection, not failure.
- **The hurricane cone** — "a point prediction with no calibrated uncertainty is a hurricane forecast with the cone deleted." Lands in one sentence.
- **The second opinion** — abstention = a radiologist saying "I'm not sure, get a second look" rather than guessing. Reframes refusal as careful practice.

---
## Judge psychology (research track ≠ dev track)
- Research/lab judges score the **question, rigor, novelty, and whether it would survive scrutiny in a real lab** — not UI polish. Do not over-invest in shine at the expense of defensibility. The scientific claim is the product.
- **What kills a strong project in Q&A:** overclaiming, and evasiveness under pressure. Being unable to name your own limitation reads as not having thought deeply. **Name the limitation before they do.**
- **The null is a feature, said out loud:** "On that out-of-distribution case the system *abstained* instead of committing — it refused to hand you a confident wrong answer. The frozen model would have returned high confidence there. A system that's confident on that null is the dangerous one."
- **Don't hide the 76.8** — present it yourself as a *self-assessed computational ceiling* (not a wet-lab claim). Volunteering it converts reflexive skepticism into a credibility signal.
- **Trust is audience-relative** — pitch the gate as letting the *scientist set their own commit threshold*. That turns abstention from a limitation into a control the expert wants.

---
## The 5 hardest questions (rehearse the exact metric names)
**Q1 — Novelty: "Isn't this just a confidence score bolted onto someone else's model?"**
> The contribution isn't a new foundation model — it's the *commit-gate abstraction + external-null validation protocol* for perturbation biology, a domain with no per-prediction confidence standard. Structural biology got pLDDT; perturbation atlases have nothing equivalent. It's a layer, not a new model — and I'll say that plainly.

**Q2 — Overfitting: "0.47→0.71 recovery of what, and did you tune recalibration on your own test set?"**
> Held-out effect-separation AUROC. Recalibration is fit on the *calibration split*; the external null is a separate out-of-distribution dataset **not used for fitting** — that's the entire point of reporting it. Anything self-assessed, I flag as self-assessed.

**Q3 — Self-scoring: "Why trust a 76.8 you gave yourself?"**
> You shouldn't trust it as an external result, and we don't present it as one. It's a self-assessment against a deliberately harsh Nobel-anchored rubric — a *computational ceiling*, the honest limit of what compute alone can claim without wet-lab confirmation. The validation protocol is the path to a number we didn't assign ourselves.

**Q4 — The null looks like failure: "Doesn't a negative external result mean it doesn't generalize?"**
> It means the gate did its job — it abstained on an OOD case instead of committing. The frozen model would have been confidently wrong there. A null the gate *catches* is the system working.

**Q5 — Transfer: "Does any of this move beyond CD4+ T cells / T1D?"**
> The gate is model-agnostic — it wraps predictions, not a biology. A new context needs two things: a labeled calibration set and one external null to test coverage. I won't claim it transfers for free; I'll claim the *recipe* transfers, and name those two as the cost.

**Bonus — "What does a biologist actually DO with a trust score?"**
> Set a commit threshold and triage: only perturbations above it go to the bench; abstained ones get a confirmatory experiment first. It's a bench-prioritization tool, not an oracle.

---
## Format assets built (this session)
1. **Live Abstention Demo** (`Live_Abstention_Demo.html`) — the highest-leverage addition. Three real genes; the gate commits (CD226), asks to verify (PRKCQ), and refuses (RASGRP1) — then the toggle shows the same RASGRP1 case going false-confident-green when frozen → correctly-grey-abstain when recalibrated. Keys 1/2/3, `t` to toggle. Zero live compute — cannot stall on stage.
2. **Experiment Workflows** (`Experiment_Workflows.html` + 5 PNGs) — one 5-second schematic per experiment, deck-ready.
3. **The deck** (`Trust_Layer_Presentation.html`) — 13-slide recovery-arc narrative with both widgets embedded.
4. **60–90s demo video** (to record): screen-record the exact 3-click + toggle sequence; 6–12 frames, ≤8s each, on-screen text echoing "0.71 / ABSTAIN / 0.47→0.71". Ships as the Devpost asset AND the live-demo backup.

*All numbers are real and hash-frozen; the null and the 76.8 stay labeled wherever they appear.*
