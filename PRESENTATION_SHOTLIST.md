# Presentation Shot-List — Trust Layer (Gladstone × Anthropic, Research track)

*Persuasion-artifact upgrade of the 11-slide deck + Trust Atlas. Honesty contract holds:
the external null stays, the 74 ceiling stays, and the recalibration recovery (AUROC 0.70)
is labelled **within-assay** everywhere it appears — it is not cross-assay transfer.*

---

## The single highest-impact change
**Re-cut the ending. The deck currently climaxes on the null and the 68/74 score and stops —
the honest-but-flat "here's where it ranks" beat.** The current build's real climax already
exists in the work but isn't in the deck: *the null had a fix, and the fix is honest about its
own scope.* Add one slide after the null — **"The null wasn't the end: recalibration recovers
it (within-assay)"** — showing AUROC **0.465 → 0.703 [0.679, 0.726]** on 2,311 held-out Zhu
genes, with the depth-confound already stripped out (Spearman 0.25 net of cell count). That
converts the arc from *idea → built → it failed → honest score* (deflating) into *idea → built
→ it failed → we diagnosed WHY and recovered it, honestly, and here's the exact ceiling*
(a recovery arc, which is the shape judges remember). One slide. Everything else is polish on top.

---

## The hero number (pick ONE, put it on three slides)
**Primary: `0.47 → 0.70`** — the recovery. Frozen gate scored at chance on external truth
(AUROC 0.465); recalibrating the base predictor + conformal layer on the new atlas recovers
real signal (0.703 on held-out genes). It is the whole thesis in two numbers: *a trust layer
that knows when it's out of distribution — and tells you exactly what it takes to fix it.*
Say it in words once: **"coin-flip to real signal, and we're honest that the fix is
within-assay."**

- **Character number (use once, in the honesty beat): `3`** — inflated results caught and
  retracted by its own adversary. This is the number that makes judges *trust the 0.70*.
- Do **not** lead with `74` — it's a ceiling, not a hook. It belongs on exactly one slide
  (the honest-ranking slide), framed as "tiers, not a percentage," never repeated.

---

## Tagline options (30-second hook / title-card)
1. **"A pLDDT for perturbation biology — the number that tells a lab which predictions to trust."**
   *(clearest one-liner; anchors to the one confidence score every judge in the room already respects.)*
2. **"Every model here predicts. We built the part that says *don't act on that one* — and proves it isn't bluffing."**
   *(the differentiator framing; strongest for a room full of prediction demos.)*
3. **"The first perturbation model honest enough to fail its own external test on the front page — then recover it in the open."**
   *(honesty-as-feature; highest-risk/highest-reward, leads with the retraction character.)*
4. **"Default-deny for the wet lab: GO, WITHHOLD, or ABSTAIN — because a wrong *go* costs a year."**
   *(payoff-first; best if the judge panel skews wet-lab / translational.)*

Recommended pairing: **#1 as the title card**, **#2 as the spoken opening line.**

---

## The 30-second hook (spoken, replaces the current cold "The problem" slide)
> "A Gladstone lab can chase maybe ten of the thousands of targets a perturbation model predicts.
> Pick a wrong one and that's a year and a freezer full of reagents gone. AlphaFold didn't just
> predict structures — pLDDT told you *which parts to believe*, and that one number changed the
> field. Perturbation biology has no pLDDT. We built one — and the interesting part is what it
> does when it's wrong."

Then straight into the refusal demo. Problem → stakes → analogy → promise, in four sentences.

---

## Shot list — slide-by-slide (ADD / CUT / REORDER / KEEP)

**KEEP & sharpen**
- S1 Problem / pLDDT analogy — keep, but compress to the 30-sec hook above; cut the paragraph prose.
- S4 Commit gate (default-deny) — keep; it's the mechanism. Tighten to the 4-check chain only.
- S5 Trust Atlas / 98-gene inspector — keep, but **promote it to the demo centerpiece** (see below).
- S9 Honest ranking / 74 — keep, ONE slide, tiers-not-percentage framing. Cut the 12-number
  score trajectory from the slide face (move to appendix; it reads as anxious, not rigorous).

**ADD (in priority order)**
1. **[NEW, highest impact] "The null wasn't the end — recalibration recovers it (within-assay)."**
   Insert immediately after the external-null slide. `0.465 → 0.703`, held-out genes, depth-stripped.
   Label within-assay in the subtitle, not a footnote. This is the new emotional peak.
2. **[NEW] "What a wet-lab scientist gets"** — one slide, one sentence of payoff + one concrete
   worked call. "Feed it your 40 candidate knockouts; it hands back GO on the 6 it can defend,
   ABSTAIN on the 20 it has no evidence for, and the receipt for every call." Priced in weeks/reagents.
3. **[NEW / reframe] "How we do it today: nothing."** — one comparison slide before the mechanism.
   Current practice = trust the point prediction blindly or eyeball effect size. That's the baseline
   you beat. Judges need the "vs what people do now" contrast to score *impact*.
4. **[NEW, optional closer] "What moves it to 80: one pre-registered CRISPRi validation, ~15–20 genes."**
   Turns the ceiling into a fundable ask, not a dead end. Ends on forward motion, not a score.

**CUT / demote to appendix**
- The attainability-certificate R²≈0.01 vs 0.98 slide → demote. It's beautiful rigor but it's the
  third consecutive "here's a limit" beat; two limits land, three deflate. Keep it for Q&A.
- The 12-round score trajectory string → appendix.
- Any slide that is a result table without a sentence telling the judge what it *means*.

**REORDER — the winning sequence**
1. Hook (problem + stakes + pLDDT analogy + promise) →
2. **Refusal demo** (the magic moment, see below) →
3. "How it's done today: nothing" (baseline) →
4. Commit gate mechanism (default-deny, 4 checks) →
5. Calibration guarantee (intervals mean what they say) →
6. External null — *executed, pre-registered, negative* (the honesty gut-punch) →
7. **Recovery: 0.47 → 0.70 within-assay** (the new peak) →
8. Honesty engine: 3 self-retractions + independent audit (why you trust all the above) →
9. Honest ceiling 74 + the one fundable next step →
10. Close on the tagline.

---

## The magic moment (make refusal the demo, not a slide)
Currently the Atlas is a browse-98-genes inspector. Turn it into a **10-second scripted beat**
that judges *watch happen*:

> "Here's RASGRP1 — one of the **largest real effects** in the external atlas. A naive tool GO's it
> on effect size. Watch our gate — **WITHHOLD** — because there's no cell-type eQTL supporting it.
> The gate judges *evidence, not effect size*. That refusal is the product."

Then flip the Atlas's frozen↔recalibrated honesty toggle live: "and here's the same layer admitting
it's out of distribution, then recovering when we recalibrate." The refusal + the toggle *is* the
wow — a system visibly declining a confident-but-unsupported call is something no prediction demo
in the room can show. Pre-record it as fallback (Murphy's Law is undefeated on live demos).

---

## Honesty as a FEATURE, not a footnote (the reframe judges reward)
Move the self-retraction story out of the credibility appendix and make it a **named asset**:
the **"adversarial honesty engine."** The pitch line:

> "Every number you just saw survived an adversary we built to destroy it. It caught us
> overclaiming three times and we retracted — on the front page, not in a footnote. In a field
> where overclaiming is the norm, a system that tells you when *not* to trust it — including when
> not to trust itself — is the entire point."

Judges in a research track are trained to distrust polished claims. A team that pre-emptively
shows its own retractions and a blind pre-registered null flips that reflex: the honesty *is* the
credibility signal that makes the 0.70 believable. Pair it with the independent-audit line
(second agent re-derived the calibration surface from the raw 16.9 GB atlas, max diff 0.0).

---

## Why this wins a RESEARCH track specifically (vs a dev track)
- Dev-track judges reward a working product and a slick build. Research-track judges reward a
  **credible claim** and reward you *more* for stating limits than for hiding them — so the null +
  retractions + certified information-limit are assets here, not liabilities. Lead with them.
- Research judges score the *problem statement* heavily and remember specificity. Spend real time
  on "a wrong T1D target = a year + reagents," name the cell type and the disease, and give the
  immunologist on the panel one non-obvious-but-correct call to nod at.
- Use the platform-native rails out loud: live genetics MCP queries, auditable per-call receipts,
  reproducible artifacts. It signals you used Claude Science as intended, which the org is judging.
- Close on the fundable next step (the pre-registered CRISPRi validation). Research judges reward
  a project that knows exactly what experiment would settle it.
