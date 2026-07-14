# Merit Ceiling — Brutal Reality Check on the +20 Goal
**Question posed:** can this project honestly gain +20 merit points (from ~52 to ~72) in-scope
(hackathon, computation-only, NO wet-lab)? What is the true ceiling?

**Anchor rubric (Nobel-scaled /100):** 50–64 = competent honest prototype / GOOD hackathon ·
65–79 = solid incremental methods, well-powered · **>65 REQUIRES validated, well-powered, NOVEL
science** · 80–89 = strong novel well-powered · 100 = Nobel, field-changing.

**Current honest position:** ~52 (trajectory 50.4→54.0→51.8). Squarely "competent honest
prototype." Reaching 72 means landing MID-WAY into the 65–79 band — i.e. crossing a qualitative
gate, not accumulating points.

---

## 1. The gate at 65 is qualitative, and the project fails it on the validation leg

">65 REQUIRES validated well-powered novel science" is a three-part AND. Score each:

| Property | Status | Evidence |
|---|---|---|
| **NOVEL (method)** | ✗ | Conformal / Mondrian / weighted conformal are established. Conformal prediction is *already* published for single-cell RNA-seq annotation with statistical guarantees (Bioinformatics 2025); calibrated uncertainty for perturbation prediction is a crowded 2025 area (PRESCRIBE, CILANTRO-SL). The integration ("calibrated trust + default-deny gate for target triage") is a fresh COMBINATION → **application-novel, not method-novel.** |
| **WELL-POWERED** | ✓ *for calibration only* | The ONE well-powered honest result is donor-blocked LODO conformal coverage (0.79/0.90/0.95). Every result that would constitute *science* (does the trust score find real targets? does the gate catch real failures?) is UNDERPOWERED: recovery n=8 (p=0.72), Check A n=24/4-pos, clinical 7/98 all-success, ClinVar 9 genes. |
| **VALIDATED (biology)** | ✗ | Calibration is validated AS calibration. The biological claim — this identifies good/bad T1D targets — is NOT. Every external validation returned null or was retracted as a surrogate: drug-target AUROC 0.829 was a surrogate (real gate MCC +0.119); mouse-KO AUROC 0.687 scored a free-text column (real gate MCC −0.15, p=0.68). |

**The straddle.** The crown result (conformal coverage) is well-powered and validated-as-
calibration but is *incremental methods on borrowed foundations*. The scientific payload (target
discovery) — the thing that would make this "novel science" — is neither validated nor powered.
Methods-quality pulls toward the low end of 65–79; science-payload-quality pins it in 50–64.

---

## 2. Lever budget (honest maxima, from ~52)

| Lever | Honest max | In-scope? |
|---|---:|:--:|
| Novelty framing ("pLDDT for perturbation biology") | +2 | ✓ |
| **Prior-art benchmark WIN vs a NAMED SOTA uncertainty method** | **+5** | ✓ |
| Powered external validation WITH FAILURES | +3 (would be +8..12 if a benchmark with real negatives existed) | ✗ |
| Methodological extension (RCPS/LTT risk control, group-conditional/adaptive coverage) | +4 | ✓ |
| Scale (universe-scale run, robustness sweeps) | +2 | ✓ |
| Reproducibility (hash-stamped receipts, frozen splits) | +1 (half already priced in) | ✓ |

**Naive additive sum = 69.** This number is a trap and must be refused. The 65 line is a
qualitative gate, not a point total. Every in-scope lever improves the SAME story — "solid
incremental methods, well-powered" — so they COMPRESS (diminishing returns), not add. Stacking
benchmark-win + methods-extension + scale yields a better *methods* contribution, not "validated
novel *science*." The one lever that could cross the gate (powered external outcomes WITH
FAILURES) is structurally out-of-scope: it is exactly what has failed repeatedly.

**Realistic honest ceiling, computation-only in-scope: ~60–64.** Best-case reach ~65 IF a clean,
powered head-to-head benchmark win AND a real methods extension both land and are presented as a
coherent methods contribution the community would credit. **72 is not reachable.**

---

## 3. What SPECIFICALLY unlocks 65+ / 72 — and whether it is in-scope

The blocker is the **validation leg**. Concretely, one of:

1. **A powered external outcome benchmark that CONTAINS FAILURES** — the committed gate must be
   *falsifiable* and shown right against real downstream success/failure with enough events per
   class for power (≥~20–30/class, not 7/98-all-success). **Out-of-scope**: no such dataset with
   real negatives exists in-hand; every attempt returned all-success/underpowered.
2. **A wet-lab / experimental collaboration** testing the withheld vs committed calls (the
   Marson-lab connection is the natural venue). **Out-of-scope** for a hackathon — it is
   post-hackathon by construction.
3. **A published-quality, powered head-to-head WIN on a recognized benchmark** (Systema /
   PerturBench / Arc Virtual Cell Challenge-style eval) where the calibrated gate beats a NAMED
   SOTA UQ method (e.g. GEARS MC-dropout, documented as poorly calibrated) on a rank/decision or
   coverage metric. **Partially in-scope** — this is the ONLY merit lever the hackathon can
   actually pull, and it lifts toward ~65, not 72.

Only option 3 is reachable now, and it does not clear 72.

---

## 4. Honest verdict

**Merit axis:** **+20 to 72 is NOT achievable computation-only in-scope.** The honest ceiling is
**~60–64** (reach ~65 only in the single best case: a real, named, *powered* benchmark win plus a
genuine methods extension, both reported without surrogates). The +20 is blocked by the
validation leg, which by construction requires falsifiable external outcomes-with-failures or
wet-lab — neither in-scope. Telling the user 72-on-merit is reachable would be exactly the
overclaim this project exists to refuse.

**Placement axis (demo / story / judging) — a different axis:** here a **+20-equivalent jump in
competitiveness is realistically achievable.** The differentiators — the honesty contract, the
self-*retraction* of two surrogate wins, the default-deny gate, the donor-blocked conformal
headline — land at a moment when the field is *publicly alarmed* about inflated perturbation
metrics (Systema: systematic variation correlates 0.91–0.95 with the standard scores; benchmarks
where the mean baseline beats scGPT/scFoundation). In a "Built with Claude: Life Sciences"
research track judged on rigor + narrative + model use, a project whose thesis IS calibrated
honesty — and which demonstrably caught its own surrogates — is an unusually credible, timely
story. That is worth the difference between mid-pack and top-3. But it is PLACEMENT, categorically
distinct from MERIT; the user must not conflate them.

**Single highest-yield move (does double duty):** land ONE real, named, powered head-to-head
benchmark win on calibrated uncertainty — the committed gate vs a documented-poorly-calibrated
SOTA (GEARS-style MC-dropout) on a recognized dataset, scored on a decision/coverage metric,
reported with the same retraction discipline. It is the most merit-productive in-scope lever
(~+5 toward the ~65 reach) AND the natural demo centerpiece (large placement gain). Everything
else (methods extension, universe-scale robustness) is supporting depth, not a gate-crosser.
