# The Honest Path to a Real +20: A Validation Proposal

**Current merit ~68/100** (Nobel-anchored; "solid incremental methods, well-powered"). The ceiling
is ~72 because the system is validated **only computationally** — it has never been tested against
experimental truth it did not build. The 80-89 band ("strong novel well-powered science") requires
exactly that. This document scopes the three real routes there, cheapest first. All costs are
ESTIMATES (not vendor quotes); all datasets/methods are cited.

---

## Route A (DO FIRST — free, weeks, biggest single jump): retrospective validation on the Zhu 2025 CD4+ T-cell Perturb-seq atlas

**The opportunity:** the Marson/Pritchard labs published **genome-scale Perturb-seq in PRIMARY
human CD4+ T cells** — the *exact* cell type the gate makes calls in — ~22M cells across resting +
stimulated conditions, ~12,748 genes targeted, MIT-licensed on the CZI Virtual Cells Platform
(Zhu, Dann, … Pritchard, Marson, bioRxiv 2025.12.23.696273 / doi 10.64898/2025.12.23.696273).

**Why it matters:** this is a **true external molecular benchmark the project did not build** — the
single thing missing from every prior merit round. Freeze the committed gate, then run it *blind*
on the Zhu perturbations and ask: are the gate's GO calls the perturbations that actually produce
the predicted transcriptional shift, and are its WITHHOLD/ABSTAIN calls genuinely the ones where
the effect is absent or context-dependent? Genuine negatives (no-effect perturbations) are present
in the atlas — those are the true-negative labels the gate must respect, which public trial data
never had (only n=1 clean failure).

- **Cost/timeline:** free (public), ~2-4 weeks of compute + analysis. In-scope for THIS project.
- **Honest merit contribution (estimate):** the **largest single jump** — converts the system from
  internally-validated to externally-validated in matched context. Plausibly lifts merit from ~68
  into the **low-to-mid 70s**, and if the gate's calls hold up at power, potentially to the 75-80
  boundary. This is the one route you could start today.
- **Risk, stated honestly:** the gate may NOT be well-calibrated out-of-distribution (prior
  cross-assay work showed frozen calibration fails on Norman K562 and needed recalibration). If the
  blind test fails, that is itself a real, publishable finding — reported, not buried.

---

## Route B (weeks-to-months, ~low-mid five figures): the minimum-viable wet-lab validation

**Design (pre-registered):** MIN-VIABLE VERSION — focused validation of ~15-20 pre-registered gate calls (ARRAYED):

SCOPE: ~15 highest-stakes gate calls (must include the committed anchors: CD226->GO, RASGRP1->ABSTAIN, PRKCQ->WITHHOLD, plus a balanced spread of GO/WITHHOLD/ABSTAIN) + 5 controls = ~20 genes. Arrayed dCas9-KRAB CRISPRi or Cas9 RNP electroporation, single-donor-resolution miniaturized format (Sci Rep 2025 s41598-025-13532-z). 6 donors, rest+stim.

READOUTS: bulk RNA-seq (molecular concordance) + flow (%FOXP3+/%CD25+ and IL-2/IFN-g/IL-17). ~20 genes x 6 donors x 2 contexts ~= 240 bulk RNA-seq samples.

COST (ESTIMATE): bulk RNA-seq ~$80/sample x 240 ~= $19k; RNP/guides + electroporation kits ~$5-10k; donor

**Power:** POWER SKETCH (two units of analysis — cells for molecular, donors for functional):

(A) MOLECULAR (per-perturbation, single-cell): detecting the gate's predicted transcriptional shift.
- 100 informative cells/perturbation -> ~25% fold-changes detectable; 300 cells -> ~10-15% (Dixit/Adamson 2016; PerturbPlan PMC13228452; 10x guidance). GATE CALLS OFTEN HINGE ON MODEST EFFECTS, so budget 300 informative cells/perturbation. With 3 guides x ~100 cells/guide that is achievable at the ~200k-recovered-

**Throughput/feasibility:** REALISTIC THROUGHPUT — pooled Perturb-seq / CRISPRi in PRIMARY human CD4+ T cells (2024-2025 state of the art):

FEASIBILITY IS ESTABLISHED. The Marson lab has now run genome-scale Perturb-seq directly in primary human CD4+ T cells (Zhu, Dann, ... Pritchard, Marson; bioRxiv 2025.12.23.696273 / SSRN 6137047, Dec 2025), targeting ~all expressed genes across rest/stim contexts with hundreds of thousands to millions of single-cell transcriptomes. This is the direct methodological precedent for the p

**Cost (ESTIMATES):** COST (all figures ESTIMATES; not vendor quotes):

Unit anchors (ESTIMATES — representative 2024-25 10x/NovaSeq consumable costs, NOT verified against PerturbPlan body text): ~$0.08-0.10 per recovered cell for 3' capture + library prep, and ~$0.35-0.40 per million reads on a NovaSeq X 25B flow cell. At 25k reads/cell that is ~$0.009/cell sequencing, so ALL-IN ~$0.09-0.10 per recovered cell for capt


- **Honest merit contribution:** a prospective wet-lab confirmation of the gate's calls is the
  cleanest route into the **80-89 band** — a novel, validated result. This is what a "+20" actually
  requires.
- **Feasibility caveat:** needs a wet-lab collaborator with primary human T-cell CRISPR capability
  (e.g. a Marson-style lab). Out of scope for a hackathon; realistic as a post-hackathon
  collaboration.

---

## Route C (the full prospective study): pre-registered pooled Perturb-seq calibration study

FULL VERSION — pre-registered pooled Perturb-seq calibration study:

SCOPE: full ~35-gene test panel (15 GO / 10 WITHHOLD / 10 ABSTAIN) + 5 positive controls + 5 non-targeting, 3-4 guides each, POOLED Perturb-seq in primary human CD4+ T cells across rest+stim, 3 donors for the molecular arm; PLUS a paired arrayed functional arm on 6-8 donors for the flow endpoints (%FOXP3/%cytokine). ~200k recovered single cells (~81k informative perturbed cells) at 300 cells/perturbation.

DESIGN: mirrors the Marson genome-scale Perturb-seq workflow (bioRxiv 2025.12.23.696273) at ~1/200th the gene count. Molecular concordance from scRNA-seq + functional confirmation from the arrayed arm feed a pre-registere


The definitive version — a pre-registered, powered, prospective test of the whole calibration +
gate on fresh primary-cell data. This is landmark-track (potential 80+, publishable), and the
natural centerpiece of a grant or a lab collaboration.

---

## The licensed-failure-data route (for completeness — NOT recommended for an academic team)

LICENSED / RESTRICTED-ACCESS SOURCES OF ADJUDICATED DRUG-TARGET FAILURES (efficacy vs safety vs business), at scale:

1) Citeline (Norstella) suite — the canonical commercial source.
   • Pharmaprojects: 90,000+ drug profiles (~20,000 active); tracks discontinuations and discontinuation reasons at drug/program level (40+ yrs curation). [FACT]
   • Trialtrove: ~400,000 trials; trial-level status, outcomes, screen-failure metrics, endpoints. [FACT]

Access: COMMERCIAL (Citeline/Clarivate/IQVIA/Evaluate/GlobalData): pricing is quote-based and NOT publicly listed; the only public anchor is that these subscriptions run 'in the tens of thousands of dollars per year' per seat/module (IntuitionLabs pipeline review). [ESTIMATE — vendors do not publish list pr


Commercial trial-outcome corpora (Citeline/Pharmaprojects, Informa) contain adjudicated
efficacy-vs-safety failures at scale — the "powered external failure set" the project lacks — but
subscriptions run tens of thousands of dollars/year and are quote-based. Not feasible for an
academic team; listed only so the option is on the record.

---

## Recommendation

1. **Start with Route A now.** It is free, in-scope, uses the exact matched cell type, and is the
   single largest honest merit move available without a wet lab. I can execute it — freeze the gate,
   pull the Zhu atlas, run the blind validation, report whatever it shows.
2. **Route B/C are the true +20.** They cross into validated-novel-science, but require a wet-lab
   collaborator and belong in a post-hackathon proposal or grant. This document is that proposal.

**The honest bottom line:** a literal +20 lives in Routes B/C (experimental validation). But Route A
is a real, free, in-scope step that could move the number meaningfully *and* is exactly the
external-truth test the merit councils have said all along is the missing piece. It is the highest
expected value action on the board.

_Citations grounded in: Zhu et al. 2025 (bioRxiv 2025.12.23.696273); Schmidt, Marson et al. 2022
(CRISPRa/i screens in primary T cells); 10x/Replogle Perturb-seq scale references; Citeline/Informa
(commercial). All cost figures are estimates, not vendor quotes._
