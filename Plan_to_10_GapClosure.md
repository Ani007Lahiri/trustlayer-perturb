# Plan to 10/10 — Gap-by-Gap Closure Map
**Built with Claude: Life Sciences (Researcher track). Rubric: Impact 25 / Claude Use 25 / Depth 20 / Demo 30.**
Generated from an 11-agent gap-closure fan-out (one agent per gap identified in the score analysis). Every plan below is web-grounded and honest about its ceiling. Dispatched Jul 12; all laptop-executable unless flagged.

## The honest headline
The earlier analysis said Impact was structurally capped at ~7.5 without new data. **That turned out to be wrong in one important way:** a genuine SECOND combinatorial dataset IS reachable in-scope. That single finding lifts the realistic ceiling on the whole submission.

**Realistic achievable ceilings after this plan (honest, not aspirational):**
- **Impact: 7.5 -> ~9** (was capped; the 2nd-dataset validation unblocks it). A clean 10 still needs wet-lab — out of scope.
- **Claude Use: 8.5 -> ~9** (reusable skill already built + architecture surfaced). A clean 10 needs a genuinely novel *capability*, not composition — honesty forbids claiming it.
- **Depth & Execution: 8 -> ~10** (pure packaging; no science gap). This one genuinely closes.

So: **not a guaranteed 10/10/10 — that would be dishonest to promise — but a credible ~9/9/10 on the non-demo axes**, up from 7.5/8.5/8. That is a large, real move, and it is mostly executable on this laptop.

---

## THE BREAKTHROUGH: Impact gap 1 — external validation (I1)
**Status: PARTIAL->YES in-scope. This is the highest-value single finding of the fan-out.**
The agent found and *live-verified on the GEO FTP*: **CaRPool-seq (Wessels/Satija/Sanjana), GEO GSE213957** — a combinatorial (single/double/triple-guide) perturbation screen, downloadable NOW, **no Figshare/Zenodo, no GPU**, fits in 8 GB. This is a real second combinatorial dataset to replicate the SCC synergy-miscoverage result on.
- Ruled out: Datlinger (single-guide only), Replogle (blocked + mostly singles).
- **Why it lifts Impact to ~9 regardless of outcome:** pre-registered, a POSITIVE replication = cross-modality generalization (tier-moving); a NULL = an honest modality boundary (still publishable, still rigorous). Either way the "does it generalize?" objection is answered with evidence.
- Needs from user: **nothing blocking** — GEO is on the allowlist. ~9-14 laptop-hours.
- Plus a no-new-data rigor upgrade (~3h) cross-validating benchmark/Lever-A claims across the 3 on-disk datasets.
   1. STEP 0 (pre-register FIRST, ~0.5h): Before downloading, freeze a hash-committed pre-registration mirroring the Norman SCC protocol: (i) synergy definition (deviation of observed double pseudobulk from additive singles mo
   2. STEP 1 (download, ~0.5h, network already open — GEO on allowlist): urlretrieve from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213957/suppl/ two files: GSE213957_THP1-CaRPool-seq.metadata.tsv.gz (2.2 MB; per-ce
   3. STEP 2 (parse + build perturbation table, ~1.5h): Untar, load the sparse count matrix (scipy.sparse — fits in 8 GB; THP-1 cells x ~20k genes stays well under RAM as CSR). Join cells to metadata; filter to the AML combina
   4. STEP 3 (pseudobulk + synergy, ~1.5h): Compute per-perturbation pseudobulk mean expression (log-normalized, same normalization as the Norman pipeline). For each double with both singles present, compute observed-vs-additi
   5. STEP 4 (run SCC vs vanilla, ~2h): Port the frozen Norman SCC pipeline onto the CaRPool-seq pseudobulk. Split into calibration/test folds by PAIR (leave-one-pair-out style, so no single's information leaks into its own do
   6. STEP 5 (inference, ~1.5h): Run the pre-registered permutation test (shuffle synergy labels, recompute slope-reduction null, n>=2000 perms -> p-value) and the powered pair-level LOPO interaction test with bootstrap CI. Co
   7. STEP 6 (report + honesty audit, ~1h): Write the result into the calibration benchmark as a second combinatorial row. State plainly whether CaRPool-seq replicated Norman (slope reduction, same direction, bar cleared?) or 

---

## IMPACT — remaining three gaps (each ~8-9 ceiling, all laptop-executable)
### I2 · Concrete biological payoff (meta -> patient-relevant) — ceiling 8-9
Build a trust-gate-endorsed **T1D target shortlist** from the on-disk Marson CD4+ data, with independent genetic support (Open Targets, GWAS Catalog, CD4 eQTL coloc), pre-registered, covariate-matched enrichment (COMMIT vs ABSTAIN vs matched-random), known-target positive controls. Honest framing: hypothesis-generating. 10 blocked by wet-lab (out of scope). ~8-12h.
   1. STEP 0 — Pre-register and freeze BEFORE touching genetics (reuse the project's existing SHA-256 pre-registration mechanism). Declare in one frozen spec: (a) the T1D-relevant CD4+ signature(s); (b) the perturbation effect
   2. STEP 1 — Define the disease signature from EXTERNAL sources only (kills circularity). PRIMARY = IL-2 / IL2RA / Treg-stability module: T1D is centrally an IL-2/Treg-axis disease — IL2RA(CD25), IL2 and CTLA4 are among the 
   3. STEP 2 — Score every perturbation for its predicted effect on the signature, using the SAME predictor + preprocessing as benchmark v2 (mean/linear/HistGBR ensemble; no torch/GEARS needed). For each perturbed target gene 
   4. STEP 3 — Apply the trust gate. Run the v2 calibrated conformal gate on each target's signature-shift prediction -> COMMIT / VERIFY / ABSTAIN + calibrated confidence. This partitions the ranked list into an ENDORSED (COMM
   5. STEP 4 — Assemble THREE orthogonal, external genetic-support tracks the model never saw (this is the credibility engine). (i) Open Targets Platform GraphQL (api.platform.opentargets.org/api/v4/graphql), disease=EFO_00013

### I3 · Positive capability, not just a caution flag — ceiling 10 (data-bounded magnitude)
**"Conformal experiment triage"**: quantify how many perturbation assays you can SKIP at a guaranteed error rate by trusting the gate; plus budget-constrained hit discovery. Converts "prevents a mistake" into "unlocks efficiency." Fully on-disk. The capability framing reaches 10; the headline *number* is measured honestly, not asserted. ~8-12h.
   1. STEP 0 (pre-register, ~0.5h): Before touching gate outputs, freeze a hash-committed protocol: (i) 'discovery/hit' label = perturbation whose measured pseudobulk effect size (L2 norm of mean expression delta vs control ov
   2. STEP 1 (reuse existing outputs, ~1h): Do NOT recompute predictions. Load, per (predictor x dataset) cell already in benchmark v2, the held-out predictions and the per-perturbation conformal nonconformity scores / calibra
   3. STEP 2 — PAYOFF A, selective experimentation / experiment savings (~2-3h): For each (predictor x dataset), rank perturbations by gate confidence and sweep a commit threshold. COMMIT = trust prediction, skip the assay; VE
   4. STEP 3 — PAYOFF B, budget-constrained discovery (~2-3h): Simulate a fixed wet-lab budget B (sweep B = 10%..50% of candidates). Goal = capture true top-quartile hits. Strategy = spend B assays on the gate's VERIFY band (u
   5. STEP 4 — significance + honesty (~1h): Permutation-test the calibration delta from STEP 2b (shuffle gate vs raw-confidence assignment, recompute skip-at-fixed-error, get p-value and CI) — mirror the existing perm-test ma

### I4 · Adoption signal — ceiling: adoption-READINESS high, demonstrated adoption out of scope
Release the calibration benchmark as a **first-of-kind open eval** (verified novel: scPerturBench/PerturBench/PertEval score accuracy, NOT calibration) — GitHub + submission format + validator + seeded leaderboard + Zenodo/HF DOI + open license. Day-1 claim = "adoptable + first-of-kind"; actual external uptake can't happen in 6 days. Needs: GitHub org + (free) Zenodo link, vendor-neutral names. ~9-18h.
   1. 1. Name and scope the package with NO model/vendor identifiers. Candidates: 'pertcal' / 'PerturbCal' / 'PCB (Perturbation Calibration Benchmark)'. Pick a GitHub org/user to host it (user decision) — the org name is publi
   2. 2. Freeze the benchmark SPEC (SPEC.md + versioned JSON schema). Define the eval contract precisely: a submission is a per-(dataset, perturbation) table with a point estimate AND an uncertainty object (prediction interval
   3. 3. Package the existing v2 results into a clean repo. Ship a model-agnostic library (the conformal gate + scorer) plus data/prepare_*.py scripts that FETCH from the public accessions (Norman GSE133344, Datlinger GSE92872
   4. 4. Build the CLI + validator so a stranger can run it in one command. `pertcal validate submission.parquet` (schema + leakage + coverage-of-required-perturbations checks) and `pertcal score --submission ... --dataset ...
   5. 5. Seed the baseline leaderboard (this IS the day-1 leaderboard). Run all 5 predictors x 3 datasets through the scorer, emit leaderboard.csv + a rendered static leaderboard (a README table plus a static HTML page via Git
   6. 6. Reproducibility harness. `pip install -e .`, a pinned environment.yml/requirements.txt, a `make reproduce` target that regenerates the baseline leaderboard numbers, and GitHub Actions CI that re-scores a small committ

---

## CLAUDE USE — three gaps
### C3 · Reusable artifact — ALREADY BUILT & DOGFOODED — ceiling 10 (highest leverage)
A portable Claude Skill **`research-claim-audit`** (SKILL.md + kernel.py + references) packaging the receipt + adversarial-council + honesty-auditor loops as claim-agnostic helpers. Built and tested this session — it **caught a real overclaim** ("doubles"=1.95x effect) on our own SCC claim and auto-rewrote it. Converts process-narrative into a shippable artifact another researcher can run. ~6h polish remaining.
   1. STEP 0 (DONE this session): scaffold the skill. Files saved as artifacts: SKILL.md (purpose/when-to-use/workflow/honesty-contract), kernel.py (preregister, verify_receipt, adversarial_council, honesty_audit, run_gauntlet
   2. STEP 1 — Verify the two execution paths (DONE, keep as a regression test). (a) Offline/deterministic: receipt freeze->verify->tamper-detect works with zero LLM; rule-mode council flags small-n and unbacked strength words
   3. STEP 2 — Package for another researcher: add README.md (install = drop folder into skills dir; 3-line quickstart), requirements (NONE beyond stdlib for receipts; host.llm OR any OpenAI/Anthropic client shim for council/a
   4. STEP 3 — Make it framework-portable so it is not Claude-Science-only: add a thin `llm_backend` adapter so `_llm_json` can target (a) host.llm inside Claude Science, (b) the Anthropic API directly, or (c) Claude Code as a
   5. STEP 4 — Dogfood on a claim OUTSIDE this project to prove generality: pick one published comp-bio abstract sentence (e.g. a GEARS or scGPT performance claim) plus a small evidence dict, run run_gauntlet, and include the 
   6. STEP 5 — Record the 30-second demo segment for the 3-min video: terminal split-screen showing `preregister(...)` printing a receipt_id, then `run_gauntlet(...)` on the team's SCC claim catching 'doubles->1.95x' and emitt

### C2 · Surface the surprise — ceiling 10 (pure packaging)
Make the invisible agent architecture VISIBLE: an architecture diagram, a replayable "flight-recorder" HTML of the council scoring + the auditor catching an error live. Lands inside both Demo (30%) and Claude Use (25%). ~6-10h.
   1. SINGLE MOST CONVINCING ARTIFACT (the answer to the gap's core question): a replayable, independently-verifiable 'watch Claude's honesty auditor catch OUR OWN overclaim and force the corrected number' moment — shown insid
   2. STEP 1 — Mine the transcripts for the hero catch (1.5-2.5h). Search this project's frames (host.frames with a regex over 'auditor'/'overclaim'/'coverage'/'retract'/'correct') for every honesty-auditor moment that changed
   3. STEP 2 — Draw the canonical architecture diagram, ONE figure (1-2h). Panel: Any Predictor -> Conformal Trust Gate (commit / verify / abstain) with a feedback loop showing [Adversarial Council: N agents] <-> [Honesty Audi
   4. STEP 3 — Build the replayable 'flight recorder' HTML (2-3h). Route A (authenticity, low effort): run simonw/claude-code-transcripts (Apache-2.0; pip/uvx install; converts Claude Code session JSON/JSONL into a browsable, 
   5. STEP 4 — Build the 'council scoreboard' panel (1h) — directly delivers the gap's literal ask ('watch the council score and the auditor catch an error'). A small timeline/table: 20 adversarial rounds down the y-axis, meri

### C1 · Novel capability — ceiling ~8-9 (honesty forbids a clean 10)
Verdict after prior-art search (CoVe, Self-Refine, Reflexion, Constitutional AI, Anthropic verification subagents): the *generic* self-critique+verifier pattern is NOT novel. The honest novel contribution = **receipt-gated auditing specialized to statistical claims + a pre-registration-adherence (p-hacking) gate**. Frame as composition, not clean-10 novelty. ~6-11h.
   1. 1. FRAME IT HONESTLY (30 min, highest leverage). In the writeup add a 'Related work & honest positioning' box: our auditor is an instance of the CoVe/Self-Refine/Reflexion/Constitutional-AI/CitationAgent/tool-receipt fam
   2. 2. FORMALIZE THE PATTERN AS A NAMED, REUSABLE SPEC (1-2 h). Write PATTERN.md defining 'Receipt-Gated Scientific-Claim Auditing (RGSCA)': (a) claim = a natural-language quantitative assertion in the writeup; (b) receipt =
   3. 3. BUILD claim_ledger.jsonl FROM THE EXISTING RESULTS (1-2 h). Enumerate every load-bearing number already in the project (cov 0.465->0.709 recalibration recovery, perm p=0.0005 Lever A; SCC miscoverage-slope reduction, 
   4. 4. BUILD auditor.py — THE MINIMAL BUILDABLE LOOP (2-3 h). For each claim in the ledger, spawn a FRESH Claude context (host.llm or a delegated sub-agent) that receives ONLY the receipt + claim text, never the generator's 
   5. 5. STAGE ONE HONEST LIVE-RETRACTION DEMO (1-2 h, this is the 30%-weighted demo money shot). Deliberately inject ONE overclaim into a draft sentence (e.g. change 'reduces the miscoverage slope' to 'eliminates miscoverage'

---

## DEPTH & EXECUTION — four gaps, this axis genuinely closes to ~10
### D1 · Narrative spine — ceiling 10 (delivered)
20 rounds mapped onto a **5-station falsification arc** (hypothesis->test->falsify->refine->honest ceiling); spine figure + narration one-pager already rendered. Reads cold in ~15s; the "1 confirmed / 1 partial / 1 dead" scoreboard is a STRENGTH. ~3.5h.
### D2 · Foreground receipts — ceiling 10
An **INTEGRITY.md ledger** as hero doc; `receipts/` + one-command `make verify`; the falsified PGW result promoted to a front-page asset. Grounded in model-cards/datasheets/ACM badging/pre-registration norms. ~6-9h.
### D3 · Product-grade repo — ceiling 10 (skeleton already written)
README + `run_all` regenerating every figure/receipt from public data + env pinning + LICENSE + clean-clone acceptance test. The one gap fully closeable to 10 (reproducibility is engineering, not science). Marson DOI 10.64898/2025.12.23.696273 verified. ~8-13h (−3-4h already done).
### D4 · Loose ends — ceiling 10
205 uncommitted paths, stale branch, 2 [VERIFY] citation stubs (flagged for live re-verification, NOT yet confirmed), 1 dangling receipt, stale README. A 10-step judge-proof checklist ending in a zero-loose-ends receipt. ~3-4.5h.

---

## Recommended execution order (by ROI, deadline Jul 13 9PM ET)
**Tier 1 — do first (unblocks the biggest points, all laptop):**
1. **I1 CaRPool-seq validation** (~9-14h) — the one move that lifts Impact from capped to ~9. Start the download immediately to de-risk.
2. **C3 finish `research-claim-audit` skill** (~6h) — already built; polishing it is the single highest-leverage Claude-Use move.
3. **D3 + D4 repo + cleanup** (~12-17h) — takes Depth to ~10 and is required by the open-source rule anyway.

**Tier 2 — high value, do if time:**
4. **I3 experiment-triage capability** (~8-12h) — turns the caution flag into a positive capability (Impact).
5. **D1 + D2 narrative spine + integrity ledger** (~10-13h) — makes the earned rigor legible.
6. **C2 surface the architecture** (~6-10h) — feeds directly into the demo video.

**Tier 3 — nice to have:** I2 (target shortlist), I4 (benchmark release), C1 (novelty framing paragraph).

## What a clean 10/10/10 would still require (stated honestly)
- **Impact 10:** wet-lab confirmation of a gate-endorsed target — impossible in a compute hackathon.
- **Claude Use 10:** a genuinely new Claude *capability*, not expert composition — we have composition.
- **Depth 10:** genuinely reachable — this is the one that closes.
So the honest target is **~9 / ~9 / 10** on the non-demo axes, plus a strong demo (the 30% swing). That is a legitimate top-3 profile — materially stronger than where we started, and every number above is grounded, not aspirational.
