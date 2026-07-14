# Integrity Ledger

> **Why this document exists.** Every predictive claim in this project was
> pre-registered — its prediction and analysis plan hashed with SHA-256 and
> frozen *before* the confirmatory test was run. This ledger is the public
> receipt. It records what we predicted, when we froze it, what actually
> happened, and whether the pre-registered bar was cleared. Nulls and
> not-cleared rows are shown with the same prominence as wins. A clean
> reported negative is a credibility signal, not a blemish — see
> [`docs/negative_result_PGW.md`](docs/negative_result_PGW.md).

## The Ledger

| Claim ID | Pre-registered prediction (verbatim) | Prereg SHA-256 (short) | Outcome | Cleared? | Receipt |
|---|---|---|---|---|---|
| **A — Cross-assay transfer** | "A conformal predictor calibrated on one assay will under-cover when applied to a second assay; recalibration on target-assay data will restore coverage toward the nominal level." | `77bcaea9` | Frozen coverage collapses to **0.465** under cross-assay transfer, recovers to **0.709** after recalibration. Permutation **p = 0.0005**. | ✅ **Cleared** — directional prediction confirmed; recalibration recovers coverage as predicted. | Permutation test (10,000 draws), pre-registered direction and estimator. |
| **SCC — Synergy-conditioned coverage** | "Vanilla conformal prediction will under-cover high-synergy genetic double perturbations; conditioning the prediction interval on a synergy prior will reduce the miscoverage-vs-synergy slope." | `7c97f9d9` | On Norman doubles (n = 131): vanilla **under-covers** high-synergy doubles; the synergy prior **reduces the miscoverage slope**. Permutation **p = 0.005**. | ⚠️ **Partially cleared** — pre-registered *direction* confirmed; frozen *composite bar* **not** fully cleared (vanilla top-50% high-syn coverage 0.864, CI [0.773, 0.939], includes 0.90). Reported as a **within-tier rigor gain, not a tier change.** | Permutation test; single dataset (Norman). |
| **SCC — Powered pair-level re-analysis** | "In a leave-one-pair-out (LOPO) design powered at the pair level, the synergy × coverage interaction coefficient will be positive and its confidence interval will exclude zero." | `77bcaea9` | Powered pair-level LOPO interaction **+1.675**, permutation **p = 0.0005**, cluster-bootstrap **95% CI [0.97, 2.99]** excludes 0. | ✅ **Cleared** — coefficient positive, CI excludes zero. | Cluster bootstrap (pair-clustered), permutation test; pre-registered LOPO estimator. |
| **SCC-EXT — Second-dataset replication (CaRPool-seq)** | "The pre-registered synergy × coverage interaction (β > 0, permutation p < 0.05) will reproduce on an independent combinatorial dataset with a different cell type and perturbation modality." | `9d48478e` (canonical, pre-download) | On CaRPool-seq (GSE213957; THP-1, Cas13 **knockdown** vs Norman's K562 CRISPRa): interaction **β = +2.20**, permutation **p = 0.0005**, bootstrap **95% CI [1.36, 4.03]** excludes 0. Leakage-safe prior ρ = 0.62 (canonical) / 0.79 (corroboration). **Two independent agent-session reconstructions converge** (canonical HVG n=142; corroboration all-gene n=158). | ⚠️ **Directional replication cleared; composite bar NOT cleared on either dataset.** The relative interaction slope reproduces robustly across a modality flip and across both reconstructions; the absolute-coverage composite bar clears on neither. **Within-tier credibility gain, not a tier change to 80.** | Pre-registered pre-download; independent blinded critique found no blocking error; sensitivity checks (below) rebut confounds. |
| **PGW — Graph-distance conformal** | "Incorporating STRING protein–protein interaction graph distance into the conformity score will improve calibration beyond an expression-only feature." | `eab893d1` | **FALSIFIED.** STRING graph adds **no** calibration information beyond one expression feature; graph-augmented score is **worse** (Wilcoxon **p = 0.009**). | ❌ **Not cleared — pre-registered NULL.** Reported as a clean honest negative. | Wilcoxon signed-rank; hash frozen *before* test. See [`docs/negative_result_PGW.md`](docs/negative_result_PGW.md). |
| **Benchmark v2 — RandomForest calibration break** | "Across a multi-predictor / multi-dataset benchmark, calibration quality will not be uniform; at least one flexible learner will break calibration on small-n datasets, and a calibration gate will be required to certify a predictor before use." | Benchmark v2 protocol (5 predictors × 3 datasets) | RandomForest calibration **breaks on small-n**: Norman coverage **0.394**, Datlinger **0.533** at nominal **0.80**. Relational predictor strongest real learner on T-cell data (**R² = 0.385, Spearman 0.607**); the gate calibrates it. | ✅ **Cleared** — non-uniform calibration confirmed; break localized to small-n; gate certifies the relational predictor. | 5 predictors × 3 public datasets; coverage at nominal 0.80. |

**Ledger rows: 6.**

---

## SCC-EXT sensitivity checks (rebut the magnitude-confound objection)

An independent blinded critique demanded these before crediting SCC-EXT. All B = 2000 permutation.

| Check | Norman β | CaRPool β | Verdict |
|---|---|---|---|
| Scale-free `synergy_frac` (removes magnitude scale) | +1.68 (p≤0.001) | +1.01 (p≤0.001) | Effect **survives** |
| `s_hat` residualized on `true_mag` (removes magnitude confound) | +2.85 | +3.95 | Effect **strengthens** → genuinely synergy-conditioned |
| `calib_frac = 1.0` (removes inert 60/40 split) | +2.15 | +2.31 | Effect **survives** |

**O3 (magnitude confound) — REBUTTED:** the interaction is not a magnitude artifact; the method name is earned.
**O2 (retained caveat):** the interaction sign is *partly mechanical* once `|residual|` correlates with the conditioner. Cross-modality replication demonstrates **statistical robustness**, not proven conserved synergy *biology*. This caveat stays disclosed.

---

## Declared Limitations & Non-Claims

We state these explicitly so they cannot be read into the results by implication.

- **SCC replicates directionally on two datasets, but the composite bar clears on neither.** The synergy × coverage interaction reproduces across a modality flip (K562 CRISPRa → THP-1 Cas13 knockdown) and across two independent reconstructions. The **absolute-coverage composite bar is not cleared on either dataset.** We claim a **cross-modality directional replication / within-tier credibility gain**, **not** achievement of nominal coverage and **not** a tier change to 80.
- **The instrument was rebuilt concurrently with the new-data analysis.** This yields ~1.5 independent validations, not a clean 2. Disclosed as a limitation on the strength of the replication.
- **PGW is a null.** Pre-registered and **falsified**. No positive claim. Retained because a frozen, reported null is evidence about the method space.
- **Merit score is internal, not the hackathon score.** The **78/100** merit figure comes from a deliberately-harsh, Nobel-anchored **internal** tracker (trajectory 50.4 → 78.0). It is **not** the hackathon scoring number and must not be presented as such. The 78 is the **independent-critique-adjudicated** number; it corrects a self-council 79.5/81 downward — integrity keeps the more conservative, better-supported figure.
- **No wet-lab validation.** Every result is computational, on public datasets. No experimental / in-vitro / in-vivo claim.

---

## Errors we caught and corrected

Integrity is not the absence of errors — it is catching them and leaving the correction on the record.

- **Correlation-label error.** An earlier analysis mislabeled a correlation quantity in the SCC pipeline. The honesty-auditor flagged the mismatch; the label was corrected and downstream numbers re-derived before freezing.
- **Overclaim on SCC coverage.** A draft stated the synergy prior "achieved nominal coverage" on high-synergy doubles. The auditor cross-checked the frozen composite bar (CI [0.773, 0.939] includes 0.90) and blocked the claim. Downgraded to a **within-tier rigor gain**.
- **Concurrent-session score reconciliation (SCC-EXT).** Two independent agent sessions ran the CaRPool-seq replication in the same repo at the same time. One self-council scored 81 (clears 80) using a preprocessing-sensitive all-gene pipeline; an independent blinded critique of the other session scored **78** (does not clear 80) using the standard HVG pipeline. On discovering the divergence we adopted the committed, tested HVG pipeline as canonical, retained the all-gene run as corroboration (it agrees on the robust finding), and **reconciled merit to 78** — keeping the conservative, better-supported number. A concurrent-write incident that briefly quarantined one session's files as "foreign" was resolved and the files reclassified as corroboration.

Both/all corrections are described further in the auditor positioning:
[`docs/RELATED_WORK_positioning.md`](docs/RELATED_WORK_positioning.md).
