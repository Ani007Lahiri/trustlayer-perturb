# T1D Target Nomination — pressure-tested shortlist

*Produced by 3 parallel research agents (T1D genetics · Treg-stability biology · druggability/novelty) → scored synthesis → a skeptical-immunologist red-team. Purpose: turn the "one prospective bet" in the plan from a placeholder into a real, expert-vetted, orderable nomination.*

> **⚠️ v2.1 UPDATE — SUPERSEDES the iteration-2 framing below (see `README.md`, `SCOPE.md`, and `revised_plans/PLAN_OF_ATTACK_v2.pdf` §3a/§3b).** A primary-source verification pass corrected three things in this memo:
> 1. **CD226 is NOT "direction-secure."** The KD-on-Treg-suppression direction is **contested**: [Ma et al., *Cell Reports* 2023;42(10):113306, PMID 37864795](https://pubmed.ncbi.nlm.nih.gov/37864795/) — verified real via NCBI — shows Treg-conditional CD226 deletion *impairs* suppression and worsens GvHD, the OPPOSITE of the Brusko/Thirawatananond *Diabetes* 2023 NOD-protection result (same Foxp3-Cre design). Reframe CD226 as the **hedged/moderate-trust calibration flagship**: the pipeline surfaces it and *openly assigns it hedged trust*, citing the Ma-vs-Thirawatananond split on camera. Keep the human-Treg KD experiment as *resolving a real, cited split* (stabilization vs blunted activation at matched activation; TIGIT-dependence) — NOT "resolving a nonexistent contradiction." The competitor is a real anti-DNAM-1 TNAX asset (describe the asset, not the "RVW101" code).
> 2. **PRKCQ is now the primary novel bet; RASGRP1 is a gated fallback.** RASGRP1's "hypermorphic risk allele → KD mimics protection" premise is **UNCONFIRMED** — the T1D lead SNP is not a significant RASGRP1 eQTL in any CD4/Treg resource (the risk↑RASGRP1 signal is a *different SNP in SLE*), and germline loss *causes* autoimmunity/lymphoma. Run a **Day-1 eQTL-direction go/no-go**; if it fails, PRKCQ (PKC-θ blockade *enhances* Treg suppression — direction verified) is the primary bet.
> 3. **RBPJ-NCOR is Sakaguchi's (Nature 2025), NOT Marson's** — only RNF20 is Marson's (Nature 2020). RBPJ stays a landmine, but not "his own paper."
>
> Everything else here (STUB1/UBASH3A analysis, other landmines, validation design, command pack) stands. The ★-ratings in the CD226 table below reflect the *pre-verification* view — read them through the corrections above.

---

## Verdict: yes — a credible, novel-enough, right-direction nomination exists

The impact claim is no longer hand-waving. There is at least one gene that is simultaneously (a) anchored to a **fine-mapped coding candidate-causal T1D variant**, (b) **risk-allele-direction secure** (though the KD-on-suppression sign is *contested* — see the v2.1 banner and the Direction row: Ma 2023 *Cell Reports* vs Brusko 2023 *Diabetes*), (c) backed by **Treg-intrinsic in-vivo causal evidence**, and (d) trivially orderable and validatable in one week. That gene is **CD226** — and per v2.1 it is deployed as the *hedged calibration flagship*, not a confident anchor. Two backups cover the angles CD226 doesn't.

The one honest caveat (the red-team was blunt about it): CD226 is *not* a blue-sky-novel gene — it's an active idea in the T1D-Treg field, **and its KD-on-suppression direction is contested (v2.1)**. So the novelty lives in the **method** (a trust-calibrated pipeline that independently re-derives a genetically-anchored target and *scores how much to believe it*) and in the **genotype-stratified framing** — not in "nobody thought of this gene." The contested direction is not a weakness to hide; it is *exactly why CD226 is the calibration flagship* — the pipeline hedges its trust and says why, on camera. Oversell it as a confident discovery and a Marson-lab judge will puncture it; present it as honestly-hedged and it becomes the strongest proof the trust layer works.

---

## 🥇 Headline nomination — CD226 (DNAM-1)

**The one-liner:** Knock down CD226, the activating costimulatory receptor whose T1D *risk* allele is a gain-of-function, and Tregs hold their FOXP3 under inflammatory stress instead of converting to pathogenic effectors — orderable in primary human CD4 next week, with a genotype-stratified twist Pritchard will love.

| Axis | Evidence | Strength |
|---|---|---|
| **T1D genetics** | GWS, **fine-mapped coding causal variant** rs763361 (Gly307Ser) — not a fragile eQTL. Ser307 adds a cytoplasmic phospho-site → ↑ERK/PI3K/pSTAT4/IFN-γ. Risk = *more* CD226 activity. | ★★★★★ |
| **Direction (KD stabilizes?)** | ⚠️ **CONTESTED — downgraded from ★★★★★ (v2.1).** Risk-allele sign is secure (GoF → KD opposes risk), and blockade **raised CD4⁺FOXP3⁺ Tregs/FOXP3 in human cells** (humanized GvHD) + NOD protection (Brusko/Thirawatananond *Diabetes* 2023). **BUT** [Ma et al., *Cell Reports* 2023;42(10):113306, PMID 37864795](https://pubmed.ncbi.nlm.nih.gov/37864795/) — same Foxp3-Cre design — shows deletion *impairs* suppression (mTOR/Myc) and worsens GvHD. So mouse is split on the KD-on-**function** sign; only the risk-allele sign "agrees." **Do NOT say "human + mouse + genetics all agree."** | ★★☆☆☆ (contested) |
| **Treg-intrinsic causal** | **Treg-conditional (Foxp3-Cre) Cd226-KO in NOD reduced diabetes** (44% vs 66.7%, p=0.042) and insulitis (p<0.0001), ↑Treg TIGIT, ↓ "ex-Treg" FOXP3 loss. In-vivo, Treg-specific — no other candidate has this. | ★★★★★ |
| **Tractability** | Surface receptor → clean CRISPRi KD, **flow-cytometry KD confirmation for free**, antibody arm for orthogonal validation. | ★★★★☆ |
| **Novelty** | Fresh *experiment* (never run as CRISPRi in primary human CD4 with a Treg-stability readout; not a T1D drug program) but *known biology* (CD226/TIGIT is an IO axis; active in the T1D-Treg field). | ★★★☆☆ |

**Caveats to state openly (don't let a judge find them first):** the NOD effect was female-restricted and modest; part of the protection is TIGIT-mediated effector suppression, not purely cell-intrinsic FOXP3 stabilization; and the biggest validation risk is that "stabilization" is really just **blunted activation** (less costim → cells look calmer). The validation design below controls for exactly that.

---

## Backups (cover the angles CD226 doesn't)

**🥈 STUB1 (CHIP) — the mechanism-perfect, genetics-free play.** Inflammatory cytokines induce STUB1, which polyubiquitinates and degrades FOXP3, flipping Tregs Th1-like; **knockdown preserves FOXP3 specifically under inflammation** — a flawless rest-vs-stim, cell-intrinsic stability mechanism that will almost certainly *score* on your readout. The catch: **not a T1D GWAS locus** (Pritchard's first question has no answer) and broad proteostasis/tox risk. Use it as the *mechanistic* hit and positive stabilizer control, not the genetically-anchored headline.

**🥉 UBASH3A — the genetics-anchored effector-conversion play.** Excellent T1D genetics (GWS, fine-mapped CD4 enhancer, direct CD4 eQTL; risk raises UBASH3A in *stimulated* CD4 → suppresses IL-2). But its effect is effector/IL-2-centric and *non-cell-autonomous*, so it belongs in an **effector-conversion** readout, not a FOXP3-stability one — and you must own the human/mouse paradox (germline NOD KO *accelerates* disease). Great as the parallel arm that gives Pritchard a fine-mapped locus.

---

## ⚠️ Landmines — do NOT nominate these (the red-team caught them)

- **RNF20** — **Marson's own published FOXP3-screen hit** ([Nature 2020](https://www.nature.com/articles/s41586-020-2246-4)). Nominating it to Marson = pitching him his own paper. **RBPJ-NCOR** — a landmine (a known published FOXP3-screen hit, [Nature 2025](https://www.nature.com/articles/s41586-025-08795-5)) but it is **Sakaguchi's, NOT Marson's** (v2.1 correction). Still don't nominate it; just don't call it "Marson's own."
- **SIRPG** — real T1D locus, but risk allele *reduces* SIRPγ, so a **knockdown mimics the pathogenic direction** (you'd need CRISPRa). Directionally disqualified for a KD screen.
- **PTPN2, IKZF2/Helios** — hot IO targets, and the tractable direction is *opposite* what T1D needs.
- **BACH2, IKZF4/Eos** — genuine T1D genes but *positive* Treg regulators; knocking them down **destabilizes** Tregs. Use them as positive-control *destabilizers*, never as therapeutic KD nominations.

---

## How this plugs into the project

This is a **pre-registered watchlist**, not a claim about your (unrun) data. The pipeline's job is to independently surface these — with a mechanism and a trust score — from the real Perturb-seq:

- **The demo anchor beat:** surface **CD226** framed as "our pipeline re-derived a genetics-anchored T1D target and — because the KD-on-suppression direction is *contested in the literature* — scored it **moderate/hedged trust**, cited the Ma-vs-Thirawatananond split, and still handed us a one-week experiment to resolve it." (This is the *hedged calibration flagship*, per v2.1 — the confident-recovery beat is carried by the non-canonical Iakovliev genes instead.)
- **Blinded validation, tier 2:** CD226 recovery is a *coding-variant-anchored* external check that complements the Iakovliev core genes (FOXP3/CTLA4/STAT1). If the model ranks CD226 high with a high trust score, that's a second, independent "the model found real biology it was never told" moment.
- **Positive controls for the screen:** FOXP3-KD and BACH2-KD (destabilizers, should score negative-direction); STUB1-KD (stabilizer, positive-direction). These make your calibration curve credible.
- **The genotype-by-perturbation angle** (test whether rs763361 risk-carriers show a bigger KD benefit) is the single most Pritchard-pleasing, hard-to-copy flourish available.

---

## Killer one-week validation for CD226

- **Perturbation:** dCas9-KRAB CRISPRi knockdown of CD226 in FACS-sorted primary human CD4⁺CD25⁺CD127ˡᵒ Tregs (≥3–4 HLA-diverse donors); run in parallel in total CD4 for the effector-conversion arm. Confirm KD by **surface CD226 flow** (free protein-level check).
- **Challenge (the ex-Treg condition):** anti-CD3 restim + **CD155-Fc** (engage the CD226/TIGIT axis) under IL-6 + low IL-2.
- **Readouts (all 1-week-feasible):** FOXP3 retention (% and gMFI) [primary]; TIGIT up-regulation [mechanism]; ex-Treg IFN-γ⁺/FOXP3-lost fraction; CTV suppression assay [function]; Th1/IFN-γ conversion in total CD4.
- **Controls:** non-targeting sgRNA; FOXP3-KD (positive destabilizer); STUB1-KD (positive stabilizer); anti-CD226 blocking-Ab (orthogonal, non-genetic); genotype donors at rs763361.
- **The risk the design must kill:** that "stabilization" is just **blunted activation**. Compare FOXP3 retention **at matched activation levels**, and test whether **TIGIT blockade abrogates the CD226-KD benefit** (if yes → the intended mechanism operates in human cells; if FOXP3 rescue persists → something cleaner).

---

## Reproducible evidence pack (gget + PubMed — run these to regenerate the log)

```bash
# --- gget: gene metadata, tractability & disease links (run in a venv: pip install gget) ---
for G in CD226 UBASH3A STUB1 IL6ST; do gget search -s human "$G" -o ${G}_search.json; done
# take the ENSG id from each search, then:
gget info      ENSG00000150637        -o CD226_info.json        # CD226 (confirm id from search)
gget opentargets ENSG00000150637 -r tractability -o CD226_tract.json
gget opentargets ENSG00000150637 -r diseases     -o CD226_dis.json   # look for EFO_0001359 (T1D)
# repeat opentargets for UBASH3A / STUB1 / IL6ST to compare tractability & T1D association scores
```

```bash
# --- PubMed E-utilities: reproducible literature log (NCBI_EMAIL in env) ---
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&email=$NCBI_EMAIL&term="
# CD226 T1D Treg causal
"$BASE" 'CD226[tiab] AND (type 1 diabetes[tiab] OR T1D[tiab]) AND (Treg[tiab] OR FOXP3[tiab])'
# rs763361 / Gly307Ser coding variant
"$BASE" '(rs763361[tiab] OR Gly307Ser[tiab] OR "CD226 307"[tiab]) AND diabetes[tiab]'
# STUB1 FOXP3 degradation
"$BASE" '(STUB1[tiab] OR CHIP[tiab]) AND FOXP3[tiab] AND (ubiquitin*[tiab] OR degradation[tiab])'
# UBASH3A T1D CD4 / IL-2
"$BASE" 'UBASH3A[tiab] AND (diabetes[tiab] OR "type 1"[tiab]) AND (IL-2[tiab] OR "CD4"[tiab])'
```
Record for each pass: exact query, date searched, result count (per the pubmed-database skill's output discipline).

---

## Honest impact read

CD226 takes your impact from "placeholder" to **a genetics-anchored, in-vivo-supported, orderable T1D nomination with a genotype-stratified validation** — stronger and more credible than the rivals' correlational drug-target rankings, and Gladstone-runnable. It is *not* a blockbuster-novel gene, and its KD-on-suppression direction is *contested* (v2.1) — so lead with the **method + honest hedged trust + the Monday experiment that resolves the split**, keep STUB1/UBASH3A as the "and it generalizes" evidence, and the impact axis is genuinely win-grade. Whether you *win* still rides mostly on Demo + Claude Use — but you no longer have a soft flank on impact, and the contested direction is converted from a liability into the demo's honesty proof.

---

## Sources
- CD226 Treg-conditional KO / T1D, coding variant, NOD **protection** (Brusko/Thirawatananond) — Diabetes 2023. https://diabetesjournals.org/diabetes/article/72/11/1629/153546/
- **[v2.1] CD226 Treg-conditional KO — OPPOSITE direction (deletion impairs suppression via mTOR/Myc, worsens GvHD)** — Ma et al., *Cell Reports* 2023;42(10):113306, PMID 37864795 (verified via NCBI). https://pubmed.ncbi.nlm.nih.gov/37864795/ — this is the *real* contradiction; the two same-design KO papers disagree on the KD-on-function sign.
- CD226 rs763361 gain-of-function signaling — Front Immunol 2022. https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.886736/full
- STUB1 degrades FOXP3 under inflammation — Immunity 2013. https://pmc.ncbi.nlm.nih.gov/articles/PMC3817295/
- UBASH3A → ↓IL-2 in CD4, T1D — Diabetes 2017. https://pmc.ncbi.nlm.nih.gov/articles/PMC5482087/
- UBASH3A germline KO accelerates NOD — Sci Rep 2020. https://www.nature.com/articles/s41598-020-68956-6
- IL6ST/ANKRD55 T1D CD4 coloc — Nat Genet 2021. https://pmc.ncbi.nlm.nih.gov/articles/PMC8273124/
- Tocilizumab (anti-IL-6R) EXTEND null in T1D — JCI Insight 2021. https://insight.jci.org/articles/view/150074
- **Marson** Treg CRISPR screen (RNF20/USP22) — Nature 2020. https://www.nature.com/articles/s41586-020-2246-4
- RBPJ-NCOR FOXP3 screen — Nature 2025 (**Sakaguchi lab, NOT Marson** — v2.1 correction). https://www.nature.com/articles/s41586-025-08795-5
- SIRPG T1D splicing/direction — Diabetes 2022. https://diabetesjournals.org/diabetes/article/71/2/350/138930/
- Iakovliev T1D omnigenic core genes — AJHG 2023. https://pubmed.ncbi.nlm.nih.gov/37164005/
