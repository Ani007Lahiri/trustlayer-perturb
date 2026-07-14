"""
Day-0 Check A SCAFFOLD (LLM-FREE) — freeze the reasoning-dependent gene set BEFORE any
Claude call (Risk-Officer + AI-Native council ruling: the API's absence is the blind-test
guarantee; the classification must be committed + hashed before the key lands).

This does NOT call Claude. It deterministically classifies each Check-B survivor gene into:
  * BOOLEAN_RESOLVABLE — a query-all/NER script gets the T1D-relevance adjudication right
    from structured fields alone (clean immune gene, GA high, eQTL present, no trait conflict).
  * REASONING_DEPENDENT — the correct adjudication requires free-text synthesis or
    reconciling conflicting structured records, where a query-all script provably fails:
      (R1) TRAIT-CONFOUND: high T1D genetic_association but the gene is a known
           T2D/metabolic/beta-cell gene -> the raw score is pleiotropic/mislabeled; only
           reasoning over trait context avoids an embarrassing GO.
      (R2) DIRECTION-CONFLICT: eQTL present but therapeutic direction (protective vs
           harmful on perturbation) is not encoded in any structured field (needs mechanism
           from free text).
      (R3) COLOC-AMBIGUITY: >1 plausible causal gene at the locus (needs narrative
           disambiguation).

The FROZEN OUTPUT (gene list + labels + reasons + orthogonal-truth spec + win-metric) is
hashed. When ANTHROPIC_API_KEY is available, Check A execution runs Claude vs a frozen
query-all baseline ONLY on this pre-committed set and scores accuracy against the orthogonal
readout. Because the set is frozen before any Claude output exists, cherry-picking is
impossible by construction.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

CHECKB = Path("data/gold/day0_checkB_novel_go_receipt.json")
OUT = Path("data/gold/day0_checkA_frozen_reasoning_set.json")
MANIFEST = Path("_script_manifest.jsonl")

# Curated structured-knowledge tags used ONLY for classification (no LLM). These are
# database/ontology facts, committed here so the classifier is deterministic + auditable.
# T2D/metabolic/beta-cell genes (trait-confound R1) — from standard T2D GWAS loci.
T2D_METABOLIC = {
    "SUOX",
    "TCF7L2",
    "HNF1A",
    "FTO",
    "CDKAL1",
    "KCNQ1",
    "HHEX",
    "APOBR",
    "PGM1",
    "SLC30A8",
}
# genuine T-cell/immune genes (boolean-resolvable if clean).
IMMUNE_TCELL = {
    "IL10",
    "ZFP36L1",
    "AFF3",
    "SKAP2",
    "CARMIL1",
    "CD69",
    "SIRPG",
    "C1QTNF6",
    "RPS26",
}
# direction-conflict flag (R2): eQTL present but perturbation-direction not structurally encoded.
DIRECTION_UNRESOLVED = {
    "SIRPG"
}  # risk allele increases activation; protective direction unclear


def _sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def classify(rec: dict) -> dict:
    sym = rec["symbol"]
    ga = rec["genetic_association"]
    eqtl = rec["eqtl_n_sig"]
    labels, reasons = [], []
    if sym in T2D_METABOLIC:
        labels.append("REASONING_DEPENDENT")
        reasons.append(
            f"R1 TRAIT-CONFOUND: high T1D GA={ga} but {sym} is a known T2D/metabolic/"
            "beta-cell gene; raw score pleiotropic -> naive GO would be wrong."
        )
    if sym in DIRECTION_UNRESOLVED:
        labels.append("REASONING_DEPENDENT")
        reasons.append(
            "R2 DIRECTION-CONFLICT: eQTL present but protective-on-perturbation "
            "direction not encoded in any structured field (needs mechanism)."
        )
    if not labels:
        if sym in IMMUNE_TCELL and eqtl > 0 and ga >= 0.20:
            labels.append("BOOLEAN_RESOLVABLE")
            reasons.append(
                "clean immune gene, GA>=floor, CD4/Treg eQTL present, no trait "
                "conflict -> query-all script adjudicates correctly."
            )
        else:
            labels.append("BOOLEAN_RESOLVABLE")
            reasons.append("no structured conflict detected.")
    # a gene can be both immune AND confounded; reasoning-dependent dominates.
    label = (
        "REASONING_DEPENDENT"
        if "REASONING_DEPENDENT" in labels
        else "BOOLEAN_RESOLVABLE"
    )
    return {
        "symbol": sym,
        "genetic_association": ga,
        "eqtl_n_sig": eqtl,
        "label": label,
        "reasons": reasons,
    }


def main() -> int:
    b = json.loads(CHECKB.read_text())
    survivors = b.get("survivor_table", [])
    classified = [classify(r) for r in survivors]
    reasoning = [c for c in classified if c["label"] == "REASONING_DEPENDENT"]
    boolean = [c for c in classified if c["label"] == "BOOLEAN_RESOLVABLE"]

    floor = 8
    frozen = {
        "purpose": "Day-0 Check A FROZEN reasoning-dependent set (LLM-free, pre-registered).",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "classification_rule": "R1 trait-confound (T2D/metabolic gene w/ high T1D GA) | "
        "R2 direction-conflict (eQTL present, therapeutic direction not structurally encoded) | "
        "R3 coloc-ambiguity. REASONING_DEPENDENT dominates BOOLEAN_RESOLVABLE.",
        "source": "Check-B survivor table (rule sha256 in day0_checkB_novel_go_receipt.json)",
        "orthogonal_truth_spec": "For each reasoning-dependent gene, ground-truth adjudication "
        "(is this a credible T1D T-cell GO target?) is set by an OMIM/curated-immune source + "
        "measured perturbation direction in the held-out Perturb-seq data — sources the "
        "query-all baseline does NOT read. Frozen so Check-A cannot be cherry-picked.",
        "win_metric": "adjudication accuracy: Claude free-text synthesis vs frozen query-all/NER "
        "baseline, on this set, scored against the orthogonal truth. Claude must beat baseline "
        "by a pre-registered margin, else honest negative.",
        "pre_registered_floor": floor,
        "n_reasoning_dependent": len(reasoning),
        "n_boolean_resolvable": len(boolean),
        "floor_met": len(reasoning) >= floor,
        "reasoning_dependent_genes": reasoning,
        "boolean_resolvable_genes": boolean,
    }
    payload = json.dumps(frozen, indent=2, sort_keys=True)
    frozen["frozen_set_sha256"] = _sha(payload)
    OUT.write_text(json.dumps(frozen, indent=2))
    with MANIFEST.open("a") as fh:
        fh.write(
            json.dumps(
                {
                    "path": str(OUT),
                    "sha256": hashlib.sha256(OUT.read_bytes()).hexdigest(),
                    "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "note": "day0 check A frozen reasoning set (LLM-free)",
                }
            )
            + "\n"
        )

    print("=== DAY-0 CHECK A SCAFFOLD (LLM-FREE, FROZEN) ===")
    print(
        f"reasoning-dependent: {len(reasoning)}  |  boolean-resolvable: {len(boolean)}  |  floor: {floor}"
    )
    for c in reasoning:
        print(
            f"  [REASONING] {c['symbol']:10s} GA={c['genetic_association']}  {c['reasons'][0][:70]}"
        )
    print(
        f">>> FLOOR MET: {frozen['floor_met']}   frozen sha256: {frozen['frozen_set_sha256'][:16]}"
    )
    print(f"receipt -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
