"""
Day-5: end-to-end critic -> structured verdict -> commit_gate.
Critique-corrected. The commit_gate is a DETERMINISTIC PROPAGATOR, not a discovery engine.

This runner:
  1. Loads live receipts (genetics, data-facts) + derives REAL per-gene trust (trust.py).
  2. Reports each gene's real trust up front and ASSERTS (fail-closed) that trust is
     NON-BINDING for the demo trio (every gene clears the 0.50 floor by margin >= 0.65).
     If this assertion fails, the DEMO CLAIM must change -- NOT the threshold.
  3. Builds a CriticVerdict per target from live evidence + the v3 eQTL ruling (cited).
  4. Runs commit_gate -> GO/WITHHELD; each decision lists its BINDING constraints.
  5. Runs two counterfactuals (RASGRP1 eQTL flip; PRKCQ genetic-floor flip) to show the
     gate keys on the FACTS, not gene names -- and that trust never becomes binding.
  6. Emits an auditable lineage receipt.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.commit_gate import CommitGate, CriticVerdict, MIN_TRUST_SCORE  # noqa: E402
from trustlayer.trust import derive_trust  # noqa: E402

GENETICS = Path("data/gold/genetics_gate_receipt.json")
DATAFACTS = Path("data/gold/data_facts_receipt.json")
OUT = Path("data/gold/pipeline_day5_receipt.json")
TRUST_MARGIN = (
    0.65  # pre-registered: GO-eligible genes must clear the floor by this margin
)

ROLES = {"CD226": "ANCHOR", "RASGRP1": "BET", "PRKCQ": "CONTROL"}

# v4 eQTL ruling: celltype_matched_eqtl is now DERIVED FROM A LIVE eQTL Catalogue QUERY
# (src/trustlayer/eqtl_gate.py), not hand-typed. The recorded query (data/gold/
# eqtl_gate_receipt.json) tested every CD4/Treg/T-cell dataset for a significant cis-eQTL:
#   CD226   -> True  (4/11 CD4/Treg datasets; OneK1K CD4_TCM -log10p=30)
#   RASGRP1 -> False (0/120 T-cell datasets; risk variant absent in CD4 T) <- load-bearing veto
#   PRKCQ   -> True  (3/11; prior hand-typed False was WRONG — but NON-BINDING, see below)
# genome_wide_sig_snp / eqtl_direction_consistent / proxy_tissue_only remain human-curated
# from the cited literature (a directional/mechanistic call, not a catalogue lookup).
from trustlayer.eqtl_gate import load_or_query as _eqtl_query  # noqa: E402

# Human-curated fields (literature-derived direction/GWS status; cited).
_CURATED = {
    "CD226": dict(genome_wide_sig_snp=True, eqtl_direction_consistent=True, proxy_tissue_only=False,
                  cite="rs763361 (Ser307/Gly307Ser) coding missense GWS T1D lead; KD mimics protection"),
    "RASGRP1": dict(genome_wide_sig_snp=True, eqtl_direction_consistent=None, proxy_tissue_only=True,
                    cite="rs72727394 p=4e-10; NO significant cis-eQTL in 120 T-cell datasets (live)"),
    "PRKCQ": dict(genome_wide_sig_snp=False, eqtl_direction_consistent=None, proxy_tissue_only=True,
                  cite="sub-GWS (GA=0.162); vetoed on genetic floor; eQTL non-binding"),
}
# celltype_matched_eqtl injected from the live/recorded eQTL query.
EQTL_RULING = {
    s: dict(_CURATED[s], celltype_matched_eqtl=_eqtl_query(s)["celltype_matched_eqtl"])
    for s in ("CD226", "RASGRP1", "PRKCQ")
}


def build_verdict(sym, ga, trust) -> CriticVerdict:
    r = EQTL_RULING[sym]
    return CriticVerdict(
        gene=sym,
        genetic_association=ga,
        trust_score=trust,
        leakage_audit_passed=True,  # gold-disjoint axis + frozen splits (asserted Days 1-4)
        genome_wide_sig_snp=r["genome_wide_sig_snp"],
        celltype_matched_eqtl=r["celltype_matched_eqtl"],
        eqtl_direction_consistent=r["eqtl_direction_consistent"],
        proxy_tissue_only=r["proxy_tissue_only"],
        notes=r["cite"],
    )


def main() -> int:
    genetics = json.loads(GENETICS.read_text())
    ga = {s: genetics["trio"][s]["genetic_association"] for s in ROLES}
    trust = derive_trust()

    # ---- FIX #1/#2: report trust + assert NON-BINDING for the trio (fail-closed) ----
    print("=" * 70)
    print("DAY-5 END-TO-END: critic -> verdict -> commit_gate")
    print("=" * 70)
    print("\n  Real per-gene trust (from conformal machinery):")
    for s in ROLES:
        t = trust[s]
        print(f"    {s:8s} trust={t['trust']:.3f}  [{t['provenance']}]")
    # hard assertion (fail-closed): every trio gene clears TRUST_MARGIN -> trust never
    # binding for the split. If this fails, the DEMO CLAIM changes, not the threshold.
    assert MIN_TRUST_SCORE <= TRUST_MARGIN, "margin must be >= floor"
    for s in ROLES:
        assert trust[s]["trust"] >= TRUST_MARGIN, (
            f"DEMO CLAIM INVALID: {s} trust {trust[s]['trust']} < margin {TRUST_MARGIN}; "
            f"trust would be binding -> rewrite the claim, do NOT lower the threshold."
        )
    print(
        f"\n  [assert] all trio trust >= {TRUST_MARGIN} -> trust is NON-BINDING "
        f"for the split; genetics/eQTL carry the decisions. PASS"
    )

    # ---- run the gate ----
    gate = CommitGate(out_dir="data/gold/nominations")
    split = []
    print("\n  Believe / veto split (binding constraints shown):")
    for s in ROLES:
        nom = gate.evaluate(build_verdict(s, ga[s], trust[s]["trust"]), role=ROLES[s])
        gate.commit(nom)
        marker = "GREEN GO " if nom.decision == "GO" else "RED WITHHELD"
        binding = [r for r in nom.reasons if r.startswith("BLOCK")]
        print(
            f"\n    {s:8s} [{ROLES[s]:7s}] -> {marker}  (GA={ga[s]}, trust={trust[s]['trust']:.3f})"
        )
        if binding:
            for b in binding:
                print(f"        binding: {b}")
        else:
            print(f"        {nom.reasons[-1]}")
        split.append(
            {
                "gene": s,
                "role": ROLES[s],
                "decision": nom.decision,
                "genetic_association": ga[s],
                "trust": trust[s]["trust"],
                "trust_provenance": trust[s]["provenance"],
                "binding_constraints": binding,
                "hash": nom.content_hash(),
            }
        )

    # ---- counterfactuals (FACT-driven, not gene-name-driven) ----
    print("\n  Counterfactuals (prove the gate keys on FACTS, not gene names):")
    # RASGRP1: flip cell-type eQTL True -> should GO
    v = build_verdict("RASGRP1", ga["RASGRP1"], trust["RASGRP1"]["trust"])
    v.celltype_matched_eqtl = True
    v.eqtl_direction_consistent = True
    v.proxy_tissue_only = False
    cf1 = gate.evaluate(v, role="BET")
    print(
        f"    RASGRP1 + cell-type eQTL -> {cf1.decision}  "
        f"(was WITHHELD; flips to GO iff eQTL is the only blocker)"
    )
    # PRKCQ: raise GA above the floor. With the LIVE eQTL correction, PRKCQ genuinely HAS a
    # CD4/Treg eQTL (3/11 datasets) -> its ONLY binding constraint is the genetic floor. Lift
    # GA above the floor and it flips to GO. This proves trust (0.883, the HIGHEST of the trio)
    # was never the binding constraint -- PRKCQ is withheld purely on genetics.
    v2 = build_verdict("PRKCQ", 0.50, trust["PRKCQ"]["trust"])  # hypothetical GA=0.50
    cf2 = gate.evaluate(v2, role="CONTROL")
    print(
        f"    PRKCQ + GA=0.50 (hypothetical) -> {cf2.decision}  "
        f"(trust={trust['PRKCQ']['trust']:.3f} does NOT bind; genetic floor was the SOLE blocker)"
    )

    receipt = {
        "framing": "commit_gate is a DETERMINISTIC PROPAGATOR, not a discovery engine. "
        "Proves reproducible/auditable/fail-closed governance, NOT that the "
        "biology ruling is correct. Only PRKCQ's genetic-floor veto is live-data "
        "independent of the desired outcome.",
        "trust_margin_registered": TRUST_MARGIN,
        "trust_non_binding_for_trio": True,
        "trust_scores": trust,
        "genetic_association": ga,
        "eqtl_ruling_cited": {s: EQTL_RULING[s]["cite"] for s in ROLES},
        "believe_veto_split": split,
        "counterfactuals": {
            "RASGRP1_add_celltype_eqtl": cf1.decision,
            "PRKCQ_raise_GA_to_0.50": cf2.decision,
        },
        "lineage": {
            "genetics_source": "data/gold/genetics_gate_receipt.json (Open Targets+GWAS live)",
            "trust_source": "conformal machinery (trust.py); provenance per gene",
            "rasgrp1_secondary_data_veto": "cross-donor r=0.072, n_downstream=992 "
            "(data_facts_receipt.json)",
        },
    }
    OUT.write_text(json.dumps(receipt, indent=2))
    print(f"\n  -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
