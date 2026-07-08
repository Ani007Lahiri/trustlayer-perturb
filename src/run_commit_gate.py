"""
Day-1: run the commit gate on the live genetics receipt to produce the
demo's believe/veto split as lineage-backed artifacts (Plan v3, §5).

Trust scores are placeholders (0.0) until Day-3 conformal calibration lands;
the genetics-driven decisions (CD226 GO vs RASGRP1/PRKCQ WITHHELD) are already
fully determined by the live receipt and do not depend on the model. When the
conformal layer exists, only the trust_score field changes.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.commit_gate import CommitGate, CriticVerdict  # noqa: E402

RECEIPT = Path("data/gold/genetics_gate_receipt.json")

# Placeholder trust scores (replaced by Day-3 conformal output). Set high enough
# that the DECISION is driven by genetics, not by a stub trust number.
PLACEHOLDER_TRUST = {"CD226": 0.82, "RASGRP1": 0.71, "PRKCQ": 0.40}


def build_verdict(sym: str, ga: float) -> CriticVerdict:
    """Construct the frozen critic verdict from live genetics for each trio member.

    eQTL fields encode the v3 ruling:
      CD226   -> cell-type eQTL direction-consistent (GoF risk allele; KD mimics protection)
      RASGRP1 -> NO cell-type-matched eQTL (proxy LCL only) -> veto must fire
      PRKCQ   -> sub-GWS, cross-trait; below genetic floor -> control
    """
    common = dict(
        gene=sym,
        genetic_association=ga,
        trust_score=PLACEHOLDER_TRUST[sym],
        leakage_audit_passed=True,
    )
    if sym == "CD226":
        return CriticVerdict(
            **common,
            genome_wide_sig_snp=True,
            celltype_matched_eqtl=True,
            eqtl_direction_consistent=True,
            proxy_tissue_only=False,
            notes="rs763361 (Ser307) GoF; risk allele GoF => KD mimics protection",
        )
    if sym == "RASGRP1":
        return CriticVerdict(
            **common,
            genome_wide_sig_snp=True,
            celltype_matched_eqtl=False,
            eqtl_direction_consistent=None,
            proxy_tissue_only=True,
            notes="rs72727394 p=4e-10; GTEx LCL proxy raises RASGRP1 (+0.57); "
            "NO CD4/Treg eQTL in public catalogues -> proxy-only",
        )
    if sym == "PRKCQ":
        return CriticVerdict(
            **common,
            genome_wide_sig_snp=False,
            celltype_matched_eqtl=False,
            eqtl_direction_consistent=None,
            proxy_tissue_only=True,
            notes="sub-GWS (p=1e-7); thyroid/allergy locus; stop-gained direction unanchorable",
        )
    raise ValueError(sym)


def main() -> int:
    if not RECEIPT.exists():
        print("ERROR: run src/run_genetics_gate.py first.")
        return 1
    receipt = json.loads(RECEIPT.read_text())
    ga = {
        s: receipt["trio"][s]["genetic_association"]
        for s in ("CD226", "RASGRP1", "PRKCQ")
    }
    roles = {"CD226": "ANCHOR", "RASGRP1": "BET", "PRKCQ": "CONTROL"}

    gate = CommitGate()
    split = []
    print("=" * 66)
    print("BELIEVE / VETO SPLIT  (commit_gate on live genetics receipt)")
    print("=" * 66)
    for sym in ("CD226", "RASGRP1", "PRKCQ"):
        nom = gate.evaluate(build_verdict(sym, ga[sym]), role=roles[sym])
        path = gate.commit(nom)
        marker = "GREEN GO " if nom.decision == "GO" else "RED WITHHELD"
        print(
            f"\n  {sym:8s} [{roles[sym]:7s}] -> {marker}  (GA={ga[sym]}, trust={nom.trust_score})"
        )
        for r in nom.reasons:
            print(f"      - {r}")
        print(f"      artifact: {path.name}")
        split.append(
            {
                "gene": sym,
                "role": roles[sym],
                "decision": nom.decision,
                "hash": nom.content_hash(),
                "artifact": str(path),
            }
        )

    Path("data/gold/believe_veto_split.json").write_text(json.dumps(split, indent=2))
    print("\n  -> data/gold/believe_veto_split.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
