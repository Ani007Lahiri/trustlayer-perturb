"""
Day-1 genetics gate runner (Plan v3, §4).

Produces a reproducible JSON receipt (data/gold/genetics_gate_receipt.json)
credentialing the CD226 / RASGRP1 / PRKCQ target hierarchy against live
Open Targets + GWAS Catalog evidence. This receipt is what the demo's
believe/veto split and the nomination memo cite -- no hand-typed numbers.

Run:  python src/run_genetics_gate.py
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.genetics import (  # noqa: E402
    ot_full_t1d_ranking,
    ot_genetic_association,
    gwas_catalog_associations,
    T1D_EFO,
)

TRIO = ["CD226", "RASGRP1", "PRKCQ"]
GOLD_CANDIDATES = [
    "SH2B3",
    "CTSH",
    "PTPN22",
    "INS",
    "IL2RA",
    "CTLA4",
    "IFIH1",
    "BACH2",
    "TYK2",
]
# v3 lead SNPs to reconfirm.
LEAD_SNPS = {"CD226": "rs763361", "RASGRP1": "rs72727394"}

OUT = Path("data/gold/genetics_gate_receipt.json")


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    print(f"[gate] Open Targets full T1D ranking ({T1D_EFO}) ...")
    ranking, total = ot_full_t1d_ranking()
    print(f"[gate] {total} targets associated with T1D.")

    print("[gate] genetic_association datatype scores for trio + gold candidates ...")
    assoc = ot_genetic_association(TRIO + GOLD_CANDIDATES)
    for a in assoc.values():
        r = ranking.get(a.symbol)
        a.rank = r[0] if r else None
        a.total_targets = total

    print("[gate] GWAS Catalog lead-SNP reconfirmation ...")
    snp_evidence = {sym: gwas_catalog_associations(rs) for sym, rs in LEAD_SNPS.items()}

    # --- v3 hierarchy assertion: CD226 must out-anchor RASGRP1 must out-anchor PRKCQ ---
    ga = {s: assoc[s].genetic_association for s in TRIO}
    hierarchy_ok = (ga["CD226"] or 0) > (ga["RASGRP1"] or 0) > (ga["PRKCQ"] or 0)

    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "disease": {"name": "type 1 diabetes mellitus", "id": T1D_EFO},
        "source": "Open Targets Platform v4 GraphQL + GWAS Catalog REST (live)",
        "total_t1d_targets": total,
        "trio": {s: assoc[s].as_dict() for s in TRIO},
        "gold_candidates": {s: assoc[s].as_dict() for s in GOLD_CANDIDATES},
        "lead_snp_evidence": snp_evidence,
        "v3_hierarchy_assertion": {
            "rule": "genetic_association(CD226) > RASGRP1 > PRKCQ",
            "values": ga,
            "passes": hierarchy_ok,
        },
        "roles": {
            "CD226": "ANCHOR (calibration anchor / GO exemplar)",
            "RASGRP1": "BET (novel bet / veto exemplar, proxy-flagged)",
            "PRKCQ": "DEMOTED (excluded cross-trait control)",
        },
    }
    OUT.write_text(json.dumps(receipt, indent=2))

    # --- human-readable summary ---
    print("\n" + "=" * 68)
    print("GENETICS GATE RECEIPT  (T1D, MONDO_0005147)")
    print("=" * 68)
    for s in TRIO:
        a = assoc[s]
        print(
            f"  {s:8s} role={receipt['roles'][s].split(' ')[0]:9s} "
            f"genetic_assoc={a.genetic_association}  overall={a.overall_score}  "
            f"rank=#{a.rank}/{total}"
        )
    print(
        f"\n  v3 hierarchy (CD226 > RASGRP1 > PRKCQ): "
        f"{'PASS' if hierarchy_ok else 'FAIL'}  {ga}"
    )
    print("\n  Gold-set candidate tier (top T1D targets, for the cited manifest):")
    for s in GOLD_CANDIDATES:
        a = assoc[s]
        print(
            f"    {s:8s} genetic_assoc={a.genetic_association}  rank=#{a.rank}/{total}"
        )
    print(f"\n  receipt -> {OUT}")
    return 0 if hierarchy_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
