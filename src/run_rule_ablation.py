"""
Item 3 — commit_gate rule ablation (flip-count-first; B5 fix).

Cycle-2-PASS methodology, honoring critique B5:

  - Report FLIP COUNT FIRST. A rule that flips zero verdicts on the decision set is labeled
    "inactive - undetermined", NOT "not pivotal".
  - The original methodology tied "pivotal" to Δ(item-1 validated precision). Item 1 was
    DROPPED (B2 pre-check: T1D has only 9 ClinVar genes, 2 in test fold -> no powered truth
    set). Therefore Δprecision is UNDEFINED here and we DO NOT fabricate it. We report the
    measurable quantity (verdict flips) and state the precision-CI test as future work
    contingent on a powered truth channel.
  - Pre-registered discount: the genetic floor will appear pivotal on the trio because
    PRKCQ is withheld on it; this is expected and, per B1, genetics-confounded -> reported
    but discounted.

Ablation = force one rule to always-pass (remove its ability to block), singly AND pairwise,
recompute the trio verdicts, count flips vs the baseline gate.

Decision set: the demo trio (CD226/RASGRP1/PRKCQ) with their live-derived verdict inputs,
reconstructed exactly as in run_pipeline_day5.py. Deterministic; local CPU; $0.
Writes: data/gold/rule_ablation_receipt.json + manifest entry.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from dataclasses import replace
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.commit_gate import (  # noqa: E402
    CommitGate,
    CriticVerdict,
    MIN_GENETIC_ASSOCIATION,
    MIN_TRUST_SCORE,
)
from trustlayer.trust import derive_trust  # noqa: E402
from trustlayer.eqtl_gate import load_or_query as _eqtl_query  # noqa: E402

GENETICS = Path("data/gold/genetics_gate_receipt.json")
OUT = Path("data/gold/rule_ablation_receipt.json")
MANIFEST = Path("_script_manifest.jsonl")

ROLES = {"CD226": "ANCHOR", "RASGRP1": "BET", "PRKCQ": "CONTROL"}
RULES = ["leakage", "genetic_floor", "trust_floor", "eqtl"]

_CURATED = {
    "CD226": dict(
        genome_wide_sig_snp=True,
        eqtl_direction_consistent=True,
        proxy_tissue_only=False,
    ),
    "RASGRP1": dict(
        genome_wide_sig_snp=True, eqtl_direction_consistent=None, proxy_tissue_only=True
    ),
    "PRKCQ": dict(
        genome_wide_sig_snp=False,
        eqtl_direction_consistent=None,
        proxy_tissue_only=True,
    ),
}


def _sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _manifest(p: Path, note: str) -> None:
    with MANIFEST.open("a") as fh:
        fh.write(
            json.dumps(
                {
                    "path": str(p),
                    "sha256": _sha256(p),
                    "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "note": note,
                }
            )
            + "\n"
        )


def build_verdicts() -> dict[str, CriticVerdict]:
    genetics = json.loads(GENETICS.read_text())
    ga = {s: genetics["trio"][s]["genetic_association"] for s in ROLES}
    trust = derive_trust()
    out = {}
    for s in ROLES:
        r = _CURATED[s]
        out[s] = CriticVerdict(
            gene=s,
            genetic_association=ga[s],
            trust_score=trust[s]["trust"],
            leakage_audit_passed=True,
            genome_wide_sig_snp=r["genome_wide_sig_snp"],
            celltype_matched_eqtl=_eqtl_query(s)["celltype_matched_eqtl"],
            eqtl_direction_consistent=r["eqtl_direction_consistent"],
            proxy_tissue_only=r["proxy_tissue_only"],
        )
    return out


def ablate(v: CriticVerdict, disabled: set[str]) -> CriticVerdict:
    """Return a verdict with the named rules neutralized (forced to pass)."""
    kw = {}
    if "leakage" in disabled:
        kw["leakage_audit_passed"] = True  # already True; neutralizing keeps it passing
    if "genetic_floor" in disabled:
        kw["genetic_association"] = max(
            v.genetic_association or 0.0, MIN_GENETIC_ASSOCIATION
        )
    if "trust_floor" in disabled:
        kw["trust_score"] = max(v.trust_score or 0.0, MIN_TRUST_SCORE)
    if "eqtl" in disabled:
        kw["celltype_matched_eqtl"] = True
        kw["eqtl_direction_consistent"] = True
        kw["proxy_tissue_only"] = False
    return replace(v, **kw)


def decisions(gate: CommitGate, verdicts: dict, disabled: set[str]) -> dict[str, str]:
    return {
        s: gate.evaluate(ablate(v, disabled), role=ROLES[s]).decision
        for s, v in verdicts.items()
    }


def main() -> int:
    gate = CommitGate(out_dir="data/gold/nominations_ablation")
    verdicts = build_verdicts()
    baseline = decisions(gate, verdicts, set())

    def flips(dis: set[str]) -> tuple[int, dict]:
        d = decisions(gate, verdicts, dis)
        changed = {
            s: {"baseline": baseline[s], "ablated": d[s]}
            for s in ROLES
            if d[s] != baseline[s]
        }
        return len(changed), changed

    single = {}
    for rule in RULES:
        n, changed = flips({rule})
        single[rule] = {
            "flips": n,
            "changed": changed,
            "label": "ACTIVE" if n > 0 else "inactive - undetermined",
        }

    pairwise = {}
    for a, b in combinations(RULES, 2):
        n, changed = flips({a, b})
        n_a = single[a]["flips"]
        n_b = single[b]["flips"]
        pairwise[f"{a}+{b}"] = {
            "flips": n,
            "changed": changed,
            "jointly_pivotal_beyond_singles": n > max(n_a, n_b),
        }

    receipt = {
        "purpose": "Item 3 rule ablation (flip-count-first, B5 fix)",
        "decision_set": "demo trio (CD226/RASGRP1/PRKCQ) — n=3",
        "baseline_verdicts": baseline,
        "single_rule_ablation": single,
        "pairwise_ablation": pairwise,
        "precision_test_status": (
            "Δ(validated precision) is UNDEFINED — item 1 was dropped (B2: T1D has only 9 "
            "ClinVar genes; no powered truth set). Flip count is the measurable quantity. A "
            "precision-CI test (rule pivotal iff Δlift CI excludes 0) is future work, "
            "contingent on a genetics-INDEPENDENT powered truth channel."
        ),
        "preregistered_discount": (
            "The genetic floor flips PRKCQ (and appears 'pivotal'); this is EXPECTED and, per "
            "B1, genetics-confounded, so its pivotalness is DISCOUNTED, not headlined."
        ),
    }
    OUT.write_text(json.dumps(receipt, indent=2))
    _manifest(OUT, "item-3 rule ablation")

    print("=== ITEM 3 — RULE ABLATION (flip-count-first) ===")
    print(f"baseline trio verdicts: {baseline}")
    print("\nSINGLE-RULE ablation (flips vs baseline):")
    for rule, r in single.items():
        print(f"  {rule:14s} flips={r['flips']}  [{r['label']}]  {r['changed'] or ''}")
    print("\nPAIRWISE ablation (jointly-pivotal beyond singles?):")
    for pair, r in pairwise.items():
        jp = "  <-- JOINTLY PIVOTAL" if r["jointly_pivotal_beyond_singles"] else ""
        print(f"  {pair:28s} flips={r['flips']}{jp}")
    print(f"\nreceipt -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
