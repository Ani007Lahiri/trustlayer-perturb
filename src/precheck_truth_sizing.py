"""
Data-sizing pre-check (B2 gate) for the NEXT-RUN item-1 lift.

Per the cycle-2-PASS methodology: BEFORE any lift computation, count TEST-fold genes that
carry a held-out truth label, split into (gold-forced) vs (random-landed), and count the
genetics-eligible stratum. If TEST-fold truth positives < 15, item 1 is REPORTED AS A CASE
STUDY (n honesty), NOT an "all-genes lift".

Truth channel (gate does NOT threshold on it): ClinVar pathogenic / likely-pathogenic genes
for Type 1 diabetes mellitus, pulled live via NCBI E-utilities and cached to a receipt.
Residual-leakage caveat is recorded (ClinVar shares generative process with OT genetic_assoc
that the gate's floor uses -> B1).

Deterministic w.r.t. the frozen split (hash 45ca2893cbe7e282). No compute cost; local CPU.
Writes: data/gold/truth_sizing_receipt.json  and appends to _script_manifest.jsonl.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import requests

SPLITS = Path("data/gold/frozen_splits.json")
GOLD = Path("data/gold/t1d_gold_set.json")
GENETICS_RECEIPT = Path("data/gold/genetics_gate_receipt.json")
OUT = Path("data/gold/truth_sizing_receipt.json")
MANIFEST = Path("_script_manifest.jsonl")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
T1D_TERM = (
    '"type 1 diabetes mellitus"[dis] AND '
    '("pathogenic"[Clinical significance] OR "likely pathogenic"[Clinical significance])'
)
CASE_STUDY_THRESHOLD = 15  # pre-registered (methodology B2)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record_manifest(path: Path, note: str) -> None:
    entry = {
        "path": str(path),
        "sha256": _sha256(path),
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "note": note,
    }
    with MANIFEST.open("a") as fh:
        fh.write(json.dumps(entry) + "\n")


def pull_clinvar_t1d_genes() -> tuple[list[str], dict]:
    """Return (sorted unique gene symbols, provenance dict) for ClinVar T1D P/LP variants."""
    s = requests.get(
        f"{EUTILS}/esearch.fcgi",
        params={"db": "clinvar", "term": T1D_TERM, "retmax": 500, "retmode": "json"},
        timeout=30,
    )
    s.raise_for_status()
    js = s.json()["esearchresult"]
    ids = js.get("idlist", [])
    total = int(js.get("count", 0))

    genes: set[str] = set()
    # esummary in batches of 200
    for i in range(0, len(ids), 200):
        batch = ids[i : i + 200]
        if not batch:
            continue
        r = requests.get(
            f"{EUTILS}/esummary.fcgi",
            params={"db": "clinvar", "id": ",".join(batch), "retmode": "json"},
            timeout=60,
        )
        r.raise_for_status()
        res = r.json().get("result", {})
        for uid in res.get("uids", []):
            rec = res.get(uid, {})
            for g in rec.get("genes", []) or []:
                sym = g.get("symbol")
                if sym:
                    genes.add(sym.upper())
        time.sleep(0.34)  # NCBI rate limit courtesy

    prov = {
        "source": "NCBI ClinVar via E-utilities (live)",
        "term": T1D_TERM,
        "n_records_reported": total,
        "n_records_fetched": len(ids),
        "n_unique_genes": len(genes),
        "caveat": (
            "RESIDUAL-LEAKAGE (B1): ClinVar P/LP for T1D shares its generative process with "
            "the OT genetic_association score the gate's floor thresholds on. This truth is "
            "genetics-flavored; item-1 lift is scoped to 'BEYOND genetics conditional on the "
            "floor-pass stratum', not 'gate beats genetics'."
        ),
    }
    return sorted(genes), prov


def main() -> None:
    splits = json.loads(SPLITS.read_text())
    test_genes = set(splits["test_genes"])
    gold_forced = set(splits["gold_in_test"])
    split_hash = splits["content_hash"]

    truth_genes, prov = pull_clinvar_t1d_genes()
    truth_set = set(truth_genes)

    # TEST-fold truth positives, split gold-forced vs random-landed
    test_truth = sorted(test_genes & truth_set)
    test_truth_gold_forced = sorted(set(test_truth) & gold_forced)
    test_truth_random_landed = sorted(set(test_truth) - gold_forced)

    # genetics-eligible stratum: genes with genetic_association >= floor (0.20).
    # We only have live GA for the trio + gold in receipts; for the stratum count at scale
    # we note this is bounded by what OT scores are available. Record what we can verify.
    ga_available = {}
    if GENETICS_RECEIPT.exists():
        gr = json.loads(GENETICS_RECEIPT.read_text())
        # genetics_gate_receipt stores per-gene association where pulled
        for g, rec in (
            (gr.get("genes") or {}).items() if isinstance(gr.get("genes"), dict) else []
        ):
            if isinstance(rec, dict) and rec.get("genetic_association") is not None:
                ga_available[g.upper()] = rec["genetic_association"]

    n_test_truth = len(test_truth)
    is_case_study = n_test_truth < CASE_STUDY_THRESHOLD

    receipt = {
        "purpose": "B2 data-sizing pre-check for item-1 lift (methodology 2026-07-09)",
        "frozen_split_hash": split_hash,
        "case_study_threshold": CASE_STUDY_THRESHOLD,
        "truth_channel": prov,
        "counts": {
            "n_test_fold_genes": len(test_genes),
            "n_truth_genes_total": len(truth_set),
            "n_test_fold_truth_positives": n_test_truth,
            "n_test_truth_gold_forced": len(test_truth_gold_forced),
            "n_test_truth_random_landed": len(test_truth_random_landed),
        },
        "test_fold_truth_positives": test_truth,
        "test_truth_gold_forced": test_truth_gold_forced,
        "test_truth_random_landed": test_truth_random_landed,
        "ga_available_in_receipts": ga_available,
        "B2_GATE_DECISION": (
            "CASE_STUDY (n<15): item 1 is reported as a test-fold decision replay with "
            "explicit n honesty, NOT an all-genes lift."
            if is_case_study
            else "LIFT_ALLOWED (n>=15): item 1 may be reported as a lift, still stratified "
            "by the genetic floor per B1."
        ),
        "is_case_study": is_case_study,
    }
    OUT.write_text(json.dumps(receipt, indent=2))
    _record_manifest(OUT, "B2 truth-sizing pre-check for item-1")

    print("=== B2 DATA-SIZING PRE-CHECK ===")
    print(f"frozen split hash:        {split_hash}")
    print(f"truth genes (ClinVar T1D):{len(truth_set)}")
    print(f"TEST-fold genes:          {len(test_genes)}")
    print(
        f"TEST-fold truth positives:{n_test_truth}  "
        f"(gold-forced={len(test_truth_gold_forced)}, "
        f"random-landed={len(test_truth_random_landed)})"
    )
    print(f"threshold:                {CASE_STUDY_THRESHOLD}")
    print(f">>> B2 GATE: {'CASE_STUDY' if is_case_study else 'LIFT_ALLOWED'}")
    print(f"truth genes in test: {test_truth}")
    print(f"receipt -> {OUT}")


if __name__ == "__main__":
    main()
