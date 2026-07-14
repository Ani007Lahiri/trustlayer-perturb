"""
Day-0 Check A — EXECUTION v4 (enriched pool; MCC-primary; single frozen win condition).

Fixes both iteration-3 BLOCKs:
  * n_GO=1 collapse -> now runs on the ENRICHED frozen set v2 (24 genes, 4 GO / 20 NO-GO),
    so the task measures PATTERN SEPARATION across non-obvious T-cell positives (SIRPG, CD5,
    STAT4, IL10) vs hard negatives (immune genes with wrong perturbation direction; metabolic
    genes that are T1D-dominant in genetics) — not a single obvious point-identification.
  * win-condition pre-registration mismatch -> the win condition is defined HERE as a single
    constant (WIN_CONDITION_TEXT) and the SAME string is stamped into the receipt; there is no
    separate frozen file carrying a stale/contradictory rule. The frozen SET (labels) and the
    win RULE are both fixed before any Claude call in this run.

Blinding unchanged: open-ended prompt from ONLY symbol + genetic_association + eqtl_n_sig
(pulled from the enriched frozen set); never reads GO labels or the ground-truth rule fields
into the prompt. Repeated trials (temperature=0 is near-deterministic, so trials mainly detect
API nondeterminism, reported honestly as such). MCC-primary scoring vs the strongest baseline.
"""

from __future__ import annotations

import json
import math
import os
import re
import sys
import time
import urllib.request
import urllib.error
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.canonical import write_hashed_json  # noqa: E402

FROZEN_SET = Path("data/gold/day0_checkA_frozen_set_v2.json")
ORIG_FROZEN = Path(
    "data/gold/day0_checkA_frozen_reasoning_set.json"
)  # for eqtl_n_sig lookup
OUT = Path("data/gold/day0_checkA_execution_receipt_v2.json")
MANIFEST = Path("_script_manifest.jsonl")

MODEL = "claude-sonnet-4-5-20250929"
API_URL = "https://api.anthropic.com/v1/messages"
N_TRIALS = 3

MCC_MARGIN = 0.15
MCC_ABS_FLOOR = 0.50
WIN_CONDITION_TEXT = (
    f"PRIMARY METRIC = Matthews Correlation Coefficient (MCC), because raw accuracy is "
    f"uninformative under class imbalance. WIN iff (Claude MCC - best_baseline_MCC) >= "
    f"{MCC_MARGIN} AND Claude MCC >= {MCC_ABS_FLOOR}. LOSS iff Claude MCC < best_baseline "
    f"MCC. Otherwise TIE_NO_WIN (honest negative). Baselines = majority-class, eqtl-only, "
    f"query-all; strongest by MCC is the bar. Reported with Wilson 95% CI on accuracy. This "
    f"is a single-model pilot on n=24 (4 positives); a WIN is a pilot signal, not a validated "
    f"general-capability claim."
)

# eQTL n_sig is not stored in the enriched set; genes outside the original 9 need a value.
# We fetch it live once (same eqtl_gate the project uses) and cache into the receipt so the
# baseline is auditable. To keep this fast + reproducible we read from a small prefetch below.
EQTL_CACHE = Path("data/gold/_checkA_eqtl_cache.json")

PROMPT_TEMPLATE = """You are evaluating whether a human gene is a credible drug/experimental target for Type 1 Diabetes (T1D) acting through CD4+ T-cell autoimmune biology.

Gene symbol: {symbol}
Open Targets T1D genetic association score: {ga}
Significant CD4/Treg cis-eQTLs detected: {eqtl}

Think about this gene's actual biology: What tissue or cell type is it primarily known to act in? What established disease(s) is it associated with? Is its genetic signal for diabetes most plausibly driven by an autoimmune CD4+ T-cell mechanism (relevant to Type 1 diabetes), or by some other mechanism (for example pancreatic beta-cell / metabolic biology relevant to Type 2 diabetes, cardiac biology, or a non-causal association at a crowded genetic locus)?

Reason briefly (2-4 sentences), then give your final answer on the LAST line in exactly this format with no other text on that line:
GO: YES   (if this is a credible T1D CD4+ T-cell target)
GO: NO    (otherwise)"""


def _call_claude(prompt: str, api_key: str, max_retries: int = 4) -> str:
    body = json.dumps(
        {
            "model": MODEL,
            "max_tokens": 400,
            "temperature": 0,
            "messages": [{"role": "user", "content": prompt}],
        }
    ).encode()
    req = urllib.request.Request(
        API_URL,
        data=body,
        headers={
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        },
        method="POST",
    )
    last = None
    for i in range(max_retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())["content"][0]["text"]
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}: {e.read().decode()[:200]}"
        except Exception as e:  # noqa
            last = str(e)
        time.sleep(2 * (i + 1))
    raise RuntimeError(f"Claude call failed: {last}")


def _parse_go(text: str):
    m = re.findall(r"GO:\s*(YES|NO)", text, flags=re.IGNORECASE)
    return None if not m else m[-1].strip().upper() == "YES"


def _mcc(tp, tn, fp, fn):
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return (tp * tn - fp * fn) / den if den else 0.0


def _conf(preds, gt):
    tp = sum(1 for s in gt if preds[s] and gt[s])
    tn = sum(1 for s in gt if not preds[s] and not gt[s])
    fp = sum(1 for s in gt if preds[s] and not gt[s])
    fn = sum(1 for s in gt if not preds[s] and gt[s])
    return tp, tn, fp, fn


def _bal_acc(tp, tn, fp, fn):
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    return (sens + spec) / 2


def _wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d
    return (round(max(0, c - h), 4), round(min(1, c + h), 4))


def _get_eqtl(symbols: list[str]) -> dict:
    """eqtl_n_sig per gene. Use cache if present; else fetch live via project eqtl_gate."""
    if EQTL_CACHE.exists():
        cache = json.loads(EQTL_CACHE.read_text())
        if all(s in cache for s in symbols):
            return cache
    from trustlayer import eqtl_gate

    orig = {
        g["symbol"]: g["eqtl_n_sig"]
        for g in json.loads(ORIG_FROZEN.read_text())["reasoning_dependent_genes"]
    }
    out = {}
    for s in symbols:
        if s in orig:
            out[s] = orig[s]
            continue
        try:
            import requests

            q = 'query S($q:String!){ search(queryString:$q, entityNames:["target"]){ hits{ id entity } } }'
            r = requests.post(
                "https://api.platform.opentargets.org/api/v4/graphql",
                json={"query": q, "variables": {"q": s}},
                timeout=30,
            ).json()
            ensg = next(
                (
                    h["id"]
                    for h in r["data"]["search"]["hits"]
                    if h.get("entity") == "target"
                ),
                None,
            )
            e = eqtl_gate.query_celltype_eqtl(ensg) if ensg else {"n_sig": 0}
            out[s] = e["n_sig"]
        except Exception:
            out[s] = 0
        time.sleep(0.2)
    EQTL_CACHE.write_text(json.dumps(out, indent=2))
    return out


def main() -> int:
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        print("ANTHROPIC_API_KEY not set — blocked.")
        return 1

    fs = json.loads(FROZEN_SET.read_text())
    genes = fs["genes"]
    gt = {g["symbol"]: g["ground_truth_credible_t1d_tcell_go"] for g in genes}
    ga = {g["symbol"]: g["ot_t1d_assoc"] for g in genes}
    syms = list(gt)
    eqtl = _get_eqtl(syms)

    print(
        f"=== DAY-0 CHECK A EXECUTION v4 (enriched n={len(syms)}, {sum(gt.values())} GO) ==="
    )
    print(f"frozen set hash: {fs['frozen_set_v2_sha256'][:16]}\n")

    claude = {}
    for s in syms:
        prompt = PROMPT_TEMPLATE.format(symbol=s, ga=ga[s], eqtl=eqtl[s])
        votes, texts = [], []
        for _ in range(N_TRIALS):
            t = _call_claude(prompt, api_key)
            votes.append(_parse_go(t))
            texts.append(t)
            time.sleep(0.3)
        valid = [v for v in votes if v is not None]
        maj = Counter(valid).most_common(1)[0][0] if valid else False
        claude[s] = {
            "go": bool(maj),
            "votes": votes,
            "stable": len(set(valid)) <= 1 and len(valid) == N_TRIALS,
            "sample": texts[0],
        }
        mark = "OK " if bool(maj) == gt[s] else "ERR"
        print(
            f"  [{mark}] {s:8s} GA={ga[s]:.3f} eqtl={eqtl[s]} votes={votes} -> {maj} (truth={gt[s]})"
        )

    cp = {s: claude[s]["go"] for s in gt}

    def score(preds):
        c = sum(1 for s in gt if preds[s] == gt[s])
        tp, tn, fp, fn = _conf(preds, gt)
        return {
            "n_correct": c,
            "accuracy": round(c / len(gt), 4),
            "balanced_accuracy": round(_bal_acc(tp, tn, fp, fn), 4),
            "mcc": round(_mcc(tp, tn, fp, fn), 4),
            "confusion_tp_tn_fp_fn": [tp, tn, fp, fn],
        }

    cs = score(cp)
    cs["accuracy_wilson95"] = _wilson(cs["n_correct"], len(gt))

    maj_label = sum(gt.values()) > len(gt) / 2
    baselines = {
        "majority_class": {s: maj_label for s in gt},
        "eqtl_only": {s: eqtl[s] > 0 for s in gt},
        "query_all": {s: (ga[s] >= 0.20 and eqtl[s] > 0) for s in gt},
    }
    bscores = {n: score(p) for n, p in baselines.items()}
    best_mcc = max(b["mcc"] for b in bscores.values())
    gap = round(cs["mcc"] - best_mcc, 4)
    outcome = (
        "WIN"
        if (gap >= MCC_MARGIN and cs["mcc"] >= MCC_ABS_FLOOR)
        else "LOSS"
        if cs["mcc"] < best_mcc
        else "TIE_NO_WIN"
    )

    receipt = {
        "purpose": "Day-0 Check A EXECUTION v4 — Claude open-ended adjudication vs baselines on "
        "the ENRICHED frozen set (pattern separation), MCC-primary, single stamped win rule.",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": MODEL,
        "n_trials_per_gene": N_TRIALS,
        "blinding_contract": "open-ended prompt from ONLY symbol + OT_T1D_assoc + eqtl_n_sig; "
        "no GO labels or ground-truth rule fields in the prompt.",
        "frozen_set_v2_sha256_referenced": fs["frozen_set_v2_sha256"],
        "n_genes": len(gt),
        "n_ground_truth_go": sum(gt.values()),
        "class_balance": fs["class_balance"],
        "primary_metric": "MCC",
        "pre_registered_win_condition": WIN_CONDITION_TEXT,
        "eqtl_n_sig_used": eqtl,
        "claude_score": cs,
        "baseline_scores": bscores,
        "best_baseline_mcc": best_mcc,
        "claude_mcc_gap_vs_best_baseline": gap,
        "outcome": outcome,
        "all_trials_stable": all(claude[s]["stable"] for s in claude),
        "per_gene": {
            s: {
                "ground_truth_go": gt[s],
                "claude_go": claude[s]["go"],
                "votes": claude[s]["votes"],
                "stable": claude[s]["stable"],
                "correct": claude[s]["go"] == gt[s],
            }
            for s in gt
        },
        "claude_sample_responses": {s: claude[s]["sample"] for s in claude},
        "honest_framing": "n=24 pilot with 4 positives — small but no longer collapsed to a "
        "single point. MCC is primary (raw acc uninformative under imbalance). temperature=0 "
        "means the 3 trials mainly detect API nondeterminism, not sampling variance; 'stable' "
        "is reported as such, not as strong robustness evidence. The ground truth is a "
        "deterministic non-LLM rule (immune-GO gate + OT T1D-vs-T2D + GWAS tags + measured "
        "direction); it can still be imperfect (e.g. the direction condition is a weak prior "
        "and the immune-GO gate is immune-broad, not strictly CD4+-T-cell). A WIN here is a "
        "pilot signal that free-text reasoning separates T-cell targets from confounders "
        "better than the function-blind structured baselines — NOT a general capability claim.",
    }
    h = write_hashed_json(OUT, receipt, "execution_receipt_v2_sha256")

    import hashlib

    with MANIFEST.open("a") as fh:
        fh.write(
            json.dumps(
                {
                    "path": str(OUT),
                    "sha256": hashlib.sha256(OUT.read_bytes()).hexdigest(),
                    "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "note": "day0 check A execution v4 (enriched pool, MCC-primary)",
                }
            )
            + "\n"
        )

    print(
        f"\n  claude:        n={cs['n_correct']}/{len(gt)}  bal_acc={cs['balanced_accuracy']}  MCC={cs['mcc']}"
    )
    for n, b in bscores.items():
        print(
            f"  {n:14s} n={b['n_correct']}/{len(gt)}  bal_acc={b['balanced_accuracy']}  MCC={b['mcc']}"
        )
    print(
        f"\n  win = MCC gap>={MCC_MARGIN} AND MCC>={MCC_ABS_FLOOR}   (best baseline MCC={best_mcc}, gap={gap})"
    )
    print(f">>> OUTCOME: {outcome}")
    print(f"hash {h[:16]}  receipt -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
