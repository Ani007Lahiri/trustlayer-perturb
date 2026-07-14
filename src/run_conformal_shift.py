"""Runner: naive-LODO vs Mondrian vs weighted covariate-shift conformal.
Writes data/gold/day0_mondrian_weighted_receipt.json (hash-stamped)."""
import json, hashlib, sys, time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.conformal_shift import run_shift_comparison

t0 = time.time()
res = run_shift_comparison()
levels = res.levels

def cov_arr(d):
    return [d[L]["row_empirical"] for L in levels]

naive = cov_arr(res.naive)
mond = cov_arr(res.mondrian)
wt = cov_arr(res.weighted)

# deviation from nominal (lower = better)
def dev(cov):
    return [round(abs(cov[i] - levels[i]), 4) for i in range(len(levels))]

receipt = {
    "check": "Principled conformal under donor covariate shift (naive-LODO vs Mondrian vs weighted)",
    "generated_utc": datetime.now(timezone.utc).isoformat(),
    "prereg_sha256": "b7743cae256bc859bbbb192f248e8bf0c90f65bd57c102ae84a4ae5309962509",
    "levels": levels,
    "coverage": {
        "naive_LODO": {str(L): res.naive[L] for L in levels},
        "mondrian_by_condition": {str(L): res.mondrian[L] for L in levels},
        "weighted_covariate_shift": {str(L): res.weighted[L] for L in levels},
    },
    "coverage_arrays": {"naive": naive, "mondrian": mond, "weighted": wt, "nominal": levels},
    "abs_deviation_from_nominal": {"naive": dev(naive), "mondrian": dev(mond), "weighted": dev(wt)},
    "domain_shift": res.domain_shift,
    "per_fold": res.per_fold,
    "guarantee_statements": {
        "naive_LODO": "NO valid conformal guarantee: donor-blocking makes calibration and test "
        "non-exchangeable (different donor composition). Coverage is empirical only.",
        "mondrian_by_condition": "Valid marginal coverage 1-alpha IF cal u test are exchangeable "
        "WITHIN each culture_condition stratum (donor does not shift residuals conditional on condition). "
        "Strictly weaker assumption than global exchangeability.",
        "weighted_covariate_shift": "Valid marginal coverage 1-alpha under PURE COVARIATE SHIFT "
        "P_test(Y|X)=P_cal(Y|X) with correct density ratio w(x) (Tibshirani et al. 2019). "
        "This is the defensible guarantee statement the naive version lacks.",
    },
    "wall_seconds": round(time.time() - t0, 1),
}

payload = json.dumps(receipt, indent=2, sort_keys=True, default=str)
receipt["_sha256_self"] = hashlib.sha256(payload.encode()).hexdigest()
out = Path(__file__).resolve().parent.parent / "data/gold/day0_mondrian_weighted_receipt.json"
with open(out, "w") as fh:
    json.dump(receipt, fh, indent=2, default=str)

print("=== COVERAGE (row empirical) ===")
print("level   naive   mondrian  weighted")
for i, L in enumerate(levels):
    print(f"{L:.2f}   {naive[i]:.4f}  {mond[i]:.4f}   {wt[i]:.4f}")
print("domain classifier AUC (mean):", res.domain_shift["domain_classifier_auc_mean"])
print("abs dev naive   :", dev(naive))
print("abs dev mondrian:", dev(mond))
print("abs dev weighted:", dev(wt))
print("sha256:", receipt["_sha256_self"])
print("wrote:", out)
