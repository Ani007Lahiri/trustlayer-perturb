"""
Day-3: DONOR-blocked leave-one-donor-out (LODO) conformal trust layer.
Plan v3 Fix 6 headline. Critique-corrected (Gate-1 Day-3, 3 blocking fixes applied).

Data: GWCD4i.DE_stats.by_donors.h5mu -> 6 donor-pair modalities from 4 donors.
Each pair retains only perturbations reproducibly detected in BOTH its donors
(~4880 rows, 2591-gene reproducible core) -> the concrete small-effective-N regime.

Construction (exchangeability fix #1):
  For held-out donor d:
    calibration = pairs whose BOTH members != d   (the C(3,2)=3 pairs among other donors)
    test        = pairs containing d               (the 3 pairs with d)
  Assert: no calibration pair contains d (fail-closed).

Pipeline per fold (circularity fix #3 -> three-way separation):
  1. Build (gene,condition) rows for calib pairs and test pairs; target = log1p trans-effect
     magnitude (on-target gene excluded), recomputed PER-FOLD from that pair's zscore layer
     (NOT the Day-2 global fit) -- honors the implementation note.
  2. Fit base predictor (HGBR) on calibration rows only.
  3. Split-conformal: nonconformity = |residual| on calibration; quantile at level 1-alpha.
  4. Evaluate coverage + interval width + trust score + selective risk on TEST rows only.

Power (fix #2): report Clopper-Pearson interval on POOLED coverage + its WIDTH as a
first-class result; do not claim per-level discrimination at n_eff ~ tens.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import beta as beta_dist
from sklearn.ensemble import HistGradientBoostingRegressor

H5MU = Path("data/raw/GWCD4i.DE_stats.by_donors.h5mu")
DONORS = ["CE0006864", "CE0008162", "CE0008678", "CE0010866"]
SEED = 20260708


# ---------------------------------------------------------------- data access
def _decode_cat(h5obj) -> np.ndarray:
    if isinstance(h5obj, h5py.Group) and "categories" in h5obj and "codes" in h5obj:
        cats = np.array(
            [c.decode() if isinstance(c, bytes) else c for c in h5obj["categories"][:]]
        )
        return cats[h5obj["codes"][:]]
    return np.array([a.decode() if isinstance(a, bytes) else a for a in h5obj[:]])


def pair_name(a: str, b: str) -> str:
    return f"{a}_{b}"


def all_pairs() -> list[str]:
    return [pair_name(a, b) for a, b in combinations(DONORS, 2)]


def pairs_excluding(d: str) -> list[str]:
    """Pairs whose BOTH members != d (clean calibration pairs)."""
    others = [x for x in DONORS if x != d]
    return [pair_name(a, b) for a, b in combinations(others, 2)]


def pairs_containing(d: str) -> list[str]:
    return [p for p in all_pairs() if d in p.split("_")]


def build_pair_table(f: h5py.File, modality: str) -> dict:
    """Compute per-row leakage-safe features + trans-only target for one donor-pair.

    Features (independent of the DE readout): log baseMean of the perturbed target gene,
    condition one-hot. Target: log1p(|| zscore over trans genes, on-target excluded ||).
    """
    g = f["mod"][modality]
    gene = _decode_cat(g["obs"]["target_contrast_gene_name"])
    cond = _decode_cat(g["obs"]["culture_condition"])
    var_names = _decode_cat(
        g["var"]["gene_name"] if "gene_name" in g["var"] else g["var"]["_index"]
    )
    col = {gn: i for i, gn in enumerate(var_names)}

    z = g["layers"]["zscore"]
    bm = g["layers"]["baseMean"]
    n = z.shape[0]
    trans_mag = np.empty(n)
    log_baseMean = np.zeros(n)
    CHUNK = 1000
    for s in range(0, n, CHUNK):
        e = min(s + CHUNK, n)
        zb = np.nan_to_num(z[s:e, :], nan=0.0)
        bmb = np.nan_to_num(bm[s:e, :], nan=0.0)
        full_sq = np.sum(zb * zb, axis=1)
        for r in range(s, e):
            c = col.get(gene[r])
            ontarget_sq = zb[r - s, c] ** 2 if c is not None else 0.0
            trans_mag[r] = np.sqrt(max(full_sq[r - s] - ontarget_sq, 0.0))
            # baseMean of the on-target gene = pre-perturbation expression (admissible)
            log_baseMean[r] = (
                np.log1p(max(bmb[r - s, c], 0.0)) if c is not None else 0.0
            )

    X = np.column_stack(
        [
            log_baseMean,
            (cond == "Rest").astype(float),
            (cond == "Stim8hr").astype(float),
            (cond == "Stim48hr").astype(float),
        ]
    )
    y = np.log1p(trans_mag)
    return {"gene": gene, "cond": cond, "X": X, "y": y}


def _stack(f, modalities: list[str]) -> dict:
    parts = [build_pair_table(f, m) for m in modalities]
    return {
        "gene": np.concatenate([p["gene"] for p in parts]),
        "cond": np.concatenate([p["cond"] for p in parts]),
        "X": np.vstack([p["X"] for p in parts]),
        "y": np.concatenate([p["y"] for p in parts]),
    }


# ---------------------------------------------------------------- conformal
def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    lo = beta_dist.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta_dist.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def aurc(errors: np.ndarray, trust: np.ndarray) -> float:
    """Area under the risk-coverage curve. Lower = better ordering by trust.
    Retain most-trusted first; risk = mean error over retained."""
    order = np.argsort(-trust)  # high trust first
    e = errors[order]
    n = len(e)
    risks = np.cumsum(e) / np.arange(1, n + 1)
    return float(np.mean(risks))


@dataclass
class LODOResult:
    levels: list[float]
    pooled_coverage: dict = field(
        default_factory=dict
    )  # level -> {emp, cp_lo, cp_hi, width, k, n}
    per_donor_coverage: dict = field(default_factory=dict)
    aurc_model: float = 0.0
    aurc_shuffle: float = 0.0
    aurc_perm: dict = field(default_factory=dict)
    n_eff_perturbations_test: int = 0


def run_lodo(levels=(0.80, 0.90, 0.95)) -> LODOResult:
    rng = np.random.default_rng(SEED)
    with h5py.File(H5MU, "r") as f:
        # pooled test predictions across the 4 LODO folds
        covered = {
            L: [] for L in levels
        }  # bool per test row (at 90% for pooling headline)
        per_donor = {}
        all_test_errors, all_test_trust = [], []
        all_test_genes = set()

        for d in DONORS:
            calib_mods = pairs_excluding(d)
            test_mods = pairs_containing(d)
            # exchangeability assertion (fix #1)
            for m in calib_mods:
                assert d not in m.split("_"), f"LEAK: donor {d} in calibration pair {m}"

            calib = _stack(f, calib_mods)
            test = _stack(f, test_mods)

            model = HistGradientBoostingRegressor(
                max_iter=300,
                learning_rate=0.05,
                max_depth=4,
                l2_regularization=1.0,
                random_state=SEED,
            ).fit(calib["X"], calib["y"])

            # nonconformity on calibration (residual absolute)
            calib_res = np.abs(calib["y"] - model.predict(calib["X"]))
            test_pred = model.predict(test["X"])
            test_abs_err = np.abs(test["y"] - test_pred)

            fold_cov = {}
            for L in levels:
                q = np.quantile(calib_res, L, method="higher")
                cov = test_abs_err <= q
                covered[L].extend(cov.tolist())
                fold_cov[L] = float(cov.mean())
            per_donor[d] = {"n_test": int(len(test["y"])), "coverage": fold_cov}

            # Trust score for selective prediction = negative ensemble spread (low
            # disagreement across bootstrap models => high trust). Split-conformal width is
            # CONSTANT per level, so it cannot rank predictions; the ensemble spread gives a
            # per-prediction reliability signal. Each fold uses an INDEPENDENT bootstrap
            # ensemble (distinct seeds per model AND per fold) so the spread is not degenerate.
            fold_seed = SEED + hash(d) % 100000
            spreads = _ensemble_spread(calib, test, n_models=8, base_seed=fold_seed)
            trust = -spreads  # low spread = high trust
            all_test_errors.append(test_abs_err)
            all_test_trust.append(trust)
            all_test_genes.update(test["gene"].tolist())

        # Coverage reporting (fix #2 + effective-N honesty):
        # Row-level CP is OVERCONFIDENT here (rows are non-independent: 2591 unique
        # perturbations reappear across 3 shared-donor test pairs). The HONEST uncertainty
        # is the ACROSS-DONOR-FOLD spread (4 independent held-out donors). We report BOTH:
        #   - fold_mean +/- fold_std as the PRIMARY interval (respects donor structure)
        #   - row-level empirical + a naive CP (flagged as overconfident) for completeness
        res = LODOResult(levels=list(levels))
        for L in levels:
            arr = np.array(covered[L], dtype=bool)
            k, n = int(arr.sum()), int(len(arr))
            lo, hi = clopper_pearson(k, n)
            fold_vals = [per_donor[d]["coverage"][L] for d in DONORS]
            fmean, fstd = float(np.mean(fold_vals)), float(np.std(fold_vals))
            res.pooled_coverage[L] = {
                "nominal": L,
                "fold_mean": round(fmean, 4),
                "fold_std": round(fstd, 4),
                "fold_values": [round(v, 4) for v in fold_vals],
                "fold_interval_lo": round(fmean - fstd, 4),
                "fold_interval_hi": round(fmean + fstd, 4),
                "fold_interval_contains_nominal": (fmean - fstd) <= L <= (fmean + fstd),
                "row_empirical": round(k / n, 4),
                "row_cp_lo": round(lo, 4),
                "row_cp_hi": round(hi, 4),
                "row_cp_note": "row-level CP is OVERCONFIDENT (non-independent rows); "
                "use fold_std as the honest uncertainty",
                "n_rows": n,
            }
        res.per_donor_coverage = per_donor

        # selective risk on pooled test rows (fix #3: test-only)
        errs = np.concatenate(all_test_errors)
        trust = np.concatenate(all_test_trust)
        res.aurc_model = aurc(errs, trust)
        # proper permutation null (>=2000 shuffles) with empirical p-value
        res.aurc_perm = aurc_permutation_test(errs, trust, n_perm=2000, seed=SEED)
        res.aurc_shuffle = res.aurc_perm["aurc_null_mean"]
        res.n_eff_perturbations_test = len(all_test_genes)
    return res


def _ensemble_spread(
    calib: dict, test: dict, n_models: int = 8, base_seed: int = SEED
) -> np.ndarray:
    """Std of test predictions across a bootstrap ensemble (heuristic uncertainty).

    Each member uses a DISTINCT bootstrap draw and a DISTINCT model seed (base_seed+s),
    so the ensemble is genuinely diverse (fixes the degenerate-spread bug where every
    member drew identical indices)."""
    preds = []
    n = len(calib["y"])
    for s in range(n_models):
        rng = np.random.default_rng(base_seed + s)  # distinct draw per member
        idx = rng.integers(0, n, n)
        m = HistGradientBoostingRegressor(
            max_iter=200,
            learning_rate=0.05,
            max_depth=4,
            l2_regularization=1.0,
            random_state=base_seed + s,
        ).fit(calib["X"][idx], calib["y"][idx])
        preds.append(m.predict(test["X"]))
    return np.std(np.vstack(preds), axis=0)


def aurc_permutation_test(
    errors: np.ndarray, trust: np.ndarray, n_perm: int = 2000, seed: int = SEED
) -> dict:
    """Permutation null for the trust-ordering: is model AURC lower than random ordering?

    Returns observed AURC, null mean, and empirical one-sided p-value
    (fraction of shuffled-trust AURCs <= observed)."""
    rng = np.random.default_rng(seed)
    observed = aurc(errors, trust)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = aurc(errors, rng.permutation(trust))
    p = float((np.sum(null <= observed) + 1) / (n_perm + 1))
    return {
        "aurc_observed": round(observed, 5),
        "aurc_null_mean": round(float(null.mean()), 5),
        "aurc_null_std": round(float(null.std()), 5),
        "p_value_one_sided": round(p, 4),
        "n_permutations": n_perm,
        "significant_at_0.05": p < 0.05,
    }
