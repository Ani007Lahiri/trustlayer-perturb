"""
Calibration & Selective-Prediction Benchmark for Perturbation Prediction (minimal public v1).

The point no existing perturbation benchmark (OP3, PerturBench, PertEval-scFM, VCBench, Arc VCC)
reports: per-prediction CALIBRATION VALIDITY and SELECTIVE-PREDICTION risk-coverage, as the primary
axis, across multiple predictors and multiple datasets — the "does the trust gate wrap ANY predictor"
claim made runnable.

Machinery is standard and cited (NOT novel): split conformal (Vovk; Angelopoulos & Bates 2021),
locally-adaptive/normalized nonconformity (Papadopoulos 2008; Lei et al. 2018).

Run:  python benchmark/calibration_benchmark.py
Wrap your own predictor: implement fit(X,y)->self / predict(X)->yhat and add it to PREDICTORS.
"""
import json, hashlib, numpy as np, pandas as pd
from pathlib import Path
from sklearn.linear_model import Ridge
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from scipy.stats import spearmanr, beta

SEED = 20260708
ROOT = Path(__file__).resolve().parents[1]
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))  # numpy 2.x renamed trapz->trapezoid
LEVELS = (0.80, 0.90, 0.95)
SWEEP = np.round(np.arange(0.05, 0.96, 0.05), 2)   # 19-level calibration sweep


# ---------- predictors (model-agnostic: anything with fit/predict works) ----------
def make_predictors():
    return {
        "mean_baseline": _MeanBaseline(),
        "linear":        Ridge(alpha=1.0, random_state=SEED) if False else Ridge(alpha=1.0),
        "hist_gbr":      HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05,
                                                       max_depth=3, random_state=SEED),
        "random_forest": RandomForestRegressor(n_estimators=200, max_depth=8,
                                               random_state=SEED, n_jobs=-1),
    }

class _MeanBaseline:
    def fit(self, X, y): self.m_ = float(np.mean(y)); return self
    def predict(self, X): return np.full(len(X), self.m_)


# ---------- conformal machinery (standard, cited) ----------
def split_conformal_q(calib_res, alpha):
    """Marginal split-conformal quantile (constant-width). Angelopoulos & Bates 2021."""
    n = len(calib_res)
    k = int(np.ceil((n + 1) * (1 - alpha)))
    k = min(k, n)
    return np.sort(calib_res)[k - 1]

def local_sigma(Xc, res_c, Xt):
    """Locally-adaptive difficulty sigma(x): auxiliary GBR fit of |resid| ~ features.
    Papadopoulos 2008 / Lei 2018 normalized nonconformity. Model-agnostic wrt base predictor."""
    if Xc.shape[1] == 0 or np.allclose(Xc, Xc[0]):
        return np.ones(len(Xt)), np.ones(len(Xc))
    aux = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.05, max_depth=3,
                                        random_state=SEED)
    aux.fit(Xc, np.abs(res_c))
    eps = 1e-6
    sc = np.clip(aux.predict(Xc), eps, None)
    st = np.clip(aux.predict(Xt), eps, None)
    return st, sc

def clopper_pearson(k, n, alpha=0.05):
    lo = 0.0 if k == 0 else beta.ppf(alpha/2, k, n-k+1)
    hi = 1.0 if k == n else beta.ppf(1-alpha/2, k+1, n-k)
    return lo, hi


# ---------- metrics ----------
def coverage_validity(y_t, pred_t, res_c):
    out = {}
    for a_level in SWEEP:
        q = split_conformal_q(res_c, 1 - a_level)  # a_level = nominal coverage
        emp = float(np.mean(np.abs(y_t - pred_t) <= q))
        out[float(a_level)] = emp
    cal_err = float(np.mean([abs(out[l] - l) for l in out]))
    headline = {}
    n = len(y_t)
    for lv in LEVELS:
        q = split_conformal_q(res_c, 1 - lv)
        cov = np.abs(y_t - pred_t) <= q
        emp = float(np.mean(cov)); lo, hi = clopper_pearson(int(cov.sum()), n)
        headline[lv] = {"empirical": round(emp, 4), "nominal": lv,
                        "gap": round(emp - lv, 4), "cp_lo": round(lo, 4), "cp_hi": round(hi, 4),
                        "within_tol": bool(abs(emp - lv) <= 0.03 or (lo <= lv <= hi))}
    return {"sweep": out, "calibration_error": round(cal_err, 4), "headline": headline}

def risk_coverage(y_t, pred_t, sigma_t, n_null=100, rng=None):
    """RMSE on most-confident retained fraction, ranked by locally-adaptive sigma(x)."""
    rng = rng or np.random.default_rng(SEED)
    err = np.abs(y_t - pred_t)
    order = np.argsort(sigma_t)                 # most confident first
    fracs = np.round(np.arange(0.1, 1.01, 0.1), 2)
    curve = {}
    for f in fracs:
        k = max(1, int(round(f * len(err))))
        curve[float(f)] = float(np.sqrt(np.mean(err[order[:k]] ** 2)))
    aurc = float(_trapz([curve[f] for f in fracs], fracs))
    nulls = []
    for _ in range(n_null):
        o = rng.permutation(len(err))
        pts = []
        for f in fracs:
            k = max(1, int(round(f * len(err))))
            pts.append(np.sqrt(np.mean(err[o[:k]] ** 2)))
        nulls.append(_trapz(pts, fracs))
    null_mean = float(np.mean(nulls)); null_std = float(np.std(nulls))
    return {"curve": curve, "aurc": round(aurc, 5),
            "random_null_aurc": round(null_mean, 5), "random_null_std": round(null_std, 5),
            "beats_null": bool(aurc < null_mean)}

def point_metrics(y_t, pred_t):
    rmse = float(np.sqrt(np.mean((y_t - pred_t) ** 2)))
    ss = float(1 - np.sum((y_t - pred_t)**2) / np.sum((y_t - np.mean(y_t))**2))
    sp = float(spearmanr(pred_t, y_t).statistic) if np.std(pred_t) > 0 else 0.0
    return {"rmse": round(rmse, 4), "r2": round(ss, 4), "spearman": round(sp, 4)}


# ---------- data loaders ----------
def load_gladstone():
    ft = pd.read_parquet(ROOT / "data/interim/day2_feature_table.parquet").copy()
    # Canonical LEAKAGE-SAFE feature set (src/trustlayer/features.py::feature_cols).
    # trans_effect_magnitude is EXCLUDED — it is a transform of the target y=log1p(trans_effect_magnitude).
    feats = ["log_baseMean", "log_n_cells",
             "cond_Rest", "cond_Stim8hr", "cond_Stim48hr", "baseMean_missing"]
    # log_baseMean is NaN for genes with undetectable baseline (log of ~0). That missingness is
    # itself informative: add an explicit indicator and fill with a below-minimum sentinel using
    # TRAIN statistics only (no leakage), rather than dropping rows.
    tr_min = float(ft.loc[ft["fold"] == "train", "log_baseMean"].min())
    sentinel = tr_min - 1.0
    ft["baseMean_missing"] = ft["log_baseMean"].isna().astype(float)
    ft["log_baseMean"] = ft["log_baseMean"].fillna(sentinel)
    d = {}
    for split in ("train", "calib", "test"):
        s = ft[ft["fold"] == split]
        d[split] = (s[feats].to_numpy(float), s["y"].to_numpy(float))
    return d, feats

def load_norman():
    X = np.load(ROOT / "handoff/norman_X.npy")[:, :1]   # only log_baseMean is real
    y = np.load(ROOT / "handoff/norman_ys.npy")
    return X, y, ["log_baseMean"]

def load_datlinger():
    """Datlinger 2017 CROP-seq Jurkat T-cell (GSE92872). PUBLIC, small-n (64 rows: 32 perts x 2 conds).
    Same target y=log1p(||z_trans||) as Gladstone/Norman."""
    df = pd.read_parquet(ROOT / "data/interim/datlinger_feature_table.parquet")
    # log_n_cells is DROPPED: at this small n it correlates -0.95 with y because the target
    # ||z_trans|| is built from Welch z-scores that scale with 1/sqrt(n_cells). That is a
    # construction artifact of the target, not perturbation signal — including it is leakage.
    feats = ["log_baseMean", "cond_unstim", "cond_stim"]
    X = df[feats].to_numpy(float)
    y = df["y"].to_numpy(float)
    return X, y, feats, df


# ---------- relational (GEARS-style) predictor ----------
class RelationalPredictor:
    """Predicts a held-out perturbation's effect as a similarity-weighted average of RELATED
    training perturbations. Relatedness = target-gene co-expression in CONTROL (unperturbed) cells
    (leakage-free: baseline data, never the perturbation response). This is GEARS' inductive bias
    ('generalize from related perturbations') — NOT GEARS itself (no GNN / GO graph / cell-level training).
    Interface takes perturbation labels + a similarity matrix, evaluated leave-one-perturbation-out."""
    def __init__(self, sim_df, k=8, temp=0.05):
        self.sim = sim_df; self.k = k; self.temp = temp
    def predict_lopo(self, genes, conds, y):
        """Leave-one-PERTURBATION-out: for each row, predict from OTHER perturbations' y,
        weighted by co-expression similarity of target genes (same-condition rows preferred)."""
        genes = np.asarray(genes); conds = np.asarray(conds); y = np.asarray(y, float)
        preds = np.full(len(y), np.nan)
        for i in range(len(y)):
            g_i, c_i = genes[i], conds[i]
            # candidate training rows: different PERTURBATION (never same gene -> no leakage)
            mask = (genes != g_i)
            if mask.sum() == 0: continue
            w = np.zeros(mask.sum()); yj = y[mask]; gj = genes[mask]; cj = conds[mask]
            for j, (gg, cc) in enumerate(zip(gj, cj)):
                s = self.sim.loc[g_i, gg] if (g_i in self.sim.index and gg in self.sim.columns) else 0.0
                cond_bonus = 1.0 if cc == c_i else 0.5   # prefer same condition
                w[j] = max(s, 0.0) * cond_bonus
            if w.sum() == 0:                              # no related pert -> fall back to global mean
                preds[i] = yj.mean(); continue
            # softmax over top-k by weight
            order = np.argsort(w)[::-1][:self.k]
            wk = w[order]; wk = np.exp(wk / self.temp); wk /= wk.sum()
            preds[i] = float(np.dot(wk, yj[order]))
        return preds


# ---------- runners ----------
def run_predictor_fixed_split(Xtr, ytr, Xc, yc, Xt, yt, model):
    model.fit(Xtr, ytr)
    pc, pt = model.predict(Xc), model.predict(Xt)
    res_c = np.abs(yc - pc)
    st, sc = local_sigma(Xc, yc - pc, Xt)
    return {"coverage": coverage_validity(yt, pt, res_c),
            "selective": risk_coverage(yt, pt, st),
            "point": point_metrics(yt, pt)}

def run_gladstone():
    d, _ = load_gladstone()
    (Xtr, ytr), (Xc, yc), (Xt, yt) = d["train"], d["calib"], d["test"]
    res = {}
    for name, model in make_predictors().items():
        res[name] = run_predictor_fixed_split(Xtr, ytr, Xc, yc, Xt, yt, model)
    return res

def run_norman(n_repeats=200):
    X, y, _ = load_norman()
    rng = np.random.default_rng(SEED)
    n = len(y); half = n // 2
    agg = {name: {"cov": {lv: [] for lv in LEVELS}, "calerr": [], "aurc": [], "null": [],
                  "rmse": [], "r2": [], "spearman": []} for name in make_predictors()}
    for _ in range(n_repeats):
        idx = rng.permutation(n); ci, ti = idx[:half], idx[half:]
        for name, model in make_predictors().items():
            model.fit(X[ci], y[ci])
            pc, pt = model.predict(X[ci]), model.predict(X[ti])
            res_c = np.abs(y[ci] - pc)
            cov = coverage_validity(y[ti], pt, res_c)
            for lv in LEVELS: agg[name]["cov"][lv].append(cov["headline"][lv]["empirical"])
            agg[name]["calerr"].append(cov["calibration_error"])
            st, _ = local_sigma(X[ci], y[ci]-pc, X[ti])
            sel = risk_coverage(y[ti], pt, st, n_null=30, rng=rng)
            agg[name]["aurc"].append(sel["aurc"]); agg[name]["null"].append(sel["random_null_aurc"])
            pm = point_metrics(y[ti], pt)
            for k in ("rmse","r2","spearman"): agg[name][k].append(pm[k])
    out = {}
    for name, a in agg.items():
        out[name] = {
            "coverage": {"headline": {lv: {"empirical": round(np.mean(a["cov"][lv]),4),
                                           "nominal": lv, "gap": round(np.mean(a["cov"][lv])-lv,4),
                                           "within_tol": bool(abs(np.mean(a["cov"][lv])-lv)<=0.05)}
                                      for lv in LEVELS},
                         "calibration_error": round(float(np.mean(a["calerr"])),4)},
            "selective": {"aurc": round(float(np.mean(a["aurc"])),5),
                          "random_null_aurc": round(float(np.mean(a["null"])),5),
                          "beats_null": bool(np.mean(a["aurc"]) < np.mean(a["null"]))},
            "point": {"rmse": round(float(np.mean(a["rmse"])),4),
                      "r2": round(float(np.mean(a["r2"])),4),
                      "spearman": round(float(np.mean(a["spearman"])),4)},
            "n_repeats": n_repeats, "n_total": int(n)}
    return out

def run_smalln(X, y, n_repeats=200, tol=0.05):
    """Generic repeated random-split runner for small-n datasets (Norman, Datlinger)."""
    rng = np.random.default_rng(SEED)
    n = len(y); half = n // 2
    agg = {name: {"cov": {lv: [] for lv in LEVELS}, "calerr": [], "aurc": [], "null": [],
                  "rmse": [], "r2": [], "spearman": []} for name in make_predictors()}
    for _ in range(n_repeats):
        idx = rng.permutation(n); ci, ti = idx[:half], idx[half:]
        for name, model in make_predictors().items():
            model.fit(X[ci], y[ci])
            pc, pt = model.predict(X[ci]), model.predict(X[ti])
            res_c = np.abs(y[ci] - pc)
            cov = coverage_validity(y[ti], pt, res_c)
            for lv in LEVELS: agg[name]["cov"][lv].append(cov["headline"][lv]["empirical"])
            agg[name]["calerr"].append(cov["calibration_error"])
            st, _ = local_sigma(X[ci], y[ci]-pc, X[ti])
            sel = risk_coverage(y[ti], pt, st, n_null=30, rng=rng)
            agg[name]["aurc"].append(sel["aurc"]); agg[name]["null"].append(sel["random_null_aurc"])
            pm = point_metrics(y[ti], pt)
            for k in ("rmse","r2","spearman"): agg[name][k].append(pm[k])
    return _agg_out(agg, n_repeats, n, tol)

def _agg_out(agg, n_repeats, n, tol):
    out = {}
    for name, a in agg.items():
        out[name] = {
            "coverage": {"headline": {lv: {"empirical": round(np.mean(a["cov"][lv]),4),
                                           "nominal": lv, "gap": round(np.mean(a["cov"][lv])-lv,4),
                                           "within_tol": bool(abs(np.mean(a["cov"][lv])-lv)<=tol)}
                                      for lv in LEVELS},
                         "calibration_error": round(float(np.mean(a["calerr"])),4)},
            "selective": {"aurc": round(float(np.mean(a["aurc"])),5),
                          "random_null_aurc": round(float(np.mean(a["null"])),5),
                          "beats_null": bool(np.mean(a["aurc"]) < np.mean(a["null"]))},
            "point": {"rmse": round(float(np.mean(a["rmse"])),4),
                      "r2": round(float(np.mean(a["r2"])),4),
                      "spearman": round(float(np.mean(a["spearman"])),4)},
            "n_repeats": n_repeats, "n_total": int(n)}
    return out

def run_datlinger(n_repeats=200):
    X, y, _, _ = load_datlinger()
    return run_smalln(X, y, n_repeats=n_repeats, tol=0.05)

def run_relational(n_repeats=200):
    """RelationalPredictor evaluated on Datlinger. LOPO predictions -> conformal via repeated
    calib/test splits of the resulting (y,pred) pairs. sigma(x) from features, same as others."""
    X, y, feats, df = load_datlinger()
    sim = pd.read_parquet(ROOT / "data/interim/datlinger_coexpr_sim.parquet")
    model = RelationalPredictor(sim, k=8, temp=0.05)
    pred = model.predict_lopo(df["gene"].to_numpy(), df["condition"].to_numpy(), y)
    ok = ~np.isnan(pred)
    Xo, yo, po = X[ok], y[ok], pred[ok]
    rng = np.random.default_rng(SEED); n = len(yo); half = n // 2
    agg = {"cov": {lv: [] for lv in LEVELS}, "calerr": [], "aurc": [], "null": [],
           "rmse": [], "r2": [], "spearman": []}
    for _ in range(n_repeats):
        idx = rng.permutation(n); ci, ti = idx[:half], idx[half:]
        res_c = np.abs(yo[ci] - po[ci])
        cov = coverage_validity(yo[ti], po[ti], res_c)
        for lv in LEVELS: agg["cov"][lv].append(cov["headline"][lv]["empirical"])
        agg["calerr"].append(cov["calibration_error"])
        st, _ = local_sigma(Xo[ci], yo[ci]-po[ci], Xo[ti])
        sel = risk_coverage(yo[ti], po[ti], st, n_null=30, rng=rng)
        agg["aurc"].append(sel["aurc"]); agg["null"].append(sel["random_null_aurc"])
        pm = point_metrics(yo[ti], po[ti])
        for k in ("rmse","r2","spearman"): agg[k].append(pm[k])
    res = _agg_out({"relational_gears_style": agg}, n_repeats, n, tol=0.05)
    res["relational_gears_style"]["note"] = "GEARS-style relational (co-expression-weighted, LOPO); NOT GEARS (no GNN)."
    return res["relational_gears_style"]

def main():
    results = {"Gladstone_CD4T_T1D": run_gladstone(),
               "Norman2019_K562": run_norman(),
               "Datlinger2017_Jurkat_Tcell": run_datlinger()}
    # 5th predictor is dataset-specific (needs the co-expression graph + LOPO) -> attach to Datlinger
    results["Datlinger2017_Jurkat_Tcell"]["relational_gears_style"] = run_relational()
    prereg = json.load(open(ROOT / "benchmark/benchmark_prereg.json"))
    receipt = {"title": prereg["title"], "prereg_sha256": prereg["prereg_sha256"],
               "seed": SEED, "levels": list(LEVELS), "results": results}
    blob = json.dumps(receipt, sort_keys=True, separators=(",", ":")).encode()
    receipt["receipt_sha256"] = hashlib.sha256(blob).hexdigest()
    json.dump(receipt, open(ROOT / "benchmark/benchmark_receipt.json", "w"), indent=1)
    return receipt

if __name__ == "__main__":
    r = main()
    print("receipt_sha256:", r["receipt_sha256"])
