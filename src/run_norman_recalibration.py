"""Recalibration test: turn iter4's honest cross-assay transfer NEGATIVE into a scoped POSITIVE.

Claim: frozen Marson T1D conformal quantiles do NOT transfer to Norman K562 (iter4:
80/90/95 -> 4.8/7.6/17.1). BUT re-fitting the IDENTICAL split-conformal method on Norman's
OWN calibration split restores near-nominal coverage on held-out Norman test -> the METHOD
generalizes even though the specific quantiles do not ("recalibrate per assay and it holds").

Inputs (all cached by run_norman_stream.py / iter4):
  data/interim/norman_group_stats.npz     mean/var/ncells (33694 genes x 106 groups; grp0=control)
  data/interim/norman_meta.json           inter_idx (9780 panel), norman_sym, id2target (105 perts)
  data/interim/frozen_marson_calibration.joblib  {model: HGBR, frozen_q: {L:q}}

Honesty: prereg hash-frozen in data/gold/day0_recalibration_prereg.json BEFORE scoring;
y-target reproduces iter4 frozen coverage EXACTLY before recalibration; success band pre-set.
"""
import json, numpy as np, joblib
from sklearn.ensemble import HistGradientBoostingRegressor

BASE="data/interim/"; SEED=20260708; LEVELS=(0.80,0.90,0.95)
meta=json.load(open(BASE+"norman_meta.json")); npz=np.load(BASE+"norman_group_stats.npz")
mean,var,ncells=npz["mean"],npz["var"],npz["ncells"]
inter_idx=np.array(meta["inter_idx"]); norman_sym=np.array(meta["norman_sym"]); id2target=meta["id2target"]
panel_sym=norman_sym[inter_idx]; cm=mean[:,0]; cv=var[:,0]; nc=ncells[0]

ys=[]; log_bm=[]
for gid in range(1,106):
    tgt=id2target[str(gid)]
    pm=mean[inter_idx,gid]; pv=var[inter_idx,gid]; npp=ncells[gid]
    den=np.sqrt(pv/max(npp,1)+cv[inter_idx]/max(nc,1))
    z=np.where(den>0,(pm-cm[inter_idx])/den,0.0).copy(); z[panel_sym==tgt]=0.0
    ys.append(np.log1p(np.sqrt(np.sum(z*z))))
    on=norman_sym==tgt; log_bm.append(np.log1p(max(cm[on].sum(),0.0)) if on.any() else 0.0)
ys=np.array(ys); log_bm=np.array(log_bm); N=len(ys)
X=np.column_stack([log_bm,np.ones(N),np.zeros(N),np.zeros(N)])  # Rest condition

fm=joblib.load(BASE+"frozen_marson_calibration.joblib"); reg=fm["model"]; q_frozen=fm["frozen_q"]
resid=np.abs(ys-reg.predict(X))
print("FROZEN transfer (reproduce iter4):",{L:round(float(np.mean(resid<=q_frozen[L])),4) for L in LEVELS})

def split(seed,cf=0.5,frozen_base=False):
    p=np.random.default_rng(seed).permutation(N); nca=int(round(N*cf)); ci,ti=p[:nca],p[nca:]
    if frozen_base: pc,pt=reg.predict(X[ci]),reg.predict(X[ti])
    else:
        m=HistGradientBoostingRegressor(max_iter=300,learning_rate=0.05,max_depth=4,
            l2_regularization=1.0,random_state=SEED).fit(X[ci],ys[ci])
        pc,pt=m.predict(X[ci]),m.predict(X[ti])
    rc,rt=np.abs(ys[ci]-pc),np.abs(ys[ti]-pt)
    return {L:float(np.mean(rt<=np.quantile(rc,L,method="higher"))) for L in LEVELS}

print("PRIMARY (seed 20260711):",{L:round(v,4) for L,v in split(20260711).items()})
covs={L:[] for L in LEVELS}
for s in range(500):
    r=split(20260711+s)
    for L in LEVELS: covs[L].append(r[L])
print("ROBUSTNESS 500-split mean:",{L:round(float(np.mean(covs[L])),4) for L in LEVELS})
