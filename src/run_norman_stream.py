"""Stream GSE133344 (Norman 2019 K562 Perturb-seq) mtx.gz -> per-group log-CP10K mean/var.
Two passes with pandas C parser. Writes data/interim/norman_group_stats.npz.
Normalization: CP10K then log1p, matching a standard scRNA workflow; DE = Welch z of
(pert mean - control mean) so the frozen Marson conformal target y=log1p(||z_trans||) can be
recomputed on Norman with the SAME definition.
"""
import numpy as np, pandas as pd, gzip, os
os.chdir("/Users/anirudhlahiri/Desktop/Hackathon_Brainstorm")
MTX="data/external/norman2019/GSE133344_filtered_matrix.mtx.gz"
group_by_col=np.load("data/interim/norman_group_by_col.npy")
ncells_group=np.load("data/interim/norman_ncells_group.npy")
ngroups=len(ncells_group); NGENE=33694; NCOL=len(group_by_col)

def reader():
    # skip 2 comment lines + 1 dims line = 3 header rows
    return pd.read_csv(MTX, compression="gzip", sep=r"\s+", skiprows=3, header=None,
                       names=["gene","col","val"], dtype={"gene":np.int32,"col":np.int32,"val":np.float32},
                       chunksize=25_000_000, engine="c")

# PASS 1: cell library sizes (sum of counts per column)
cell_tot=np.zeros(NCOL, dtype=np.float64)
n=0
for ch in reader():
    np.add.at(cell_tot, ch["col"].values-1, ch["val"].values)
    n+=len(ch); print("pass1 rows",n,flush=True)
np.save("data/interim/norman_cell_tot.npy", cell_tot)
print("PASS1 done. cells with zero total:", int((cell_tot==0).sum()),flush=True)

# PASS 2: per-(gene,group) sum & sumsq of log1p(CP10K)
S=np.zeros(NGENE*ngroups, dtype=np.float64)
SS=np.zeros(NGENE*ngroups, dtype=np.float64)
safe_tot=np.where(cell_tot>0,cell_tot,1.0)
n=0
for ch in reader():
    g=ch["gene"].values-1; c=ch["col"].values-1; v=ch["val"].values
    grp=group_by_col[c]
    keep=grp>=0
    if not keep.any(): 
        n+=len(ch); continue
    g=g[keep]; c=c[keep]; v=v[keep]; grp=grp[keep]
    x=np.log1p(1e4*v/safe_tot[c])          # log1p CP10K
    idx=g*ngroups+grp
    np.add.at(S, idx, x)
    np.add.at(SS, idx, x*x)
    n+=len(ch); print("pass2 rows",n,flush=True)
S=S.reshape(NGENE,ngroups); SS=SS.reshape(NGENE,ngroups)
mean=S/np.maximum(ncells_group,1)          # zeros included via ncells_group denom
ex2=SS/np.maximum(ncells_group,1)
var=np.maximum(ex2-mean**2, 0.0)
np.savez_compressed("data/interim/norman_group_stats.npz",
                    mean=mean.astype(np.float32), var=var.astype(np.float32),
                    ncells=ncells_group)
print("PASS2 done -> norman_group_stats.npz", mean.shape, flush=True)
