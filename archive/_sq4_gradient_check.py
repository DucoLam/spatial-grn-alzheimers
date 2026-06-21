"""Global-activity-gradient check (L2/3 IT, SQ2 clusters 0/1).

Question: are the top spatial regulators distinctively up in cluster 0, or is cluster 0
just a globally louder GRN domain that ALL regulons ride?

Two measures:
 (A) Direct global: total |GRN weight| per cell (sum over ALL edges), mean per cluster.
 (B) Per-regulator: for every regulator, Sig|w| in c0 vs c1 (signed mean over cluster
     cells, then |.|.sum()), directional ratio c0/c1. Distribution over all regulators,
     and where the 20 robust hits land in it.
"""
import numpy as np
import pandas as pd
import analysis_common as ac

ST = "L2/3 IT"
OUT = ac.ANALYSIS / "SQ4" / "sq4_gradient_check.csv"
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"

X, cell_ids, meta, _ = ac.load_cache()
idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
edge_src = idx["edge_src"].astype(str)

rows = ac.subtype_rows(meta, ST)
Xr = X[rows].tocsr()
cids = cell_ids[rows]
lab_df = pd.read_csv(LAB_CSV)
lab_map = dict(zip(lab_df.cell_id.astype(str), lab_df.grn_cluster.astype(int)))
lab = np.array([lab_map.get(c, -1) for c in cids])
valid = lab >= 0
m0mask = valid & (lab == 0)
m1mask = valid & (lab == 1)
print(f"L2/3 IT: {valid.sum()} labelled cells | c0={m0mask.sum()} c1={m1mask.sum()}")

# (A) direct global activity gradient: total |w| per cell over ALL edges
tot = np.asarray(np.abs(Xr).sum(axis=1)).ravel()
g0, g1 = tot[m0mask].mean(), tot[m1mask].mean()
print("\n=== (A) GLOBAL total |GRN weight| per cell ===")
print(f"  cluster 0 mean total |w| = {g0:.4g}")
print(f"  cluster 1 mean total |w| = {g1:.4g}")
print(f"  ratio c0/c1 = {g0/g1:.2f}   (>>1 => cluster 0 is a globally louder GRN domain)")

# (B) per-regulator magnitude by cluster, over ALL regulators
master = pd.read_csv(ac.ANALYSIS / "SQ4" / "sq4_regulator_morans.csv")
mm = master[master.subtype == ST]
robust = set(mm[(mm.q_fdr < 0.05) & (mm.q_eqsec < 0.05)]
             .sort_values("I_R", ascending=False).head(20)["regulator"])

regs = pd.unique(edge_src)
recs = []
for R in regs:
    cols = np.flatnonzero(edge_src == R)
    sub = Xr[:, cols]
    a0 = np.asarray(sub[m0mask].mean(axis=0)).ravel()
    a1 = np.asarray(sub[m1mask].mean(axis=0)).ravel()
    mag0, mag1 = np.abs(a0).sum(), np.abs(a1).sum()
    recs.append({"regulator": R, "mag_c0": mag0, "mag_c1": mag1,
                 "ratio_c0_over_c1": mag0 / (mag1 + 1e-12),
                 "louder": "c0" if mag0 > mag1 else "c1",
                 "is_top20": R in robust})
df = pd.DataFrame(recs)
df.to_csv(OUT, index=False)

frac_c0 = (df.louder == "c0").mean()
r = df.ratio_c0_over_c1.replace([np.inf, -np.inf], np.nan).dropna()
print("\n=== (B) per-regulator c0/c1 magnitude ratio, ALL regulators ===")
print(f"  n regulators = {len(df)}")
print(f"  fraction louder in c0 = {frac_c0:.3f}")
print(f"  ratio c0/c1  median={r.median():.2f}  "
      f"25%={r.quantile(.25):.2f}  75%={r.quantile(.75):.2f}  "
      f"mean={r.mean():.2f}")
print(f"  log2(ratio)  median={np.log2(r).median():.2f}  "
      f"(0 = no gradient; >0 = c0 louder)")

top = df[df.is_top20]
print("\n=== top-20 robust hits vs the full distribution ===")
print(f"  top20 fraction louder in c0 = {(top.louder=='c0').mean():.3f}")
print(f"  top20 ratio c0/c1 median = {top.ratio_c0_over_c1.median():.2f}")
# percentile of each top20 gene's ratio within the full distribution
allr = df.ratio_c0_over_c1.values
for _, row in top.sort_values("ratio_c0_over_c1", ascending=False).iterrows():
    pct = (allr < row.ratio_c0_over_c1).mean() * 100
    print(f"    {row.regulator:<13} ratio={row.ratio_c0_over_c1:6.2f}  "
          f"(pctile {pct:5.1f} of all regulators)")
print("\nsaved", OUT)
