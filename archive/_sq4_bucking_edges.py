"""Edge-level bucking scan (L2/3 IT, SQ2 clusters 0/1).

Within a single TF's regulon, most target edges follow the regulon's dominant spatial
direction (e.g. louder in c0). A "bucking" edge is a specific TF->target relation that
goes the OTHER way -- stronger in the minority cluster. Those are genuine spatial
rewiring of an individual edge, which is exactly what SQ4 is after.

Weights = signed mean over ACTIVE cells per cluster (isolates wiring shape from the
on/off domain effect). To exclude the near-zero noise swaps we caught before, a bucking
edge must be SUBSTANTIAL: max(|w_dom|,|w_sub|) in the top quartile of the regulon's
targets, AND the minority-cluster weight |w_sub| >= ABS_FLOOR.
"""
import numpy as np
import pandas as pd
import analysis_common as ac

ST = "L2/3 IT"
ABS_FLOOR = 1e-4          # noise floor on the minority-cluster weight
TOPQ = 0.75               # "substantial" = top quartile by max(|w0|,|w1|)
OUTDIR = ac.ANALYSIS / "SQ4"
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"

X, cell_ids, meta, _ = ac.load_cache()
idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
edge_src = idx["edge_src"].astype(str)
edge_dst = idx["edge_dst"].astype(str)
rows = ac.subtype_rows(meta, ST)
Xr = X[rows].tocsr()
cids = cell_ids[rows]
lab_df = pd.read_csv(LAB_CSV)
lab_map = dict(zip(lab_df.cell_id.astype(str), lab_df.grn_cluster.astype(int)))
lab = np.array([lab_map.get(c, -1) for c in cids])
valid = lab >= 0

# top-5 most autocorrelated robust TFs + LINC02248 (the standout whole-regulon reverser)
master = pd.read_csv(OUTDIR / "sq4_regulator_morans.csv")
mm = master[master.subtype == ST]
top5 = (mm[(mm.q_fdr < 0.05) & (mm.q_eqsec < 0.05)]
        .sort_values("I_R", ascending=False).head(5)["regulator"].tolist())
GENES = top5 + (["LINC02248"] if "LINC02248" not in top5 else [])
print(f"scanning TFs: {GENES}")

def mean_active(M, sel):
    sub = M[sel]; act = np.abs(sub).sum(1) > 0
    return (sub[act].mean(0) if act.sum() > 0 else np.zeros(M.shape[1])), int(act.sum())

all_recs = []
for R in GENES:
    cols = np.flatnonzero(edge_src == R)
    targets = edge_dst[cols]
    M = Xr[:, cols].toarray().astype(np.float64)
    a0, n0 = mean_active(M, valid & (lab == 0))
    a1, n1 = mean_active(M, valid & (lab == 1))
    mag0, mag1 = np.abs(a0).sum(), np.abs(a1).sum()
    dom, sub = (0, 1) if mag0 >= mag1 else (1, 0)
    wdom = a0 if dom == 0 else a1
    wsub = a1 if dom == 0 else a0
    irow = mm[mm.regulator == R].iloc[0]

    mx = np.maximum(np.abs(a0), np.abs(a1))
    thr = np.quantile(mx[mx > 0], TOPQ) if (mx > 0).any() else np.inf
    substantial = mx >= thr
    buck = substantial & (np.abs(wsub) > np.abs(wdom)) & (np.abs(wsub) >= ABS_FLOOR)

    print(f"\n=== {R}  I_R={float(irow.I_R):.3f}  (d={len(cols)} targets, active n0={n0} n1={n1}) ===")
    print(f"  regulon louder in cluster {dom}  (mag c0={mag0:.4g}, c1={mag1:.4g})")
    print(f"  substantial targets (top quartile by |w|): {int(substantial.sum())}")
    print(f"  of those, BUCKING (stronger in minority c{sub}, |w|>={ABS_FLOOR}): {int(buck.sum())}")
    if buck.any():
        bi = np.flatnonzero(buck)
        bi = bi[np.argsort((np.abs(wsub) - np.abs(wdom))[bi])[::-1]]
        print(f"  top bucking edges (target: w_c0, w_c1, ratio |w_minority|/|w_dom|):")
        for t in bi[:12]:
            ratio = abs(wsub[t]) / (abs(wdom[t]) + 1e-12)
            print(f"    {targets[t]:<14} w0={a0[t]:+.5f}  w1={a1[t]:+.5f}  ratio={ratio:5.1f}x")
    for t in range(len(cols)):
        all_recs.append({"regulator": R, "target": targets[t],
                         "w_c0": a0[t], "w_c1": a1[t],
                         "dom_cluster": dom, "substantial": bool(substantial[t]),
                         "bucking": bool(buck[t])})

df = pd.DataFrame(all_recs)
df.to_csv(OUTDIR / "sq4_bucking_edges.csv", index=False)
print("\nsaved", OUTDIR / "sq4_bucking_edges.csv")
