"""Diagnose the shape-cosine gap: weight concentration, cosine vs top-K, and a noise
floor (within-cluster random-split cosine). Also redraw the regulon with more edges.

cosine is scale-invariant, so shape cosine = cosine(mean0, mean1) over the chosen targets.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import analysis_common as ac

GENES = ["AC109492.1", "AC106845.1"]
ST = "L2/3 IT"
DRAW_K = 40
OUTDIR = ac.ANALYSIS / "SQ4"
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"

X, cell_ids, meta, _ = ac.load_cache()
idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
edge_src = idx["edge_src"].astype(str)
edge_dst = idx["edge_dst"].astype(str)
rows = ac.subtype_rows(meta, ST)
cids = cell_ids[rows]
lab_df = pd.read_csv(LAB_CSV)
lab_map = dict(zip(lab_df.cell_id.astype(str), lab_df.grn_cluster.astype(int)))
lab = np.array([lab_map.get(c, -1) for c in cids])
valid = lab >= 0
rng = np.random.default_rng(ac.SEED)

def cos(a, b):
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else np.nan

def mean_active(M, sel):
    sub = M[sel]; act = np.abs(sub).sum(1) > 0
    return sub[act].mean(0) if act.sum() > 0 else np.zeros(M.shape[1]), int(act.sum())

for R in GENES:
    cols = np.flatnonzero(edge_src == R)
    targets = edge_dst[cols]
    M = X[rows][:, cols].toarray().astype(np.float64)
    d = len(cols)
    m0, n0 = mean_active(M, valid & (lab == 0))
    m1, n1 = mean_active(M, valid & (lab == 1))

    # concentration of the overall regulon
    imp = np.abs(M[valid & (np.abs(M).sum(1) > 0)].mean(0))
    order = np.argsort(imp)[::-1]
    cum = np.cumsum(imp[order]) / imp.sum()
    def kfor(frac):
        return int(np.searchsorted(cum, frac) + 1)
    print(f"\n=== {R}  (d={d} targets, active n0={n0} n1={n1}) ===")
    print(f"  weight in top16={cum[min(15,d-1)]:.2f} top32={cum[min(31,d-1)]:.2f} "
          f"top64={cum[min(63,d-1)]:.2f} | targets for 50%={kfor(.5)} 80%={kfor(.8)} 90%={kfor(.9)}")

    # cosine restricted to the top-K most important targets
    print("  cosine(cluster0, cluster1) by top-K targets:")
    for K in [8, 16, 32, 64, 128, d]:
        ix = order[:K]
        print(f"    top{K:>4}: cos={cos(m0[ix], m1[ix]):.3f}")

    # noise floor: split each cluster's active cells in half, cosine between halves
    floors = {}
    for c, lc in ((0, n0), (1, n1)):
        sub = M[valid & (lab == c)]
        sub = sub[np.abs(sub).sum(1) > 0]
        cs = []
        for _ in range(50):
            perm = rng.permutation(sub.shape[0]); h = sub.shape[0] // 2
            cs.append(cos(sub[perm[:h]].mean(0), sub[perm[h:]].mean(0)))
        floors[c] = float(np.mean(cs))
    print(f"  noise-floor cosine (within-cluster random split): "
          f"c0={floors[0]:.3f}  c1={floors[1]:.3f}   <-- compare to between-cluster cos")

    # redraw with more edges
    top = order[:min(DRAW_K, d)]
    ang = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(top), endpoint=False)
    pos = np.c_[np.cos(ang), np.sin(ang)]
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.6))
    for ax, mv, ti, nc in ((axes[0], m0, f"cluster 0 (active n={n0})", "#7fb8d6"),
                           (axes[1], m1, f"cluster 1 (active n={n1})", "#08519c")):
        w = mv[top]; mx = np.abs(w).max() or 1.0
        for (px, py), wt in zip(pos, w):
            col = "#c0392b" if wt < 0 else "#2c7fb8"
            ax.plot([0, px], [0, py], color=col, lw=0.3 + 4 * abs(wt) / mx, alpha=0.8)
        ax.scatter(pos[:, 0], pos[:, 1], s=34, color="0.88", edgecolors="0.5", lw=0.4)
        for (px, py), t in zip(pos, top):
            ax.text(px * 1.16, py * 1.16, targets[t], fontsize=6, ha="center", va="center")
        ax.scatter([0], [0], s=300, color=nc)
        ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.4); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(ti, fontsize=10)
    fig.suptitle(f"{R} — top {len(top)} edges per SQ2 cluster (L2/3 IT)", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUTDIR / f"sq4_regulon_top{DRAW_K}_{R.replace('.', '_').lower()}.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)
print("\ndone")
