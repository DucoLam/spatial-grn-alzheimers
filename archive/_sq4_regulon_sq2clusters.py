"""Compare AC106845.1's mean regulon between SQ2's GRN clusters (L2/3 IT, k=2).

SQ2's clusters are defined on the WHOLE GRN (independent of this regulator) and were
shown spatially organized (neighbour purity p=0.001). So "AC106845.1's regulon differs
between SQ2 clusters" is a clean, non-circular statement. Map = SQ2 result (cells coloured
by cluster); networks = AC106845.1's mean regulon in each cluster.
Mild caveat: AC106845.1's ~499 edges are among the 246,890 that defined SQ2's clusters,
so not 100% independent -- but negligible (1 of 500 regulators in a 20-PC clustering).
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import analysis_common as ac

R = sys.argv[1] if len(sys.argv) > 1 else "AC106845.1"
ST = "L2/3 IT"
TOPK = 14
PURITY_P = 0.001   # from SQ2 clustering_summary.csv (L2/3 IT)
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"
OUT = ac.ANALYSIS / "SQ4" / f"sq4_regulon_sq2clusters_{R.replace('.', '_').lower()}.png"

X, cell_ids, meta, _ = ac.load_cache()
idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
edge_src = idx["edge_src"].astype(str)
edge_dst = idx["edge_dst"].astype(str)
cols = np.flatnonzero(edge_src == R)
targets = edge_dst[cols]
print(f"{R}: {len(cols)} target edges")

rows = ac.subtype_rows(meta, ST)
cids = cell_ids[rows]
sect = meta["section"].values[rows]
xy = meta[["x", "y"]].values[rows]

lab_df = pd.read_csv(LAB_CSV)
lab_map = dict(zip(lab_df.cell_id.astype(str), lab_df.grn_cluster.astype(int)))
lab = np.array([lab_map.get(c, -1) for c in cids])
valid = lab >= 0
print(f"matched {valid.sum()}/{len(rows)} L2/3 cells to SQ2 labels; "
      f"cluster sizes {np.bincount(lab[valid])}")

M = X[rows][:, cols].toarray().astype(np.float64)        # signed regulon (n_l23, d)
clusters = sorted(np.unique(lab[valid]))
means = {c: M[(lab == c)].mean(0) for c in clusters}
for c in clusters:
    sub = M[lab == c]
    nact = int((np.abs(sub).sum(1) > 0).sum())
    print(f"  cluster {c}: n={sub.shape[0]} n_active={nact} "
          f"mean|regulon|={np.abs(means[c]).sum():.4g}")

# honesty: magnitude vs shape difference between the two clusters
if len(clusters) == 2:
    m0, m1 = means[clusters[0]], means[clusters[1]]
    mag0, mag1 = np.abs(m0).sum(), np.abs(m1).sum()
    s0, s1 = m0 / (mag0 or 1), m1 / (mag1 or 1)
    cos = float(s0 @ s1 / (np.linalg.norm(s0) * np.linalg.norm(s1) + 1e-12))
    print(f"magnitude ratio c{clusters[1]}/c{clusters[0]} = {mag1/mag0:.2f} | "
          f"shape cosine = {cos:.3f}")
    dT = np.argsort(np.abs(m1 - m0))[::-1][:10]
    print("top differential targets:",
          [(targets[t], round(m0[t], 6), round(m1[t], 6)) for t in dT])

# shared top-k targets + circular layout
score = np.maximum.reduce([np.abs(means[c]) for c in clusters])
top = np.argsort(score)[::-1][:TOPK]
ang = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(top), endpoint=False)
pos = np.c_[np.cos(ang), np.sin(ang)]

# Two side-by-side GRN panels (one per spatial cluster), identical layout + shared
# edge-width scale: the on/off magnitude contrast IS the story. Straight edges,
# thick lines and large labels so it stays readable when printed small.
from matplotlib.patches import FancyBboxPatch
clr = {0: "#3b6fb6", 1: "#e8703a", 2: "#3aa86b", 3: "#9467bd", 4: "#d6a800"}
gmax = max(np.abs(means[c][top]).max() for c in clusters) or 1.0
SCALE = 13.0   # max edge linewidth (doubled for small-size readability)

ncl = len(clusters)
fig, axes = plt.subplots(1, ncl, figsize=(4.9 * ncl, 5.4))
axes = np.atleast_1d(axes)

def draw_net(ax, mvec, col):
    for (px, py), wt in zip(pos, mvec[top]):
        ax.plot([0, px], [0, py], color=col, lw=0.8 + SCALE * abs(wt) / gmax,
                alpha=0.7, solid_capstyle="round", zorder=1)
    ax.scatter(pos[:, 0], pos[:, 1], s=150, color="white", edgecolors=col,
               linewidths=1.8, zorder=2)
    for (px, py), t in zip(pos, top):
        ha = "left" if px > 0.05 else ("right" if px < -0.05 else "center")
        ax.text(px * 1.18, py * 1.18, targets[t], fontsize=12, ha=ha, va="center",
                color="0.12")
    hub = FancyBboxPatch((-0.48, -0.15), 0.96, 0.30,
                         boxstyle="round,pad=0.02,rounding_size=0.15",
                         linewidth=0, facecolor=col, zorder=3)
    ax.add_patch(hub)
    ax.text(0, 0, R, color="white", fontsize=12, ha="center", va="center",
            fontweight="bold", zorder=4)
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal"); ax.axis("off")

for ax, c in zip(axes, clusters):
    draw_net(ax, means[c], clr[c % 5])
    ax.set_title(f"Spatial cluster {c}\n$\\Sigma|w|$ = {np.abs(means[c]).sum():.1e}",
                 fontsize=15, color=clr[c % 5], pad=8)

# extra horizontal gap so the inner labels of the two panels don't collide
fig.subplots_adjust(wspace=0.28, left=0.01, right=0.99, top=0.88, bottom=0.02)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("saved", OUT)
