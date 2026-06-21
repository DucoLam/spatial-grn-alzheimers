"""Four NON-network ways to show one regulator's weights across SQ2 clusters.
Renders a 2x2 contact sheet so we can pick a direction.
A = dumbbell | B = slopegraph | C = grouped bars | D = heatmap
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import analysis_common as ac

R = sys.argv[1] if len(sys.argv) > 1 else "AC106845.1"
ST = "L2/3 IT"
TOPK = 14
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"
OUT = ac.ANALYSIS / "SQ4" / f"sq4_regulon_candidates_{R.replace('.', '_').lower()}.png"

X, cell_ids, meta, _ = ac.load_cache()
idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
edge_src = idx["edge_src"].astype(str)
edge_dst = idx["edge_dst"].astype(str)
cols = np.flatnonzero(edge_src == R)
targets = edge_dst[cols]

rows = ac.subtype_rows(meta, ST)
cids = cell_ids[rows]
lab_df = pd.read_csv(LAB_CSV)
lab_map = dict(zip(lab_df.cell_id.astype(str), lab_df.grn_cluster.astype(int)))
lab = np.array([lab_map.get(c, -1) for c in cids])
valid = lab >= 0
M = X[rows][:, cols].toarray().astype(np.float64)
clusters = sorted(np.unique(lab[valid]))
means = {c: M[lab == c].mean(0) for c in clusters}
c0, c1 = clusters[0], clusters[1]

top = np.argsort(np.abs(means[c0]))[::-1][:TOPK]
genes = [targets[t] for t in top]
SCL = 1e5
w0 = np.abs(means[c0][top]) * SCL
w1 = np.abs(means[c1][top]) * SCL
tot0, tot1 = np.abs(means[c0]).sum(), np.abs(means[c1]).sum()
C0, C1 = "#3b6fb6", "#e8703a"
XLAB = r"mean edge weight $|w|$  ($\times10^{-5}$)"

fig, axs = plt.subplots(2, 2, figsize=(13, 11))
(axA, axB), (axC, axD) = axs

# --- A: dumbbell ---------------------------------------------------------
y = np.arange(len(top))[::-1]
for yi, a, b in zip(y, w1, w0):
    axA.plot([a, b], [yi, yi], color="0.82", lw=2.4, zorder=1, solid_capstyle="round")
axA.scatter(w1, y, s=70, color=C1, edgecolors="white", linewidths=0.7, zorder=3)
axA.scatter(w0, y, s=90, color=C0, edgecolors="white", linewidths=0.7, zorder=4)
axA.set_yticks(y); axA.set_yticklabels(genes, fontsize=9)
axA.set_xlim(left=0); axA.set_xlabel(XLAB, fontsize=10)
axA.tick_params(axis="y", length=0)
for s in ("top", "right", "left"):
    axA.spines[s].set_visible(False)
axA.grid(axis="x", ls=":", lw=0.7, color="0.86"); axA.set_axisbelow(True)
axA.set_title("A   Dumbbell", fontsize=13, loc="left", fontweight="bold")

# --- B: slopegraph -------------------------------------------------------
for g, a, b in zip(genes, w0, w1):
    axB.plot([0, 1], [a, b], color="0.7", lw=1.1, zorder=1)
axB.scatter(np.zeros(len(top)), w0, s=55, color=C0, zorder=3)
axB.scatter(np.ones(len(top)), w1, s=55, color=C1, zorder=3)
for g, a in zip(genes, w0):
    axB.text(-0.04, a, g, ha="right", va="center", fontsize=8)
axB.set_xlim(-0.5, 1.25); axB.set_xticks([0, 1])
axB.set_xticklabels([f"Cluster {c0}", f"Cluster {c1}"], fontsize=10)
axB.set_ylabel(XLAB, fontsize=10)
for s in ("top", "right", "bottom"):
    axB.spines[s].set_visible(False)
axB.tick_params(axis="x", length=0)
axB.set_title("B   Slopegraph", fontsize=13, loc="left", fontweight="bold")

# --- C: grouped horizontal bars -----------------------------------------
yb = np.arange(len(top))[::-1]
axC.barh(yb + 0.2, w0, height=0.38, color=C0, label=f"Cluster {c0}")
axC.barh(yb - 0.2, w1, height=0.38, color=C1, label=f"Cluster {c1}")
axC.set_yticks(yb); axC.set_yticklabels(genes, fontsize=9)
axC.set_xlabel(XLAB, fontsize=10); axC.tick_params(axis="y", length=0)
for s in ("top", "right", "left"):
    axC.spines[s].set_visible(False)
axC.grid(axis="x", ls=":", lw=0.7, color="0.86"); axC.set_axisbelow(True)
axC.legend(frameon=False, fontsize=10, loc="lower right")
axC.set_title("C   Grouped bars", fontsize=13, loc="left", fontweight="bold")

# --- D: heatmap ----------------------------------------------------------
mat = np.vstack([w0, w1]).T                       # genes x 2
im = axD.imshow(mat, aspect="auto", cmap="rocket_r" if False else "magma")
axD.set_yticks(range(len(top))); axD.set_yticklabels(genes, fontsize=9)
axD.set_xticks([0, 1]); axD.set_xticklabels([f"Cluster {c0}", f"Cluster {c1}"], fontsize=10)
for i in range(len(top)):
    for j, v in enumerate(mat[i]):
        axD.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7,
                 color="white" if v < mat.max() * 0.6 else "black")
cb = fig.colorbar(im, ax=axD, fraction=0.046, pad=0.04)
cb.set_label(XLAB, fontsize=9)
axD.set_title("D   Heatmap", fontsize=13, loc="left", fontweight="bold")

fig.suptitle(f"{R} target weights across SQ2 spatial clusters  "
             f"(Σ|w|: c{c0}={tot0:.1e}, c{c1}={tot1:.1e})", fontsize=14, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.98))
fig.savefig(OUT, dpi=150)
print("saved", OUT)
