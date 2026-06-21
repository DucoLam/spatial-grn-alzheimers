"""Route B figure: cluster cells by AC106845.1's regulon (L2/3 IT, largest section),
show the clusters are spatially separated (the finding), and draw each cluster's mean
regulon as a star network (what differs). Map = evidence; networks = illustration."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import analysis_common as ac

R = "AC106845.1"
ST = "L2/3 IT"
TOPK = 18
OUT = ac.ANALYSIS / "SQ4" / "sq4_regulon_clusters_ac106845.png"

X, cell_ids, meta, _ = ac.load_cache()
idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
edge_src = idx["edge_src"].astype(str)
edge_dst = idx["edge_dst"].astype(str)
cols = np.flatnonzero(edge_src == R)
targets = edge_dst[cols]
print(f"{R}: {len(cols)} target edges")

rows = ac.subtype_rows(meta, ST)
sub_sect = meta["section"].values[rows]
Mabs_all = np.abs(X[rows][:, cols].toarray())
active_all = Mabs_all.sum(1) > 0

# largest section by # active cells (best-estimated network)
best_sec, best_n = None, 0
for s in pd.unique(sub_sect):
    n = int(((sub_sect == s) & active_all).sum())
    if n > best_n:
        best_n, best_sec = n, s
print(f"largest L2/3 section = {best_sec} with {best_n} active cells")

mask = (sub_sect == best_sec) & active_all
sel = rows[np.flatnonzero(mask)]                 # global indices, one section
M = X[sel][:, cols].toarray().astype(np.float64) # signed regulon (n, d)
xy = meta[["x", "y"]].values[sel]

# cluster cells by their regulon
km = KMeans(n_clusters=2, n_init=10, random_state=ac.SEED).fit(M)
lab = km.labels_
# order clusters so 0 = lower mean activity (for consistent colour)
act = np.array([np.abs(M[lab == c]).sum(1).mean() for c in (0, 1)])
if act[0] > act[1]:
    lab = 1 - lab

# spatial separation: 6-NN neighbour purity vs label-permutation null
W = ac.within_section_knn_W(meta, sel).tocoo()
aa, bb = W.row, W.col
obs = float(np.mean(lab[aa] == lab[bb]))
rng = np.random.default_rng(ac.SEED)
null = np.empty(999)
for i in range(999):
    lp = lab[rng.permutation(len(lab))]
    null[i] = np.mean(lp[aa] == lp[bb])
p_sep = (np.sum(null >= obs) + 1) / 1000.0

mean0 = M[lab == 0].mean(0)
mean1 = M[lab == 1].mean(0)
# honesty diagnostics: magnitude vs shape difference
mag0, mag1 = np.abs(mean0).sum(), np.abs(mean1).sum()
s0, s1 = mean0 / (mag0 or 1), mean1 / (mag1 or 1)
shape_cos = float(s0 @ s1 / (np.linalg.norm(s0) * np.linalg.norm(s1) + 1e-12))
print(f"cluster sizes: {np.bincount(lab)} | purity obs={obs:.3f} p={p_sep:.4f}")
print(f"mean |regulon| (magnitude): c0={mag0:.4g} c1={mag1:.4g} ratio={mag1/mag0:.2f}")
print(f"shape cosine between cluster mean regulons = {shape_cos:.3f}")
dT = np.argsort(np.abs(mean1 - mean0))[::-1][:8]
print("top differential targets:", [(targets[t], round(mean0[t], 5), round(mean1[t], 5))
                                     for t in dT])

# shared top-k targets + circular layout
score = np.maximum(np.abs(mean0), np.abs(mean1))
top = np.argsort(score)[::-1][:TOPK]
ang = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(top), endpoint=False)
pos = np.c_[np.cos(ang), np.sin(ang)]

fig = plt.figure(figsize=(16.5, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1])

# left: spatial map of the section, coloured by regulon cluster
axm = fig.add_subplot(gs[0, 0])
cmap = {0: "#7fb8d6", 1: "#08519c"}
for c in (0, 1):
    m = lab == c
    axm.scatter(xy[m, 0], xy[m, 1], s=14, color=cmap[c], alpha=0.8,
                label=f"cluster {c} (n={int(m.sum())})")
axm.set_aspect("equal"); axm.set_xticks([]); axm.set_yticks([])
axm.set_title(f"L2/3 IT section {best_sec}\ncells coloured by {R} regulon cluster\n"
              f"neighbour purity p={p_sep:.3f}", fontsize=11)
axm.legend(fontsize=9, loc="best")

def draw_net(ax, means, title):
    w = means[top]
    mx = np.abs(w).max() or 1.0
    for (px, py), wt in zip(pos, w):
        col = "#c0392b" if wt < 0 else "#2c7fb8"
        ax.plot([0, px], [0, py], color=col, lw=0.4 + 5 * abs(wt) / mx,
                alpha=0.85, zorder=1, solid_capstyle="round")
    ax.scatter(pos[:, 0], pos[:, 1], s=70, color="0.88", edgecolors="0.4",
               linewidths=0.6, zorder=2)
    for (px, py), t in zip(pos, top):
        ax.text(px * 1.18, py * 1.18, targets[t], fontsize=8, ha="center",
                va="center", rotation=0)
    ax.scatter([0], [0], s=420, color="#333", zorder=3)
    ax.text(0, 0, R, color="white", fontsize=8.5, ha="center", va="center",
            fontweight="bold", zorder=4)
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=11)

axA = fig.add_subplot(gs[0, 1])
axB = fig.add_subplot(gs[0, 2])
draw_net(axA, mean0, f"cluster 0 — mean regulon\n(lower activity, |w|={mag0:.3g})")
draw_net(axB, mean1, f"cluster 1 — mean regulon\n(higher activity, |w|={mag1:.3g})")
fig.suptitle(f"{R} regulon varies across space in L2/3 IT  "
             f"(blue = activation, red = repression; line width ∝ |weight|)",
             fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print("saved", OUT)
