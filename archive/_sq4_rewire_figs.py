"""Phase 2: draw the rewiring winners chosen by Phase 1.

Reads sq4_rewire_scan.csv, picks the top-3 rewiring candidates (highest rewiring_score
among genes active in both SQ2 clusters), draws each gene's regulon network per SQ2
cluster, and for the single best candidate draws it across the 3 largest L2/3 sections
(within-section SQ2 split). If nothing genuinely rewires, it says so and still draws the
top-scoring gene for the record.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import analysis_common as ac

ST = "L2/3 IT"
TOPK = 16
OUTDIR = ac.ANALYSIS / "SQ4"
SCAN = OUTDIR / "sq4_rewire_scan.csv"
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"

scan = pd.read_csv(SCAN)
cand = scan[scan.both_active == True].sort_values("rewiring_score", ascending=False)
if len(cand) == 0:
    print("NO genes active in both clusters -> no rewiring assessable; using top I_R gene")
    winners = scan.sort_values("I_R", ascending=False).head(1)["regulator"].tolist()
else:
    winners = cand.head(3)["regulator"].tolist()
best = winners[0]
print(f"rewiring winners: {winners} | best={best}")
print(scan[["regulator", "shape_cosine", "rewiring_score", "both_active"]].to_string(index=False))

X, cell_ids, meta, _ = ac.load_cache()
idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
edge_src = idx["edge_src"].astype(str)
edge_dst = idx["edge_dst"].astype(str)

rows = ac.subtype_rows(meta, ST)
cids = cell_ids[rows]
sect = meta["section"].values[rows]
lab_df = pd.read_csv(LAB_CSV)
lab_map = dict(zip(lab_df.cell_id.astype(str), lab_df.grn_cluster.astype(int)))
lab = np.array([lab_map.get(c, -1) for c in cids])
valid = lab >= 0
PAL = ["#7fb8d6", "#08519c"]


def mean_active(M, sel):
    sub = M[sel]
    act = np.abs(sub).sum(1) > 0
    return (sub[act].mean(0) if act.sum() > 0 else np.zeros(M.shape[1])), int(act.sum())


def draw_net(ax, mvec, top, targets, title, node_color):
    w = mvec[top]; mx = np.abs(w).max() or 1.0
    ang = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(top), endpoint=False)
    pos = np.c_[np.cos(ang), np.sin(ang)]
    for (px, py), wt in zip(pos, w):
        col = "#c0392b" if wt < 0 else "#2c7fb8"
        ax.plot([0, px], [0, py], color=col, lw=0.4 + 5 * abs(wt) / mx,
                alpha=0.85, zorder=1, solid_capstyle="round")
    ax.scatter(pos[:, 0], pos[:, 1], s=55, color="0.88", edgecolors="0.4",
               linewidths=0.5, zorder=2)
    for (px, py), t in zip(pos, top):
        ax.text(px * 1.2, py * 1.2, targets[t], fontsize=7, ha="center", va="center")
    ax.scatter([0], [0], s=360, color=node_color, zorder=3)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=10)


# ---- per-gene SQ2-cluster figures for the winners ----
for R in winners:
    cols = np.flatnonzero(edge_src == R)
    targets = edge_dst[cols]
    M = X[rows][:, cols].toarray().astype(np.float64)
    m0, n0 = mean_active(M, valid & (lab == 0))
    m1, n1 = mean_active(M, valid & (lab == 1))
    s0 = m0 / (np.abs(m0).sum() + 1e-12); s1 = m1 / (np.abs(m1).sum() + 1e-12)
    cos = float(s0 @ s1 / (np.linalg.norm(s0) * np.linalg.norm(s1) + 1e-12))
    top = np.argsort(np.maximum(np.abs(m0), np.abs(m1)))[::-1][:TOPK]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4))
    draw_net(axes[0], m0, top, targets, f"SQ2 cluster 0 (active n={n0})", PAL[0])
    draw_net(axes[1], m1, top, targets, f"SQ2 cluster 1 (active n={n1})", PAL[1])
    fig.suptitle(f"{R} — regulon across SQ2 spatial clusters (L2/3 IT)\n"
                 f"shape cosine={cos:.3f}  (low = rewiring; high = on/off)  ·  "
                 f"blue=activation red=repression, width∝|w|", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUTDIR / f"sq4_rewire_fig_{R.replace('.', '_').lower()}.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  drew {R} (cosine {cos:.3f})")

# ---- best gene across the 3 largest L2/3 sections ----
big = pd.Series(sect[valid]).value_counts().head(3).index.tolist()
cols = np.flatnonzero(edge_src == best)
targets = edge_dst[cols]
M = X[rows][:, cols].toarray().astype(np.float64)
top = np.argsort(np.abs(M[valid]).mean(0))[::-1][:TOPK]
fig, axes = plt.subplots(len(big), 2, figsize=(11, 5.0 * len(big)), squeeze=False)
for r, s in enumerate(big):
    for c in (0, 1):
        mv, na = mean_active(M, valid & (sect == s) & (lab == c))
        draw_net(axes[r][c], mv, top, targets,
                 f"section {s[:18]}…\ncluster {c} (active n={na})", PAL[c])
fig.suptitle(f"{best} — regulon by SQ2 cluster across the 3 largest L2/3 sections\n"
             f"(rows = sections; note: different sections are different donors)",
             fontsize=11)
fig.tight_layout()
fig.savefig(OUTDIR / f"sq4_rewire_persection_{best.replace('.', '_').lower()}.png",
            dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  drew per-section panel for {best}")
print("done")
