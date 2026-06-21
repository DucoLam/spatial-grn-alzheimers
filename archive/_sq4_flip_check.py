"""Is it a flip or scaled reweighting? Per-target weight scatter (cluster0 vs cluster1),
Pearson/Spearman, biggest rank-swappers, and a SHARED-SCALE network redraw (honest
magnitude). cluster means = mean over active cells, signed."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import analysis_common as ac

GENES = ["AC109492.1", "AC106845.1"]
ST = "L2/3 IT"
DRAW_K = 40
OUTDIR = ac.ANALYSIS / "SQ4"
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"

X, cell_ids, meta, _ = ac.load_cache()
idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
edge_src = idx["edge_src"].astype(str); edge_dst = idx["edge_dst"].astype(str)
rows = ac.subtype_rows(meta, ST)
cids = cell_ids[rows]
lab_df = pd.read_csv(LAB_CSV)
lab_map = dict(zip(lab_df.cell_id.astype(str), lab_df.grn_cluster.astype(int)))
lab = np.array([lab_map.get(c, -1) for c in cids]); valid = lab >= 0

def mean_active(M, sel):
    sub = M[sel]; act = np.abs(sub).sum(1) > 0
    return sub[act].mean(0) if act.sum() > 0 else np.zeros(M.shape[1])

for R in GENES:
    cols = np.flatnonzero(edge_src == R); targets = edge_dst[cols]
    M = X[rows][:, cols].toarray().astype(np.float64)
    m0 = mean_active(M, valid & (lab == 0))
    m1 = mean_active(M, valid & (lab == 1))

    pear = pearsonr(m0, m1)[0]
    spear_signed = spearmanr(m0, m1).correlation
    spear_mag = spearmanr(np.abs(m0), np.abs(m1)).correlation   # does loud stay loud?
    print(f"\n=== {R} ===")
    print(f"  Pearson(signed w)        = {pear:.3f}")
    print(f"  Spearman(signed w)       = {spear_signed:.3f}")
    print(f"  Spearman(|w|) loud-stays-loud = {spear_mag:.3f}   "
          f"(near 1 = no flip; <=0 = flip)")
    # biggest magnitude-rank swappers
    r0 = pd.Series(np.abs(m0)).rank(ascending=False)
    r1 = pd.Series(np.abs(m1)).rank(ascending=False)
    dchg = (r0 - r1)
    big = np.argsort(np.abs(dchg.values))[::-1][:10]
    print("  biggest |w|-rank swaps (target: rank_c0 -> rank_c1):")
    for t in big:
        print(f"    {targets[t]:<14} {int(r0[t]):>4} -> {int(r1[t]):<4} "
              f"(w0={m0[t]:.5f}, w1={m1[t]:.5f})")

    # figure: scatter + shared-scale nets
    fig = plt.figure(figsize=(16, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1])
    ax = fig.add_subplot(gs[0, 0])
    lim = max(np.abs(m0).max(), np.abs(m1).max()) * 1.05
    ax.scatter(m0, m1, s=12, alpha=0.5, color="#08519c")
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8, label="y = x (same shape)")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("target weight in cluster 0"); ax.set_ylabel("target weight in cluster 1")
    ax.set_title(f"per-target weights\nPearson={pear:.2f}  Spearman|w|={spear_mag:.2f}",
                 fontsize=10)
    ax.legend(fontsize=8)

    order = np.argsort(np.maximum(np.abs(m0), np.abs(m1)))[::-1][:DRAW_K]
    ang = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(order), endpoint=False)
    pos = np.c_[np.cos(ang), np.sin(ang)]
    gmax = max(np.abs(m0[order]).max(), np.abs(m1[order]).max()) or 1.0   # SHARED scale
    for j, (mv, ti, nc) in enumerate(((m0, "cluster 0", "#7fb8d6"),
                                      (m1, "cluster 1", "#08519c"))):
        axn = fig.add_subplot(gs[0, j + 1])
        for (px, py), wt in zip(pos, mv[order]):
            col = "#c0392b" if wt < 0 else "#2c7fb8"
            axn.plot([0, px], [0, py], color=col, lw=0.2 + 5 * abs(wt) / gmax, alpha=0.8)
        axn.scatter(pos[:, 0], pos[:, 1], s=26, color="0.88", edgecolors="0.5", lw=0.4)
        for (px, py), t in zip(pos, order):
            axn.text(px * 1.16, py * 1.16, targets[t], fontsize=6, ha="center", va="center")
        axn.scatter([0], [0], s=260, color=nc)
        axn.set_xlim(-1.4, 1.4); axn.set_ylim(-1.4, 1.4); axn.set_aspect("equal"); axn.axis("off")
        axn.set_title(f"{ti} (SHARED width scale)", fontsize=10)
    fig.suptitle(f"{R} — flip check: per-target scatter + shared-scale regulon (L2/3 IT)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(OUTDIR / f"sq4_flipcheck_{R.replace('.', '_').lower()}.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)
print("\ndone")
