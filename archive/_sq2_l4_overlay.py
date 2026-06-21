"""Test the hypothesis: in one section, SQ2 cluster 1 (of L2/3 IT) forms a band along
the L4 IT boundary (i.e. the deep side of L2/3), while cluster 0 is the rest.
Overlay L4 IT cells (grey) under the L2/3 clusters, and quantify each L2/3 cell's
distance to the nearest L4 IT cell, c0 vs c1."""
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
import analysis_common as ac

PURITY_P = 0.001
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"
palette = ["#3b6fb6", "#e8703a"]   # c0 blue, c1 orange
SEC = sys.argv[1] if len(sys.argv) > 1 else "H21.33.019.Cx30.MTG.02.007.5.1"

meta = pd.read_parquet(ac.CACHE / "spatial_meta.parquet")
meta["section"] = meta["section"].astype(str)
sec = next((s for s in meta["section"].unique() if s.startswith(SEC) or SEC in s), SEC)

l23 = meta[(meta.subtype == "L2/3 IT") & (meta.section == sec)].copy()
l4 = meta[(meta.subtype == "L4 IT") & (meta.section == sec)].copy()
lab = pd.read_csv(LAB_CSV)
lmap = dict(zip(lab.cell_id.astype(str), lab.grn_cluster.astype(int)))
l23["cluster"] = l23.index.astype(str).map(lmap)
l23 = l23[l23.cluster.notna()].copy(); l23["cluster"] = l23.cluster.astype(int)
print(f"section {sec}: L2/3={len(l23)} (c0={int((l23.cluster==0).sum())}, "
      f"c1={int((l23.cluster==1).sum())}), L4 IT={len(l4)}")

# distance from each L2/3 cell to nearest L4 IT cell
if len(l4) > 0:
    tree = cKDTree(l4[["x", "y"]].values)
    d, _ = tree.query(l23[["x", "y"]].values, k=1)
    l23["d_to_L4"] = d
    d0 = l23.loc[l23.cluster == 0, "d_to_L4"]; d1 = l23.loc[l23.cluster == 1, "d_to_L4"]
    U, p = mannwhitneyu(d1, d0, alternative="less")  # H1: c1 closer to L4 than c0
    print(f"median dist to nearest L4 IT:  c0={d0.median():.0f}  c1={d1.median():.0f} (µm)")
    print(f"Mann-Whitney (c1 < c0): U={U:.0f}, p={p:.2e}")
else:
    print("no L4 IT cells in this section")

fig, ax = plt.subplots(figsize=(7.5, 7.5))
if len(l4):
    ax.scatter(l4.x, l4.y, s=14, color="0.7", alpha=0.35, edgecolors="none",
               label=f"L4 IT  (n={len(l4)})", zorder=1)
for c in (0, 1):
    s = l23[l23.cluster == c]
    ax.scatter(s.x, s.y, s=22, color=palette[c], alpha=0.55, edgecolors="none",
               label=f"L2/3 SQ2 cluster {c}  (n={len(s)})", zorder=2)
ax.set_aspect("equal"); ax.set_axisbelow(True); ax.minorticks_on()
ax.grid(True, which="major", ls=":", lw=0.7, color="0.7", alpha=0.85)
ax.grid(True, which="minor", ls=":", lw=0.35, color="0.85", alpha=0.5)
ax.tick_params(which="major", direction="out", length=5, labelsize=9)
ax.set_xlabel(r"x position ($\mu$m)", fontsize=11)
ax.set_ylabel(r"y position ($\mu$m)", fontsize=11)
ax.set_title(f"L2/3 IT GRN clusters vs L4 IT — section {sec}\n"
             f"(does cluster 1 line the L4 boundary?)", fontsize=11)
ax.legend(fontsize=10, loc="best", framealpha=0.9)
fig.tight_layout()
out = ac.ANALYSIS / "SQ2" / f"sq2_l4_overlay_{sec.replace('.', '_')}.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved", out)
