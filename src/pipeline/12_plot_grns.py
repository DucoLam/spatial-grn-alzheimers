"""
12_plot_grns.py -- visualize a few nearby L2/3 IT GRNs (SQ2 sanity preview).

Picks the section with the most located L2/3 IT cells, takes the 5 cells closest
to that section's centroid, loads their wScReNI networks, and draws them side by
side over a COMMON set of top edges (fixed layout) so they're directly comparable.
Network diagram if networkx is available, else heatmaps. No scipy dependency.
"""
import os
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TG   = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/data/tangram/Tangram_cell_locations_l23_it.csv")
GRN  = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/data/wscreni_networks_pool4/wScReNI")
OUT  = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/plots/grn_examples")
OUT.mkdir(parents=True, exist_ok=True)
N, TOPK = 5, 40

df = pd.read_csv(TG)
df["sec"] = df["tangram_mapped_section"].astype(str)
s = df["sec"].value_counts().idxmax()
sub = df[df["sec"] == s].reset_index(drop=True)
xy = sub[["tangram_x", "tangram_y"]].values
cen = xy.mean(0)
seed = int(np.argmin(((xy - cen) ** 2).sum(1)))
order = np.argsort(np.sqrt(((xy - xy[seed]) ** 2).sum(1)))[:N]
chosen = sub.iloc[order].reset_index(drop=True)
print(f"section {s}  (L2/3 IT here: {len(sub)});  picked {len(chosen)} nearby cells")

# barcode -> path index (one scan of the big dir)
bc2path = {}
for fn in os.listdir(GRN):
    if fn.endswith(".network.txt"):
        core = fn[:-len(".network.txt")]
        bc = core.split(".", 1)[1] if "." in core else core
        bc2path[bc] = GRN / fn

mats, names, coords = [], [], []
for _, r in chosen.iterrows():
    f = bc2path.get(r["cell_id"])
    if f is None:
        print("  no GRN for", r["cell_id"]); continue
    mats.append(pd.read_csv(f, sep="\t", index_col=0))
    names.append(r["cell_id"][:16]); coords.append((r["tangram_x"], r["tangram_y"]))
print("loaded", len(mats), "GRNs")

genes = mats[0].index.values
A = np.stack([m.values for m in mats])          # (k, 500, 500)
meanabs = np.abs(A).mean(0)
ii, jj = np.unravel_index(np.argsort(meanabs.flatten())[-TOPK:], meanabs.shape)
involved = sorted(set(ii.tolist()) | set(jj.tolist()))
labels = [genes[k] for k in involved]
gi = {g: i for i, g in enumerate(genes)}

try:
    import networkx as nx
    Gu = nx.DiGraph()
    for a, b in zip(ii, jj):
        Gu.add_edge(genes[a], genes[b])
    pos = nx.spring_layout(Gu, seed=1, k=0.6)
    fig, axes = plt.subplots(1, len(mats), figsize=(5 * len(mats), 5.5))
    axes = np.atleast_1d(axes)
    for ax, m, nm, xyc in zip(axes, mats, names, coords):
        M = m.values
        ews = [abs(M[gi[a], gi[b]]) for a, b in Gu.edges]
        mx = max(ews) if ews else 1
        nx.draw_networkx_nodes(Gu, pos, ax=ax, node_size=70, node_color="#e0e0e0")
        nx.draw_networkx_labels(Gu, pos, ax=ax, font_size=5)
        nx.draw_networkx_edges(Gu, pos, ax=ax, width=[0.3 + 3 * w / mx for w in ews],
                               edge_color="steelblue", alpha=0.6, arrows=False)
        ax.set_title(f"{nm}\n({xyc[0]:.0f}, {xyc[1]:.0f})", fontsize=7); ax.axis("off")
    mode = "network"
except Exception as e:
    print("networkx unavailable -> heatmaps:", e)
    fig, axes = plt.subplots(1, len(mats), figsize=(4.5 * len(mats), 5))
    axes = np.atleast_1d(axes)
    vmax = np.abs(A[:, involved][:, :, involved]).max()
    for ax, m, nm, xyc in zip(axes, mats, names, coords):
        subM = m.values[np.ix_(involved, involved)]
        ax.imshow(subM, cmap="viridis", vmin=0, vmax=vmax, aspect="auto")
        ax.set_title(f"{nm}\n({xyc[0]:.0f}, {xyc[1]:.0f})", fontsize=7)
        ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=4, rotation=90)
        ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=4)
    mode = "heatmap"

fig.suptitle(f"{len(mats)} nearby L2/3 IT GRNs — section {'.'.join(s.split('.')[-3:])} "
             f"(top {TOPK} edges, {mode})", fontsize=10)
out = OUT / "grn_examples_l23it.png"
fig.savefig(out, dpi=140, bbox_inches="tight")
print("wrote", out, "mode", mode, "cells", len(mats))
