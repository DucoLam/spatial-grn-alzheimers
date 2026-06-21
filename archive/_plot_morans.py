"""Paper figures for SQ2 Moran's I, built from the morans_*.csv result files.
  (1) morans_summary_box.png  -- GRN vs expression Moran's I across the 20 PCs,
                                 per subtype (the headline + circularity baseline)
  (2) morans_perpc.png        -- per-PC paired bars (GRN vs expression), 2x2 subtypes
"""
import os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.chdir("/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ2/morans_i")
ORDER = [("Astrocyte", "morans_astrocyte.csv"), ("L2/3 IT", "morans_l23_it.csv"),
         ("L4 IT", "morans_l4_it.csv"), ("Oligodendrocyte", "morans_oligo.csv")]
data = {name: pd.read_csv(f) for name, f in ORDER if os.path.exists(f)}

# (1) summary boxplot: GRN vs expression Moran's I distribution over PCs ------
fig, ax = plt.subplots(figsize=(8, 4.8))
pos, ticks, ticklab = 0, [], []
for name, df in data.items():
    g = df[df.space == "GRN"]["morans_I"].values
    e = df[df.space == "expression"]["morans_I"].values
    bp = ax.boxplot([g, e], positions=[pos, pos + 0.7], widths=0.55,
                    patch_artist=True, showfliers=False)
    for patch, col in zip(bp["boxes"], ["#1f77b4", "#bbbbbb"]):
        patch.set_facecolor(col)
    for med in bp["medians"]:
        med.set_color("black")
    ticks.append(pos + 0.35); ticklab.append(name); pos += 2
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(ticks); ax.set_xticklabels(ticklab)
ax.set_ylabel("Moran's I  (spatial autocorrelation)")
ax.set_title("SQ2: spatial structure of GRN components vs expression baseline\n"
             "(each box = 20 principal components, within-section)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor="#1f77b4", label="GRN"),
                   Patch(facecolor="#bbbbbb", label="expression baseline")],
          frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig("morans_summary_box.png", dpi=160)
plt.close(fig)

# (2) per-PC paired bars, 2x2 -------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
for ax, (name, df) in zip(axes.ravel(), data.items()):
    g = df[df.space == "GRN"].sort_values("pc")
    e = df[df.space == "expression"].sort_values("pc")
    x = g["pc"].values
    ax.bar(x - 0.2, g["morans_I"], width=0.4, color="#1f77b4", label="GRN")
    ax.bar(x + 0.2, e["morans_I"].reindex(range(len(x))).fillna(0).values
           if len(e) != len(x) else e["morans_I"].values,
           width=0.4, color="#bbbbbb", label="expression")
    # star significant GRN PCs
    for xi, (val, p) in zip(x, zip(g["morans_I"], g["p_perm"])):
        if p < 0.05:
            ax.text(xi - 0.2, val + 0.01, "*", ha="center", fontsize=9, color="#1f77b4")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title(name); ax.set_xticks(x[::2])
    ax.set_ylabel("Moran's I")
axes.ravel()[0].legend(frameon=False)
for ax in axes[1]:
    ax.set_xlabel("principal component")
fig.suptitle("SQ2: Moran's I per GRN principal component (* = p<0.05, permutation)",
             y=1.0)
fig.tight_layout()
fig.savefig("morans_perpc.png", dpi=160)
plt.close(fig)
print("wrote morans_summary_box.png and morans_perpc.png")
