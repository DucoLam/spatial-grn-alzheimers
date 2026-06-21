"""GRN-only version of morans_perpc.png (no expression baseline).
Per-PC Moran's I for the GRN, 4x1 stacked subtypes, significant PCs highlighted.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.chdir("/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ2/morans_i")
ORDER = [("Astrocyte", "morans_astrocyte.csv"), ("L2/3 IT", "morans_l23_it.csv"),
         ("L4 IT", "morans_l4_it.csv"), ("Oligodendrocyte", "morans_oligo.csv")]
data = {name: pd.read_csv(f) for name, f in ORDER if os.path.exists(f)}

fig, axes = plt.subplots(4, 1, figsize=(7, 11), sharex=True)
for ax, (name, df) in zip(axes, data.items()):
    g = df[df.space == "GRN"].sort_values("pc")
    x = g["pc"].values
    sig = g["p_perm"].values < 0.05
    colors = ["#1f77b4" if s else "#c9c9c9" for s in sig]
    ax.bar(x, g["morans_I"].values, width=0.7, color=colors)
    for xi, val, s in zip(x, g["morans_I"].values, sig):
        if s:
            ax.text(xi, val + 0.01, "*", ha="center", fontsize=9, color="#1f77b4")
    ax.axhline(0, color="k", lw=0.8)
    nsig = int(sig.sum())
    ax.set_title(f"{name}  ({nsig}/{len(x)} PCs significant)", fontsize=11)
    ax.set_ylabel("Moran's I")
axes[-1].set_xticks(x)
axes[-1].set_xlabel("GRN principal component")
fig.suptitle("SQ2: Moran's I per GRN principal component "
             "(blue = p<0.05 permutation; grey = n.s.)", y=0.995, fontsize=12)
fig.tight_layout()
fig.savefig("morans_perpc_grnonly.png", dpi=160)
plt.close(fig)
print("wrote morans_perpc_grnonly.png")
