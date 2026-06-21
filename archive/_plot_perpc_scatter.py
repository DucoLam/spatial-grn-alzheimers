"""Per-PC scatter plots for SQ3: for each subtype, the top-5 spatial PCs,
each panel = per-donor Moran's I (y) vs CPS (x), with Spearman r/p/q in title.
Reads sq3_perpc_morans.csv + sq3_perpc_corr.csv produced by sq3_spatial_severity.py.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ3"
SLUG = {"Astrocyte": "astrocyte", "L2/3 IT": "l23_it",
        "L4 IT": "l4_it", "Oligodendrocyte": "oligo"}
SUBTYPES = ["Astrocyte", "L2/3 IT", "L4 IT", "Oligodendrocyte"]

pp = pd.read_csv(os.path.join(D, "sq3_perpc_morans.csv"))
corr = pd.read_csv(os.path.join(D, "sq3_perpc_corr.csv"))

for st in SUBTYPES:
    cs = corr[(corr.subtype == st) & (corr.in_topk)].copy()
    if cs.empty:
        continue
    # order panels by spatial strength (pooled_I), descending
    cs = cs.sort_values("pooled_I", ascending=False)
    pcs = cs["pc"].tolist()
    n = len(pcs)
    fig, axs = plt.subplots(1, n, figsize=(3.4 * n, 3.4), squeeze=False)
    axs = axs[0]
    d_st = pp[pp.subtype == st]
    for ax, pc in zip(axs, pcs):
        sub = d_st[d_st.pc == pc].dropna(subset=["morans_I", "CPS"])
        ax.scatter(sub["CPS"], sub["morans_I"], s=55, zorder=3, color="#c0392b")
        for _, row in sub.iterrows():
            ax.annotate(str(row["donor"])[-6:], (row["CPS"], row["morans_I"]),
                        fontsize=6, xytext=(3, 3), textcoords="offset points")
        c = cs[cs.pc == pc].iloc[0]
        ax.set_title(f"PC{pc} (pooled I={c['pooled_I']:.2f})\n"
                     f"r={c['spearman_r_vs_CPS']:.2f} p={c['p']:.3f} q={c['fdr_q']:.2f}",
                     fontsize=9)
        ax.set_xlabel("CPS")
        ax.axhline(0, color="0.6", lw=0.7)
    axs[0].set_ylabel("Per-donor Moran's I")
    fig.suptitle(f"SQ3 per-PC — {st}  (top-5 spatial PCs, N={int(cs.iloc[0]['N'])})",
                 fontsize=11)
    fig.tight_layout()
    out = os.path.join(D, f"sq3_perpc_scatter_{SLUG[st]}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)
print("done")
