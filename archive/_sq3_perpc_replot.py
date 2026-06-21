"""Replot sq3_perpc_<subtype>.png from sq3_perpc_corr.csv (no recompute).
red = top-3 spatial PCs (by pooled Moran's I); star = p<0.05 (any PC)."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import analysis_common as ac

D = ac.ANALYSIS / "SQ3"
RED = "#c0392b"
ppc = pd.read_csv(D / "sq3_perpc_corr.csv")

for st in ppc.subtype.unique():
    sub = ppc[ppc.subtype == st].sort_values("pc")
    if sub.empty:
        continue
    colors = [RED if r <= 5 else "0.8" for r in sub["spatial_rank"]]

    fig, ax = plt.subplots(figsize=(8.5, 4))
    ax.bar(sub["pc"], sub["spearman_r_vs_CPS"], color=colors,
           edgecolor="0.3", linewidth=0.4, zorder=3)

    # data-driven y-limits: fit to the actual bar range (always include 0),
    # leaving headroom for the bar tops and the significance stars
    vals = sub["spearman_r_vs_CPS"].dropna()
    lo = min(0.0, float(vals.min()))
    hi = max(0.0, float(vals.max()))
    span = (hi - lo) or 1.0
    y_lo, y_hi = lo - 0.14 * span, hi + 0.12 * span

    # star ANY significant PC (p < 0.05)
    star_off = 0.03 * span
    for _, rr in sub.iterrows():
        if not np.isnan(rr["p"]) and rr["p"] < 0.05:
            y = rr["spearman_r_vs_CPS"]
            ax.text(rr["pc"], y + (star_off if y >= 0 else -3 * star_off), "*",
                    ha="center", va="bottom", fontsize=14, zorder=4)

    ax.axhline(0, color="k", lw=0.8, zorder=2)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xticks(sorted(sub["pc"].unique()))
    ax.tick_params(axis="x", labelsize=9)
    ax.tick_params(axis="y", labelsize=9)
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", ls=":", lw=0.7, color="0.6", alpha=0.85, zorder=0)
    ax.grid(True, axis="x", ls=":", lw=0.4, color="0.8", alpha=0.5, zorder=0)
    ax.set_xlabel("GRN principal component", fontsize=11)
    ax.set_ylabel("Spearman r (per-donor Moran's I vs CPS)", fontsize=11)
    ax.set_title(f"Spearman Correlation between {st} GRN PCs and CPS", fontsize=12)

    legend_handles = [
        Patch(facecolor=RED, edgecolor="0.3", label="top 5 spatial PCs"),
        Line2D([0], [0], marker="*", color="none", markerfacecolor="k",
               markeredgecolor="k", markersize=12, linestyle="None", label="p < 0.05"),
    ]
    ax.legend(handles=legend_handles, fontsize=10, frameon=True, framealpha=0.9,
              loc="best")
    fig.tight_layout()
    out = D / f"sq3_perpc_{ac.SLUG[st]}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("saved", out)
