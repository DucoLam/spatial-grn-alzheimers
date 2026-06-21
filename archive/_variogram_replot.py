"""Replot variogram_all.png from the saved per-subtype curves (no recompute)."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import analysis_common as ac

VDIR = ac.ANALYSIS / "SQ2" / "variogram"
# slug -> (display name, colour)  -- same per-subtype colours as the SQ4 scatter
subs = [("l23_it", "L2/3 IT", "#1f77b4"),
        ("l4_it", "L4 IT", "#2ca02c"),
        ("astrocyte", "Astrocyte", "#d62728"),
        ("oligo", "Oligodendrocyte", "#ff7f0e")]

fig, ax = plt.subplots(figsize=(7.2, 3.3))   # shorter
for slug, name, col in subs:
    df = pd.read_csv(VDIR / f"variogram_{slug}.csv").dropna(subset=["grn_mean_dissim"])
    ax.errorbar(df.dist_center, df.grn_mean_dissim, yerr=df.grn_sem,
                fmt="-o", ms=3.5, lw=1.6, color=col, label=name,
                elinewidth=0.7, capsize=2, alpha=0.95)

ax.set_axisbelow(True)
ax.minorticks_on()
ax.grid(True, which="major", ls=":", lw=0.7, color="0.55", alpha=0.9)
ax.grid(True, which="minor", ls=":", lw=0.4, color="0.72", alpha=0.7)
ax.tick_params(which="major", direction="out", length=5, labelsize=9)
ax.tick_params(which="minor", direction="out", length=2.5)
ax.set_xlabel(r"Within-section Spatial Distance ($\mu$m)", fontsize=11)
ax.set_ylabel("Mean GRN Dissimilarity", fontsize=11)
ax.set_title("Spatial GRN Variogram for 4 Selected Subtypes", fontsize=12)
ax.legend(frameon=False, fontsize=9.5, ncol=2)
fig.tight_layout()
fig.savefig(VDIR / "variogram_all.png", dpi=150)
print("saved", VDIR / "variogram_all.png")
