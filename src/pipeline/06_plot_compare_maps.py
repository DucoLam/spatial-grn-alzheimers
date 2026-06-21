"""
06_plot_compare_maps.py -- per-section INFERRED vs ORIGINAL maps, side by side.

For each overlapping tissue section, two panels in the same coordinate frame:
  LEFT  = INFERRED: our RNA cells placed by Tangram (argmax), coloured by subtype.
  RIGHT = ORIGINAL: the real MERFISH cells of those same 4 subtypes, coloured by
          their true Subclass.
Grey = all MERFISH cells of the section (tissue outline). One figure per donor.
"""
import time
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

TG_DIR  = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/data/tangram")
MERFISH = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
OUT_DIR = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/plots/compare_maps")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUB = {"astrocyte": "Astrocyte", "l23_it": "L2/3 IT", "l4_it": "L4 IT", "oligo": "Oligodendrocyte"}
COL = {"Astrocyte": "#d62728", "L2/3 IT": "#1f77b4", "L4 IT": "#2ca02c", "Oligodendrocyte": "#9467bd"}
BG = 8000

print("start", time.strftime("%H:%M:%S"))
dfs = []
for p, l in SUB.items():
    f = TG_DIR / f"Tangram_cell_locations_{p}.csv"
    if f.exists():
        d = pd.read_csv(f); d["subtype"] = l; dfs.append(d)
inf = pd.concat(dfs, ignore_index=True)
inf["sec"] = inf["tangram_mapped_section"].astype(str)
inf["donor"] = inf["sec"].map(lambda s: ".".join(s.split(".")[:3]))
print("inferred located cells:", len(inf))

sp = ad.read_h5ad(MERFISH, backed="r")
ckey = next(k for k in ("X_spatial_raw", "spatial", "X_spatial") if k in sp.obsm)
coords = np.asarray(sp.obsm[ckey])
sec_arr = sp.obs["Section"].astype(str).values
subc = sp.obs["Subclass"].astype(str).values
rng = np.random.default_rng(0)

def style(ax, title):
    ax.set_title(title, fontsize=8); ax.set_aspect("equal")
    ax.grid(True, ls=":", lw=0.4, color="grey", alpha=0.5)
    ax.tick_params(labelsize=6); ax.set_xlabel("x (um)", fontsize=6); ax.set_ylabel("y (um)", fontsize=6)
    ax.legend(fontsize=5, markerscale=2, loc="upper right")

for d in sorted(inf["donor"].unique()):
    sub = inf[inf["donor"] == d]
    secs = sorted(sub["sec"].unique())
    fig, axes = plt.subplots(len(secs), 2, figsize=(11, 5*len(secs)), squeeze=False)
    for r, s in enumerate(secs):
        bgm = sec_arr == s
        bg = coords[bgm]
        if len(bg) > BG:
            bg = bg[rng.choice(len(bg), BG, replace=False)]
        short = ".".join(s.split(".")[-3:])

        # LEFT: inferred
        ax = axes[r][0]
        ax.scatter(bg[:, 0], bg[:, 1], s=1, c="lightgrey", linewidths=0, rasterized=True)
        ss = sub[sub["sec"] == s]
        for l in SUB.values():
            q = ss[ss["subtype"] == l]
            if len(q):
                ax.scatter(q["tangram_x"], q["tangram_y"], s=7, c=COL[l], linewidths=0, alpha=0.85, label=f"{l} ({len(q)})")
        style(ax, f"INFERRED  {short}")

        # RIGHT: original (true MERFISH subclass)
        ax = axes[r][1]
        ax.scatter(bg[:, 0], bg[:, 1], s=1, c="lightgrey", linewidths=0, rasterized=True)
        for l in SUB.values():
            m = bgm & (subc == l)
            c2 = coords[m]
            if len(c2):
                ax.scatter(c2[:, 0], c2[:, 1], s=7, c=COL[l], linewidths=0, alpha=0.6, label=f"{l} ({len(c2)})")
        style(ax, f"ORIGINAL (true MERFISH)  {short}")

    fig.suptitle(f"Inferred vs Original — donor {d}", fontsize=12)
    fig.tight_layout()
    out = OUT_DIR / f"compare_{d}.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out.name}  ({len(secs)} sections)")
print("done", time.strftime("%H:%M:%S"))
