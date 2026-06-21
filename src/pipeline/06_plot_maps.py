"""
06_plot_maps.py -- spatial maps of inferred cell locations, per tissue section.

For each MERFISH section (= one tissue sample), plots:
  - grey background = that section's real MERFISH cells (tissue outline)
  - coloured points = our RNA cells placed there by Tangram (argmax coords),
    coloured by cell subtype.

One figure per donor (its sections side by side). Uses the production
Tangram_cell_locations_<prefix>.csv files (argmax = tangram_x/tangram_y).
"""
import time
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

TG_DIR   = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/data/tangram")
MERFISH  = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
OUT_DIR  = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/plots/maps")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUBTYPES = {"astrocyte": "Astrocyte", "l23_it": "L2/3 IT", "l4_it": "L4 IT", "oligo": "Oligodendrocyte"}
COLORS   = {"Astrocyte": "#d62728", "L2/3 IT": "#1f77b4", "L4 IT": "#2ca02c", "Oligodendrocyte": "#9467bd"}
BG_MAX   = 8000   # background cells plotted per section

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

print("06_plot_maps start:", now())

# ---- load + combine inferred locations -------------------------------------
dfs = []
for prefix, label in SUBTYPES.items():
    f = TG_DIR / f"Tangram_cell_locations_{prefix}.csv"
    if f.exists():
        d = pd.read_csv(f)
        d["subtype"] = label
        dfs.append(d)
        print(f"  loaded {f.name}: {len(d):,} cells")
    else:
        print(f"  MISSING {f.name}")
df = pd.concat(dfs, ignore_index=True)
if "tangram_mapped_section" not in df.columns:
    raise SystemExit("No tangram_mapped_section column in CSVs.")
df["tangram_mapped_section"] = df["tangram_mapped_section"].astype(str)
print(f"  total inferred-located cells: {len(df):,}")

# ---- MERFISH background coords ---------------------------------------------
print("Loading MERFISH (backed) for tissue background ...")
sp = ad.read_h5ad(MERFISH, backed="r")
ckey = next(k for k in ("X_spatial_raw", "spatial", "X_spatial") if k in sp.obsm)
bg_coords = np.asarray(sp.obsm[ckey])
bg_sec = sp.obs["Section"].astype(str).values
print(f"  background coords {bg_coords.shape}, key '{ckey}'")

def donor_of(section):
    return ".".join(section.split(".")[:3])

df["donor"] = df["tangram_mapped_section"].map(donor_of)
donors = sorted(df["donor"].unique())
print(f"Donors with mapped cells: {len(donors)}")

rng = np.random.default_rng(0)

for donor in donors:
    sub = df[df["donor"] == donor]
    secs = sorted(sub["tangram_mapped_section"].unique())
    ncol = min(len(secs), 4)
    nrow = int(np.ceil(len(secs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5*ncol, 5*nrow), squeeze=False)
    axes = axes.ravel()
    for ax, sec in zip(axes, secs):
        bgm = bg_sec == sec
        bg = bg_coords[bgm]
        if len(bg) > BG_MAX:
            bg = bg[rng.choice(len(bg), BG_MAX, replace=False)]
        ax.scatter(bg[:, 0], bg[:, 1], s=1, c="lightgrey", linewidths=0, rasterized=True)
        s = sub[sub["tangram_mapped_section"] == sec]
        for label in SUBTYPES.values():
            ss = s[s["subtype"] == label]
            if len(ss):
                ax.scatter(ss["tangram_x"], ss["tangram_y"], s=6, c=COLORS[label],
                           label=f"{label} ({len(ss)})", linewidths=0, alpha=0.8)
        ax.set_title(".".join(sec.split(".")[-3:]), fontsize=8)
        ax.set_aspect("equal")
        # coordinate grid + ticks so positions are readable (um)
        ax.grid(True, which="both", ls=":", lw=0.5, color="grey", alpha=0.6)
        ax.minorticks_on()
        ax.tick_params(axis="both", which="major", labelsize=6, length=4)
        ax.tick_params(axis="both", which="minor", length=2)
        ax.set_xlabel("x (um)", fontsize=7)
        ax.set_ylabel("y (um)", fontsize=7)
        ax.legend(fontsize=6, markerscale=2, loc="upper right")
    for ax in axes[len(secs):]:
        ax.axis("off")
    fig.suptitle(f"Inferred cell locations — donor {donor}", fontsize=12)
    fig.tight_layout()
    out = OUT_DIR / f"map_{donor}.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"  wrote {out.name}  ({len(secs)} sections, {len(sub):,} cells)")

print("06_plot_maps done:", now())
