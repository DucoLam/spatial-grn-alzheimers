import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import os
import time

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

# ── Config ────────────────────────────────────────────────────────────────────

MAPPING_PATH   = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/Tangram_rna_with_coords_matched_140.h5ad"
GENE           = "GFAP"
OUT_DIR        = "/tudelft.net/staff-umbrella/ScReNI/dflam/plots/07b_inferred_locations_140"
DPI            = 300
EXPR_THRESHOLD = 0.0

# Tangram params used — shown in plot titles
TANGRAM_GENES  = 140
TANGRAM_TOPK   = "0.5%"
TANGRAM_EPOCHS = 1000

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("07b_plot_inferred_locations")
print(f"Gene:    {GENE}")
print(f"Tangram: {TANGRAM_GENES} genes (140 native MERFISH genes), top {TANGRAM_TOPK}, {TANGRAM_EPOCHS} epochs")
print(f"Out:     {OUT_DIR}")
print("Start:", now())
print("=" * 60)

# ── Load ──────────────────────────────────────────────────────────────────────

print(f"\nLoading: {MAPPING_PATH}")
adata = ad.read_h5ad(MAPPING_PATH)
print(f"  Shape:   {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"  obs cols: {list(adata.obs.columns)}")

if "tangram_x" not in adata.obs.columns or "tangram_y" not in adata.obs.columns:
    raise ValueError("tangram_x / tangram_y not found in obs. Run 06b first.")

if GENE not in adata.var_names:
    raise ValueError(f"Gene '{GENE}' not found. Available (first 20): {list(adata.var_names[:20])}")

# ── Extract gene expression ───────────────────────────────────────────────────

gene_idx = adata.var_names.get_loc(GENE)
expr = adata.X[:, gene_idx]
if hasattr(expr, 'toarray'):
    expr = expr.toarray().flatten()
expr = np.array(expr, dtype=np.float32)

n_expressing   = (expr > EXPR_THRESHOLD).sum()
pct_expressing = n_expressing / len(expr) * 100
vmin = float(expr[expr > EXPR_THRESHOLD].min()) if n_expressing > 0 else 0
vmax = float(expr[expr > EXPR_THRESHOLD].max()) if n_expressing > 0 else 0

print(f"\n{GENE} expression stats:")
print(f"  Total cells:      {len(expr):,}")
print(f"  Expressing:       {n_expressing:,}  ({pct_expressing:.1f}%)")
print(f"  Range (expr):     [{vmin:.4f}, {vmax:.4f}]")

# ── Subtype colormap ──────────────────────────────────────────────────────────

col_type   = next((c for c in ("Supertype", "cell_type", "celltype") if c in adata.obs.columns), None)
subtypes   = sorted(adata.obs[col_type].astype(str).unique()) if col_type else []
cmap_sub   = plt.get_cmap("tab20", max(len(subtypes), 1))
subtype_colors = {st: cmap_sub(i) for i, st in enumerate(subtypes)}

# ── Plot per donor ────────────────────────────────────────────────────────────

donors = sorted(adata.obs["Donor ID"].astype(str).unique())
print(f"\nPlotting {len(donors)} donors to '{OUT_DIR}/'...")

for i, donor in enumerate(donors):
    mask      = adata.obs["Donor ID"].astype(str) == donor
    adata_d   = adata[mask]
    expr_d    = expr[mask]
    n_cells   = mask.sum()

    x = adata_d.obs["tangram_x"].values.astype(float)
    y = adata_d.obs["tangram_y"].values.astype(float)

    obs0     = adata_d.obs.iloc[0]
    braak    = obs0.get("Braak",    "N/A")
    apoe4    = obs0.get("APOE4 Status", "N/A")
    cps      = obs0.get("Continuous Pseudo-progression Score", None)
    sex      = obs0.get("Sex", "N/A")
    age      = obs0.get("Age at Death", "N/A")
    cps_str  = f"{cps:.3f}" if cps is not None and not pd.isna(cps) else "N/A"

    n_expr   = (expr_d > EXPR_THRESHOLD).sum()
    pct_expr = n_expr / n_cells * 100

    print(f"  [{i+1}/{len(donors)}] {donor}:  {n_cells:,} cells,  {n_expr:,} expressing {GENE} ({pct_expr:.1f}%)")

    fig, axes = plt.subplots(1, 2, figsize=(24, 11))

    # ── Left: gene expression ─────────────────────────────────────────────────
    ax = axes[0]
    expressed     = expr_d >  EXPR_THRESHOLD
    not_expressed = expr_d <= EXPR_THRESHOLD

    if not_expressed.sum() > 0:
        ax.scatter(x[not_expressed], y[not_expressed],
                   c='#cccccc', s=8, linewidths=0, alpha=0.4, rasterized=True)
    if expressed.sum() > 0:
        sc = ax.scatter(x[expressed], y[expressed],
                        c=expr_d[expressed], cmap='Reds',
                        norm=mcolors.Normalize(vmin=vmin, vmax=vmax),
                        s=8, linewidths=0, alpha=0.8, rasterized=True)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.01)
        cbar.set_label(f"{GENE} expression (log1p normalized)", fontsize=11)

    ax.set_title(f"{GENE} expression — Inferred locations\n"
                 f"Tangram: {TANGRAM_GENES} genes (140 native MERFISH genes), "
                 f"top {TANGRAM_TOPK} weighted centroid, {TANGRAM_EPOCHS} epochs",
                 fontsize=10, pad=8)
    ax.set_xlabel("Inferred x (µm)", fontsize=11)
    ax.set_ylabel("Inferred y (µm)", fontsize=11)
    ax.set_aspect('equal')
    ax.legend(handles=[Patch(facecolor='#cccccc', label='Not expressed')],
              fontsize=9, loc='lower right')

    # ── Right: subtype coloring ───────────────────────────────────────────────
    ax2 = axes[1]
    if col_type:
        cell_subtypes = adata_d.obs[col_type].astype(str).values
        for st in subtypes:
            st_mask = cell_subtypes == st
            if st_mask.sum() > 0:
                ax2.scatter(x[st_mask], y[st_mask],
                            c=[subtype_colors[st]], s=8, linewidths=0,
                            alpha=0.7, rasterized=True, label=st)
        legend_handles = [Line2D([0], [0], marker='o', color='w',
                                 markerfacecolor=subtype_colors[st],
                                 markersize=7, label=st) for st in subtypes]
        ax2.legend(handles=legend_handles, fontsize=7, loc='lower right',
                   ncol=2, framealpha=0.8)
    ax2.set_title(f"Astrocyte subtypes — Inferred locations\n"
                  f"Tangram: {TANGRAM_GENES} genes (140 native MERFISH genes), "
                  f"top {TANGRAM_TOPK} weighted centroid, {TANGRAM_EPOCHS} epochs",
                  fontsize=10, pad=8)
    ax2.set_xlabel("Inferred x (µm)", fontsize=11)
    ax2.set_ylabel("Inferred y (µm)", fontsize=11)
    ax2.set_aspect('equal')

    # ── Shared donor info ─────────────────────────────────────────────────────
    apoe4_labels = {'N': 'Non-carrier', 'Y': 'Carrier', 'homo': 'Homozygous',
                    '0': 'Non-carrier', '1': 'Heterozygous', '2': 'Homozygous'}
    apoe4_str = apoe4_labels.get(str(apoe4), str(apoe4))

    fig.suptitle(
        f"Donor: {donor}  |  Sex: {sex}  |  Age: {age}  |  "
        f"Braak: {braak}  |  APOE4: {apoe4_str}  |  AD Score: {cps_str}\n"
        f"Cells: {n_cells:,}  |  {GENE} expressing: {n_expr:,} ({pct_expr:.1f}%)",
        fontsize=12, y=1.01
    )

    plt.tight_layout()
    safe_donor = donor.replace("/", "-").replace(" ", "_")
    out_path = os.path.join(OUT_DIR, f"{safe_donor}_{GENE}_inferred.png")
    plt.savefig(out_path, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {out_path}")

print("\n" + "=" * 60)
print("SUMMARY")
print(f"  Gene:             {GENE}")
print(f"  Tangram params:   {TANGRAM_GENES} genes, top {TANGRAM_TOPK}, {TANGRAM_EPOCHS} epochs")
print(f"  Donors plotted:   {len(donors)}")
print(f"  Output:           {OUT_DIR}/")
print(f"End: {now()}")
print("=" * 60)
