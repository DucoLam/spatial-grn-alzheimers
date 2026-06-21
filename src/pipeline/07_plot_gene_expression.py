import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import os
import time

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

# ── Config ────────────────────────────────────────────────────────────────────

PATH           = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad_spatial/SpaGE_500_HVGs_FULL.h5ad"
GENE           = "GFAP"
OBSM_KEY       = "X_spatial_raw"
OUT_DIR        = "section_plots"
DPI            = 300
EXPR_THRESHOLD = 0.0

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("07_plot_gene_expression")
print(f"Gene: {GENE}  |  DPI: {DPI}  |  Threshold: {EXPR_THRESHOLD}")
print("Start:", now())
print("=" * 60)

# ── Load ──────────────────────────────────────────────────────────────────────

print(f"\nLoading data: {PATH}")
adata = ad.read_h5ad(PATH)
print(f"  Shape:    {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"  Sections: {adata.obs['Section'].nunique():,}")
print(f"  Donors:   {adata.obs['Donor ID'].nunique() if 'Donor ID' in adata.obs else 'N/A':,}")

if GENE not in adata.var_names:
    raise ValueError(f"Gene '{GENE}' not found in dataset. Available genes (first 20): {list(adata.var_names[:20])}")

# ── Gene expression stats ─────────────────────────────────────────────────────

gene_idx = adata.var_names.get_loc(GENE)
expr = adata.X[:, gene_idx]
if hasattr(expr, 'toarray'):
    expr = expr.toarray().flatten()
expr = np.array(expr, dtype=np.float32)

n_expressing    = (expr > EXPR_THRESHOLD).sum()
pct_expressing  = n_expressing / len(expr) * 100
vmin = float(expr[expr > EXPR_THRESHOLD].min()) if n_expressing > 0 else 0
vmax = float(expr[expr > EXPR_THRESHOLD].max()) if n_expressing > 0 else 0

print(f"\n{GENE} expression stats (across all cells):")
print(f"  Total cells:         {len(expr):,}")
print(f"  Expressing cells:    {n_expressing:,}  ({pct_expressing:.1f}%)")
print(f"  Not expressing:      {(len(expr) - n_expressing):,}  ({100 - pct_expressing:.1f}%)")
print(f"  Min (expressed):     {vmin:.4f}")
print(f"  Max:                 {vmax:.4f}")
print(f"  Mean (all cells):    {expr.mean():.4f}")
print(f"  Mean (expressing):   {expr[expr > EXPR_THRESHOLD].mean():.4f}" if n_expressing > 0 else "  Mean (expressing):   N/A")
print(f"  Colorbar range:      [{vmin:.4f}, {vmax:.4f}]")

# ── Plot per section ──────────────────────────────────────────────────────────

sections = sorted(adata.obs['Section'].unique())
print(f"\nPlotting {len(sections)} sections to '{OUT_DIR}/'...")

for i, section in enumerate(sections):
    mask      = adata.obs['Section'] == section
    adata_sec = adata[mask]
    expr_sec  = expr[mask]

    n_sec_total      = mask.sum()
    n_sec_expressing = (expr_sec > EXPR_THRESHOLD).sum()
    pct_sec          = n_sec_expressing / n_sec_total * 100

    print(f"  [{i+1:3d}/{len(sections)}] {section}:  {n_sec_total:,} cells, "
          f"{n_sec_expressing:,} expressing ({pct_sec:.1f}%)")

    coords = adata_sec.obsm[OBSM_KEY]
    x, y   = coords[:, 0], coords[:, 1]

    obs      = adata_sec.obs.iloc[0]
    donor_id = obs['Donor ID']
    braak    = obs['Braak']
    apoe4    = obs['APOE4 Status']
    cps      = obs['Continuous Pseudo-progression Score']
    sex      = obs['Sex']
    age      = obs['Age at Death']

    apoe4_labels = {'N': 'Non-carrier', 'Y': 'Carrier', 'homo': 'Homozygous',
                    '0': 'Non-carrier (0)', '1': 'Heterozygous (1)', '2': 'Homozygous (2)'}
    apoe4_str = apoe4_labels.get(str(apoe4), str(apoe4))

    expressed     = expr_sec >  EXPR_THRESHOLD
    not_expressed = expr_sec <= EXPR_THRESHOLD

    fig, ax = plt.subplots(figsize=(14, 14))

    if not_expressed.sum() > 0:
        ax.scatter(x[not_expressed], y[not_expressed],
                   c='#cccccc', s=5.0, linewidths=0, alpha=0.4, rasterized=True)

    if expressed.sum() > 0:
        sc = ax.scatter(x[expressed], y[expressed],
                        c=expr_sec[expressed], cmap='Reds',
                        norm=mcolors.Normalize(vmin=vmin, vmax=vmax),
                        s=5.0, linewidths=0, alpha=0.7, rasterized=True)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.01)
        cbar.set_label(f"{GENE} expression (log1p normalized)", fontsize=12)

    title = (
        f"{GENE} expression — Section: {section}\n"
        f"Donor: {donor_id}  |  Sex: {sex}  |  Age: {age}\n"
        f"Braak: {braak} (tau spread, 0–6)  |  APOE4: {apoe4_str}  |  AD Progression Score: {cps:.3f}\n"
        f"Cells: {n_sec_total:,}  |  Expressing: {n_sec_expressing:,} ({pct_sec:.1f}%)"
    )
    ax.set_title(title, fontsize=11, pad=12)
    ax.set_xlabel("x (µm)", fontsize=11)
    ax.set_ylabel("y (µm)", fontsize=11)
    ax.set_aspect('equal')
    ax.legend(handles=[Patch(facecolor='#cccccc', label='Not expressed')],
              fontsize=10, loc='lower right')

    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, f"{section}_{GENE}.png")
    plt.savefig(out_path, dpi=DPI, bbox_inches='tight')
    plt.close()

print("\n" + "=" * 60)
print("SUMMARY")
print(f"  Gene:                 {GENE}")
print(f"  Total cells:          {len(expr):,}")
print(f"  Expressing cells:     {n_expressing:,}  ({pct_expressing:.1f}%)")
print(f"  Expression range:     [{vmin:.4f}, {vmax:.4f}]")
print(f"  Sections plotted:     {len(sections)}")
print(f"  Output directory:     {OUT_DIR}/")
print(f"End: {now()}")
print("=" * 60)
