import anndata as ad
import pandas as pd
import time

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

PREDICTED_PATH = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad_spatial/SpaGE_predicted_HVGs_FULL_50PV.h5ad"
SPATIAL_PATH   = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
HVG_PATH       = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/hvg_names_sub.csv"
OUTPUT_PATH    = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad_spatial/SpaGE_500_HVGs_FULL.h5ad"

print("=" * 60)
print("04_merge_spage_predicted")
print("Start:", now())
print("=" * 60)

# ── Load predicted HVGs ───────────────────────────────────────────────────────

print(f"\nLoading predicted HVGs: {PREDICTED_PATH}")
adata_predicted = ad.read_h5ad(PREDICTED_PATH)
print(f"  Shape: {adata_predicted.n_obs:,} cells x {adata_predicted.n_vars:,} genes")

# ── Load original MERFISH ─────────────────────────────────────────────────────

print(f"\nLoading original MERFISH: {SPATIAL_PATH}")
adata_spatial = ad.read_h5ad(SPATIAL_PATH)
print(f"  Shape:    {adata_spatial.n_obs:,} cells x {adata_spatial.n_vars:,} genes")
print(f"  Sections: {adata_spatial.obs['Section'].nunique() if 'Section' in adata_spatial.obs else 'N/A'}")
print(f"  Donors:   {adata_spatial.obs['Donor ID'].nunique() if 'Donor ID' in adata_spatial.obs else 'N/A'}")

# ── Remove Blank genes ────────────────────────────────────────────────────────

n_before_blank = adata_spatial.n_vars
blank_mask = ~adata_spatial.var_names.str.startswith('Blank')
n_blanks = (~blank_mask).sum()
adata_spatial = adata_spatial[:, blank_mask]
print(f"\nBlank gene removal:")
print(f"  Genes before: {n_before_blank:,}")
print(f"  Blank genes removed: {n_blanks:,}")
print(f"  Genes after:  {adata_spatial.n_vars:,}")

# ── HVG accounting ────────────────────────────────────────────────────────────

hvg_names = pd.read_csv(HVG_PATH)['Gene Name'].tolist()
print(f"\nHVG list: {HVG_PATH}")
print(f"  Total HVGs: {len(hvg_names):,}")

already_measured = [g for g in hvg_names if g in adata_spatial.var_names]
predicted_genes  = [g for g in hvg_names if g in adata_predicted.var_names]
in_neither       = [g for g in hvg_names if g not in adata_spatial.var_names and g not in adata_predicted.var_names]
overlap          = set(already_measured) & set(predicted_genes)

print(f"\nHVG breakdown:")
print(f"  Already measured in MERFISH: {len(already_measured):,}")
print(f"  Predicted by SpaGE:          {len(predicted_genes):,}")
print(f"  In neither (missing):        {len(in_neither):,}")
print(f"  Overlap (measured+predicted): {len(overlap):,}  (should be 0)")
if in_neither:
    print(f"  Missing genes: {in_neither}")

adata_measured = adata_spatial[:, already_measured].copy()
print(f"\nMeasured subset shape: {adata_measured.n_obs:,} x {adata_measured.n_vars:,}")

# ── Cell count consistency check ──────────────────────────────────────────────

print(f"\nCell count check:")
print(f"  MERFISH cells:    {adata_spatial.n_obs:,}")
print(f"  Predicted cells:  {adata_predicted.n_obs:,}")
if adata_spatial.n_obs != adata_predicted.n_obs:
    print("  WARNING: cell counts do not match!")
else:
    print("  OK: cell counts match")

# ── Merge ─────────────────────────────────────────────────────────────────────

print("\nMerging measured + predicted HVGs...")
adata_full = ad.concat([adata_measured, adata_predicted], axis=1, merge='same')
adata_full.obs  = adata_spatial.obs
adata_full.obsm = adata_spatial.obsm

print(f"  Merged shape: {adata_full.n_obs:,} cells x {adata_full.n_vars:,} genes")

# ── Save ──────────────────────────────────────────────────────────────────────

print(f"\nSaving to: {OUTPUT_PATH}")
adata_full.write_h5ad(OUTPUT_PATH)

print("\n" + "=" * 60)
print("SUMMARY")
print(f"  MERFISH input cells:          {adata_spatial.n_obs:,}")
print(f"  MERFISH genes (post-blank):   {adata_spatial.n_vars:,}")
print(f"  Blank genes removed:          {n_blanks:,}")
print(f"  HVGs already measured:        {len(already_measured):,}")
print(f"  HVGs from SpaGE prediction:   {len(predicted_genes):,}")
print(f"  HVGs missing entirely:        {len(in_neither):,}")
print(f"  Final shape:                  {adata_full.n_obs:,} x {adata_full.n_vars:,}")
print(f"  Output path:                  {OUTPUT_PATH}")
print(f"End: {now()}")
print("=" * 60)
