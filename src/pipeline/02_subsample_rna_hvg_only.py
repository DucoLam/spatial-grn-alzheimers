import anndata as ad
import pandas as pd
import time

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

HVG_PATH         = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/hvg_names_sub.csv"
SPATIAL_HVG_PATH = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/spatial_hvg_names_sub.csv"
RNA_PATH         = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
OUTPUT_PATH      = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/SEAAD_MTG_RNAseq_HVG_SPATIAL_ONLY.h5ad"

print("=" * 60)
print("02_subsample_rna_hvg_only")
print("Start:", now())
print("=" * 60)

# ── Load HVG lists ────────────────────────────────────────────────────────────

print(f"\nLoading HVG list from:         {HVG_PATH}")
hvg_file = pd.read_csv(HVG_PATH)
print(f"  HVGs in list:                {len(hvg_file):,}")

print(f"Loading spatial HVG list from: {SPATIAL_HVG_PATH}")
hvg_file_spatial = pd.read_csv(SPATIAL_HVG_PATH)
print(f"  Spatial HVGs in list:        {len(hvg_file_spatial):,}")

genes_to_keep_raw = hvg_file['Gene Name'].tolist() + hvg_file_spatial['Gene Name'].tolist()
genes_to_keep = list(dict.fromkeys(genes_to_keep_raw))
print(f"\nCombined unique genes to keep: {len(genes_to_keep):,}")
print(f"  (deduplicated from {len(genes_to_keep_raw):,} total)")

# ── Load RNA ──────────────────────────────────────────────────────────────────

print(f"\nLoading RNA from: {RNA_PATH}")
adata = ad.read_h5ad(RNA_PATH)
print(f"  Input shape:   {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# ── Filter genes ──────────────────────────────────────────────────────────────

genes_present = [g for g in genes_to_keep if g in adata.var_names]
missing_genes = set(genes_to_keep) - set(genes_present)

print(f"\nGene filtering:")
print(f"  Requested:   {len(genes_to_keep):,}")
print(f"  Found:       {len(genes_present):,}")
print(f"  Missing:     {len(missing_genes):,}")
if missing_genes:
    print(f"  Missing list: {sorted(missing_genes)}")

# ── Subset ────────────────────────────────────────────────────────────────────

adata_subset = adata[:, genes_present].copy()
print(f"\nOutput shape:  {adata_subset.n_obs:,} cells x {adata_subset.n_vars:,} genes")

# ── Save ──────────────────────────────────────────────────────────────────────

print(f"\nSaving to: {OUTPUT_PATH}")
adata_subset.write_h5ad(OUTPUT_PATH)

print("\n" + "=" * 60)
print("SUMMARY")
print(f"  Input cells:      {adata.n_obs:,}")
print(f"  Input genes:      {adata.n_vars:,}")
print(f"  Output cells:     {adata_subset.n_obs:,}")
print(f"  Output genes:     {adata_subset.n_vars:,}")
print(f"  Missing HVGs:     {len(missing_genes):,}")
print(f"  Output path:      {OUTPUT_PATH}")
print(f"End: {now()}")
print("=" * 60)
