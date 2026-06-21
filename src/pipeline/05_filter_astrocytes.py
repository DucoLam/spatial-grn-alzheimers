import anndata as ad
import time

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

RNA_PATH = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/seaad_paired_rna.h5ad"
RNA_OUT  = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/astrocytes_seaad_paired_rna.h5ad"

print("=" * 60)
print("05_filter_astrocytes")
print("Start:", now())
print("=" * 60)

def find_column(adata, candidates):
    for col in candidates:
        if col in adata.obs.columns:
            return col
    raise KeyError(f"No matching column found among {candidates}. Available: {list(adata.obs.columns)}")

# ── Load RNA ──────────────────────────────────────────────────────────────────

print(f"\nLoading scRNA-seq: {RNA_PATH}")
adata_rna = ad.read_h5ad(RNA_PATH)
print(f"  Shape:         {adata_rna.n_obs:,} cells x {adata_rna.n_vars:,} genes")
print(f"  Donors:        {adata_rna.obs['Donor ID'].nunique() if 'Donor ID' in adata_rna.obs else 'N/A'}")

col_type     = find_column(adata_rna, ("cell_type", "celltype", "cell_type_alias_label", "Supertype", "cluster"))
all_types    = adata_rna.obs[col_type].value_counts()
astro_labels = [v for v in adata_rna.obs[col_type].unique() if "astro" in str(v).lower()]

print(f"  Cell type col: {col_type}")
print(f"  Unique types:  {len(all_types):,}")
print(f"  All types:     {all_types.to_dict()}")
print(f"\n  Astrocyte labels found: {astro_labels}")
print(f"  Astrocyte cells:        {adata_rna.obs[col_type].isin(astro_labels).sum():,}")

# ── Filter ────────────────────────────────────────────────────────────────────

adata_astro = adata_rna[adata_rna.obs[col_type].isin(astro_labels)].copy()
print(f"\nAfter filter: {adata_astro.n_obs:,} cells ({adata_astro.n_obs/adata_rna.n_obs*100:.1f}% of input)")

if 'Donor ID' in adata_astro.obs:
    print(f"  Donors retained: {adata_astro.obs['Donor ID'].nunique():,}")

# ── Save ──────────────────────────────────────────────────────────────────────

print(f"\nSaving: {RNA_OUT}")
adata_astro.write_h5ad(RNA_OUT)

print("\n" + "=" * 60)
print("SUMMARY")
print(f"  Input:   {adata_rna.n_obs:,} cells x {adata_rna.n_vars:,} genes")
print(f"  Output:  {adata_astro.n_obs:,} astrocyte cells x {adata_astro.n_vars:,} genes")
print(f"  Labels:  {astro_labels}")
print(f"  Path:    {RNA_OUT}")
print(f"  → Feed into 06_run_tangram.py as RNA_PATH")
print(f"End: {now()}")
print("=" * 60)
