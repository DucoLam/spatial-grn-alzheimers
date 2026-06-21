"""Subsample all astrocytes from SEA-AD paired HVG/HVP files.

Reads the shared HVG/HVP files, filters to Astrocyte cell type (all 6,777 cells),
computes a k=20 KNN on the joint WNN-input embedding (obsm['X_pca']), and writes
matched RNA/ATAC files to dflam's data directory.

Outputs:
    dflam/data/seaad_paired_rna_astrocytes.h5ad
    dflam/data/seaad_paired_atac_astrocytes.h5ad
"""

import anndata as ad
import numpy as np
from sklearn.neighbors import NearestNeighbors
import time

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

RNA_IN   = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/seaad_paired_rna_hvg.h5ad"
ATAC_IN  = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/seaad_paired_atac_hvp.h5ad"
RNA_OUT  = "/tudelft.net/staff-umbrella/ScReNI/dflam/data/seaad_paired_rna_astrocytes.h5ad"
ATAC_OUT = "/tudelft.net/staff-umbrella/ScReNI/dflam/data/seaad_paired_atac_astrocytes.h5ad"
K        = 20

print("=" * 60)
print("08_subsample_astrocytes")
print(f"Selecting all Astrocyte cells, k={K} KNN")
print("Start:", now())
print("=" * 60)

# Load RNA
print(f"\nLoading RNA: {RNA_IN}")
rna = ad.read_h5ad(RNA_IN)
print(f"  Full shape: {rna.shape}")

# Verify barcodes align before loading ATAC
print(f"\nLoading ATAC: {ATAC_IN}")
atac = ad.read_h5ad(ATAC_IN)
print(f"  Full shape: {atac.shape}")

if not (rna.obs_names == atac.obs_names).all():
    raise RuntimeError("RNA and ATAC obs_names are not aligned — cannot subsample safely.")
print("  Barcode alignment: OK")

# Filter to astrocytes
mask = rna.obs["cell_type"] == "Astrocyte"
rna_astro  = rna[mask].copy()
atac_astro = atac[mask].copy()
print(f"\nAstrocytes selected: {rna_astro.n_obs:,}")
print(f"  Donors: {rna_astro.obs['Donor ID'].nunique()}")
print(f"  Subtypes: {rna_astro.obs['Subclass'].value_counts().to_dict() if 'Subclass' in rna_astro.obs.columns else 'N/A'}")

assert (rna_astro.obs_names == atac_astro.obs_names).all(), "Post-filter barcode mismatch"

# Compute KNN on joint embedding
print(f"\nComputing k={K} KNN on obsm['X_pca'] ({rna_astro.obsm['X_pca'].shape})...")
t0 = time.time()
nn = NearestNeighbors(n_neighbors=K, n_jobs=-1).fit(rna_astro.obsm["X_pca"])
_, knn_indices = nn.kneighbors(rna_astro.obsm["X_pca"])
knn_indices = knn_indices.astype(np.int64)
print(f"  Done in {time.time()-t0:.1f}s  |  shape: {knn_indices.shape}")

# Store KNN in both files
rna_astro.uns["knn_indices"]  = knn_indices
atac_astro.uns["knn_indices"] = knn_indices

# Remove any stale full-dataset KNN
for key in ["wnn_neighbor_indices", "neighbors"]:
    rna_astro.uns.pop(key, None)
    atac_astro.uns.pop(key, None)

# Save
print(f"\nSaving RNA:  {RNA_OUT}")
rna_astro.write_h5ad(RNA_OUT)
print(f"Saving ATAC: {ATAC_OUT}")
atac_astro.write_h5ad(ATAC_OUT)

print(f"\n{'='*60}")
print("SUMMARY")
print(f"  Cells:   {rna_astro.n_obs:,}  (all Astrocyte cells)")
print(f"  Genes:   {rna_astro.n_vars:,}  HVGs")
print(f"  Peaks:   {atac_astro.n_vars:,}  HVPs")
print(f"  KNN k:   {K}")
print(f"  RNA out: {RNA_OUT}")
print(f"  ATAC out:{ATAC_OUT}")
print("End:", now())
print("=" * 60)
