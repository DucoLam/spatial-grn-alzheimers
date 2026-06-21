"""
08_subsample.py  — parameterized cell type subsample
Usage:
    python 08_subsample.py --cell_type "Astrocyte"    --prefix astrocytes
    python 08_subsample.py --cell_type "L2/3 IT"      --prefix l23_it
    python 08_subsample.py --cell_type "L4 IT"        --prefix l4_it
    python 08_subsample.py --cell_type "Oligodendrocyte" --prefix oligodendrocyte
"""
import sys, argparse, time
import anndata as ad
import numpy as np
from sklearn.neighbors import NearestNeighbors
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--cell_type", required=True, help="Subclass label in seaad_paired_rna_hvg.h5ad")
parser.add_argument("--prefix",    required=True, help="Output file prefix (e.g. l23_it)")
parser.add_argument("--k",         type=int, default=20, help="KNN k (default 20)")
args = parser.parse_args()

BSC      = Path("/tudelft.net/staff-umbrella/ScReNI/bsc-screni")
DFLAM    = Path("/tudelft.net/staff-umbrella/ScReNI/dflam")
RNA_IN   = BSC  / "data/processed/seaad/seaad_paired_rna_hvg.h5ad"
ATAC_IN  = BSC  / "data/processed/seaad/seaad_paired_atac_hvp.h5ad"
RNA_OUT  = DFLAM / f"data/seaad_paired_rna_{args.prefix}.h5ad"
ATAC_OUT = DFLAM / f"data/seaad_paired_atac_{args.prefix}.h5ad"
(DFLAM / "data").mkdir(parents=True, exist_ok=True)

print("=" * 60)
print(f"08_subsample  —  cell_type={args.cell_type}  prefix={args.prefix}")
print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)

print(f"\nLoading RNA: {RNA_IN}")
rna = ad.read_h5ad(RNA_IN)
print(f"  Full shape: {rna.shape}")

print(f"\nLoading ATAC: {ATAC_IN}")
atac = ad.read_h5ad(ATAC_IN)
print(f"  Full shape: {atac.shape}")

# Filter by Subclass
col = "Subclass" if "Subclass" in rna.obs.columns else "cell_type"
mask = rna.obs[col] == args.cell_type
print(f"\nFiltering on obs[{col}] == {args.cell_type}")
print(f"  Cells selected: {mask.sum():,} / {rna.n_obs:,}")
if mask.sum() == 0:
    available = sorted(rna.obs[col].unique())
    print(f"  ERROR: no cells found. Available: {available}")
    sys.exit(1)

rna_sub  = rna[mask].copy()
atac_sub = atac[mask].copy()

# Verify barcode alignment
assert list(rna_sub.obs_names) == list(atac_sub.obs_names), "Barcode mismatch!"
print(f"  Barcode alignment: OK")

# KNN on PCA
print(f"\nComputing k={args.k} KNN on obsm[X_pca]...")
t0 = time.time()
nn = NearestNeighbors(n_neighbors=args.k, n_jobs=-1).fit(rna_sub.obsm["X_pca"])
_, knn_indices = nn.kneighbors(rna_sub.obsm["X_pca"])
rna_sub.uns["knn_indices"] = knn_indices
print(f"  Done in {time.time()-t0:.1f}s  |  shape: {knn_indices.shape}")

print(f"\nSaving RNA:  {RNA_OUT}")
rna_sub.write_h5ad(RNA_OUT)
print(f"Saving ATAC: {ATAC_OUT}")
atac_sub.write_h5ad(ATAC_OUT)

print("\n" + "=" * 60)
print(f"SUMMARY")
print(f"  Cell type: {args.cell_type}")
print(f"  Cells:     {rna_sub.n_obs:,}")
print(f"  Genes:     {rna_sub.n_vars:,} HVGs")
print(f"  Peaks:     {atac_sub.n_vars:,} HVPs")
print(f"  KNN k:     {args.k}")
print(f"  RNA out:   {RNA_OUT}")
print(f"  ATAC out:  {ATAC_OUT}")
print(f"End: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
