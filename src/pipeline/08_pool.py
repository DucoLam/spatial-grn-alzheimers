"""
08_pool.py  -- POOLED multi-subtype subsample (ScReNI paper-faithful)

Pools ALL cells of the given cell subtypes into one dataset (no count cap),
then computes KNN ACROSS all pooled cells on the joint WNN embedding
(obsm['X_pca']). This matches the original ScReNI design: cells of all types
are clustered together and neighbours are found over the whole dataset.

Default dataset = our 4 subtypes: Astrocyte, L2/3 IT, L4 IT, Oligodendrocyte.

Outputs:
  dflam/data/seaad_paired_rna_<prefix>.h5ad   (pooled cells x 500 HVG, knn in .uns)
  dflam/data/seaad_paired_atac_<prefix>.h5ad  (pooled cells x 10k HVP)
"""
import sys, argparse, time
import anndata as ad
import numpy as np
from sklearn.neighbors import NearestNeighbors
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--cell_types", nargs="+",
                default=["Astrocyte", "L2/3 IT", "L4 IT", "Oligodendrocyte"],
                help="Subclasses to pool (space-separated)")
ap.add_argument("--prefix", default="pool4")
ap.add_argument("--k", type=int, default=20, help="KNN k across the pooled set")
args = ap.parse_args()

BSC      = Path("/tudelft.net/staff-umbrella/ScReNI/bsc-screni")
DFLAM    = Path("/tudelft.net/staff-umbrella/ScReNI/dflam")
RNA_IN   = BSC  / "data/processed/seaad/seaad_paired_rna_hvg.h5ad"
ATAC_IN  = BSC  / "data/processed/seaad/seaad_paired_atac_hvp.h5ad"
RNA_OUT  = DFLAM / f"data/seaad_paired_rna_{args.prefix}.h5ad"
ATAC_OUT = DFLAM / f"data/seaad_paired_atac_{args.prefix}.h5ad"
(DFLAM / "data").mkdir(parents=True, exist_ok=True)

print("=" * 64)
print(f"08_pool  (POOLED, paper-faithful)")
print(f"  cell_types : {args.cell_types}")
print(f"  prefix     : {args.prefix}   k={args.k}")
print(f"  Start      : {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 64)

print(f"\nLoading RNA  : {RNA_IN}")
rna = ad.read_h5ad(RNA_IN)
print(f"  full: {rna.shape}")
print(f"Loading ATAC : {ATAC_IN}")
atac = ad.read_h5ad(ATAC_IN)
print(f"  full: {atac.shape}")

col = "Subclass" if "Subclass" in rna.obs.columns else "cell_type"
mask = rna.obs[col].isin(args.cell_types)
print(f"\nPooling on obs['{col}']:  {int(mask.sum()):,} / {rna.n_obs:,} cells")
print(rna.obs.loc[mask, col].value_counts().to_string())
if mask.sum() == 0:
    sys.exit("ERROR: no cells matched the requested subtypes.")

rna_sub  = rna[mask.values].copy()
atac_sub = atac[mask.values].copy()
assert list(rna_sub.obs_names) == list(atac_sub.obs_names), "Barcode mismatch!"
print(f"  Barcode alignment OK  ({rna_sub.n_obs:,} pooled cells)")

if "X_pca" not in rna_sub.obsm:
    sys.exit("ERROR: obsm['X_pca'] (joint WNN embedding) missing from HVG file.")

print(f"\nComputing k={args.k} KNN ACROSS all {rna_sub.n_obs:,} pooled cells on obsm['X_pca'] ...")
t0 = time.time()
nn = NearestNeighbors(n_neighbors=args.k, n_jobs=-1).fit(rna_sub.obsm["X_pca"])
_, knn_indices = nn.kneighbors(rna_sub.obsm["X_pca"])
rna_sub.uns["knn_indices"] = knn_indices
print(f"  done in {time.time()-t0:.1f}s  |  shape {knn_indices.shape}")

# sanity: how within-subtype are the neighbourhoods in practice?
sub_arr = rna_sub.obs[col].astype(str).values
same = np.mean([np.mean(sub_arr[knn_indices[i]] == sub_arr[i]) for i in range(rna_sub.n_obs)])
print(f"  mean fraction of neighbours sharing the cell's subtype: {same*100:.1f}%")

print(f"\nSaving RNA  : {RNA_OUT}")
rna_sub.write_h5ad(RNA_OUT)
print(f"Saving ATAC : {ATAC_OUT}")
atac_sub.write_h5ad(ATAC_OUT)

print("\n" + "=" * 64)
print("SUMMARY")
print(f"  Pooled cells : {rna_sub.n_obs:,}")
print(f"  Genes (HVG)  : {rna_sub.n_vars:,}")
print(f"  Peaks (HVP)  : {atac_sub.n_vars:,}")
print(f"  KNN k        : {args.k}  (across pooled set)")
print(f"  End          : {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 64)
