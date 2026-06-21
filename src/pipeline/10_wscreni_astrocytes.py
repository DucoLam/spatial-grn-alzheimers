"""
Step 3 – wScReNI inference for astrocytes
Infers a cell-specific GRN for each of the 6,777 astrocytes using
paired RNA + ATAC data and the gene-peak-TF triplets from step 2.
"""

import sys
import time
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/src")
from screni.inference.wscreni import infer_wscreni_networks
from screni.data.gene_peak_relations import GenePeakOverlapLabs

# ── paths ─────────────────────────────────────────────────────────────────────
DFLAM    = Path("/tudelft.net/staff-umbrella/ScReNI/dflam")
DATA_DIR = DFLAM / "data"
OUT_DIR  = DFLAM / "data" / "wscreni_networks_astrocytes"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PREFIX   = "astrocytes"
N_JOBS   = 8
SEED     = 42

# ── load RNA + ATAC ───────────────────────────────────────────────────────────
rna_path  = DATA_DIR / "seaad_paired_rna_astrocytes.h5ad"
atac_path = DATA_DIR / "seaad_paired_atac_astrocytes.h5ad"

for p in [rna_path, atac_path]:
    if not p.exists():
        sys.exit(f"ERROR: {p} not found. Run step 1 (08_subsample_astrocytes) first.")

print(f"Loading RNA  : {rna_path}")
rna  = ad.read_h5ad(str(rna_path))
print(f"Loading ATAC : {atac_path}")
atac = ad.read_h5ad(str(atac_path))

if "knn_indices" not in rna.uns:
    sys.exit("ERROR: rna.uns['knn_indices'] missing — re-run step 1.")
knn_indices = rna.uns["knn_indices"]

# ── load phase 3 outputs (gene-peak triplets) ─────────────────────────────────
triplets_path  = DATA_DIR / f"{PREFIX}_triplets.csv"
peak_mat_path  = DATA_DIR / f"{PREFIX}_peak_overlap_matrix.npz"
peak_info_path = DATA_DIR / f"{PREFIX}_peak_info.csv"

for p in [triplets_path, peak_mat_path, peak_info_path]:
    if not p.exists():
        sys.exit(f"ERROR: {p} not found. Run step 2 (09_gene_peak_astrocytes) first.")

triplets   = pd.read_csv(str(triplets_path))
peak_info  = pd.read_csv(str(peak_info_path), index_col=0)
peak_mat   = np.load(str(peak_mat_path))["peak_matrix"].astype(np.float64)
peak_names = list(peak_info.index)

# Gaussian noise for RF numerical stability (mirrors upstream pipeline)
rng = np.random.default_rng(seed=SEED)
peak_mat += rng.normal(0, 1e-5, peak_mat.shape)

# Build GenePeakOverlapLabs from triplets
labs = GenePeakOverlapLabs.from_dataframe(
    triplets.rename(columns={"target_gene": "gene.name", "peak": "peak.name"})
)

# ── summary ───────────────────────────────────────────────────────────────────
print("=" * 72)
print(f"wScReNI inference — {PREFIX}")
print("=" * 72)
print(f"  RNA  : {rna.shape[0]} cells x {rna.shape[1]} genes")
print(f"  ATAC : {atac.shape[0]} cells x {atac.shape[1]} peaks")
print(f"  Peak matrix (triplet peaks only): {peak_mat.shape}")
print(f"  Triplets: {len(triplets)}")
print(f"  KNN k: {knn_indices.shape[1]}")
print(f"  n_jobs: {N_JOBS}  |  n_trees: 100  |  seed: {SEED}")
print(f"  Network output: {OUT_DIR}")
print()

# ── infer wScReNI ─────────────────────────────────────────────────────────────
print("Inferring wScReNI networks (paired RNA+ATAC) ...")
t0 = time.time()

w_nets = infer_wscreni_networks(
    expr=rna,
    peak_mat=peak_mat,
    peak_names=peak_names,
    labs=labs,
    nearest_neighbors_idx=knn_indices,
    network_path=str(OUT_DIR),
    n_jobs=N_JOBS,
    n_trees=100,
    seed=SEED,
)

elapsed = time.time() - t0
print(f"  Done: {len(w_nets)} cell networks in {elapsed/60:.1f} min")
print()
print("=" * 72)
print("Step 3 complete.")
print(f"Networks written to: {OUT_DIR}")
print("Next: run step 4 (11_spatial_link.py) to join with Tangram coordinates.")
print("=" * 72)
