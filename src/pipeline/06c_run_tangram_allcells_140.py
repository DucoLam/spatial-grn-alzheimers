"""
06c_run_tangram_allcells_140
Per-donor matched Tangram mapping using ALL spatial cell types as reference.
Queries: RNA-seq astrocytes (~186/donor)
Reference: ALL MERFISH cells for that donor (~70k/donor)
Genes: 140 native MERFISH genes, top 0.5% weighted centroid
"""

import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

sys.path.insert(0, "/tudelft.net/staff-umbrella/ScReNI/dflam/Tangram")
import tangram as tg

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

RNA_PATH     = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/astrocytes_seaad_paired_rna.h5ad"
MERFISH_PATH = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
MAPPING_OUT  = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/Tangram_rna_with_coords_allcells_140.h5ad"
CSV_OUT      = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/Tangram_cell_locations_allcells_140.csv"

TOP_K_FRAC   = 0.005
NUM_EPOCHS   = 1000
DEVICE       = "cpu"
RANDOM_SEED  = 42

print("=" * 60)
print("06c_run_tangram_allcells_140")
print(f"Mode:   per-donor matched, ALL cell types as spatial reference")
print(f"Genes:  140 native MERFISH genes")
print(f"Coords: top {TOP_K_FRAC*100:.1f}% weighted centroid")
print(f"Config: NUM_EPOCHS={NUM_EPOCHS}, DEVICE={DEVICE}, SEED={RANDOM_SEED}")
print("Start:", now())
print("=" * 60)

# Load MERFISH once to get gene list
print(f"\nLoading MERFISH: {MERFISH_PATH}")
adata_sp_full = ad.read_h5ad(MERFISH_PATH)
blank_mask = ~adata_sp_full.var_names.str.startswith("Blank")
use_genes  = list(adata_sp_full.var_names[blank_mask])
adata_sp_full = adata_sp_full[:, use_genes].copy()
print(f"  Total cells: {adata_sp_full.n_obs:,}  |  Genes: {len(use_genes)}")

col_sp_type = next((c for c in ("Supertype", "cell_type", "celltype") if c in adata_sp_full.obs.columns), None)
n_astro = adata_sp_full.obs[col_sp_type].str.lower().str.contains("astro", na=False).sum()
print(f"  Astrocytes: {n_astro:,}  ({n_astro/adata_sp_full.n_obs*100:.1f}%)")

coords_key = next((k for k in ("X_spatial_raw", "spatial", "X_spatial") if k in adata_sp_full.obsm), None)
print(f"  Coords key: '{coords_key}'")

# Load RNA astrocytes
print(f"\nLoading RNA (astrocytes): {RNA_PATH}")
adata_rna = ad.read_h5ad(RNA_PATH)
print(f"  Shape: {adata_rna.n_obs:,} x {adata_rna.n_vars:,}  |  Donors: {adata_rna.obs['Donor ID'].nunique()}")

# Overlapping donors
rna_donors = set(adata_rna.obs["Donor ID"].astype(str).unique())
sp_donors  = set(adata_sp_full.obs["Donor ID"].astype(str).unique())
shared     = sorted(rna_donors & sp_donors)
print(f"\nOverlapping donors: {len(shared)}  {shared}")

rna_matched = adata_rna[adata_rna.obs["Donor ID"].astype(str).isin(shared)].copy()
print(f"RNA cells in overlapping donors: {rna_matched.n_obs:,}")

# Per-donor mapping
results  = []
t_total  = time.time()

for i, donor in enumerate(shared):
    print(f"\n{'='*60}")
    print(f"Donor {i+1}/{len(shared)}: {donor}")

    sc_d = adata_rna[adata_rna.obs["Donor ID"].astype(str) == donor].copy()
    sp_d = adata_sp_full[adata_sp_full.obs["Donor ID"].astype(str) == donor].copy()

    n_rna = sc_d.n_obs
    n_sp  = sp_d.n_obs
    n_sp_astro = sp_d.obs[col_sp_type].str.lower().str.contains("astro", na=False).sum()
    mb = n_rna * n_sp * 4 / 1e6
    print(f"  RNA astrocytes:     {n_rna:,}")
    print(f"  Spatial total:      {n_sp:,}  (astro: {n_sp_astro:,}, other: {n_sp - n_sp_astro:,})")
    print(f"  Matrix size:        {n_rna:,} x {n_sp:,}  ({mb:.1f} MB)")

    if n_rna == 0 or n_sp == 0:
        print(f"  Skipping: no cells")
        continue

    sc_d_full = sc_d.copy()

    donor_genes = [g for g in use_genes if g in sc_d.var_names and g in sp_d.var_names]
    print(f"  Shared genes: {len(donor_genes)} / {len(use_genes)}")
    tg.pp_adatas(sc_d, sp_d, genes=donor_genes)

    t0 = time.time()
    adata_map = tg.map_cells_to_space(
        sc_d, sp_d,
        mode="cells",
        density_prior="uniform",
        num_epochs=NUM_EPOCHS,
        device=DEVICE,
        random_state=RANDOM_SEED,
        verbose=False,
    )
    elapsed = time.time() - t0
    print(f"  Mapping done in {elapsed/60:.1f} min")

    M = adata_map.X
    if sp.issparse(M):
        M = M.toarray()

    sp_coords  = sp_d.obsm[coords_key]
    best_idx   = M.argmax(axis=1)
    best_prob  = M.max(axis=1)
    best_coords = sp_coords[best_idx]

    # Top 0.5% weighted centroid
    k = max(1, int(np.ceil(TOP_K_FRAC * n_sp)))
    top_k_idx  = np.argsort(M, axis=1)[:, -k:]
    M_topk     = np.zeros_like(M)
    rows       = np.arange(M.shape[0])[:, None]
    M_topk[rows, top_k_idx] = M[rows, top_k_idx]
    row_sums   = M_topk.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    M_topk    /= row_sums
    weighted_x = (M_topk * sp_coords[:, 0]).sum(axis=1)
    weighted_y = (M_topk * sp_coords[:, 1]).sum(axis=1)
    print(f"  Top-k: k={k} ({TOP_K_FRAC*100:.1f}% of {n_sp:,} spatial cells)")

    # Check what cell type argmax landed on
    mapped_types = sp_d.obs[col_sp_type].iloc[best_idx].values
    astro_hits   = pd.Series(mapped_types).str.lower().str.contains("astro", na=False).mean()
    print(f"  Argmax landed on astrocyte: {astro_hits*100:.1f}%")

    sp_obs_best = sp_d.obs.iloc[best_idx].reset_index(drop=True)
    sp_obs_best.index = sc_d_full.obs_names

    row = pd.DataFrame({
        "tangram_x":        weighted_x,
        "tangram_y":        weighted_y,
        "tangram_x_argmax": best_coords[:, 0],
        "tangram_y_argmax": best_coords[:, 1],
        "tangram_prob_max": best_prob,
        "tangram_mapped_cell_type": mapped_types,
    }, index=sc_d_full.obs_names)

    for col in ("Section", "Donor ID", "Supertype", "Braak", "APOE4 Status"):
        if col in sp_obs_best.columns:
            row["tangram_mapped_" + col.lower().replace(" ", "_")] = sp_obs_best[col].values

    results.append(row)
    print(f"  Done: {donor}")

print(f"\nCombining {len(results)} donors...")
all_results = pd.concat(results)
print(f"  Mapped: {len(all_results):,} cells")

for col in all_results.columns:
    rna_matched.obs[col] = all_results.loc[rna_matched.obs_names, col]

print(f"Saving h5ad: {MAPPING_OUT}")
rna_matched.write_h5ad(MAPPING_OUT)

print(f"Saving CSV:  {CSV_OUT}")
csv_df = all_results[["tangram_x", "tangram_y", "tangram_prob_max", "tangram_mapped_cell_type"]].copy()
csv_df.index.name = "cell_id"
csv_df.to_csv(CSV_OUT)

total_elapsed = time.time() - t_total
print("\n" + "=" * 60)
print("SUMMARY")
print(f"  Mode:   per-donor matched, ALL spatial cell types")
print(f"  Genes:  {len(use_genes)} native MERFISH")
print(f"  Donors: {len(results)}/{len(shared)}")
print(f"  Cells:  {len(all_results):,}")
print(f"  Time:   {total_elapsed/60:.1f} min")
print(f"  h5ad:   {MAPPING_OUT}")
print(f"  csv:    {CSV_OUT}")
print(f"End: {now()}")
print("=" * 60)
