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
MAPPING_OUT  = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/Tangram_rna_with_coords_matched_140.h5ad"
CSV_OUT      = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/Tangram_cell_locations_140.csv"

TOP_K_FRAC   = 0.005  # top 0.5% — same validated method as 100-gene run
NUM_EPOCHS   = 1000
DEVICE       = "cpu"
RANDOM_SEED  = 42

print("=" * 60)
print("06b_run_tangram_matched_140")
print(f"Mode:   per-donor matched mapping")
print(f"Genes:  all 140 raw MERFISH genes (no HVG filter, no SpaGE)")
print(f"Coords: top {TOP_K_FRAC*100:.1f}% weighted centroid")
print(f"Config: NUM_EPOCHS={NUM_EPOCHS}, DEVICE={DEVICE}, SEED={RANDOM_SEED}")
print("Start:", now())
print("=" * 60)

# Load raw MERFISH to get the 140 measured genes
print(f"\nLoading MERFISH for gene list: {MERFISH_PATH}")
adata_merfish = ad.read_h5ad(MERFISH_PATH)
blank_mask    = ~adata_merfish.var_names.str.startswith("Blank")
use_genes     = list(adata_merfish.var_names[blank_mask])
print(f"  Measured genes (non-blank): {len(use_genes)}")
del adata_merfish

# Load RNA astrocytes
print(f"\nLoading RNA (astrocytes): {RNA_PATH}")
adata_rna = ad.read_h5ad(RNA_PATH)
print(f"  Shape:   {adata_rna.n_obs:,} cells x {adata_rna.n_vars:,} genes")
print(f"  Donors:  {adata_rna.obs['Donor ID'].nunique():,}")

# Load raw MERFISH spatial, filter to astrocytes
print(f"\nLoading raw MERFISH spatial: {MERFISH_PATH}")
adata_sp_full = ad.read_h5ad(MERFISH_PATH)
print(f"  Full shape: {adata_sp_full.n_obs:,} cells x {adata_sp_full.n_vars:,} genes")

# Remove blank genes
adata_sp_full = adata_sp_full[:, use_genes].copy()
print(f"  After removing blanks: {adata_sp_full.n_obs:,} cells x {adata_sp_full.n_vars:,} genes")

col_sp_type = next((c for c in ("Supertype", "cell_type", "celltype") if c in adata_sp_full.obs.columns), None)
astro_mask  = adata_sp_full.obs[col_sp_type].str.lower().str.contains("astro", na=False)
adata_sp    = adata_sp_full[astro_mask].copy()
print(f"  Astrocytes: {adata_sp.n_obs:,} cells  (col='{col_sp_type}')")

# Find overlapping donors
rna_donors = set(adata_rna.obs["Donor ID"].astype(str).unique())
sp_donors  = set(adata_sp.obs["Donor ID"].astype(str).unique())
shared     = sorted(rna_donors & sp_donors)
rna_only   = sorted(rna_donors - sp_donors)

print(f"\nDonor overlap:")
print(f"  RNA donors:        {len(rna_donors):,}")
print(f"  Spatial donors:    {len(sp_donors):,}")
print(f"  Overlapping:       {len(shared):,}  {shared}")
print(f"  RNA only (skip):   {len(rna_only):,}")

rna_matched = adata_rna[adata_rna.obs["Donor ID"].astype(str).isin(shared)].copy()
print(f"\nRNA cells in overlapping donors: {rna_matched.n_obs:,}")

coords_key = next((k for k in ("X_spatial_raw", "spatial", "X_spatial") if k in adata_sp.obsm), None)
if coords_key is None:
    raise KeyError(f"No spatial coords found. obsm keys: {list(adata_sp.obsm.keys())}")
print(f"Spatial coords from obsm['{coords_key}']")

# Per-donor Tangram mapping
results = []
t_total = time.time()

for i, donor in enumerate(shared):
    print(f"\n{'='*60}")
    print(f"Donor {i+1}/{len(shared)}: {donor}")

    sc_d = adata_rna[adata_rna.obs["Donor ID"].astype(str) == donor].copy()
    sp_d = adata_sp[adata_sp.obs["Donor ID"].astype(str) == donor].copy()

    n_rna = sc_d.n_obs
    n_sp  = sp_d.n_obs
    mb    = n_rna * n_sp * 4 / 1e6
    print(f"  RNA astrocytes:     {n_rna:,}")
    print(f"  Spatial astrocytes: {n_sp:,}")
    print(f"  Matrix size:        {n_rna:,} x {n_sp:,}  ({mb:.1f} MB)")

    if n_rna == 0 or n_sp == 0:
        print(f"  Skipping: no cells for one modality")
        continue

    sc_d_full = sc_d.copy()

    # Use 140 MERFISH genes intersected with what's in RNA
    donor_genes = [g for g in use_genes if g in sc_d.var_names and g in sp_d.var_names]
    print(f"  Genes available in both modalities: {len(donor_genes)} / {len(use_genes)}")
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

    sp_coords = sp_d.obsm[coords_key]

    # Top 0.5% weighted centroid
    k = max(1, int(np.ceil(TOP_K_FRAC * sp_d.n_obs)))
    top_k_idx  = np.argsort(M, axis=1)[:, -k:]
    M_topk     = np.zeros_like(M)
    rows       = np.arange(M.shape[0])[:, None]
    M_topk[rows, top_k_idx] = M[rows, top_k_idx]
    row_sums   = M_topk.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    M_topk    /= row_sums
    weighted_x = (M_topk * sp_coords[:, 0]).sum(axis=1)
    weighted_y = (M_topk * sp_coords[:, 1]).sum(axis=1)
    print(f"  Top-k: k={k} ({TOP_K_FRAC*100:.1f}% of {sp_d.n_obs} spatial cells)")

    # Argmax for reference
    best_idx    = M.argmax(axis=1)
    best_coords = sp_coords[best_idx]
    best_prob   = M.max(axis=1)

    print(f"  Coord range (top-k): x=[{weighted_x.min():.0f}, {weighted_x.max():.0f}]  "
          f"y=[{weighted_y.min():.0f}, {weighted_y.max():.0f}]")
    print(f"  Prob range:  [{best_prob.min():.5f}, {best_prob.max():.5f}]")

    sp_obs_best = sp_d.obs.iloc[best_idx].reset_index(drop=True)
    sp_obs_best.index = sc_d_full.obs_names

    row = pd.DataFrame({
        "tangram_x":        weighted_x,
        "tangram_y":        weighted_y,
        "tangram_x_argmax": best_coords[:, 0],
        "tangram_y_argmax": best_coords[:, 1],
        "tangram_prob_max": best_prob,
    }, index=sc_d_full.obs_names)

    for col in ("Section", "Donor ID", "Supertype", "Braak", "APOE4 Status"):
        if col in sp_obs_best.columns:
            col_key = "tangram_mapped_" + col.lower().replace(" ", "_")
            row[col_key] = sp_obs_best[col].values

    results.append(row)
    print(f"  Done: {donor}")

# Combine
print(f"\nCombining results from {len(results)} donors...")
all_results = pd.concat(results)
print(f"  Mapped RNA cells: {len(all_results):,} / {rna_matched.n_obs:,}")

for col in all_results.columns:
    rna_matched.obs[col] = all_results.loc[rna_matched.obs_names, col]

# Save h5ad
print(f"\nSaving h5ad: {MAPPING_OUT}")
rna_matched.write_h5ad(MAPPING_OUT)

# Save CSV
print(f"Saving CSV:  {CSV_OUT}")
csv_df = all_results[["tangram_x", "tangram_y", "tangram_prob_max"]].copy()
csv_df.index.name = "cell_id"
csv_df.to_csv(CSV_OUT)
print(f"  Rows: {len(csv_df):,}  |  Columns: cell_id, tangram_x, tangram_y, tangram_prob_max")

# Comparison stats vs 100-gene run
print(f"\nComparison metric (higher prob_max = more confident mapping):")
print(f"  Mean prob_max (140 genes): {all_results['tangram_prob_max'].mean():.5f}")
print(f"  Median prob_max:           {all_results['tangram_prob_max'].median():.5f}")

total_elapsed = time.time() - t_total
print("\n" + "=" * 60)
print("SUMMARY")
print(f"  Mode:              per-donor matched mapping")
print(f"  Genes:             {len(use_genes)} raw MERFISH (no HVG filter)")
print(f"  Coordinate method: top {TOP_K_FRAC*100:.1f}% weighted centroid")
print(f"  Donors mapped:     {len(results):,} / {len(shared):,}")
print(f"  RNA cells mapped:  {len(all_results):,}")
print(f"  Total time:        {total_elapsed/3600:.2f}h ({total_elapsed/60:.1f} min)")
print(f"  Outputs:")
print(f"    h5ad: {MAPPING_OUT}")
print(f"    csv:  {CSV_OUT}")
print(f"End: {now()}")
print("=" * 60)
