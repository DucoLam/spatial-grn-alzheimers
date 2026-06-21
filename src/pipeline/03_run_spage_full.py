import sys
import time
import os
import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, "SpaGE")
from SpaGE.main import SpaGE

# ── Paths ─────────────────────────────────────────────────────────────────────

SPATIAL_PATH = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
RNA_PATH     = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/SEAAD_MTG_RNAseq_HVG_SPATIAL_ONLY.h5ad"
HVG_PATH     = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/hvg_names_sub.csv"
OUTPUT_PATH  = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad_spatial/SpaGE_predicted_HVGs_FULL_50PV.h5ad"

# ── Config ────────────────────────────────────────────────────────────────────

N_PV               = 50
SPATIAL_SUBSAMPLE_N = 1887729
RANDOM_SEED        = 42

# ── Helpers ───────────────────────────────────────────────────────────────────

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

def to_df(adata):
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    return pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)

# ── Main ──────────────────────────────────────────────────────────────────────

print("=" * 60)
print("03_run_spage_full")
print("Start:", now())
print(f"Config: N_PV={N_PV}, SPATIAL_SUBSAMPLE_N={SPATIAL_SUBSAMPLE_N:,}, RANDOM_SEED={RANDOM_SEED}")
print("=" * 60)

print(f"\nLoading spatial data: {SPATIAL_PATH}")
adata_spatial = ad.read_h5ad(SPATIAL_PATH)
print(f"  Spatial shape:  {adata_spatial.n_obs:,} cells x {adata_spatial.n_vars:,} genes")
print(f"  Sections:       {adata_spatial.obs['Section'].nunique() if 'Section' in adata_spatial.obs else 'N/A'}")
print(f"  Donors:         {adata_spatial.obs['Donor ID'].nunique() if 'Donor ID' in adata_spatial.obs else 'N/A'}")

if SPATIAL_SUBSAMPLE_N is not None:
    print(f"\nSubsampling spatial to {SPATIAL_SUBSAMPLE_N:,} cells...")
    rng = np.random.default_rng(RANDOM_SEED)
    idx = rng.choice(adata_spatial.n_obs, size=min(SPATIAL_SUBSAMPLE_N, adata_spatial.n_obs), replace=False)
    adata_spatial = adata_spatial[idx].copy()
    print(f"  After subsample: {adata_spatial.n_obs:,} cells x {adata_spatial.n_vars:,} genes")

print(f"\nLoading RNA data: {RNA_PATH}")
adata_rna = ad.read_h5ad(RNA_PATH)
print(f"  RNA shape: {adata_rna.n_obs:,} cells x {adata_rna.n_vars:,} genes")

print("\nConverting to DataFrames...")
spatial_df = to_df(adata_spatial)
rna_df     = to_df(adata_rna)
print(f"  Spatial df: {spatial_df.shape[0]:,} x {spatial_df.shape[1]:,}")
print(f"  RNA df:     {rna_df.shape[0]:,} x {rna_df.shape[1]:,}")

shared_genes = spatial_df.columns.intersection(rna_df.columns)
print(f"\nShared genes (used as PVs): {len(shared_genes):,}")
print(f"  First 10: {list(shared_genes[:10])}")

hvg_names = pd.read_csv(HVG_PATH)["Gene Name"].tolist()
hvg_names = list(dict.fromkeys(hvg_names))
genes_to_predict = [g for g in hvg_names if g in rna_df.columns and g not in spatial_df.columns]
n_hvg_already    = len([g for g in hvg_names if g in spatial_df.columns])
n_hvg_missing    = len([g for g in hvg_names if g not in rna_df.columns])

print(f"\nHVG accounting:")
print(f"  Total HVGs in list:           {len(hvg_names):,}")
print(f"  Already measured in MERFISH:  {n_hvg_already:,}")
print(f"  To predict with SpaGE:        {len(genes_to_predict):,}")
print(f"  Not in RNA (skipped):         {n_hvg_missing:,}")
print(f"  First 10 to predict: {genes_to_predict[:10]}")

assert len(shared_genes) == 140, f"Expected 140 shared genes, got {len(shared_genes)}"
assert len(hvg_names) == 500, f"Expected 500 HVGs, got {len(hvg_names)}"
assert n_hvg_already + len(genes_to_predict) == len(hvg_names), "HVG accounting mismatch"
assert N_PV <= len(shared_genes), f"N_PV={N_PV} > shared genes={len(shared_genes)}"

os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

print(f"\nStarting SpaGE: {now()}")
start = time.time()

predicted = SpaGE(
    Spatial_data=spatial_df,
    RNA_data=rna_df,
    n_pv=N_PV,
    genes_to_predict=genes_to_predict
)

elapsed = time.time() - start
print(f"SpaGE done: {now()}")
print(f"  Elapsed: {elapsed / 3600:.2f} h ({elapsed / 60:.1f} min)")
print(f"  Predicted shape: {predicted.shape[0]:,} cells x {predicted.shape[1]:,} genes")

assert predicted.shape[0] == adata_spatial.n_obs
assert predicted.shape[1] == len(genes_to_predict)

predicted.index = adata_spatial.obs_names
adata_pred = ad.AnnData(
    X=predicted.values,
    obs=adata_spatial.obs.copy(),
    var=pd.DataFrame(index=predicted.columns)
)
for key in adata_spatial.obsm.keys():
    adata_pred.obsm[key] = adata_spatial.obsm[key].copy()

adata_pred.uns["SpaGE"] = {
    "spatial_input": SPATIAL_PATH,
    "rna_input": RNA_PATH,
    "hvg_input": HVG_PATH,
    "n_pv": N_PV,
    "n_shared_genes": int(len(shared_genes)),
    "n_genes_predicted": int(len(genes_to_predict)),
    "spatial_subsample_n": SPATIAL_SUBSAMPLE_N,
    "random_seed": RANDOM_SEED,
    "elapsed_hours": float(elapsed / 3600),
    "created_at": now()
}
adata_pred.uns["SpaGE_shared_genes"]      = list(shared_genes)
adata_pred.uns["SpaGE_genes_to_predict"]  = list(genes_to_predict)

print(f"\nSaving to: {OUTPUT_PATH}")
adata_pred.write_h5ad(OUTPUT_PATH)

print("\n" + "=" * 60)
print("SUMMARY")
print(f"  Spatial input cells:    {adata_spatial.n_obs:,}")
print(f"  RNA input cells:        {adata_rna.n_obs:,}")
print(f"  RNA input genes:        {adata_rna.n_vars:,}")
print(f"  Shared genes (PVs):     {len(shared_genes):,}")
print(f"  N_PV used:              {N_PV}")
print(f"  HVGs predicted:         {len(genes_to_predict):,}")
print(f"  Output shape:           {adata_pred.n_obs:,} x {adata_pred.n_vars:,}")
print(f"  Elapsed:                {elapsed / 3600:.2f} h")
print(f"  Output path:            {OUTPUT_PATH}")
print(f"End: {now()}")
print("=" * 60)
