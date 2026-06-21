"""
06_validate_crossdonor_panel.py -- per-subtype CROSS-DONOR negative control.

Mirrors the donor-matched production config (native 140 genes, argmax) but maps
each donor's cells against a DIFFERENT donor's cells (leave-one-donor-out),
instead of the same donor/section. This is the negative control showing why
donor matching is needed. Reference is SUBSAMPLED (cap) for tractability.

For each donor d (query): map a sample of d's cells against a capped sample of
the same subtype's cells from ALL OTHER donors; error = Euclidean distance from
the argmax reference cell's coords to the query cell's true coords. Aggregate
over donors.

Output: dflam/data/validation/crossdonor_<prefix>.csv
Usage:  python 06_validate_crossdonor_panel.py --cell_type "L2/3 IT" --prefix l23_it
"""
import sys, time, argparse
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from pathlib import Path

sys.path.insert(0, "/tudelft.net/staff-umbrella/ScReNI/dflam/Tangram")
import tangram as tg

ap = argparse.ArgumentParser()
ap.add_argument("--cell_type", required=True)
ap.add_argument("--prefix", required=True)
ap.add_argument("--ref_cap", type=int, default=4000)
ap.add_argument("--query_cap", type=int, default=1000)
ap.add_argument("--num_epochs", type=int, default=1000)
ap.add_argument("--seed", type=int, default=42)
args = ap.parse_args()

BSC  = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni"
MERF = f"{BSC}/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
OUT  = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/data/validation")
OUT.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT / f"crossdonor_{args.prefix}.csv"

print(f"cross-donor (subsampled, ref_cap={args.ref_cap}) {args.cell_type}  start {time.strftime('%H:%M:%S')}")
adata = ad.read_h5ad(MERF, backed="r")
genes = [g for g in adata.var_names if not g.startswith("Blank")]
mask = (adata.obs["Subclass"].astype(str) == args.cell_type).values
sub = adata[mask].to_memory(); del adata
sub = sub[:, genes].copy()
coords_key = next((k for k in ("X_spatial_raw", "spatial", "X_spatial") if k in sub.obsm), None)
sub.obs["Donor ID"] = sub.obs["Donor ID"].astype(str)
donors = sorted(sub.obs["Donor ID"].unique())
print(f"  {args.cell_type}: {sub.n_obs:,} cells, {len(donors)} donors")

rng = np.random.default_rng(args.seed)
errs = []
for d in donors:
    q_idx = np.flatnonzero((sub.obs["Donor ID"] == d).values)
    r_idx = np.flatnonzero((sub.obs["Donor ID"] != d).values)
    if len(q_idx) < 40 or len(r_idx) < 100:
        continue
    if len(q_idx) > args.query_cap:
        q_idx = np.sort(rng.choice(q_idx, args.query_cap, replace=False))
    if len(r_idx) > args.ref_cap:
        r_idx = np.sort(rng.choice(r_idx, args.ref_cap, replace=False))
    sc, spx = sub[q_idx].copy(), sub[r_idx].copy()
    tg.pp_adatas(sc, spx, genes=genes)
    amap = tg.map_cells_to_space(sc, spx, mode="cells", density_prior="uniform",
                                 num_epochs=args.num_epochs, device="cpu",
                                 random_state=args.seed, verbose=False)
    M = amap.X.toarray() if sp.issparse(amap.X) else np.asarray(amap.X)
    true_c = sub[q_idx].obsm[coords_key]
    ref_c = sub[r_idx].obsm[coords_key]
    e = np.sqrt(((ref_c[M.argmax(1)] - true_c) ** 2).sum(1))
    errs.extend(e.tolist())
    print(f"  donor {d}: q{len(q_idx)} ref{len(r_idx)} median {np.median(e):.0f}um "
          f"within500 {(e<500).mean()*100:.1f}%")

e = np.array(errs)
row = {"cell_type": args.cell_type, "config": "cross-donor", "method": "argmax",
       "n_cells": int(len(e)), "median_um": float(np.median(e)), "mean_um": float(e.mean()),
       "within_500_pct": float((e < 500).mean() * 100),
       "within_1000_pct": float((e < 1000).mean() * 100)}
pd.DataFrame([row]).to_csv(OUT_CSV, index=False)
print(row)
print(f"saved {OUT_CSV}  done {time.strftime('%H:%M:%S')}")
