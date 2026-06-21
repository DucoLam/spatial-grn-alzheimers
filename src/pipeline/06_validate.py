"""
06_validate.py  -- parameterized per-subtype Tangram holdout validation

Mirrors per-subtype matched production (same-subtype reference, donor-matched).
For each donor: hold out HOLDOUT_FRAC of this subtype's MERFISH cells, map them
back against the remaining SAME-subtype cells, measure Euclidean error (um).
Sweeps argmax + several top-k weighted-centroid thresholds (validated method = top 0.5%).

Output: dflam/data/validation/<prefix>_validation.csv  (one row per method)

Usage:
    python 06_validate.py --cell_type "L2/3 IT" --prefix l23_it
"""
import sys, time, argparse
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from pathlib import Path

sys.path.insert(0, "/tudelft.net/staff-umbrella/ScReNI/dflam/Tangram")
import tangram as tg

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

ap = argparse.ArgumentParser()
ap.add_argument("--cell_type", required=True)
ap.add_argument("--prefix",    required=True)
ap.add_argument("--holdout",   type=float, default=0.05)
ap.add_argument("--num_epochs", type=int, default=1000)
ap.add_argument("--seed",      type=int, default=42)
args = ap.parse_args()

# top-k sweep — 0.5% (0.005) is the validated method; neighbours for context
TOPK_FRACS = [0.001, 0.005, 0.01, 0.05, 0.1]

BSC          = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni"
MERFISH_PATH = f"{BSC}/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
OUT_DIR      = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/data/validation")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV      = OUT_DIR / f"{args.prefix}_validation_persection.csv"

print("=" * 64)
print("06_validate (per-subtype, same-subtype matched reference)")
print(f"  cell_type: {args.cell_type}   prefix: {args.prefix}")
print(f"  holdout  : {args.holdout*100:g}%   epochs: {args.num_epochs}")
print(f"  top-k    : {[f'{f*100:g}%' for f in TOPK_FRACS]}")
print(f"  Start    : {now()}")
print("=" * 64)

print(f"\nLoading MERFISH: {MERFISH_PATH}")
adata = ad.read_h5ad(MERFISH_PATH)
use_genes = list(adata.var_names[~adata.var_names.str.startswith("Blank")])
adata = adata[:, use_genes].copy()
if "Subclass" not in adata.obs.columns:
    sys.exit("ERROR: MERFISH has no 'Subclass' column.")
adata = adata[adata.obs["Subclass"].astype(str) == args.cell_type].copy()
print(f"  {args.cell_type} cells: {adata.n_obs:,}")
if adata.n_obs < 40:
    sys.exit(f"ERROR: too few cells ({adata.n_obs}) for validation.")

coords_key = next((k for k in ("X_spatial_raw", "spatial", "X_spatial") if k in adata.obsm), None)
# SECTION-AWARE: each MERFISH Section is its own coordinate frame. We validate
# per (donor x section) so holdout/reference/centroid all live in ONE frame.
adata.obs["Section"] = adata.obs["Section"].astype(str)
sections = sorted(adata.obs["Section"].unique())
print(f"  Sections: {len(sections)}  (validating per-section, not per-donor)")

rng = np.random.default_rng(args.seed)
methods = ["argmax"] + [f"top_{f*100:g}%" for f in TOPK_FRACS]
errs_by_method = {m: [] for m in methods}
n_used = 0

for i, sec in enumerate(sections):
    sub = adata[adata.obs["Section"] == sec].copy()
    if sub.n_obs < 40:
        print(f"\nSection {i+1}/{len(sections)} {sec}: only {sub.n_obs} cells, skipping")
        continue
    n_hold = max(1, int(sub.n_obs * args.holdout))
    held = rng.choice(sub.n_obs, size=n_hold, replace=False)
    ref  = np.setdiff1d(np.arange(sub.n_obs), held)
    sp_held = sub[held].copy()
    sp_ref  = sub[ref].copy()
    print(f"\nSection {i+1}/{len(sections)} {sec}: {sub.n_obs:,} cells -> ref {len(ref):,} / held {len(held):,}")

    genes = [g for g in use_genes if g in sp_ref.var_names and g in sp_held.var_names]
    sc_tmp, sp_tmp = sp_held.copy(), sp_ref.copy()
    tg.pp_adatas(sc_tmp, sp_tmp, genes=genes)
    t0 = time.time()
    amap = tg.map_cells_to_space(
        sc_tmp, sp_tmp, mode="cells", density_prior="uniform",
        num_epochs=args.num_epochs, device="cpu", random_state=args.seed, verbose=False,
    )
    print(f"  mapped in {(time.time()-t0)/60:.1f} min")

    M = amap.X.toarray() if sp.issparse(amap.X) else amap.X
    true_c = sp_held.obsm[coords_key]
    ref_c  = sp_tmp.obsm[coords_key]

    # argmax
    e = np.linalg.norm(ref_c[M.argmax(axis=1)] - true_c, axis=1)
    errs_by_method["argmax"].extend(e.tolist())

    for f in TOPK_FRACS:
        k = max(1, int(np.ceil(f * len(ref))))
        topk = np.argsort(M, axis=1)[:, -k:]
        Mk = np.zeros_like(M); rows = np.arange(M.shape[0])[:, None]
        Mk[rows, topk] = M[rows, topk]
        rs = Mk.sum(axis=1, keepdims=True); rs[rs == 0] = 1; Mk /= rs
        pred = np.stack([(Mk*ref_c[:,0]).sum(1), (Mk*ref_c[:,1]).sum(1)], axis=1)
        e = np.linalg.norm(pred - true_c, axis=1)
        errs_by_method[f"top_{f*100:g}%"].extend(e.tolist())
    n_used += 1

# ---- summary + write -------------------------------------------------------
print(f"\n{'='*64}")
print(f"RESULTS  -- {args.cell_type}  ({n_used} sections, "
      f"{len(errs_by_method['argmax']):,} held-out cells)")
print(f"{'method':<14}{'median_um':>12}{'mean_um':>12}{'within_500':>12}{'within_1000':>13}")
print("-" * 64)
rows = []
for m in methods:
    e = np.array(errs_by_method[m])
    med, mean = np.median(e), np.mean(e)
    w5, w10 = (e < 500).mean()*100, (e < 1000).mean()*100
    print(f"{m:<14}{med:>11.0f}u{mean:>11.0f}u{w5:>11.1f}%{w10:>12.1f}%")
    rows.append({"cell_type": args.cell_type, "method": m, "n_cells": len(e),
                 "median_um": med, "mean_um": mean,
                 "within_500_pct": w5, "within_1000_pct": w10})

pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f"\nWrote {OUT_CSV}")

# Save per-cell error distances (for cumulative-distribution / CDF plots)
ERR_NPZ = OUT_DIR / f"{args.prefix}_errors.npz"
np.savez(ERR_NPZ, **{m: np.asarray(errs_by_method[m]) for m in methods})
print(f"Wrote {ERR_NPZ}  (per-cell errors, {len(methods)} methods)")
print(f"End: {now()}")
print("=" * 64)
