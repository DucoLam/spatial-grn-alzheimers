"""
06_validate_spage_panel.py -- per-subtype, per-section Tangram holdout validation
for a SpaGE-imputed gene panel (100 or 500 genes). Mirrors 06_validate.py exactly
(section-aware, same-subtype reference, 5% holdout, argmax + top-k sweep), but the
spatial features come from the SpaGE-imputed h5ad instead of the 140 native genes.

Panel construction (matches the original SpaGE sweep):
  measured  = HVGs also present natively in MERFISH (34)
  predicted = HVGs NOT in MERFISH (SpaGE-imputed, 466), ordered by HVG rank
  100-gene  = measured + predicted[:66]
  500-gene  = measured + predicted   (all)

Output: dflam/data/validation/panel<P>_<prefix>_persection.csv
Usage:  python 06_validate_spage_panel.py --panel 500 --cell_type "L2/3 IT" --prefix l23_it
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
ap.add_argument("--panel", required=True, choices=["100", "500"])
ap.add_argument("--cell_type", required=True)
ap.add_argument("--prefix", required=True)
ap.add_argument("--holdout", type=float, default=0.05)
ap.add_argument("--num_epochs", type=int, default=1000)
ap.add_argument("--seed", type=int, default=42)
args = ap.parse_args()
TOPK_FRACS = [0.001, 0.005, 0.01, 0.05, 0.1]

BSC   = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni"
SPAGE = f"{BSC}/data/processed/seaad_spatial/SpaGE_500_HVGs_FULL.h5ad"
MERF  = f"{BSC}/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
HVG   = f"{BSC}/data/processed/seaad/hvg_names_sub.csv"
OUT   = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/data/validation")
OUT.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT / f"panel{args.panel}_{args.prefix}_persection.csv"

print("=" * 64)
print(f"06_validate_spage_panel  panel={args.panel}  cell_type={args.cell_type}")
print(f"  holdout {args.holdout*100:g}%  epochs {args.num_epochs}  seed {args.seed}")
print(f"  Start {now()}")
print("=" * 64)

# ---- build panel gene list ------------------------------------------------
mf = ad.read_h5ad(MERF, backed="r")
merf_genes = set(g for g in mf.var_names if not g.startswith("Blank"))
del mf
hvg = pd.read_csv(HVG)["Gene Name"].tolist()
measured  = [g for g in hvg if g in merf_genes]
predicted = [g for g in hvg if g not in merf_genes]
panel_genes = measured + (predicted[:66] if args.panel == "100" else predicted)
print(f"panel {args.panel}: measured {len(measured)} + imputed "
      f"{len(panel_genes)-len(measured)} = {len(panel_genes)}")

# ---- load SpaGE spatial (backed -> subtype subset -> memory) --------------
sp_b = ad.read_h5ad(SPAGE, backed="r")
mask = (sp_b.obs["Subclass"].astype(str) == args.cell_type).values
adata = sp_b[mask].to_memory()
del sp_b
panel_genes = [g for g in panel_genes if g in adata.var_names]
adata = adata[:, panel_genes].copy()
print(f"  {args.cell_type}: {adata.n_obs:,} cells x {len(panel_genes)} genes")
if adata.n_obs < 40:
    sys.exit(f"ERROR: too few cells ({adata.n_obs}).")

coords_key = next((k for k in ("X_spatial_raw", "spatial", "X_spatial") if k in adata.obsm), None)
adata.obs["Section"] = adata.obs["Section"].astype(str)
sections = sorted(adata.obs["Section"].unique())
print(f"  Sections: {len(sections)} (per-section holdout)")

rng = np.random.default_rng(args.seed)
methods = ["argmax"] + [f"top_{f*100:g}%" for f in TOPK_FRACS]
errs = {m: [] for m in methods}

for i, sec in enumerate(sections):
    sub = adata[adata.obs["Section"] == sec].copy()
    if sub.n_obs < 40:
        continue
    n_hold = max(1, int(sub.n_obs * args.holdout))
    held = rng.choice(sub.n_obs, size=n_hold, replace=False)
    ref  = np.setdiff1d(np.arange(sub.n_obs), held)
    sp_held, sp_ref = sub[held].copy(), sub[ref].copy()
    genes = [g for g in panel_genes if g in sp_ref.var_names and g in sp_held.var_names]
    sc_tmp, sp_tmp = sp_held.copy(), sp_ref.copy()
    tg.pp_adatas(sc_tmp, sp_tmp, genes=genes)
    t0 = time.time()
    amap = tg.map_cells_to_space(
        sc_tmp, sp_tmp, mode="cells", density_prior="uniform",
        num_epochs=args.num_epochs, device="cpu", random_state=args.seed, verbose=False)
    M = amap.X.toarray() if sp.issparse(amap.X) else np.asarray(amap.X)
    true_c = sp_held.obsm[coords_key]
    ref_c  = sp_ref.obsm[coords_key]
    bi = M.argmax(axis=1)
    errs["argmax"].extend(np.sqrt(((ref_c[bi] - true_c) ** 2).sum(1)).tolist())
    for f in TOPK_FRACS:
        k = max(1, int(np.ceil(f * M.shape[1])))
        topk = np.argsort(M, axis=1)[:, -k:]
        Mk = np.zeros_like(M); r = np.arange(M.shape[0])[:, None]; Mk[r, topk] = M[r, topk]
        rs = Mk.sum(1, keepdims=True); rs[rs == 0] = 1; Mk /= rs
        pw = np.column_stack([(Mk * ref_c[:, 0]).sum(1), (Mk * ref_c[:, 1]).sum(1)])
        errs[f"top_{f*100:g}%"].extend(np.sqrt(((pw - true_c) ** 2).sum(1)).tolist())
    print(f"  [{i+1}/{len(sections)}] {sec}: ref {len(ref)} / held {len(held)} "
          f"({(time.time()-t0)/60:.1f} min)")

rows = []
for m in methods:
    e = np.array(errs[m])
    rows.append({"cell_type": args.cell_type, "panel": args.panel, "method": m,
                 "n_cells": len(e), "median_um": float(np.median(e)),
                 "mean_um": float(e.mean()),
                 "within_500_pct": float((e < 500).mean() * 100),
                 "within_1000_pct": float((e < 1000).mean() * 100)})
df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)
print(df.to_string(index=False))
print(f"saved {OUT_CSV}  |  done {now()}")
