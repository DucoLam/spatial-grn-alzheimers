"""
06_run_tangram.py  -- parameterized per-subtype matched Tangram mapping

Replicates the validated 06b baseline (140 native MERFISH genes, top-k weighted
centroid, donor-matched) but for ANY cell subtype via --cell_type.

Query     : RNA cells of the given Subclass (full gene set, intersected w/ 140 MERFISH genes)
Reference : SAME-subtype MERFISH cells (donor-matched)
Output    : dflam/data/tangram/Tangram_rna_with_coords_<prefix>.h5ad
            dflam/data/tangram/Tangram_cell_locations_<prefix>.csv

Usage:
    python 06_run_tangram.py --cell_type "L2/3 IT" --prefix l23_it
    python 06_run_tangram.py --cell_type "Oligodendrocyte" --prefix oligo --top_k_frac 0.005
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
ap.add_argument("--cell_type",  required=True, help="Subclass label, e.g. 'L2/3 IT'")
ap.add_argument("--prefix",     required=True, help="Output prefix, e.g. l23_it")
ap.add_argument("--top_k_frac", type=float, default=0.005, help="Top-k fraction for weighted centroid (default 0.005 = 0.5%)")
ap.add_argument("--num_epochs", type=int, default=1000)
ap.add_argument("--seed",       type=int, default=42)
args = ap.parse_args()

BSC          = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni"
RNA_PATH     = f"{BSC}/data/processed/seaad/seaad_paired_rna.h5ad"
MERFISH_PATH = f"{BSC}/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
OUT_DIR      = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/data/tangram")
OUT_DIR.mkdir(parents=True, exist_ok=True)
MAPPING_OUT  = OUT_DIR / f"Tangram_rna_with_coords_{args.prefix}.h5ad"
CSV_OUT      = OUT_DIR / f"Tangram_cell_locations_{args.prefix}.csv"

DEVICE = "cpu"

print("=" * 64)
print("06_run_tangram (per-subtype matched)")
print(f"  cell_type : {args.cell_type}")
print(f"  prefix    : {args.prefix}")
print(f"  top_k_frac: {args.top_k_frac}  ({args.top_k_frac*100:g}%)")
print(f"  epochs    : {args.num_epochs}   seed: {args.seed}")
print(f"  Start     : {now()}")
print("=" * 64)

# ---- MERFISH: non-blank genes + filter to subtype --------------------------
print(f"\nLoading MERFISH: {MERFISH_PATH}")
adata_sp_full = ad.read_h5ad(MERFISH_PATH)
use_genes = list(adata_sp_full.var_names[~adata_sp_full.var_names.str.startswith("Blank")])
adata_sp_full = adata_sp_full[:, use_genes].copy()
print(f"  Cells: {adata_sp_full.n_obs:,}  |  non-blank genes: {len(use_genes)}")

if "Subclass" not in adata_sp_full.obs.columns:
    sys.exit("ERROR: MERFISH has no 'Subclass' column.")
adata_sp = adata_sp_full[adata_sp_full.obs["Subclass"].astype(str) == args.cell_type].copy()
print(f"  {args.cell_type} spatial cells: {adata_sp.n_obs:,}")
if adata_sp.n_obs == 0:
    sys.exit(f"ERROR: no MERFISH cells with Subclass == '{args.cell_type}'.")
del adata_sp_full

coords_key = next((k for k in ("X_spatial_raw", "spatial", "X_spatial") if k in adata_sp.obsm), None)
if coords_key is None:
    sys.exit(f"ERROR: no spatial coords in obsm: {list(adata_sp.obsm.keys())}")
print(f"  Coords key: '{coords_key}'")

# ---- RNA: backed read, filter to subtype + overlapping donors --------------
print(f"\nLoading RNA (backed): {RNA_PATH}")
rna_b = ad.read_h5ad(RNA_PATH, backed="r")
rna_donors_all = set(rna_b.obs.loc[rna_b.obs["Subclass"].astype(str) == args.cell_type, "Donor ID"].astype(str))
sp_donors      = set(adata_sp.obs["Donor ID"].astype(str).unique())
shared = sorted(rna_donors_all & sp_donors)
print(f"  RNA donors w/ {args.cell_type}: {len(rna_donors_all)}  |  spatial donors: {len(sp_donors)}")
print(f"  Overlapping donors: {len(shared)}  {shared}")
if not shared:
    sys.exit("ERROR: no overlapping donors for this subtype.")

mask = (rna_b.obs["Subclass"].astype(str) == args.cell_type) & \
       (rna_b.obs["Donor ID"].astype(str).isin(shared))
adata_rna = rna_b[mask.values].to_memory()
del rna_b
print(f"  RNA {args.cell_type} cells (overlapping donors): {adata_rna.n_obs:,}")

# ---- per-donor Tangram mapping ---------------------------------------------
results = []
t_total = time.time()

for i, donor in enumerate(shared):
    print(f"\n{'-'*64}\nDonor {i+1}/{len(shared)}: {donor}")
    sc_d = adata_rna[adata_rna.obs["Donor ID"].astype(str) == donor].copy()
    sp_d = adata_sp [adata_sp.obs["Donor ID"].astype(str)  == donor].copy()
    n_rna, n_sp = sc_d.n_obs, sp_d.n_obs
    print(f"  RNA: {n_rna:,}  Spatial: {n_sp:,}  Matrix: {n_rna*n_sp*4/1e6:.0f} MB")
    if n_rna == 0 or n_sp == 0:
        print("  Skipping (empty modality)")
        continue

    sc_full = sc_d.copy()
    donor_genes = [g for g in use_genes if g in sc_d.var_names and g in sp_d.var_names]
    tg.pp_adatas(sc_d, sp_d, genes=donor_genes)

    t0 = time.time()
    amap = tg.map_cells_to_space(
        sc_d, sp_d, mode="cells", density_prior="uniform",
        num_epochs=args.num_epochs, device=DEVICE,
        random_state=args.seed, verbose=False,
    )
    print(f"  Mapped {len(donor_genes)} genes in {(time.time()-t0)/60:.1f} min")

    M = amap.X.toarray() if sp.issparse(amap.X) else amap.X
    coords  = sp_d.obsm[coords_key]
    sec_arr = sp_d.obs["Section"].astype(str).values   # section (= coord frame) per spatial cell

    # argmax (reference) + the section each query cell maps into
    bi = M.argmax(axis=1)
    bc = coords[bi]
    bp = M.max(axis=1)
    best_sec = sec_arr[bi]

    # SECTION-AWARE top-k weighted centroid: each Section is an independent
    # coordinate frame, so confine the top-k to the argmax cell's section.
    # Never average (x,y) across sections.
    nq = M.shape[0]
    wx = np.empty(nq); wy = np.empty(nq)
    for s in np.unique(best_sec):
        qmask = best_sec == s          # query cells assigned to section s
        cmask = sec_arr == s           # spatial cols belonging to section s
        Msub  = M[np.ix_(qmask, cmask)]
        csub  = coords[cmask]
        k = max(1, int(np.ceil(args.top_k_frac * cmask.sum())))
        topk = np.argsort(Msub, axis=1)[:, -k:]
        Mk = np.zeros_like(Msub); rows = np.arange(Msub.shape[0])[:, None]
        Mk[rows, topk] = Msub[rows, topk]
        rs = Mk.sum(axis=1, keepdims=True); rs[rs == 0] = 1; Mk /= rs
        wx[qmask] = (Mk * csub[:, 0]).sum(axis=1)
        wy[qmask] = (Mk * csub[:, 1]).sum(axis=1)
    print(f"  sections={len(np.unique(sec_arr))}  top_k={args.top_k_frac*100:g}%  "
          f"prob_max[{bp.min():.4f},{bp.max():.4f}]")

    row = pd.DataFrame({
        "tangram_x": wx, "tangram_y": wy,
        "tangram_x_argmax": bc[:, 0], "tangram_y_argmax": bc[:, 1],
        "tangram_prob_max": bp,
    }, index=sc_full.obs_names)
    sp_best = sp_d.obs.iloc[bi].reset_index(drop=True)
    for c in ("Section", "Donor ID", "Subclass", "Supertype", "Braak", "APOE4 Status"):
        if c in sp_best.columns:
            row["tangram_mapped_" + c.lower().replace(" ", "_")] = sp_best[c].values
    results.append(row)

# ---- combine + save --------------------------------------------------------
all_res = pd.concat(results)
print(f"\nMapped {len(all_res):,} / {adata_rna.n_obs:,} cells across {len(results)} donors")
for c in all_res.columns:
    adata_rna.obs[c] = all_res.loc[adata_rna.obs_names, c]

print(f"Saving h5ad: {MAPPING_OUT}")
adata_rna.write_h5ad(MAPPING_OUT)
print(f"Saving CSV : {CSV_OUT}")
# Validation shows argmax > section-confined centroid on both median error and
# within-500um, so the primary tangram_x/tangram_y = argmax; centroid kept too.
csv = pd.DataFrame({
    "cell_subtype":          args.cell_type,
    "tangram_x":             all_res["tangram_x_argmax"],
    "tangram_y":             all_res["tangram_y_argmax"],
    "tangram_x_centroid":    all_res["tangram_x"],
    "tangram_y_centroid":    all_res["tangram_y"],
    "tangram_prob_max":      all_res["tangram_prob_max"],
}, index=all_res.index)
if "tangram_mapped_section" in all_res.columns:
    csv["tangram_mapped_section"] = all_res["tangram_mapped_section"]
csv.index.name = "cell_id"
csv.to_csv(CSV_OUT)

print("\n" + "=" * 64)
print("SUMMARY")
print(f"  cell_type : {args.cell_type}")
print(f"  cells     : {len(all_res):,}")
print(f"  donors    : {len(results)}/{len(shared)}")
print(f"  mean prob_max : {all_res['tangram_prob_max'].mean():.5f}")
print(f"  time      : {(time.time()-t_total)/60:.1f} min")
print(f"  h5ad      : {MAPPING_OUT}")
print(f"  csv       : {CSV_OUT}")
print(f"  End       : {now()}")
print("=" * 64)
