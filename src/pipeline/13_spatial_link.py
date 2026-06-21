"""
13_spatial_link.py -- STEP 4: join located cells to their GRNs + metadata,
build & cache the analysis-ready matrices for SQ2-SQ4.

Outputs (dflam/data/spatial_link/):
  spatial_grn_matrix.npz  -- sparse cells x edges weight matrix
  spatial_grn_index.npz   -- cell_ids, edge gene pairs (src/dst), flat edge idx
  spatial_meta.parquet    -- per-cell: subtype, x/y (argmax) + x_centroid/y_centroid (top-0.5%),
                             section, donor, Braak, progression, APOE, prob_max
  spatial_expr.npz        -- cells x 500 HVG expression (SQ2 expression-vs-space baseline)
"""
import os, time
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from pathlib import Path
from joblib import Parallel, delayed

DFLAM    = Path("/tudelft.net/staff-umbrella/ScReNI/dflam")
TG_DIR   = DFLAM / "data/tangram"
GRN_DIR  = DFLAM / "data/wscreni_networks_pool4/wScReNI"
POOL_RNA = DFLAM / "data/seaad_paired_rna_pool4.h5ad"
OUT      = DFLAM / "data/spatial_link"; OUT.mkdir(parents=True, exist_ok=True)
SUB = {"astrocyte": "Astrocyte", "l23_it": "L2/3 IT", "l4_it": "L4 IT", "oligo": "Oligodendrocyte"}
NJOBS = 8

def now(): return time.strftime("%H:%M:%S")
print("STEP 4 start", now())

# 1. located cells (join key = barcode) -----------------------------------
dfs = []
for p, l in SUB.items():
    f = TG_DIR / f"Tangram_cell_locations_{p}.csv"
    d = pd.read_csv(f); d["subtype"] = l; dfs.append(d)
loc = pd.concat(dfs, ignore_index=True).set_index("cell_id")
print(f"located cells: {len(loc):,}")

# 2. barcode -> GRN file path (one dir scan) ------------------------------
bc2path = {}
for fn in os.listdir(GRN_DIR):
    if fn.endswith(".network.txt"):
        bc = fn[:-len(".network.txt")].split(".", 1)[1]
        bc2path[bc] = GRN_DIR / fn
print(f"GRN files indexed: {len(bc2path):,}")

cell_ids = [c for c in loc.index if c in bc2path]
print(f"located cells WITH a GRN file: {len(cell_ids):,}  (missing: {len(loc)-len(cell_ids)})")
loc = loc.loc[cell_ids]

# 3. metadata from pool4 RNA obs (keep BOTH coordinate readouts) ----------
rna = ad.read_h5ad(POOL_RNA, backed="r")
obs = rna.obs
stage_cols = [c for c in ("Donor ID", "Braak", "Overall AD neuropathological Change",
                          "Continuous Pseudo-progression Score", "APOE4 Status") if c in obs.columns]
coord_cols = ["subtype", "tangram_x", "tangram_y", "tangram_mapped_section", "tangram_prob_max"]
for c in ("tangram_x_centroid", "tangram_y_centroid"):
    if c in loc.columns:
        coord_cols.append(c)
meta = loc[coord_cols].copy().join(obs.loc[cell_ids, stage_cols])
meta.rename(columns={"tangram_x": "x", "tangram_y": "y",
                     "tangram_x_centroid": "x_centroid", "tangram_y_centroid": "y_centroid",
                     "tangram_mapped_section": "section"}, inplace=True)
print("metadata cols:", list(meta.columns))

# 4. read each GRN -> nonzero (flattened idx, weight), in parallel --------
def read_nz(cid):
    M = pd.read_csv(bc2path[cid], sep="\t", index_col=0)
    v = M.values.ravel()
    nz = np.flatnonzero(v)
    return cid, nz.astype(np.int32), v[nz].astype(np.float32), M.shape[0], list(M.index)

print("reading GRNs ...", now())
t0 = time.time()
results = Parallel(n_jobs=NJOBS, prefer="threads")(delayed(read_nz)(c) for c in cell_ids)
print(f"  read {len(results)} GRNs in {(time.time()-t0)/60:.1f} min")

genes = results[0][4]; ngenes = results[0][3]
order = [r[0] for r in results]
rows = np.concatenate([np.full(len(r[1]), i, np.int32) for i, r in enumerate(results)])
cols = np.concatenate([r[1] for r in results])
vals = np.concatenate([r[2] for r in results])
X = sp.csr_matrix((vals, (rows, cols)), shape=(len(order), ngenes * ngenes))
del rows, cols, vals, results
print(f"  full sparse matrix: {X.shape}  nnz={X.nnz:,}")

# drop always-zero columns (keep the active edge scaffold) ----------------
col_nnz = np.asarray((X != 0).sum(axis=0)).ravel()
keep = np.flatnonzero(col_nnz > 0)
Xk = X[:, keep].tocsr()
gi, gj = np.unravel_index(keep, (ngenes, ngenes))
edge_src = np.array(genes)[gi]
edge_dst = np.array(genes)[gj]
print(f"  active-edge matrix: {Xk.shape}  ({len(keep):,} edges)")

# 5. expression for the same cells (HVG) ----------------------------------
expr = rna[order].to_memory()
E = expr.X.toarray() if sp.issparse(expr.X) else np.asarray(expr.X)
E = E.astype(np.float32)
print(f"  expression matrix: {E.shape}")

# 6. save -----------------------------------------------------------------
sp.save_npz(OUT / "spatial_grn_matrix.npz", Xk)
np.savez(OUT / "spatial_grn_index.npz",
         cell_ids=np.array(order), edge_src=edge_src, edge_dst=edge_dst,
         edge_flat_idx=keep.astype(np.int64))
np.savez(OUT / "spatial_expr.npz", cell_ids=np.array(order),
         expr=E, gene_names=np.array(list(expr.var_names)))
meta.loc[order].to_parquet(OUT / "spatial_meta.parquet")

print("\nSUMMARY")
print(f"  cells: {len(order):,}  | edges: {Xk.shape[1]:,}  | HVG: {E.shape[1]}")
print(f"  coords kept: argmax (x,y)" + (" + centroid (x_centroid,y_centroid)" if "x_centroid" in meta.columns else ""))
print(f"  per-subtype: {meta['subtype'].value_counts().to_dict()}")
print(f"  outputs in {OUT}")
print("STEP 4 done", now())
