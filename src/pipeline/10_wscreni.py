"""
10_wscreni.py  -- parameterized wScReNI inference (with optional cell-sharding)

Sharding: pass --n_shards N --shard_id i to process only this shard's slice of
cells. Uses the function's built-in `cell_index` parameter, so the FULL
expression matrix + FULL global KNN are always passed (neighbour lookups
unchanged) and each cell's result is identical regardless of shard. Cells whose
output file already exists are skipped (no recompute / no overwrite).

Usage:
    python 10_wscreni.py --prefix pool4                       # all cells (1 job)
    python 10_wscreni.py --prefix pool4 --n_shards 16 --shard_id 3   # one shard
"""
import sys, argparse, time
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--prefix",   required=True)
parser.add_argument("--n_jobs",   type=int, default=8)
parser.add_argument("--n_trees",  type=int, default=100)
parser.add_argument("--seed",     type=int, default=42)
parser.add_argument("--n_shards", type=int, default=1, help="total number of shards")
parser.add_argument("--shard_id", type=int, default=0, help="this shard's 0-based id")
args = parser.parse_args()

sys.path.insert(0, "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/src")
from screni.data.inference import infer_wscreni_networks, GenePeakOverlapLabs

DFLAM    = Path("/tudelft.net/staff-umbrella/ScReNI/dflam")
DATA_DIR = DFLAM / "data"
OUT_DIR  = DFLAM / "data" / f"wscreni_networks_{args.prefix}"
OUT_DIR.mkdir(parents=True, exist_ok=True)

for p in [DATA_DIR / f"seaad_paired_rna_{args.prefix}.h5ad",
          DATA_DIR / f"seaad_paired_atac_{args.prefix}.h5ad",
          DATA_DIR / f"{args.prefix}_triplets.csv",
          DATA_DIR / f"{args.prefix}_peak_overlap_matrix.npz",
          DATA_DIR / f"{args.prefix}_peak_info.csv"]:
    if not p.exists():
        sys.exit(f"ERROR: {p} not found.")

rna        = ad.read_h5ad(DATA_DIR / f"seaad_paired_rna_{args.prefix}.h5ad")
atac       = ad.read_h5ad(DATA_DIR / f"seaad_paired_atac_{args.prefix}.h5ad")
triplets   = pd.read_csv(DATA_DIR / f"{args.prefix}_triplets.csv")
peak_info  = pd.read_csv(DATA_DIR / f"{args.prefix}_peak_info.csv", index_col=0)
peak_mat   = np.load(DATA_DIR / f"{args.prefix}_peak_overlap_matrix.npz")["peak_matrix"].astype(np.float64)
peak_names = list(peak_info.index)
knn_indices= rna.uns["knn_indices"]

# Gaussian noise for RF numerical stability (fixed seed -> identical across shards)
rng = np.random.default_rng(seed=args.seed)
peak_mat += rng.normal(0, 1e-5, peak_mat.shape)

labs = GenePeakOverlapLabs.from_dataframe(
    triplets.rename(columns={"target_gene": "gene.name", "peak": "peak.name"})
)

# ---- determine this shard's cells (skip ones already written) --------------
n_cells = rna.n_obs
shard_cells = np.array_split(np.arange(n_cells), args.n_shards)[args.shard_id]

wscreni_subdir = OUT_DIR / "wScReNI"
done = set()
if wscreni_subdir.exists():
    for f in wscreni_subdir.glob("*.network.txt"):
        try:
            done.add(int(f.name.split(".")[0]) - 1)   # file_num = cell_i + 1
        except ValueError:
            pass

todo = [int(i) for i in shard_cells if int(i) not in done]
print(f"wScReNI prefix={args.prefix}  shard {args.shard_id+1}/{args.n_shards}  "
      f"cells={rna.n_obs}  peaks={peak_mat.shape[1]}  triplets={len(triplets)}")
print(f"  shard cells: {len(shard_cells)}  | already done: {len(shard_cells)-len(todo)}  | to compute: {len(todo)}")
if not todo:
    print("  nothing to do for this shard — exiting cleanly.")
    sys.exit(0)

t0 = time.time()
infer_wscreni_networks(
    expr=rna, peak_mat=peak_mat, peak_names=peak_names, labs=labs,
    nearest_neighbors_idx=knn_indices, network_path=str(OUT_DIR),
    cell_index=todo,                       # <-- only this shard's cells
    n_jobs=args.n_jobs, n_trees=args.n_trees, seed=args.seed,
)
print(f"Done shard {args.shard_id+1}/{args.n_shards}: {len(todo)} networks in "
      f"{(time.time()-t0)/60:.1f} min  ->  {wscreni_subdir}")
