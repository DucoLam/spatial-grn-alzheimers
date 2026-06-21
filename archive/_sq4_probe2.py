"""Cheap SQ4 probe: TF count, value type, meta cols (index-only, no 3.7GB load)."""
import os, numpy as np, scipy.sparse as sp, pandas as pd
os.chdir("/tudelft.net/staff-umbrella/ScReNI/dflam")

z = np.load("data/spatial_link/spatial_grn_index.npz", allow_pickle=True)
src = z["edge_src"]; dst = z["edge_dst"]
print("n_edges:", len(src))
print("n_unique_TF (edge_src):", len(np.unique(src)))
print("n_unique_target (edge_dst):", len(np.unique(dst)))
u, c = np.unique(src, return_counts=True)
o = np.argsort(c)[::-1]
print("top-10 TFs by #targets:", list(zip(u[o][:10], c[o][:10])))
print("regulon size: median", int(np.median(c)), "min", c.min(), "max", c.max())

# peek at matrix values (load just data array header cheaply via npz)
X = sp.load_npz("data/spatial_link/spatial_grn_matrix.npz").tocsr()
print("matrix shape:", X.shape, "dtype:", X.dtype, "nnz/cell:", round(X.nnz/X.shape[0],1))
print("value sample (first 12 nz):", X.data[:12])
print("value min/max/mean:", X.data.min(), X.data.max(), round(float(X.data.mean()),4))

m = pd.read_parquet("data/spatial_link/spatial_meta.parquet")
print("meta cols:", list(m.columns))
print("subtype counts:\n", m["subtype"].value_counts())
