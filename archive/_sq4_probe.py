"""SQ4 feasibility probe: confirm GRN edge->TF mapping and TF-activity availability."""
import os, numpy as np, scipy.sparse as sp, pandas as pd
os.chdir("/tudelft.net/staff-umbrella/ScReNI/dflam")

print("=== spatial_link dir ===")
for f in sorted(os.listdir("data/spatial_link")):
    p = os.path.join("data/spatial_link", f)
    print(f"  {f}  {os.path.getsize(p)/1e6:.1f} MB")

print("\n=== grn index npz keys ===")
z = np.load("data/spatial_link/spatial_grn_index.npz", allow_pickle=True)
for k in z.keys():
    a = z[k]
    print(f"  {k}: shape={getattr(a,'shape',None)} dtype={getattr(a,'dtype',None)}")
    try:
        print("     head:", a.ravel()[:6])
    except Exception as e:
        print("     (head err)", e)

print("\n=== GRN matrix ===")
X = sp.load_npz("data/spatial_link/spatial_grn_matrix.npz")
print("  shape:", X.shape, "nnz/cell:", X.nnz / X.shape[0])

# look for a feature/edge name file anywhere obvious
print("\n=== candidate feature-name files ===")
for root in ["data/spatial_link", "data", "analysis"]:
    for f in os.listdir(root):
        if any(t in f.lower() for t in ["feature", "edge", "tf", "regul", "gene", "col", "var"]):
            print("  ", os.path.join(root, f))

print("\n=== meta columns ===")
m = pd.read_parquet("data/spatial_link/spatial_meta.parquet")
print("  ", list(m.columns))
