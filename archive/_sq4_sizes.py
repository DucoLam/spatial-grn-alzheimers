"""Cheap sizes probe for SQ4 cost planning (meta + index only, no 3.7GB matrix)."""
import os, numpy as np, pandas as pd
os.chdir("/tudelft.net/staff-umbrella/ScReNI/dflam")

m = pd.read_parquet("data/spatial_link/spatial_meta.parquet")
m["section"] = m["section"].astype(str)
print("meta cols:", list(m.columns))
print("\n=== cells per subtype ===")
print(m["subtype"].value_counts())

print("\n=== per (subtype, section) cell counts ===")
for st in ["Astrocyte", "L2/3 IT", "L4 IT", "Oligodendrocyte"]:
    d = m[m.subtype == st]
    sc = d.groupby("section").size()
    sc = sc[sc >= 30]                       # sections used (>=30 cells)
    print(f"\n{st}: total={len(d)}, sections>=30={len(sc)}, "
          f"cells_in_used_sections={int(sc.sum())}")
    print("  section sizes:", sorted(sc.values.tolist(), reverse=True))

z = np.load("data/spatial_link/spatial_grn_index.npz", allow_pickle=True)
print("\nn regulators (unique edge_src):", len(np.unique(z["edge_src"])))
