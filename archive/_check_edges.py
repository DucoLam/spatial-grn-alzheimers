import numpy as np
import analysis_common as ac

idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
print("keys:", list(idx.keys()))
src = idx["edge_src"].astype(str)
dst = idx["edge_dst"].astype(str)
usrc = np.unique(src)
udst = np.unique(dst)
print("n edges:", len(src))
print("n unique regulators (edge_src):", len(usrc))
print("n unique targets   (edge_dst):", len(udst))
print("src subset of dst? (gene-gene):", set(usrc).issubset(set(udst)))
print("n src that are also targets:", len(set(usrc) & set(udst)))
print("\nfirst 15 regulators:", list(usrc[:15]))
# median targets per regulator
import collections
cnt = collections.Counter(src)
vals = np.array(list(cnt.values()))
print("targets per regulator: min/median/max =", vals.min(), int(np.median(vals)), vals.max())
