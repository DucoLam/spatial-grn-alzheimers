import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
import analysis_common as ac

ST = "L2/3 IT"
SECTION = "H21.33.019.Cx30.MTG.02.007.5.1"
KS = [2, 4, 6, 8]
OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "grn_clustering")

X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
rows = ac.subtype_rows(meta, ST)          # sorted ascending
sect = meta["section"].values
xy = meta[["x", "y"]].values

in_sec = rows[sect[rows] == SECTION]
if len(in_sec) == 0:
    print("NO cells for section", SECTION)
    print("available L2/3 sections containing 'H21.33.019':")
    for s in sorted(set(sect[rows])):
        if "H21.33.019" in s:
            print("   ", s, int((sect[rows] == s).sum()))
    raise SystemExit

print(f"{ST} in {SECTION}: {len(in_sec)} cells")
G = TruncatedSVD(n_components=min(ac.N_PCS, len(rows) - 1),
                 random_state=ac.SEED).fit_transform(X[rows])
row_index = {r: i for i, r in enumerate(rows)}
gi = np.array([row_index[r] for r in in_sec])

fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for ax, k in zip(axes.ravel(), KS):
    km = KMeans(n_clusters=k, n_init=10, random_state=ac.SEED).fit(G)
    labels = km.labels_[gi]
    ax.scatter(xy[in_sec, 0], xy[in_sec, 1], c=labels, cmap="tab10",
               s=12, vmin=0, vmax=9)
    ax.set_title(f"k = {k}")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(f"GRN clusters in space — {ST}, section {SECTION} (n={len(in_sec)})")
fig.tight_layout()
fig.savefig(OUT / "map_l23_H21019_k2468.png", dpi=160)
print("saved map_l23_H21019_k2468.png")
