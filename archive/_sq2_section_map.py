"""Standalone SQ2 spatial cluster map for ONE section: L2/3 IT cells plotted at their
(x,y) coordinates, coloured by SQ2 GRN cluster. Optional argv[1] = section name."""
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import analysis_common as ac

ST = "L2/3 IT"
PURITY_P = 0.001
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"
palette = ["#3b6fb6", "#e8703a", "#3aa86b", "#9467bd", "#d6a800"]  # high-contrast

meta = pd.read_parquet(ac.CACHE / "spatial_meta.parquet")
meta["section"] = meta["section"].astype(str)
m = meta[meta["subtype"] == ST].copy()

lab_df = pd.read_csv(LAB_CSV)
lab_map = dict(zip(lab_df.cell_id.astype(str), lab_df.grn_cluster.astype(int)))
m["cluster"] = m.index.astype(str).map(lab_map)
m = m[m["cluster"].notna()].copy()
m["cluster"] = m["cluster"].astype(int)

avail = m["section"].value_counts()
if len(sys.argv) > 1:
    q = sys.argv[1]
    if q in set(avail.index):
        sec = q
    else:  # robust match: prefix / substring (handles partial section names)
        hits = [s for s in avail.index if s.startswith(q) or q in s]
        if not hits:
            print(f"no section matching '{q}'. available L2/3 sections:")
            for s, n in avail.items():
                print(f"  {s}  (n={n})")
            sys.exit(1)
        sec = max(hits, key=lambda s: avail[s])
        print(f"matched '{q}' -> {sec}")
else:
    sec = avail.index[0]
g = m[m["section"] == sec]
print(f"section {sec}: {len(g)} L2/3 cells | cluster sizes {np.bincount(g.cluster)}")

fig, ax = plt.subplots(figsize=(7.5, 7.5))
for c in sorted(g["cluster"].unique()):
    sub = g[g["cluster"] == c]
    ax.scatter(sub["x"], sub["y"], s=22, color=palette[c % len(palette)],
               alpha=0.5, edgecolors="none", label=f"SQ2 cluster {c}  (n={len(sub)})")
ax.set_aspect("equal")
ax.set_axisbelow(True)
ax.minorticks_on()
ax.grid(True, which="major", ls=":", lw=0.7, color="0.7", alpha=0.85)
ax.grid(True, which="minor", ls=":", lw=0.35, color="0.85", alpha=0.5)
ax.tick_params(which="major", direction="out", length=5, labelsize=9)
ax.tick_params(which="minor", direction="out", length=2.5)
ax.set_xlabel(r"x position ($\mu$m)", fontsize=11)
ax.set_ylabel(r"y position ($\mu$m)", fontsize=11)
ax.set_title(f"{ST} — section {sec}\ncells coloured by SQ2 GRN cluster "
             f"(neighbour-purity p={PURITY_P})", fontsize=12)
ax.legend(fontsize=11, loc="best", framealpha=0.9)
fig.tight_layout()
out = ac.ANALYSIS / "SQ2" / f"sq2_section_map_{sec.replace('.', '_')}.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved", out)
