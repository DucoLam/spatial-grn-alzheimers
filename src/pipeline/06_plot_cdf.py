"""
06_plot_cdf.py -- cumulative error-distribution (CDF) plots for Tangram validation.

x = distance threshold (um), y = % of held-out cells with error <= x.
Reads dflam/data/validation/<prefix>_errors.npz (per-cell error arrays per method).

Produces:
  1) CDF per subtype for the chosen method (default argmax) -- one curve per subtype
  2) CDF per method for each subtype -- argmax vs top-k comparison
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

VAL_DIR  = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/data/validation")
OUT_DIR  = Path("/tudelft.net/staff-umbrella/ScReNI/dflam/plots/cdf")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# prefix -> nice label
SUBTYPES = {
    "astrocyte": "Astrocyte",
    "l23_it":    "L2/3 IT",
    "l4_it":     "L4 IT",
    "oligo":     "Oligodendrocyte",
}
PRIMARY_METHOD = "argmax"
XMAX = 6000          # um
REF_LINES = [500, 1000, 2000]

def cdf_xy(errs, xmax=XMAX, n=600):
    """Return (x grid, % of cells <= x)."""
    errs = np.sort(np.asarray(errs))
    xs = np.linspace(0, xmax, n)
    ys = np.searchsorted(errs, xs, side="right") / len(errs) * 100.0
    return xs, ys

# ---- load ------------------------------------------------------------------
data = {}
for prefix, label in SUBTYPES.items():
    f = VAL_DIR / f"{prefix}_errors.npz"
    if f.exists():
        data[prefix] = dict(np.load(f, allow_pickle=True))
        print(f"loaded {f.name}: methods={list(data[prefix].keys())}, "
              f"n={len(next(iter(data[prefix].values())))}")
    else:
        print(f"MISSING: {f}")

if not data:
    sys.exit("No *_errors.npz files found. Run validation (with error saving) first.")

# ---- Figure 1: one curve per subtype, primary method -----------------------
fig, ax = plt.subplots(figsize=(8, 6))
colors = {"astrocyte": "#d62728", "l23_it": "#1f77b4", "l4_it": "#2ca02c", "oligo": "#9467bd"}
for prefix, label in SUBTYPES.items():
    if prefix not in data or PRIMARY_METHOD not in data[prefix]:
        continue
    errs = data[prefix][PRIMARY_METHOD]
    xs, ys = cdf_xy(errs)
    med = np.median(errs)
    ax.plot(xs, ys, lw=2.2, color=colors.get(prefix), label=f"{label}  (median {med:.0f} um, n={len(errs):,})")

for r in REF_LINES:
    ax.axvline(r, ls=":", color="grey", lw=1)
    ax.text(r, 2, f"{r}", rotation=90, va="bottom", ha="right", fontsize=8, color="grey")

ax.set_xlabel("Distance threshold (um)")
ax.set_ylabel("% of held-out cells within distance")
ax.set_title(f"Tangram spatial-mapping accuracy (CDF) -- method: {PRIMARY_METHOD}\n"
             f"per-section matched, 140 MERFISH genes, 5% holdout")
ax.set_xlim(0, XMAX); ax.set_ylim(0, 100)
ax.grid(alpha=0.3); ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
out1 = OUT_DIR / "cdf_by_subtype_argmax.png"
fig.savefig(out1, dpi=150); print(f"wrote {out1}")

# ---- Figure 2: per-subtype method comparison (small multiples) -------------
present = [p for p in SUBTYPES if p in data]
fig, axes = plt.subplots(1, len(present), figsize=(5*len(present), 4.5), sharey=True)
if len(present) == 1:
    axes = [axes]
for ax, prefix in zip(axes, present):
    for m in data[prefix].keys():
        errs = data[prefix][m]
        xs, ys = cdf_xy(errs)
        lw = 2.5 if m == PRIMARY_METHOD else 1.2
        ax.plot(xs, ys, lw=lw, label=m)
    for r in REF_LINES:
        ax.axvline(r, ls=":", color="grey", lw=0.8)
    ax.set_title(SUBTYPES[prefix]); ax.set_xlabel("Distance (um)")
    ax.set_xlim(0, XMAX); ax.set_ylim(0, 100); ax.grid(alpha=0.3)
axes[0].set_ylabel("% within distance")
axes[-1].legend(fontsize=7, loc="lower right")
fig.suptitle("CDF by method (argmax bold) -- per-section matched", y=1.02)
fig.tight_layout()
out2 = OUT_DIR / "cdf_by_method.png"
fig.savefig(out2, dpi=150, bbox_inches="tight"); print(f"wrote {out2}")
print("done")
