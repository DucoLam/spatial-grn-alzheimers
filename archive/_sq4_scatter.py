"""SQ4 Stage-1 scatter: n_active (log) vs spatial effect size I_R, per subtype.
Dots SOLID = significant under both poolings (q_fdr<.05 AND q_eqsec<.05, "robust"),
FADED = not. Labelled = top-10 robust by I_R + LINC02248. A star on each labelled gene
shows its spatial DIRECTION (blue = regulon louder in SQ2 cluster 0, red = louder in
cluster 1) -- so LINC02248, the lone cluster-1 gene, stands out as different."""
import csv, math, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ4"
P = f"{BASE}/sq4_regulator_morans.csv"
G = f"{BASE}/sq4_gradient_check.csv"          # per-regulator mag in each SQ2 cluster
OUT = f"{BASE}/sq4_regulator_scatter.png"
EXTRA = "LINC02248"

def f(x):
    try:
        v = float(x); return None if math.isnan(v) else v
    except (ValueError, TypeError): return None

subs = ["L2/3 IT", "L4 IT", "Astrocyte", "Oligodendrocyte"]
colors = {"L2/3 IT": "#1f77b4", "L4 IT": "#2ca02c",
          "Astrocyte": "#d62728", "Oligodendrocyte": "#ff7f0e"}
C0, C1 = "#08519c", "#d62728"                  # direction colours

# louder-cluster per regulator (from the gradient check, computed on L2/3)
direction = {}
for r in csv.DictReader(open(G)):
    m0, m1 = f(r["mag_c0"]), f(r["mag_c1"])
    if m0 is not None and m1 is not None:
        direction[r["regulator"]] = "c0" if m0 >= m1 else "c1"

rows = list(csv.DictReader(open(P)))
pts = {st: {"rob": [], "non": []} for st in subs}
robpts = []
for r in rows:
    st = r["subtype"]; ir = f(r["I_R"]); na = f(r["n_active"])
    qf = f(r["q_fdr"]); qe = f(r["q_eqsec"])
    if st not in subs or ir is None or na is None or na < 1:
        continue
    robust = (qf is not None and qf < 0.05) and (qe is not None and qe < 0.05)
    pts[st]["rob" if robust else "non"].append((na, ir))
    if robust:
        robpts.append((na, ir, r["regulator"], st, qf if qf is not None else 1.0))

fig, ax = plt.subplots(figsize=(12, 7))
for st in subs:
    non = pts[st]["non"]; rob = pts[st]["rob"]
    if non:
        ax.scatter([d[0] for d in non], [d[1] for d in non], s=34, alpha=0.13,
                   color=colors[st], edgecolors="none")
    if rob:
        ax.scatter([d[0] for d in rob], [d[1] for d in rob], s=54, alpha=0.55,
                   color=colors[st], edgecolors="black", linewidths=0.3)
    n_tot = len(non) + len(rob)
    ax.scatter([], [], s=60, alpha=0.7, color=colors[st], edgecolors="black",
               linewidths=0.3, label=f"{st}  (robust {len(rob)}/{n_tot})")

ax.set_xscale("log")
ax.set_xlabel("Number of Active Cells  (Log Scale)", fontsize=17)
ax.set_ylabel("Spatial Autocorrelation  $I_R$  (Effect Size)", fontsize=17)
ax.set_title("Spatial Autocorrelation of Each Regulator's Target Sub-network",
             fontsize=17, pad=10)
from matplotlib.ticker import AutoMinorLocator
ax.tick_params(labelsize=14)
ax.set_axisbelow(True)
ax.yaxis.set_minor_locator(AutoMinorLocator())          # finer y ticks
ax.minorticks_on()
ax.tick_params(which="major", direction="out", length=7, width=1.1)
ax.tick_params(which="minor", direction="out", length=3.5, width=0.8)
ax.grid(True, which="major", ls=":", lw=0.6, color="0.78", alpha=0.8)
ax.grid(True, which="minor", axis="x", ls=":", lw=0.4, color="0.88", alpha=0.6)
ax.grid(True, which="minor", axis="y", ls=":", lw=0.4, color="0.9", alpha=0.55)
leg1 = ax.legend(fontsize=13, frameon=False, loc="lower right")
ax.add_artist(leg1)

# per-subtype MEDIAN I_R lines (neuronal only; median: I_R is right-skewed)
med_handles, med_labels = [], []
for st in ["L2/3 IT", "L4 IT"]:
    allir = [d[1] for d in pts[st]["rob"]] + [d[1] for d in pts[st]["non"]]
    if not allir:
        continue
    med = statistics.median(allir)
    line = ax.axhline(med, color=colors[st], ls="--", lw=1.4, alpha=0.75, zorder=1)
    med_handles.append(line)
    med_labels.append(f"{st}  median $I_R$={med:.3f}")
leg2 = ax.legend(med_handles, med_labels, fontsize=12, frameon=False,
                 loc="upper right")
ax.add_artist(leg2)

# labelled set: top-5 robust by I_R + LINC02248
top = sorted(robpts, key=lambda p: -p[1])[:5]
extra = [p for p in robpts if p[2] == EXTRA and p[2] not in {q[2] for q in top}]
labelled = top + extra

# zoom out a little for headroom around the labels
ymin, ymax = ax.get_ylim(); xmin, xmax = ax.get_xlim()
yr = ymax - ymin
ax.set_ylim(ymin - 0.04 * yr, ymax + 0.20 * yr)
xr = math.log10(xmax) - math.log10(xmin)
ax.set_xlim(10 ** (math.log10(xmin) - 0.06 * xr), 10 ** (math.log10(xmax) + 0.06 * xr))

# leader-line labels (single column, left)
ymin, ymax = ax.get_ylim(); xmin, xmax = ax.get_xlim()
def xlog(frac): return 10 ** (math.log10(xmin) + frac * (math.log10(xmax) - math.log10(xmin)))
x_text = xlog(0.40)
y_top = ymax - 0.03 * (ymax - ymin)
step = 0.083 * (ymax - ymin)
for i, (na, ir, name, st, qf) in enumerate(labelled):
    d = direction.get(name, "c0")
    is_x = name == EXTRA
    yl = y_top - i * step
    label = f"{name}   $I_R$={ir:.3f}, q={qf:.3f}  (↑{d})"
    ax.annotate(label, xy=(na, ir), xytext=(x_text, yl),
                fontsize=13.5, color=C1 if is_x else colors[st],
                fontweight="bold", va="center", ha="right",
                arrowprops=dict(arrowstyle="-", lw=0.8,
                                color=C1 if is_x else "0.55"))
fig.tight_layout()
fig.savefig(OUT, dpi=150)
print("saved", OUT)
for (na, ir, name, st, qf) in labelled:
    print(f"  {name:<13} {st:<16} I_R={ir:.3f} q={qf:.3f} louder={direction.get(name,'?')}")
