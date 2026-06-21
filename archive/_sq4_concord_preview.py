"""Free preview of sensitivity #1: does variance-weighted I_R agree with an equal-vote
proxy (median per-section I_s) already stored in the master CSV? No re-run needed."""
import csv, math
P = "/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ4/sq4_regulator_morans.csv"

def f(x):
    try:
        v = float(x); return None if math.isnan(v) else v
    except (ValueError, TypeError): return None

def rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v); i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]: j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1): r[order[k]] = avg
        i = j + 1
    return r

def spearman(a, b):
    ra, rb = rank(a), rank(b); n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    va = sum((x - ma) ** 2 for x in ra); vb = sum((y - mb) ** 2 for y in rb)
    return cov / math.sqrt(va * vb) if va > 0 and vb > 0 else float("nan")

rows = list(csv.DictReader(open(P)))
for st in ["Astrocyte", "L2/3 IT", "L4 IT", "Oligodendrocyte"]:
    g = [(f(r["I_R"]), f(r["I_s_median"]), r["regulator"])
         for r in rows if r["subtype"] == st]
    g = [x for x in g if x[0] is not None and x[1] is not None]
    IR = [x[0] for x in g]; IM = [x[1] for x in g]
    rho = spearman(IR, IM)
    topw = set(x[2] for x in sorted(g, key=lambda x: -x[0])[:20])
    tope = set(x[2] for x in sorted(g, key=lambda x: -x[1])[:20])
    print(f"{st:<16} n={len(g):>3}  Spearman(I_R, median I_s)={rho:.3f}  "
          f"top20 overlap={len(topw & tope)}/20")
