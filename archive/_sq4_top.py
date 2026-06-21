"""Read the aggregated SQ4 master CSV: top regulators per subtype + cross-subtype recurrence."""
import csv, math
from collections import defaultdict
P = "/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ4/sq4_regulator_morans.csv"

rows = list(csv.DictReader(open(P)))
def f(x):
    try:
        v = float(x); return None if math.isnan(v) else v
    except (ValueError, TypeError): return None

subs = ["Astrocyte", "L2/3 IT", "L4 IT", "Oligodendrocyte"]
topsets = {}
for st in subs:
    sr = [r for r in rows if r["subtype"] == st and f(r["I_R"]) is not None]
    sr.sort(key=lambda r: f(r["I_R"]), reverse=True)
    topsets[st] = [r["regulator"] for r in sr[:20]]
    nsig = sum(1 for r in sr if f(r["q_fdr"]) is not None and f(r["q_fdr"]) < 0.05)
    print(f"=== {st}  testable={len(sr)}  q<0.05={nsig} ===")
    print(f"  {'rank regulator':<16}{'I_R':>8}{'p':>8}{'q':>8}{'n_act':>7}{'n_sec':>6}{'Is_med':>8}")
    for i, r in enumerate(sr[:12], 1):
        q = f(r["q_fdr"])
        print(f"  {i:>2} {r['regulator']:<13}{f(r['I_R']):>8.4f}{f(r['p_perm']):>8.4f}"
              f"{(q if q is not None else float('nan')):>8.4f}{r['n_active']:>7}"
              f"{r['n_sections']:>6}{f(r['I_s_median']):>8.4f}")
    print()

# cross-subtype recurrence among each subtype's top 20
cnt = defaultdict(list)
for st, regs in topsets.items():
    for g in regs:
        cnt[g].append(st)
rec = {g: sts for g, sts in cnt.items() if len(sts) >= 2}
print(f"=== regulators in the top-20 of >=2 subtypes ({len(rec)}) ===")
for g, sts in sorted(rec.items(), key=lambda kv: -len(kv[1])):
    print(f"  {g:<14} {len(sts)}x: {', '.join(sts)}")
