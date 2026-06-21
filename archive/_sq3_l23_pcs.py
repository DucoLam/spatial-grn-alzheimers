"""Print L2/3 IT per-PC: CPS correlation (SQ3) next to spatial strength (SQ2 Moran's I)."""
import csv, math
CORR = "/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ3/sq3_perpc_corr.csv"
MOR  = "/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ2/morans_i/morans_l23_it.csv"

def f(x):
    try: return float(x)
    except (ValueError, TypeError): return None

# SQ2 spatial strength per GRN PC
sq2 = {}
for r in csv.DictReader(open(MOR)):
    if r.get("space") == "GRN":
        sq2[int(float(r["pc"]))] = f(r["morans_I"])

rows = [r for r in csv.DictReader(open(CORR)) if r.get("subtype") == "L2/3 IT"]
print("cols in sq3_perpc_corr.csv:", list(rows[0].keys()) if rows else "EMPTY")
print()
print(f"{'PC':>3} {'CPS_r':>8} {'CPS_p':>8} {'q':>8} {'SQ2_Moran':>10}  selected?")
# detect a selected/spatial flag column if present
flag = next((c for c in (rows[0].keys() if rows else [])
             if c.lower() in ("selected", "is_spatial", "top5", "kept")), None)
for r in sorted(rows, key=lambda r: int(float(r["pc"]))):
    pc = int(float(r["pc"]))
    rr = f(r.get("r") or r.get("spearman_r") or r.get("rho"))
    pp = f(r.get("p") or r.get("p_value") or r.get("pval"))
    qq = f(r.get("q") or r.get("q_fdr") or r.get("bh"))
    sel = r.get(flag) if flag else ""
    rs = f"{rr:>8.3f}" if rr is not None else f"{'NA':>8}"
    ps = f"{pp:>8.3f}" if pp is not None else f"{'NA':>8}"
    qs = f"{qq:>8.3f}" if qq is not None else f"{'NA':>8}"
    ms = f"{sq2.get(pc, float('nan')):>10.3f}"
    print(f"{pc:>3} {rs} {ps} {qs} {ms}  {sel}")
