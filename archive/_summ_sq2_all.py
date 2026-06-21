import pandas as pd, os
pd.set_option("display.width", 200, "display.max_columns", 30)
base = "/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ2"
slugs = {"Astrocyte": "astrocyte", "L2/3 IT": "l23_it", "L4 IT": "l4_it", "Oligodendrocyte": "oligo"}

print("================ MORANّS I (raw GRN vs expression) ================")
for name, s in slugs.items():
    f = f"{base}/morans_i/morans_{s}.csv"
    if not os.path.exists(f):
        continue
    d = pd.read_csv(f)
    for sp in ["GRN", "expression"]:
        g = d[d.space == sp]
        print(f"{name:16s} {sp:11s} meanI={g.morans_I.mean():.3f} maxI={g.morans_I.max():.3f} "
              f"sig(p<.05)={int((g.p_perm<.05).sum())}/{len(g)}")

print("\n================ CONDITIONAL: residual Moran's I ================")
for name, s in slugs.items():
    f = f"{base}/conditional/residual_morans_{s}.csv"
    if not os.path.exists(f):
        continue
    d = pd.read_csv(f)
    print(f"{name:16s} rawMeanI={d.morans_raw.mean():.3f} residMeanI={d.morans_resid.mean():.3f} "
          f"sigRaw={int((d.p_raw<.05).sum())}/{len(d)} sigResid={int((d.p_resid<.05).sum())}/{len(d)} "
          f"meanExprR2={d.expr_R2.mean():.3f} maxExprR2={d.expr_R2.max():.3f}")

print("\n================ PARTIAL MANTEL ================")
f = f"{base}/conditional/partial_mantel.csv"
if os.path.exists(f):
    print(pd.read_csv(f).to_string(index=False))
