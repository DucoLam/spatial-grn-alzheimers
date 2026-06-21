import pandas as pd, glob, os
os.chdir("/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ2/morans_i")
for f in sorted(glob.glob("morans_*.csv")):
    d = pd.read_csv(f)
    g = d.groupby("space")
    agg = g["morans_I"].agg(["mean", "max"]).round(3)
    sig = d.assign(s=d.p_perm < 0.05).groupby("space")["s"].sum().to_dict()
    print(f"\n=== {f} ===")
    print(agg)
    print("  sig PCs (p<0.05):", sig)
