import pandas as pd, numpy as np, anndata as ad, glob
pd.set_option("display.width", 220, "display.max_columns", 40)
CACHE = "/tudelft.net/staff-umbrella/ScReNI/dflam/data/spatial_link"
meta = pd.read_parquet(CACHE + "/spatial_meta.parquet")
meta["Donor ID"] = meta["Donor ID"].astype(str)
meta["section"] = meta["section"].astype(str)

print("===== Q1: tissue sections in the analysis set =====")
print("unique sections total:", meta["section"].nunique())
print("(donor, section) pairs:", meta.groupby(["Donor ID", "section"], observed=True).ngroups)
print("sections per donor:")
print(meta.groupby("Donor ID", observed=True)["section"].nunique().sort_values(ascending=False))

print("\n===== Q2: severity metadata across the located donors =====")
loc = set(meta["Donor ID"].unique())
rna = ad.read_h5ad("/tudelft.net/staff-umbrella/ScReNI/dflam/data/seaad_paired_rna_pool4.h5ad", backed="r")
obs = rna.obs.copy()
obs["Donor ID"] = obs["Donor ID"].astype(str)
sev = ["Braak", "Thal", "CERAD score", "Overall AD neuropathological Change",
       "Overall CAA Score", "Continuous Pseudo-progression Score",
       "Cognitive Status", "Last MMSE Score", "APOE Genotype", "Age at Death", "Sex"]
sev = [c for c in sev if c in obs.columns]
dd = obs[obs["Donor ID"].isin(loc)].groupby("Donor ID", observed=True)[sev].first()
print(dd.to_string())

print("\n===== Q3: variance explained by the 20 GRN PCs (from SQ2 morans csvs) =====")
for f in sorted(glob.glob("/tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ2/morans_i/morans_*.csv")):
    d = pd.read_csv(f)
    g = d[d.space == "GRN"].sort_values("pc")
    if "var_explained" in g:
        print(f.split('/')[-1], "| 20 GRN PCs cumulative var explained:",
              round(float(g["var_explained"].sum()), 4),
              "| PC1..3:", [round(x, 3) for x in g["var_explained"].head(3)])
