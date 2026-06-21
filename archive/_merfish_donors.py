"""Count sections/donors in the FULL MERFISH dataset, and build a per-donor
demographic/pathology table for the 10 overlapping donors (sex, age, Braak, tau, CPS)."""
import anndata as ad
import pandas as pd
import numpy as np
import analysis_common as ac

H5 = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
adata = ad.read_h5ad(H5, backed="r")
obs = adata.obs.copy()
print("OBS COLS:", list(obs.columns))
print("FULL MERFISH n_cells:", adata.n_obs)

def find(*keys):
    # exact (case-insensitive) match first, then substring
    for c in obs.columns:
        if c.lower() in keys:
            return c
    for c in obs.columns:
        cl = c.lower()
        if all(k in cl for k in keys):
            return c
    return None

sec = "Section" if "Section" in obs.columns else find("section")
bc = "Specimen Barcode" if "Specimen Barcode" in obs.columns else None
don = "Donor ID" if "Donor ID" in obs.columns else find("donor")
print("section col:", sec, "| barcode col:", bc, "| donor col:", don)
print("FULL MERFISH n unique Section:", obs[sec].astype(str).nunique() if sec else "??")
if bc:
    print("FULL MERFISH n unique Specimen Barcode:", obs[bc].astype(str).nunique())
print("FULL MERFISH n_donors:", obs[don].nunique() if don else "??")
if sec and don:
    spd = obs.groupby(don, observed=True)[sec].apply(lambda s: s.astype(str).nunique())
    print("FULL MERFISH sections/donor: min/median/max =",
          int(spd.min()), float(spd.median()), int(spd.max()),
          "| total =", int(spd.sum()))

# locate demographic / pathology columns
cols = {
    "sex":   find("sex") or find("gender"),
    "age":   find("age", "death") or find("age"),
    "braak": find("braak"),
    "thal":  find("thal"),
    "cerad": find("cerad"),
    "adnc":  find("neuropathological", "change"),
    "tau":   find("tau"),
    "cps":   find("pseudo") or find("cps"),
    "apoe":  find("apoe"),
    "cog":   find("cognitive", "status"),
}
print("\nresolved metadata columns:")
for k, v in cols.items():
    print(f"  {k:6} -> {v}")

# the 10 overlap donors (from our analysis cache)
m = pd.read_parquet(ac.CACHE / "spatial_meta.parquet")
overlap = sorted(m["Donor ID"].astype(str).unique())
print(f"\n10 overlap donors: {overlap}")

obs[don] = obs[don].astype(str)
sub = obs[obs[don].isin(overlap)]
rows = []
for d, g in sub.groupby(don, observed=True):
    rec = {"donor": d, "merfish_cells": len(g),
           "merfish_sections": g[sec].astype(str).nunique() if sec else None}
    for k, v in cols.items():
        if v is not None:
            vals = g[v].dropna().unique()
            rec[k] = vals[0] if len(vals) else None
    rows.append(rec)
tab = pd.DataFrame(rows).sort_values("donor")
pd.set_option("display.width", 200, "display.max_columns", 30)
print("\n=== per-donor table (10 overlap donors) ===")
print(tab.to_string(index=False))
tab.to_csv(ac.ANALYSIS / "overlap_donor_metadata.csv", index=False)

# quick cohort summary
print("\n=== cohort summary ===")
if cols["sex"]:
    print("sex:", sub.groupby(don, observed=True)[cols["sex"]].first().value_counts().to_dict())
if cols["age"]:
    raw = sub.groupby(don, observed=True)[cols["age"]].first().astype(str)
    ages = pd.to_numeric(raw.str.replace("+", "", regex=False), errors="coerce")
    n90 = (raw.str.contains(r"90\+")).sum()
    print(f"age: {ages.min():.0f}-{ages.max():.0f} (median {ages.median():.0f}); "
          f"{n90} donor(s) coded 90+")
print("braak:", sub.groupby(don, observed=True)[cols['braak']].first().value_counts().to_dict())
print("cog:", sub.groupby(don, observed=True)[cols['cog']].first().value_counts().to_dict())
print("CPS range:", round(sub.groupby(don, observed=True)[cols['cps']].first().min(),3),
      "-", round(sub.groupby(don, observed=True)[cols['cps']].first().max(),3))
print("saved", ac.ANALYSIS / "overlap_donor_metadata.csv")
