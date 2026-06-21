import pandas as pd
pd.set_option("display.width", 200, "display.max_columns", 30)
meta = pd.read_parquet("/tudelft.net/staff-umbrella/ScReNI/dflam/data/spatial_link/spatial_meta.parquet")
print("n cells:", len(meta), "| columns:", list(meta.columns))
print("n donors:", meta["Donor ID"].nunique())

sev = ["Braak", "Overall AD neuropathological Change", "Continuous Pseudo-progression Score"]
sev = [c for c in sev if c in meta.columns]

print("\n--- severity per donor (first value) + cell count ---")
dd = meta.groupby("Donor ID").agg(n_cells=("subtype", "size"),
                                  **{c: (c, "first") for c in sev})
print(dd)

print("\n--- cells per donor x subtype ---")
print(pd.crosstab(meta["Donor ID"], meta["subtype"]))

print("\n--- donors that have >=30 located cells, per subtype ---")
ct = pd.crosstab(meta["Donor ID"], meta["subtype"])
print((ct >= 30).sum().to_dict())

for c in sev:
    print(f"\n--- value counts: {c} ---")
    print(meta.groupby("Donor ID")[c].first().value_counts(dropna=False))
