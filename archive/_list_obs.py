import anndata as ad
rna = ad.read_h5ad("/tudelft.net/staff-umbrella/ScReNI/dflam/data/seaad_paired_rna_pool4.h5ad", backed="r")
cols = list(rna.obs.columns)
print("n obs columns:", len(cols))
print("\n--- ALL obs columns ---")
for c in cols:
    print("  ", c)

kw = ("braak", "thal", "cerad", "amyloid", "apoe", "patholog", "ad ", "adnc",
      "cognit", "demen", "lewy", "late", "molecular", "neuropath", "progression",
      "plaque", "tau", "mmse", "age", "sex")
print("\n--- columns matching pathology/clinical keywords ---")
for c in cols:
    if any(k in c.lower() for k in kw):
        try:
            nun = rna.obs[c].nunique()
            print(f"  {c}   ({nun} unique)")
        except Exception:
            print(f"  {c}")
