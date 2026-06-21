import anndata as ad
import numpy as np

SPAGE = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad_spatial/SpaGE_500_HVGs_FULL.h5ad"
MERF  = "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad"

sp = ad.read_h5ad(SPAGE, backed="r")
print("SpaGE shape:", sp.shape)
print("SpaGE obs cols:", list(sp.obs.columns))
print("SpaGE obsm keys:", list(sp.obsm.keys()))
for c in ("Subclass", "Section", "Donor ID", "Supertype", "cell_type"):
    if c in sp.obs.columns:
        print(f"  has {c}: {sp.obs[c].nunique()} unique")

mf = ad.read_h5ad(MERF, backed="r")
print("\nMERFISH shape:", mf.shape)
# do obs_names align?
same_order = bool(np.array_equal(sp.obs_names.values[:10000], mf.obs_names.values[:10000]))
same_set = None
try:
    same_set = bool(set(sp.obs_names[:50000]) <= set(mf.obs_names))
except Exception as e:
    same_set = f"err {e}"
print("obs_names first-10k identical order:", same_order)
print("SpaGE first-50k subset of MERFISH names:", same_set)

# subtype labels present in SpaGE for our 4?
for col in ("Subclass", "Supertype", "cell_type"):
    if col in sp.obs.columns:
        vals = sp.obs[col].astype(str).unique()
        hits = [v for v in vals if any(t in v for t in ("Astro", "L2/3 IT", "L4 IT", "Oligo"))]
        print(f"  {col} subtype-like values:", hits[:12])
