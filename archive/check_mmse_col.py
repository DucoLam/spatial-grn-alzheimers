import anndata as ad
r = ad.read_h5ad("/tudelft.net/staff-umbrella/ScReNI/dflam/data/seaad_paired_rna_pool4.h5ad", backed="r")
hits = [c for c in r.obs.columns if any(x in c.lower() for x in ["mmse", "cogn", "mental", "score"])]
print(hits)
