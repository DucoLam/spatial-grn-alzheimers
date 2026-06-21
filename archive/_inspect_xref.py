import anndata as ad
sp = ad.read_h5ad('/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad', backed='r')
rna = ad.read_h5ad('/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/seaad_paired_rna.h5ad', backed='r')

print('=== MERFISH Section examples ===')
print(list(sp.obs['Section'].astype(str).unique())[:4])
print('MERFISH sections total:', sp.obs['Section'].nunique())

print()
print('=== RNA candidate sample/section columns ===')
for c in ('sample_id','sample_name','Section','Specimen Barcode','Donor ID'):
    if c in rna.obs.columns:
        vals = list(rna.obs[c].astype(str).unique())
        print(f'  {c}: nunique={len(vals)}  e.g. {vals[:3]}')

# How many RNA samples per donor (overlapping donors)?
donors = ['H20.33.001','H20.33.004','H20.33.012','H20.33.015','H20.33.025','H20.33.040','H20.33.044','H21.33.001','H21.33.006','H21.33.019']
print()
print('=== RNA samples per donor ===')
scol = 'sample_id' if 'sample_id' in rna.obs.columns else 'sample_name'
for d in donors:
    sub = rna.obs[rna.obs['Donor ID'].astype(str)==d]
    print(f'  {d}: {len(sub):,} cells, {sub[scol].nunique()} {scol}(s): {list(sub[scol].astype(str).unique())[:4]}')

# direct intersection of any RNA id set with MERFISH sections
mer_secs = set(sp.obs['Section'].astype(str).unique())
for c in ('sample_id','sample_name','Specimen Barcode'):
    if c in rna.obs.columns:
        inter = mer_secs & set(rna.obs[c].astype(str).unique())
        print(f'Intersection MERFISH Section ∩ RNA {c}: {len(inter)}')
