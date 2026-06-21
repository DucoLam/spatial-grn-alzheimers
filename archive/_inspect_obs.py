import anndata as ad
sp = ad.read_h5ad('/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad', backed='r')
print('MERFISH n_obs:', sp.n_obs)
print('MERFISH obs cols:', list(sp.obs.columns))
for c in ('Class','Subclass','Supertype','cell_type'):
    if c in sp.obs.columns:
        print('  MERFISH', c, 'nunique:', sp.obs[c].nunique())
if 'Subclass' in sp.obs.columns:
    print('  MERFISH Subclass cats:', sorted(sp.obs['Subclass'].astype(str).unique())[:40])
print('MERFISH obsm:', list(sp.obsm.keys()))
print()
rna = ad.read_h5ad('/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/seaad_paired_rna.h5ad', backed='r')
print('RNA shape:', rna.shape)
print('RNA obs cols:', list(rna.obs.columns))
if 'Subclass' in rna.obs.columns:
    print('  RNA Subclass nunique:', rna.obs['Subclass'].nunique())
    print('  RNA Subclass cats:', sorted(rna.obs['Subclass'].astype(str).unique()))
print('RNA Donor ID nunique:', rna.obs['Donor ID'].nunique() if 'Donor ID' in rna.obs.columns else 'n/a')
