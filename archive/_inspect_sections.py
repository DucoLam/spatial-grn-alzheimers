import anndata as ad
import numpy as np
sp = ad.read_h5ad('/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/seaad_spatial/SEAAD_MTG_MERFISH.2024-12-11.h5ad', backed='r')
obs = sp.obs
print('Section col present:', 'Section' in obs.columns)
print('Total sections:', obs['Section'].nunique())
print()
# sections per donor (overlapping donors)
donors = ['H20.33.001','H20.33.004','H20.33.012','H20.33.015','H20.33.025','H20.33.040','H20.33.044','H21.33.001','H21.33.006','H21.33.019']
for d in donors:
    sub = obs[obs['Donor ID'].astype(str)==d]
    secs = sub['Section'].astype(str).unique()
    print(f'{d}: {len(sub):,} cells, {len(secs)} section(s): {list(secs)[:6]}')
print()
# coordinate frames: do sections share a coordinate range or are they separate?
coords = sp.obsm['X_spatial_raw']
d = 'H20.33.001'
sub_mask = (obs['Donor ID'].astype(str)==d).values
sub_obs = obs[obs['Donor ID'].astype(str)==d]
sub_coords = coords[sub_mask]
print(f'--- {d} coordinate ranges per section (X_spatial_raw) ---')
for s in sub_obs['Section'].astype(str).unique()[:6]:
    m = (sub_obs['Section'].astype(str)==s).values
    c = sub_coords[m]
    print(f'  {s}: x[{c[:,0].min():.0f},{c[:,0].max():.0f}] y[{c[:,1].min():.0f},{c[:,1].max():.0f}]  n={m.sum():,}')
