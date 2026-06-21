import anndata as ad, warnings
warnings.filterwarnings('ignore')
a = ad.read_h5ad('/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/processed/seaad/Tangram_rna_with_coords_matched.h5ad')
cols = ['Donor ID','Sex','Age at Death','Braak','APOE Genotype','Continuous Pseudo-progression Score','Cognitive Status']
df = a.obs[cols].drop_duplicates(subset=['Donor ID']).sort_values('Donor ID')
print(df.to_string())
