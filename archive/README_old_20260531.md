# Spatially Annotated Gene Regulatory Networks in Alzheimer's Disease Astrocytes

**Author**: D. Lam — BSc Thesis, TU Delft  
**Cluster**: DAIC HPC (`login.daic.tudelft.nl`)  
**Project root**: `/tudelft.net/staff-umbrella/ScReNI/dflam/`  
**Data root**: `/tudelft.net/staff-umbrella/ScReNI/bsc-screni/data/`

---

## Research Question

> Can changes in the spatial location of astrocytes in human Alzheimer's disease tissue be related to changes in gene regulatory network structure and the activity of key transcription factors?

---

## Project Overview

This project combines two complementary data modalities from the **SEA-AD** (Seattle Alzheimer's Disease Brain Cell Atlas) dataset on the **Middle Temporal Gyrus (MTG)**:

1. **Paired multiome (RNA + ATAC-seq)** — 6,777 astrocytes from 28 donors, enabling per-cell gene regulatory network (GRN) inference via **ScReNI**
2. **MERFISH spatial transcriptomics** — 1,887,729 cells from 27 donors, providing ground-truth spatial coordinates for 140 genes

Since only 10 donors have both modalities, we use **Tangram** to infer the spatial location of each RNA cell by mapping its gene expression onto the MERFISH reference. The result is a table of **(cell_id → inferred x, y)** that, combined with ScReNI GRNs, produces **spatially annotated per-cell GRNs** — enabling analysis of where in the tissue GRN disruptions occur in AD.

### Why Astrocytes?
Astrocytes are the most abundant glial cell type and play a central role in AD pathology through reactive astrogliosis. They show substantial heterogeneity in both gene expression and spatial distribution across cortical layers, making them a strong candidate for spatially-resolved GRN analysis.

---

## Methodology

```
Paired RNA+ATAC (28 donors)
        │
        ▼
   [ScReNI]  ──────────────────────────────────────────────────────────────┐
   Per-cell GRNs (TF activity, peak-gene links, triplets)                  │
                                                                            │
MERFISH spatial (27 donors, 140 genes)                                     │
        │                                                                   │
        ├──[SpaGE]── Impute 466 HVGs not in MERFISH panel                  │
        │                                                                   │
        ▼                                                                   │
   SpaGE_500_HVGs_FULL.h5ad (1.88M cells × 500 genes)                     │
        │                                                                   │
        │  ← [Tangram, per-donor matched, 100 genes, top-0.5%]             │
        │                                                                   │
        ▼                                                                   │
   Inferred (x, y) per RNA cell  ◄─────────────────────────────────────────┘
        │
        ▼
   Spatially Annotated GRNs
   → Do GRN disruptions in AD cluster spatially?
   → Which TFs show location-dependent activity changes?
```

### Key Design Choices

| Decision | Choice | Reason |
|----------|--------|--------|
| Cell type | Astrocytes only | Central AD pathology role; manageable dataset size |
| Mapping mode | Per-donor matched (06b) | Avoids cross-donor coordinate misalignment |
| Gene set | 100 genes (34 measured MERFISH HVGs + top 66 SpaGE by HVG rank) | Best balance of real signal vs noise; validated |
| Coordinate method | Top 0.5% weighted centroid | Validated optimal: 2,412 µm median error, 7× better than random |
| Epochs | 1,000 | Standard Tangram default; sufficient convergence |

### Spatial Accuracy (Validation Results)

Validated using within-donor holdout (5% of spatial astrocytes per donor, infer location from remaining 95%):

| Gene set | Median error | Within 500 µm | vs Random (0.6%) |
|----------|-------------|---------------|-----------------|
| 34 measured HVGs | 3,101 µm | 6.6% | 11× better |
| 100 genes (production) | 2,412 µm | 4.1% | 7× better |
| 140 raw MERFISH | 2,826 µm | 16.8% | 28× better |
| 500 genes (34+466 SpaGE) | 3,142 µm | 5.1% | 8× better |

**Tissue dimensions**: ~11,000 × 12,000 µm. Median error of 2,412 µm ≈ 21% of tissue width — sufficient for broad regional analysis but not fine-grained laminar or plaque-proximity analysis.

---

## Data

### Input
| File | Description |
|------|-------------|
| `seaad_paired_rna.h5ad` | Full paired RNA (138,118 cells × 36,601 genes) |
| `seaad_paired_atac.h5ad` | Full paired ATAC |
| `SEAAD_MTG_MERFISH.2024-12-11.h5ad` | Raw MERFISH (1.88M cells × 140+ genes) |
| `hvg_names_sub.csv` | 500 HVG list ordered by HVG rank |

### Intermediate
| File | Description |
|------|-------------|
| `astrocytes_seaad_paired_rna.h5ad` | Filtered RNA — astrocytes only (6,777 × 36,601) |
| `SpaGE_predicted_HVGs_FULL_50PV.h5ad` | SpaGE-imputed 466 HVGs (1.88M × 466) |
| `SpaGE_500_HVGs_FULL.h5ad` | Merged: 34 measured + 466 SpaGE HVGs (1.88M × 500) |

### ScReNI Outputs
| File | Description |
|------|-------------|
| `seaad_paired_sq1_triplets.csv` | GRN triplets (TF → peak → gene) |
| `seaad_paired_sq1_peak_gene_pairs.csv` | Peak-gene regulatory links |
| `seaad_paired_sq1_motif_peak_pairs.csv` | TF motif → ATAC peak associations |
| `seaad_paired_sq1_gene_labels.csv` | Per-cell gene activity labels |
| *(sub42 variants)* | Subsample-42 GRN outputs |

### Final Outputs
| File | Description |
|------|-------------|
| `Tangram_rna_with_coords_matched.h5ad` | RNA astrocytes + inferred coordinates in obs |
| `Tangram_cell_locations.csv` | `cell_id, tangram_x, tangram_y, tangram_prob_max` |

---

## Pipeline — Run Order

> All scripts are submitted via SLURM. Run from `/tudelft.net/staff-umbrella/ScReNI/dflam/`.

```bash
# 1. Select 500 HVGs from RNA data
sbatch 01_hvg_selection.sh

# 2. Subsample RNA to HVG genes only (for SpaGE input)
sbatch 02_subsample_rna_hvg_only.sh

# 3. Run SpaGE to impute 466 HVGs not measured in MERFISH (~12h)
sbatch 03_run_spage_full.sh

# 4. Merge SpaGE predictions with measured MERFISH genes → 500-gene spatial
sbatch 04_merge_spage_predicted.sh

# 5. Filter RNA to astrocytes only → 6,777 cells
sbatch 05_filter_astrocytes.sh

# 6b. Per-donor matched Tangram mapping (PRODUCTION — ~10 min)
#     100 genes, top 0.5% weighted centroid, saves CSV + h5ad
sbatch 06b_run_tangram_matched.sh

# 7b. Plot inferred locations per donor (runs after 06b via dependency)
sbatch --dependency=afterok:<06b_jobid> 07b_plot_inferred_locations.sh

# Optional: plot actual MERFISH spatial gene expression
sbatch 07_plot_gene_expression.sh
```

> **Note**: Step 06 (`06_run_tangram.sh`) maps all 6,777 RNA cells to all 1.88M spatial cells (cross-donor, 400 GB RAM, ~5–6h). This is provided for reference but **06b is the production pipeline** as it uses per-donor matched mapping which is geometrically valid.

---

## Validation Scripts (`debug/`)

Run these to reproduce accuracy benchmarks:

```bash
# Within-donor holdout validation — different gene sets
sbatch debug/06_1_validate_tangram.sh           # 500 genes
sbatch debug/06_1_validate_tangram_measured.sh  # 34 measured only
sbatch debug/06_1_validate_tangram_top100.sh    # 100 genes
sbatch debug/06_1_validate_tangram_merfish140.sh # 140 raw MERFISH

# Top-k coordinate threshold sweep
sbatch debug/06_1_validate_topk_sweep.sh        # 100 genes, 8 thresholds
sbatch debug/06_1_validate_topk_sweep_500.sh    # 500 genes, 8 thresholds

# Cross-donor validation (mirrors cross-donor mapping)
sbatch debug/06_2_validate_tangram_crossdonor.sh
```

---

## Environment

All steps run inside an Apptainer container:

```
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
```

Python packages managed by `pixi` at `/opt/app/pixi.toml`. Additional packages installed via `pip install --user` (torch, squidpy, zarr<3) accessible via `--bind $HOME:$HOME`.

Custom Tangram fork: `/tudelft.net/staff-umbrella/ScReNI/dflam/Tangram/`

---

## SLURM Resource Requirements

| Job | Partition | QOS | CPUs | Memory | Time |
|-----|-----------|-----|------|--------|------|
| 01–05 | general | short | 4–8 | 32–64 GB | <1h |
| 03 (SpaGE) | general | medium | 8 | 128 GB | ~13h |
| 06 (full Tangram) | general | medium | 8 | 400 GB | ~6h |
| **06b (production)** | general | short | 8 | 64 GB | ~15 min |
| 07b (plots) | general | short | 2 | 32 GB | ~5 min |
| Validation | general | short | 8 | 64 GB | 15–40 min |

---

## Key Findings

- **SpaGE-predicted genes hurt localization**: adding 466 SpaGE genes on top of 34 measured ones does not improve spatial accuracy — the predictions are too correlated with the reference spatial data to add independent signal
- **140 raw MERFISH genes give best localization** (16.8% within 500 µm) but can't be used in production since the RNA data lacks those non-HVG genes
- **100-gene set is optimal for production**: best balance between RNA–spatial gene overlap and spatial signal quality
- **Top 0.5% weighted centroid is optimal**: avoids noise of argmax (single cell) and the long-tail pull of full weighted centroid
- **Spatial accuracy is sufficient for broad regional analysis**: 2,412 µm median error on an 11–12 mm tissue = ~21% of tissue width, ~7× better than random

---

## Limitations

- Spatial coordinates are inferred at region-level precision (~2.4 mm), not cell-level — plaque-proximity and laminar analyses are not feasible
- Only 10 of 28 RNA donors have matched spatial data; the other 18 are mapped cross-donor (coordinate alignment not guaranteed)
- MERFISH panel covers only 140 genes vs 36,601 in RNA — the gene overlap bottleneck is a fundamental constraint regardless of mapping tool
- Tangram assumes conserved gene expression patterns across astrocyte subtypes between modalities

---

## References

- **Tangram**: Biancalani et al. (2021) *Nature Methods* — deep learning spatial cell mapping
- **SpaGE**: Abdelaal et al. (2020) *Nucleic Acids Research* — spatial gene expression prediction
- **ScReNI**: [ScReNI paper/tool] — per-cell gene regulatory network inference from multiome data
- **SEA-AD**: Seattle Alzheimer's Disease Brain Cell Atlas — `seaad.alleninstitute.org`
