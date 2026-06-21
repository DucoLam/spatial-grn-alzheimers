# Spatially Annotated Gene Regulatory Networks in Alzheimer's Disease

**Author:** D. Lam — BSc Thesis, TU Delft
**Cluster:** DAIC HPC (`login.daic.tudelft.nl`, SLURM)
**Project (my) root:** `/tudelft.net/staff-umbrella/ScReNI/dflam/`
**Shared base repo:** `/tudelft.net/staff-umbrella/ScReNI/bsc-screni/` (read-only to us)

> Last updated mid-pipeline: spatial inference complete + validated; GRN
> inference (wScReNI) running; spatial analysis (RQ1) not yet started.

---

## Research question

**Main:** *Can changes in the spatial location of cells in human Alzheimer's disease
tissue be related to changes in gene regulatory network structure and the activity
of key transcription factors?*

**Sub-questions:**
- **SQ1 (methodological):** Can cell-specific GRNs be reliably inferred for cells in the selected AD samples using the ScReNI workflow?
- **SQ2:** Can differences in GRN structure between cells be related to differences in their spatial location?
- **SQ3:** Can spatial differences in GRN structure between cells be related to AD **severity** (Braak / neuropathological progression / APOE)?
- **SQ4:** Can spatial differences in GRN structure be related to the activity of specific TFs, and which other factors correlate with spatial differentiation when no significant structural differences are found?

Questions are framed **generally (cells)**; the analysis is run and reported **per
subtype** over the full dataset **{Astrocyte, L2/3 IT, L4 IT, Oligodendrocyte}**
(2 glia + 2 excitatory neurons) — subtypes appear in Methods/Results, not the questions.

**No control samples** — only AD across a **severity continuum** (Braak, AD
neuropathological progression, APOE). So SQ3 is a continuous/ordinal relationship,
and disease-stage varies at the **donor level** (effective N ≈ donors, not cells).

### SQ → pipeline mapping
| SQ | Needs | Status |
|---|---|---|
| SQ1 | ScReNI GRNs + a quality check (precision-recall vs ChIP atlas, `bsc-screni/evaluation.py`) | wScReNI running; quality check TODO |
| SQ2 | spatial link + GRN-vs-location (see `RQ1_grn_spatial_plan.md`) | not started |
| SQ3 | join Braak/progression/APOE to spatially-linked GRNs | not started |
| SQ4 | TF activity / cell-enriched regulators vs location | not started |

---

## What this project is (and isn't)

- **GRN inference = the established ScReNI pipeline** (`bsc-screni`, Xu et al.
  2025, GPB). We **run it unchanged** — never edit `bsc-screni/`; we only call
  its functions and write outputs to `dflam/`.
- **Spatial inference (Tangram) = our own additive layer.** bsc-screni has no
  spatial component; this is the thesis's novel contribution.
- The two branches are independent and meet only at the **spatial link** step.

```
Branch 1 (molecular → GRN):   RNA+ATAC → pool → KNN → gene_peak → wScReNI → per-cell GRNs
Branch 2 (molecular → space): RNA (140 shared genes) → Tangram → inferred (x,y)/section
                              ↓ join on cell_id (step 4) ↓
                         Spatially annotated per-cell GRNs → RQ1 analysis
```

---

## Data (SEA-AD, Middle Temporal Gyrus)

| Dataset | Content | Size | Role |
|---|---|---|---|
| Paired snMultiome (RNA+ATAC) | full transcriptome + chromatin, **no location** | 28 donors, 138,118 cells | GRN inference |
| MERFISH spatial | 140 genes + true (x,y), **no ATAC** | 1,887,729 cells, 27 donors, **69 sections** | spatial reference |
| Overlap (both modalities) | — | **10 donors** | enables mapping |

- 4-subtype paired-RNA counts: L2/3 IT 38,413 · L4 IT 19,054 · Oligodendrocyte 9,712 · Astrocyte 6,777 → **~73,956 pooled cells**.
- MERFISH is dominated by **oligodendrocytes (578,890 cells)** — why oligo validation is slow.
- **Each donor has 2–5 MERFISH sections**, each with its **own coordinate frame** (not co-registered). RNA samples ≠ MERFISH sections (0 ID overlap) — RNA is 1 dissociated library/donor.

---

## Pipeline (post-SpaGE; SpaGE is NOT used anymore)

All scripts live in `dflam/`. SLURM logs in `dflam/slurm/logs/`. Outputs in `dflam/data/`.

### Branch 2 — Spatial inference (Tangram), per subtype

| Script | Does | Run |
|---|---|---|
| `06_run_tangram.py` / `.sh` | per-subtype **matched** Tangram (140 native genes, **argmax** primary, section-aware), → coords | `sbatch --job-name=tg_<p> 06_run_tangram.sh '<CellType>' <prefix>` |
| `06_validate.py` / `.sh` | per-**section** 5% holdout validation; argmax + top-k sweep; saves `*_errors.npz` | `sbatch 06_validate.sh '<CellType>' <prefix>` |
| `06_plot_cdf.py` / `.sh` | CDF of error (% cells within distance), per subtype + per method | `sbatch 06_plot_cdf.sh` |
| `06_plot_maps.py` / `.sh` | tissue maps: grey = real MERFISH cells, colored = inferred cells by subtype, per section | `sbatch 06_plot_maps.sh` |

Outputs: `data/tangram/Tangram_cell_locations_<prefix>.csv` (cols: cell_id,
cell_subtype, **tangram_x/y = argmax**, tangram_x/y_centroid, tangram_prob_max,
tangram_mapped_section) + `..._with_coords_<prefix>.h5ad`;
`data/validation/<prefix>_validation_persection.csv` + `_errors.npz`;
`plots/cdf/`, `plots/maps/`.

### Branch 1 — GRN inference (ScReNI), POOLED across the 4 subtypes (paper-faithful)

| Script | Does | Run |
|---|---|---|
| `08_pool.py` / `.sh` | pool **all** cells of the 4 subtypes (no subsample), KNN **across the pooled set** (joint WNN embedding `obsm['X_pca']`, k=20) | `sbatch 08_pool.sh` |
| `09_gene_peak.py` / `09_gene_peak_pool.sh` | gene–peak Spearman (\|r\|>0.1) + **MOODS** motif scan → TF→peak→gene triplets | `sbatch 09_gene_peak_pool.sh pool4` |
| `10_wscreni.py` / `10_wscreni_shard.sh` | per-cell wScReNI GRNs; **sharded** via `cell_index` (16-job array) | `sbatch 10_wscreni_shard.sh pool4` |

Outputs: `data/seaad_paired_{rna,atac}_pool4.h5ad`, `data/pool4_triplets.csv`
(+ peak_overlap_matrix.npz, peak_info.csv, gene_labels.csv),
`data/wscreni_networks_pool4/wScReNI/<idx>.<cell>.network.txt` (one per cell).

*(Per-subtype `08_subsample.py` exists but is superseded by the pooled `08_pool.py`.)*

---

## Key methodological decisions

| Decision | Choice | Why |
|---|---|---|
| Spatial gene panel | **140 native MERFISH genes** | beat 100 (SpaGE) and 500 (HVG) on validation |
| Coordinate readout | **argmax** (single best cell) | beats top-k centroid on within-distance precision; no section ambiguity |
| Tangram matching | **per-subtype, donor-matched** | same-subtype reference; donor matching beats cross-donor (3,451 vs 2,357 µm) |
| Section handling | **section-aware** | each MERFISH section is its own frame; pooling them corrupts centroids |
| GRN feature selection | **500 HVG / 10k HVP** (from bsc-screni) | paper uses 2,000 HVG / 10k HVP — bsc-screni uses 500 HVG (inherited deviation, disclose) |
| GRN cells/KNN/correlation | **pooled across 4 subtypes, all cells** | ScReNI paper pools all cell types; KNN over whole set (neighbours are 99.4% same-subtype anyway) |
| Container | **`container_0-1-3.sif`** for gene_peak/wscreni | has MOODS compiled; the small `container.sif` does NOT (silent numpy fallback → timeouts) |

---

## Results so far

**Spatial accuracy (per-section validation, argmax):**

| Subtype | n | median | within 500 µm | within 1000 µm |
|---|---|---|---|---|
| L2/3 IT | 8,527 | 1,435 µm | 28.9% | 40.3% |
| L4 IT | 4,831 | 1,585 µm | 27.3% | 38.7% |
| Oligodendrocyte | 28,911 | 1,787 µm | 24.1% | 33.9% |
| Astrocyte | 10,541 | 1,947 µm | 22.5% | 31.0% |

(top-0.5% centroid is worse on within-distance everywhere; see `*_validation_persection.csv` for both.)

- **Neurons map better than glia** (layered expression → "self-forgiving": Tangram is most accurate for the most spatially-structured cells).
- **Section fix** cut medians ~400–800 µm vs the earlier (broken) per-donor pooling.
- **~24,461 cells** have both a location and (will have) a GRN — the analysis set (located ⊂ all 73,956 GRN'd cells; only the 10 overlapping donors get coords).
- Plots: `plots/cdf/` (CDF), `plots/maps/` (per-donor tissue maps with gridlines/ticks). Downloaded locally to `~/Downloads/screni_plots/`.

**GRN:** `08_pool` ✅ (73,956 cells), `gene_peak` ✅ (MOODS, **5,122 triplets**, 206 peaks). wScReNI **running** (16 shards).

---

## Caveats to disclose in the thesis

- **Resolution** ≈ cortical-depth-gradient (~1.4–1.9 mm median), **not** single-layer. Use depth-from-pia / coarse bins, not per-cell layer calls.
- **Centroid disclaimer:** the top-k centroid must be confined to one section; we pick that section by argmax (a hybrid). argmax used throughout.
- **Circularity (SQ2):** Tangram places cells by expression similarity → nearby cells are expression-similar by construction → GRNs (expression-derived) may look spatially structured artifactually. *(ScReNI itself assumes cells close in molecular space have similar GRNs.)* Must compare to an **expression-vs-space baseline** + biological validation.
- Inferred locations are **never** an input to GRN inference (the branches are independent).

---

## Next steps

1. Finish wScReNI (sharded, ETA hours) → GRN per cell for all 73,956.
2. **SQ1 quality check:** assess GRN reliability (precision-recall vs ChIP atlas).
3. **Step 4 — spatial link:** inner-join located cells (~24k) to their GRNs on `cell_id`.
4. **SQ2 analysis** (see `~/Documents/Personal/RQ1_grn_spatial_plan.md`): per-subtype GRN-dissimilarity vs cortical depth — cosine/Jaccard similarity, variogram + Mantel + Moran's I, GRN-clustering-vs-depth, with expression-baseline control. (Light compute: minutes after a one-time matrix build.)
5. **SQ3:** join Braak/progression/APOE → does the spatial-GRN signal scale with severity? (donor-level N).
6. **SQ4 (later):** cell-enriched regulators / TF activity vs location.

---

## Cluster cheat-sheet

```bash
sshpass -p '<pw>' ssh -o StrictHostKeyChecking=no dflam@login.daic.tudelft.nl "<cmd>"
squeue -u dflam                 # jobs
sacct -j <id> --format=JobID,State,Elapsed
sbatch <script.sh> [args]       # submit
sbatch --dependency=afterok:<id> ...   # chain
```
- Containers run via `apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net <SIF> pixi run --manifest-path /opt/app/pixi.toml python -u <script>`.
- QOS limits: short 4h, medium 36h, long 7d.
- **Never** write to `bsc-screni/` (shared). All our outputs → `dflam/`.
