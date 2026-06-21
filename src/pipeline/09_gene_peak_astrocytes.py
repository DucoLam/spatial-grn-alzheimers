"""Phase 3 gene-peak-TF relationships for astrocytes only.

Runs run_phase3() on the 6,777 astrocyte cells (dflam/data/),
using reference files from bsc-screni. Outputs go to dflam/data/
with prefix "astrocytes".

Outputs:
    dflam/data/astrocytes_triplets.csv
    dflam/data/astrocytes_gene_labels.csv
    dflam/data/astrocytes_peak_overlap_matrix.npz
    dflam/data/astrocytes_peak_gene_pairs.csv
    dflam/data/astrocytes_peak_info.csv
    dflam/data/astrocytes_motif_peak_pairs.csv
"""

import sys
import logging
import time
from pathlib import Path

# screni lives in bsc-screni/src
sys.path.insert(0, "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/src")

import anndata as ad
from screni.data.gene_peak_relations import run_phase3, load_transfac_motifs
from screni.data.utils import load_gene_annotations

logging.basicConfig(level=logging.INFO, format="%(message)s")

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

# ── Paths ──────────────────────────────────────────────────────────────────
DFLAM     = Path("/tudelft.net/staff-umbrella/ScReNI/dflam")
BSC       = Path("/tudelft.net/staff-umbrella/ScReNI/bsc-screni")
REF_DIR   = BSC / "data/reference"
PAPER_REF = BSC / "data/paper/reference"

RNA_IN    = BSC / "data/processed/seaad/astrocytes_seaad_paired_rna.h5ad"
ATAC_IN   = DFLAM / "data/seaad_paired_atac_astrocytes.h5ad"
OUT_DIR   = DFLAM / "data"
PREFIX    = "astrocytes"

print("=" * 60)
print("09_gene_peak_astrocytes")
print(f"Input RNA:  {RNA_IN}")
print(f"Input ATAC: {ATAC_IN}")
print(f"Output dir: {OUT_DIR}")
print(f"Prefix:     {PREFIX}")
print("Start:", now())
print("=" * 60)

# ── Load reference files ───────────────────────────────────────────────────
print("\nLoading TRANSFAC motif data...")
pwm_dict, motif_db = load_transfac_motifs(
    PAPER_REF / "all_motif_pwm.rds",
    PAPER_REF / "Tranfac201803_Hs_MotifTFsFinal",
)

print("Loading gene annotations (hg38 Ensembl 98)...")
gene_ann = load_gene_annotations(REF_DIR / "hg38.ensembl98.gtf.gz")

# ── Load astrocyte data ────────────────────────────────────────────────────
print(f"\nLoading RNA:  {RNA_IN}")
rna = ad.read_h5ad(RNA_IN)
print(f"  Shape: {rna.shape}")

print(f"Loading ATAC: {ATAC_IN}")
atac = ad.read_h5ad(ATAC_IN)
print(f"  Shape: {atac.shape}")

# ── Run Phase 3 ────────────────────────────────────────────────────────────
t0 = time.time()
results = run_phase3(
    rna_adata      = rna,
    atac_adata     = atac,
    gene_annotations = gene_ann,
    genome_fasta   = REF_DIR / "hg38.fa",
    pwm_dict       = pwm_dict,
    motif_db       = motif_db,
    output_dir     = OUT_DIR,
    prefix         = PREFIX,
)
elapsed = time.time() - t0

# ── Summary ────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("SUMMARY")
print(f"  Cells:          {rna.n_obs:,}")
print(f"  Genes (HVGs):   {rna.n_vars:,}")
print(f"  Peaks (HVPs):   {atac.n_vars:,}")
print(f"  Triplets:       {len(results[triplets]):,}")
print(f"  Unique TFs:     {results[triplets][TF].nunique()}")
print(f"  Target genes:   {results[triplets][target_gene].nunique()}")
print(f"  Peaks kept:     {results[peak_info].shape[0]}")
print(f"  Peak matrix:    {results[peak_matrix].shape}")
print(f"  Runtime:        {elapsed:.1f}s ({elapsed/60:.1f} min)")
print(f"  Output prefix:  {OUT_DIR}/{PREFIX}_*")
print("End:", now())
print("=" * 60)
