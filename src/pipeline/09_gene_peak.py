"""
09_gene_peak.py  — parameterized gene-peak-TF triplets
Usage:
    python 09_gene_peak.py --prefix astrocytes
    python 09_gene_peak.py --prefix l23_it
"""
import sys, argparse, logging, time
from pathlib import Path
import anndata as ad

parser = argparse.ArgumentParser()
parser.add_argument("--prefix", required=True, help="Cell type prefix (matches 08_subsample output)")
args = parser.parse_args()

sys.path.insert(0, "/tudelft.net/staff-umbrella/ScReNI/bsc-screni/src")
from screni.data.gene_peak_relations import run_phase3, load_transfac_motifs
from screni.data.utils import load_gene_annotations

logging.basicConfig(level=logging.INFO, format="%(levelname)-8s %(message)s")

DFLAM     = Path("/tudelft.net/staff-umbrella/ScReNI/dflam")
BSC       = Path("/tudelft.net/staff-umbrella/ScReNI/bsc-screni")
REF_DIR   = BSC / "data/reference"
PAPER_REF = BSC / "data/paper/reference"
RNA_IN    = DFLAM / f"data/seaad_paired_rna_{args.prefix}.h5ad"
ATAC_IN   = DFLAM / f"data/seaad_paired_atac_{args.prefix}.h5ad"
OUT_DIR   = DFLAM / "data"

print("=" * 60)
print(f"09_gene_peak  —  prefix={args.prefix}")
print(f"Input RNA:  {RNA_IN}")
print(f"Input ATAC: {ATAC_IN}")
print(f"Output dir: {OUT_DIR}")
print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)

for p in [RNA_IN, ATAC_IN]:
    if not p.exists():
        sys.exit(f"ERROR: {p} not found. Run 08_subsample.py --prefix {args.prefix} first.")

print("\nLoading TRANSFAC motif data...")
pwm_dict, motif_db = load_transfac_motifs(
    PAPER_REF / "all_motif_pwm.rds",
    PAPER_REF / "Tranfac201803_Hs_MotifTFsFinal",
)
print("Loading gene annotations (hg38 Ensembl 98)...")
gene_ann = load_gene_annotations(REF_DIR / "hg38.ensembl98.gtf.gz")

print(f"\nLoading RNA:  {RNA_IN}")
rna = ad.read_h5ad(RNA_IN)
print(f"  Shape: {rna.shape}")
print(f"Loading ATAC: {ATAC_IN}")
atac = ad.read_h5ad(ATAC_IN)
print(f"  Shape: {atac.shape}")

results = run_phase3(
    rna_adata       = rna,
    atac_adata      = atac,
    gene_annotations= gene_ann,
    genome_fasta    = REF_DIR / "hg38.fa",
    pwm_dict        = pwm_dict,
    motif_db        = motif_db,
    output_dir      = OUT_DIR,
    prefix          = args.prefix,
)

print("\n" + "=" * 60)
print(f"09_gene_peak complete  —  prefix={args.prefix}")
print(f"  Output prefix:  {OUT_DIR}/{args.prefix}_*")
print(f"End: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
