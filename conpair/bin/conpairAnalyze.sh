#!/bin/bash
#SBATCH -J conpairAnalyze
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -t 02:00:00
#SBATCH --partition cpushort,cmobic_short
##
## Concordance and contamination for every sample in OUTDIR/manifest.tsv
## (conpairCohort.py, in the container), then the report (conpairReport.R,
## host R). Submitted by runConpairCohort.sh after the pileups finish; can
## also be resubmitted by hand.
##

set -euo pipefail

if [ -n "${SBATCH_SCRIPT_DIR:-}" ]; then
    SDIR="${SBATCH_SCRIPT_DIR}"
else
    SDIR=$(dirname "$(readlink -f "$0")")
fi
REPO=$(dirname "$(dirname "$SDIR")")

if [ "$#" != "3" ]; then
    echo -e "\n   usage: sbatch conpairAnalyze.sh OUTDIR GENOME CONTAINER\n"
    exit 1
fi

OUTDIR=$1
GENOME=$2
CONTAINER=$3

MARKERS_DIR=/usr/bin/conpair/data/markers
case "$GENOME" in
    GRCh37) MARKERS_TXT=$MARKERS_DIR/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt ;;
    GRCh38) MARKERS_TXT=$MARKERS_DIR/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt ;;
    *) echo "ERROR: unknown genome [$GENOME]" >&2; exit 1 ;;
esac

command -v Rscript > /dev/null || { echo "ERROR: Rscript not on PATH" >&2; exit 1; }

echo "Start: conpairAnalyze $(date)"
echo "  output  $OUTDIR"
echo "  bin     $SDIR"

singularity exec --cleanenv \
    --bind "$OUTDIR" \
    --bind "$REPO" \
    --pwd "$OUTDIR" \
    "$CONTAINER" \
    python2.7 "$SDIR/conpairCohort.py" \
        --manifest "$OUTDIR/manifest.tsv" \
        --pileup-dir "$OUTDIR/pileup" \
        --markers "$MARKERS_TXT" \
        --out "$OUTDIR" \
        --threads "${SLURM_CPUS_PER_TASK:-1}"

Rscript "$SDIR/conpairReport.R" "$OUTDIR"

echo "End: conpairAnalyze $(date)"
