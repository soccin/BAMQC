#!/bin/bash
#SBATCH -J conpairPileup
#SBATCH -c 2
#SBATCH --mem 8G
#SBATCH -t 01:00:00
#SBATCH --partition cpushort,cmobic_short
##
## One SLURM array task: GATK 3 Pileup at the Conpair markers for row
## SLURM_ARRAY_TASK_ID of OUTDIR/manifest.tsv, run through Conpair's own
## run_gatk_pileup_for_sample.py in the container. Submitted by
## runConpairCohort.sh.
##

set -euo pipefail

if [ -n "${SBATCH_SCRIPT_DIR:-}" ]; then
    SDIR="${SBATCH_SCRIPT_DIR}"
else
    SDIR=$(dirname "$(readlink -f "$0")")
fi

if [ "$#" != "4" ]; then
    echo -e "\n   usage: sbatch --array 1-N conpairPileup.sh OUTDIR GENOME FASTA CONTAINER\n"
    exit 1
fi

OUTDIR=$1
GENOME=$2
FASTA=$3
CONTAINER=$4
TASK=${SLURM_ARRAY_TASK_ID:?must run as a SLURM array task}

ROW=$(awk -F'\t' -v n=$((TASK + 1)) 'NR==n' "$OUTDIR/manifest.tsv")
SAMPLE=$(cut -f1 <<< "$ROW")
BAM=$(cut -f2 <<< "$ROW")
[ -n "$SAMPLE" ] || { echo "ERROR: no manifest row for task $TASK" >&2; exit 1; }

MARKERS_DIR=/usr/bin/conpair/data/markers
case "$GENOME" in
    GRCh37) MARKERS_BED=$MARKERS_DIR/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed ;;
    GRCh38) MARKERS_BED=$MARKERS_DIR/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed ;;
    *) echo "ERROR: unknown genome [$GENOME]" >&2; exit 1 ;;
esac

PILEUP=$OUTDIR/pileup/$SAMPLE.pileup
echo "Start: conpairPileup $SAMPLE $(date)"
echo "  bam     $BAM"
echo "  output  $PILEUP"
echo "  bin     $SDIR"

if [ -s "$PILEUP" ] && [ "$PILEUP" -nt "$BAM" ]; then
    echo "Pileup exists and is newer than the BAM; reusing it."
    echo "End: conpairPileup $SAMPLE $(date)"
    exit 0
fi

TMPOUT=$OUTDIR/pileup/.$SAMPLE.pileup.$SLURM_JOB_ID
JAVATMP=$OUTDIR/tmp/$SAMPLE.$SLURM_JOB_ID
mkdir -p "$JAVATMP"

singularity exec --cleanenv \
    --bind "$(dirname "$BAM")" \
    --bind "$(dirname "$FASTA")" \
    --bind "$OUTDIR" \
    --pwd "$OUTDIR" \
    "$CONTAINER" \
    python2.7 /usr/bin/conpair/scripts/run_gatk_pileup_for_sample.py \
        --gatk=/usr/bin/GenomeAnalysisTK.jar \
        --bam="$BAM" \
        --markers="$MARKERS_BED" \
        --reference="$FASTA" \
        --xmx_java=6g \
        --temp_dir_java="$JAVATMP" \
        --outfile="$TMPOUT"

rmdir --ignore-fail-on-non-empty "$JAVATMP"

##
## run_gatk_pileup_for_sample.py ignores GATK's exit status, so check the
## output: GATK writes a closing [REDUCE RESULT] line only on success.
##
if [ ! -s "$TMPOUT" ] || ! tail -1 "$TMPOUT" | grep -q '^\[REDUCE RESULT\]'; then
    echo "ERROR: GATK Pileup did not complete for $SAMPLE" >&2
    rm -f "$TMPOUT"
    exit 1
fi
mv "$TMPOUT" "$PILEUP"

echo "Sites: $(grep -vc '^\[REDUCE RESULT\]' "$PILEUP")"
echo "End: conpairPileup $SAMPLE $(date)"
