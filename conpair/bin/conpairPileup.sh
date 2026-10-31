#!/bin/bash
#SBATCH -J conpairPileup
#SBATCH -c 2
#SBATCH --mem 8G
#SBATCH -t 01:00:00
#SBATCH --partition cpushort,cmobic_short
##
## One SLURM array task: GATK 3 Pileup at the Conpair markers for row
## SLURM_ARRAY_TASK_ID of OUTDIR/manifest.tsv, run through Conpair's own
## run_gatk_pileup_for_sample.py in the container. The BAM's build and
## reference FASTA come from the manifest, written by runConpairCohort.sh,
## which submits this script.
##

set -euo pipefail

## Exit status line at the end of the log, for slurmTools' checkRun.sh.
## SLURMTOOLS_DIR is exported by runConpairCohort.sh.
if [ -n "${SLURMTOOLS_DIR:-}" ]; then
    # shellcheck source-path=SCRIPTDIR source=../../../slurmTools/lib/jobTrailer.sh
    . "$SLURMTOOLS_DIR/lib/jobTrailer.sh"
fi

if [ -n "${SBATCH_SCRIPT_DIR:-}" ]; then
    SDIR="${SBATCH_SCRIPT_DIR}"
else
    SDIR=$(dirname "$(readlink -f "$0")")
fi
# shellcheck source-path=SCRIPTDIR source=conpairMarkers.sh
. "$SDIR/conpairMarkers.sh"

if [ "$#" != "2" ]; then
    echo -e "\n   usage: sbatch --array 1-N conpairPileup.sh OUTDIR CONTAINER\n"
    exit 1
fi

OUTDIR=$1
CONTAINER=$2
TASK=${SLURM_ARRAY_TASK_ID:?must run as a SLURM array task}

ROW=$(awk -F'\t' -v n=$((TASK + 1)) 'NR==n' "$OUTDIR/manifest.tsv")
SAMPLE=$(cut -f1 <<< "$ROW")
BAM=$(cut -f2 <<< "$ROW")
BUILD=$(cut -f4 <<< "$ROW")
FASTA=$(cut -f5 <<< "$ROW")
[ -n "$SAMPLE" ] || { echo "ERROR: no manifest row for task $TASK" >&2; exit 1; }

MARKERS=$(conpairMarkers "$BUILD") \
    || { echo "ERROR: Conpair has no markers for build [$BUILD]" >&2; exit 1; }
MARKERS_BED=$MARKERS.bed

##
## Fingerprint of a BAM: its size and the md5 of its first and last 4 MB.
## The first 4 MB hold the header (read groups, @PG command lines), the
## last the final alignments and the EOF block. Reads 8 MB rather than
## the whole BAM, and does not depend on the path, so a moved BAM keeps
## its pileup. Not a content hash: it only tells one BAM from another.
##
bamFingerprint() {
    local size md5
    size=$(stat -L -c %s "$1")
    md5=$( { head -c 4M "$1"; tail -c 4M "$1"; } | md5sum | cut -d' ' -f1)
    echo "$size $md5"
}

PILEUP=$OUTDIR/pileup/$SAMPLE.pileup
FINGERPRINT_FILE=$OUTDIR/pileup/$SAMPLE.fingerprint
FINGERPRINT=$(bamFingerprint "$BAM")
echo "Start: conpairPileup $SAMPLE $(date)"
echo "  bam          $BAM"
echo "  fingerprint  $FINGERPRINT"
echo "  build        $BUILD"
echo "  fasta        $FASTA"
echo "  output       $PILEUP"
echo "  bin          $SDIR"

if [ -s "$PILEUP" ] && [ -f "$FINGERPRINT_FILE" ] \
        && [ "$(cat "$FINGERPRINT_FILE")" == "$FINGERPRINT" ]; then
    echo "Pileup exists and the BAM fingerprint is unchanged; reusing it."
    echo "End: conpairPileup $SAMPLE $(date)"
    exit 0
fi
## Removed before the pileup is replaced, so a task that dies part way
## leaves no fingerprint, and the next run piles up again.
rm -f "$FINGERPRINT_FILE"

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
echo "$FINGERPRINT" > "$FINGERPRINT_FILE"

echo "Sites: $(grep -vc '^\[REDUCE RESULT\]' "$PILEUP")"
echo "End: conpairPileup $SAMPLE $(date)"
