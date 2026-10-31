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
# shellcheck source-path=SCRIPTDIR source=conpairMarkers.sh
. "$SDIR/conpairMarkers.sh"

if [ "$#" != "2" ]; then
    echo -e "\n   usage: sbatch conpairAnalyze.sh OUTDIR CONTAINER\n"
    exit 1
fi

OUTDIR=$1
CONTAINER=$2

## The markers follow from the builds in the manifest; runConpairCohort.sh
## has already checked that they all share one marker set.
MARKER_SETS=$(tail -n +2 "$OUTDIR/manifest.tsv" | cut -f4 | sort -u \
    | while read -r BUILD; do conpairMarkers "$BUILD" || echo "none:$BUILD"; done | sort -u)
if [ "$(wc -l <<< "$MARKER_SETS")" -ne 1 ] || [[ "$MARKER_SETS" == none:* ]]; then
    echo "ERROR: the manifest builds do not share one Conpair marker set:" >&2
    echo "$MARKER_SETS" >&2
    exit 1
fi
MARKERS_TXT=$MARKER_SETS.txt

command -v Rscript > /dev/null || { echo "ERROR: Rscript not on PATH" >&2; exit 1; }

echo "Start: conpairAnalyze $(date)"
echo "  output   $OUTDIR"
echo "  markers  $MARKERS_TXT"
echo "  bin      $SDIR"

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
