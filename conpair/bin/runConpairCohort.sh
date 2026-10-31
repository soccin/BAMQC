#!/bin/bash
##
## Cohort-wide Conpair. Checks the manifest, reads each BAM's genome build
## from its header and picks the reference FASTA and Conpair markers from
## the build, then submits one SLURM array job for the per-BAM pileups and
## one job, dependent on it, for the concordance, contamination and report.
## Runs nothing heavy itself, so it is safe on a login node.
##

set -euo pipefail

SDIR=$(dirname "$(readlink -f "$0")")
REPO=$(dirname "$(dirname "$SDIR")")
VERSION=$(cat "$SDIR/../VERSION")
# shellcheck source-path=SCRIPTDIR source=conpairMarkers.sh
. "$SDIR/conpairMarkers.sh"

CONTAINER="${CONPAIR_SIF:-/home/soccin/Code/ConPair/conpair_v0.3.3.sif}"
PARTITION=""
MANIFEST=""
OUTDIR=""

usage() {
    cat <<EOF

runConpairCohort.sh $VERSION -- all-vs-all Conpair for a cohort of BAMs

Usage:
  runConpairCohort.sh --manifest FILE --out DIR [options]

Required:
  -m, --manifest FILE   TSV with a header. Columns:
                          sample      sample ID, [A-Za-z0-9._-] only, unique
                          bam         path to an indexed BAM
                          individual  optional; expected identity, e.g.
                                      cmoPatientId. Without it samples are
                                      grouped by genotype but not flagged.
  -o, --out DIR         Output directory, on shared storage.

Options:
  -c, --container FILE  Conpair container. Default: \$CONPAIR_SIF, else
                        /home/soccin/Code/ConPair/conpair_v0.3.3.sif
  -p, --partition LIST  SLURM partition(s) for both jobs. Default: the
                        #SBATCH value in each job script.
  -h, --help            Show this message and exit.

The genome build of each BAM is read from its header
(bin/getGenomeBuildBAM.sh); the build selects the reference FASTA
(bin/getGenomeFasta.sh) and the Conpair markers. BAMs on different builds
can be mixed only when they share the markers (e.g. b37 and GRCh37-lite).

Pileups already in DIR/pileup and newer than their BAM are reused, so a
rerun after adding samples only piles up the new ones.

EOF
    exit 1
}

ARGS=$(getopt -o 'm:o:c:p:h' \
    --long 'manifest:,out:,container:,partition:,help' \
    -n 'runConpairCohort.sh' -- "$@") || usage
eval set -- "$ARGS"
while true; do
    case "$1" in
        -m|--manifest)  MANIFEST="$2"; shift 2 ;;
        -o|--out)       OUTDIR="$2"; shift 2 ;;
        -c|--container) CONTAINER="$2"; shift 2 ;;
        -p|--partition) PARTITION="$2"; shift 2 ;;
        -h|--help)      usage ;;
        --)             shift; break ;;
        *)              usage ;;
    esac
done

die() {
    echo -e "\n    ERROR: $*\n" >&2
    exit 1
}

[ -n "$MANIFEST" ] || usage
[ -n "$OUTDIR" ] || usage
[ -f "$MANIFEST" ] || die "manifest not found: $MANIFEST"
[ -f "$CONTAINER" ] || die "container not found: $CONTAINER"
command -v samtools > /dev/null || die "samtools not on PATH (needed to read BAM headers)"
CONTAINER=$(readlink -f "$CONTAINER")

##
## Check the manifest and write a normalized copy into the output
## directory: absolute BAM paths, always the columns sample, bam,
## individual, build, fasta. Every job reads that copy, never the original.
##
colIndex() {
    head -1 "$MANIFEST" | tr -d '\r' \
        | awk -F'\t' -v name="$1" '{for(i=1;i<=NF;i++) if($i==name){print i; exit}}'
}
SCOL=$(colIndex sample)
BCOL=$(colIndex bam)
ICOL=$(colIndex individual)
[ -n "$SCOL" ] || die "manifest has no 'sample' column"
[ -n "$BCOL" ] || die "manifest has no 'bam' column"

mkdir -p "$OUTDIR"
OUTDIR=$(readlink -f "$OUTDIR")
NEWMANIFEST=$(mktemp "$OUTDIR/.manifest.XXXXXX")
trap 'rm -f "$NEWMANIFEST"' EXIT

printf 'sample\tbam\tindividual\tbuild\tfasta\n' > "$NEWMANIFEST"
declare -A SEEN FASTA_OF MARKERS_OF BUILD_COUNT
ERRORS=0
LINENO_=1
while IFS= read -r LINE || [ -n "$LINE" ]; do
    LINENO_=$((LINENO_ + 1))
    LINE=${LINE%$'\r'}
    [ -n "${LINE//[[:space:]]/}" ] || continue
    IFS=$'\t' read -r -a F <<< "$LINE"
    SAMPLE=${F[$((SCOL - 1))]:-}
    BAM=${F[$((BCOL - 1))]:-}
    INDIVIDUAL=""
    [ -n "$ICOL" ] && INDIVIDUAL=${F[$((ICOL - 1))]:-}

    if [[ ! "$SAMPLE" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "    line $LINENO_: bad sample ID [$SAMPLE]" >&2; ERRORS=$((ERRORS + 1)); continue
    fi
    if [ -n "${SEEN[$SAMPLE]:-}" ]; then
        echo "    line $LINENO_: duplicate sample ID [$SAMPLE]" >&2; ERRORS=$((ERRORS + 1)); continue
    fi
    SEEN[$SAMPLE]=1
    REALBAM=$(readlink -f "$BAM" || true)
    if [ ! -f "$REALBAM" ] || [[ "$REALBAM" != *.bam ]]; then
        echo "    line $LINENO_: BAM not found or not .bam [$BAM]" >&2; ERRORS=$((ERRORS + 1)); continue
    fi
    if [ ! -f "$REALBAM.bai" ] && [ ! -f "${REALBAM%.bam}.bai" ]; then
        echo "    line $LINENO_: no index beside [$REALBAM]" >&2; ERRORS=$((ERRORS + 1)); continue
    fi
    if ! BUILD=$("$REPO/bin/getGenomeBuildBAM.sh" "$REALBAM" 2>&1); then
        echo "    line $LINENO_: genome build not recognized [$BUILD] for [$REALBAM]" >&2; ERRORS=$((ERRORS + 1)); continue
    fi
    ## Looked up once per build; getGenomeFasta.sh reports its own errors.
    if [ -z "${FASTA_OF[$BUILD]+set}" ]; then
        FASTA_OF[$BUILD]=$("$REPO/bin/getGenomeFasta.sh" "$BUILD") || FASTA_OF[$BUILD]=""
        MARKERS_OF[$BUILD]=$(conpairMarkers "$BUILD") || MARKERS_OF[$BUILD]=""
    fi
    if [ -z "${FASTA_OF[$BUILD]}" ]; then
        echo "    line $LINENO_: no reference FASTA for build [$BUILD] of [$REALBAM]" >&2; ERRORS=$((ERRORS + 1)); continue
    fi
    if [ -z "${MARKERS_OF[$BUILD]}" ]; then
        echo "    line $LINENO_: Conpair has no markers for build [$BUILD] of [$REALBAM]" >&2; ERRORS=$((ERRORS + 1)); continue
    fi
    BUILD_COUNT[$BUILD]=$(( ${BUILD_COUNT[$BUILD]:-0} + 1 ))
    printf '%s\t%s\t%s\t%s\t%s\n' "$SAMPLE" "$REALBAM" "$INDIVIDUAL" "$BUILD" "${FASTA_OF[$BUILD]}" >> "$NEWMANIFEST"
done < <(tail -n +2 "$MANIFEST")

[ "$ERRORS" -eq 0 ] || die "$ERRORS manifest error(s), nothing submitted"
N=$(( $(wc -l < "$NEWMANIFEST") - 1 ))
[ "$N" -ge 2 ] || die "need at least 2 samples, found $N"

MARKER_SETS=$(for B in "${!BUILD_COUNT[@]}"; do echo "${MARKERS_OF[$B]}"; done | sort -u)
[ "$(wc -l <<< "$MARKER_SETS")" -eq 1 ] \
    || die "the BAMs are on builds with different Conpair markers (${!BUILD_COUNT[*]}); run them separately"

## Copied with a redirect, not mv: mktemp files are mode 600, and the
## output directory is usually shared with a group.
cat "$NEWMANIFEST" > "$OUTDIR/manifest.tsv"
mkdir -p "$OUTDIR/pileup" "$OUTDIR/SLM"

PARTITION_ARGS=()
[ -n "$PARTITION" ] && PARTITION_ARGS=(--partition "$PARTITION")

PILEUP_JOB=$(sbatch --parsable \
    --array "1-$N" \
    --output "$OUTDIR/SLM/conpairPileup.%A_%a.out" \
    --export "ALL,SBATCH_SCRIPT_DIR=$SDIR" \
    "${PARTITION_ARGS[@]}" \
    "$SDIR/conpairPileup.sh" "$OUTDIR" "$CONTAINER")

ANALYZE_JOB=$(sbatch --parsable \
    --dependency "afterok:$PILEUP_JOB" \
    --output "$OUTDIR/SLM/conpairAnalyze.%j.out" \
    --export "ALL,SBATCH_SCRIPT_DIR=$SDIR" \
    "${PARTITION_ARGS[@]}" \
    "$SDIR/conpairAnalyze.sh" "$OUTDIR" "$CONTAINER")

echo
echo "runConpairCohort.sh $VERSION"
echo "  samples    $N"
for B in "${!BUILD_COUNT[@]}"; do
    echo "  build      $B (${BUILD_COUNT[$B]} BAMs): ${FASTA_OF[$B]}"
done
cat <<EOF
  markers    $(basename "$MARKER_SETS")
  container  $CONTAINER
  output     $OUTDIR
  pileups    job $PILEUP_JOB (array 1-$N)
  analysis   job $ANALYZE_JOB (after the pileups succeed)
  logs       $OUTDIR/SLM/

EOF
