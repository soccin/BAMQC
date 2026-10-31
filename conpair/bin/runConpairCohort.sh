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
  -m, --manifest FILE   TSV or CSV with a header; a tab in the header
                        means TSV, else CSV. Columns:
                          sample      sample ID, [A-Za-z0-9._-] only, unique
                          bam         path to an indexed BAM
                          individual  optional; expected identity, e.g.
                                      cmoPatientId. Without it samples are
                                      grouped by genotype but not flagged.
                        Fields may be in double quotes.
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
## Read the manifest, TSV or CSV: a tab in the header means TSV, else a
## comma means CSV. Fields may be in double quotes ("" inside is a literal
## quote); surrounding spaces, a UTF-8 BOM, CR line ends and rows with no
## values are dropped. Prints one line per row: line number, parse error,
## sample, bam, individual. They are separated by US (0x1f), not tab,
## because bash read merges consecutive tabs and so loses empty fields.
## On a header problem prints the message and exits 1.
##
US=$'\x1f'
readManifest() {
    LC_ALL=C awk -v US="$US" -v q="'" '
        function trim(s) {
            gsub(/^[ \t]+|[ \t]+$/, "", s)
            return s
        }
        ## Split line on delim into f[1..n]. Returns n, or -1 when a
        ## double quote is not closed.
        function splitFields(line, f,    n, i, c, field, inQuote) {
            delete f
            n = 0; field = ""; inQuote = 0
            for (i = 1; i <= length(line); i++) {
                c = substr(line, i, 1)
                if (inQuote) {
                    if (c != "\"") field = field c
                    else if (substr(line, i + 1, 1) == "\"") { field = field c; i++ }
                    else inQuote = 0
                } else if (c == "\"") inQuote = 1
                else if (c == delim) { f[++n] = trim(field); field = "" }
                else field = field c
            }
            f[++n] = trim(field)
            return inQuote ? -1 : n
        }
        { sub(/\r$/, "") }
        NR == 1 {
            sub(/^\357\273\277/, "")
            if (index($0, "\t")) delim = "\t"
            else if (index($0, ",")) delim = ","
            else { print "manifest header has no tab or comma; need TSV or CSV"; exit 1 }
            n = splitFields($0, h)
            if (n < 0) { print "unclosed double quote in the manifest header"; exit 1 }
            for (i = n; i >= 1; i--) col[h[i]] = i
            if (!("sample" in col)) { print "manifest has no " q "sample" q " column"; exit 1 }
            if (!("bam" in col)) { print "manifest has no " q "bam" q " column"; exit 1 }
            next
        }
        {
            n = splitFields($0, f)
            if (n < 0) { print NR US "unclosed double quote" US US US; next }
            values = ""
            for (i = 1; i <= n; i++) values = values f[i]
            if (values == "") next
            s = f[col["sample"]]; b = f[col["bam"]]
            ind = ("individual" in col) ? f[col["individual"]] : ""
            if (index(s b ind, "\t")) { print NR US "tab inside a field" US US US; next }
            print NR US "" US s US b US ind
        }
        END { if (NR == 0) { print "manifest is empty"; exit 1 } }
    ' "$1"
}

##
## Check the manifest and write a normalized copy into the output
## directory: absolute BAM paths, always the columns sample, bam,
## individual, build, fasta. Every job reads that copy, never the original.
##
ROWS=$(readManifest "$MANIFEST") || die "$ROWS"

mkdir -p "$OUTDIR"
OUTDIR=$(readlink -f "$OUTDIR")
NEWMANIFEST=$(mktemp "$OUTDIR/.manifest.XXXXXX")
trap 'rm -f "$NEWMANIFEST"' EXIT

printf 'sample\tbam\tindividual\tbuild\tfasta\n' > "$NEWMANIFEST"
declare -A SEEN FASTA_OF MARKERS_OF BUILD_COUNT
ERRORS=0
while IFS=$US read -r LINENO_ PARSE_ERROR SAMPLE BAM INDIVIDUAL; do
    ## A manifest with no data rows gives one empty line.
    [ -n "$LINENO_" ] || continue
    if [ -n "$PARSE_ERROR" ]; then
        echo "    line $LINENO_: $PARSE_ERROR" >&2; ERRORS=$((ERRORS + 1)); continue
    fi
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
done <<< "$ROWS"

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
