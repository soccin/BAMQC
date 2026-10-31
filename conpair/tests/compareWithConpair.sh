#!/bin/bash
#SBATCH -J compareWithConpair
#SBATCH -o compareWithConpair.%j.out
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -t 00:30:00
#SBATCH --partition cpushort,cmobic_short
##
## Cross-check conpairCohort.py against Conpair's own scripts for one pair
## of samples from a finished run. Runs the container's
## verify_concordances.py (with and without --normal_homozygous_markers_only)
## and estimate_tumor_normal_contaminations.py with REFERENCE as the normal
## and QUERY as the tumor, and prints Conpair's values beside the ones in
## OUTDIR. They should be identical. Exits 1 if any differ.
##
## The paired contamination estimate is only compared when REFERENCE is
## QUERY's partner in OUTDIR/conpair_contamination.tsv. The markers are
## the ones the run used, from OUTDIR/conpair_run.tsv.
##
## Takes a few minutes (Conpair's contamination estimate). Submit it with
## sbatch; do not use srun from inside an interactive allocation, which
## runs the step on the interactive node.
##

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo -e "\n   usage: sbatch compareWithConpair.sh OUTDIR REFERENCE QUERY [CONTAINER]\n"
    exit 1
fi

OUTDIR=$(readlink -f "$1")
REF=$2
QUERY=$3
CONTAINER=${4:-${CONPAIR_SIF:-/home/soccin/Code/ConPair/conpair_v0.3.3.sif}}

REF_PILEUP=$OUTDIR/pileup/$REF.pileup
QUERY_PILEUP=$OUTDIR/pileup/$QUERY.pileup
for F in "$REF_PILEUP" "$QUERY_PILEUP" "$OUTDIR/conpair_concordance.tsv" \
         "$OUTDIR/conpair_contamination.tsv" "$OUTDIR/conpair_run.tsv"; do
    [ -f "$F" ] || { echo "ERROR: not found: $F" >&2; exit 1; }
done
MARKERS_TXT=$(awk -F'\t' '$1=="markers"{print $2}' "$OUTDIR/conpair_run.tsv")
[ -n "$MARKERS_TXT" ] || { echo "ERROR: no markers in $OUTDIR/conpair_run.tsv" >&2; exit 1; }

##
## Conpair's scripts write into the working directory and call Rscript for
## a plot; the empty .Rprofile keeps the user's own out of the container's
## R. Only files this script creates are removed.
##
WORK=$(mktemp -d "$OUTDIR/.compare.XXXXXX")
# shellcheck disable=SC2329  # invoked by the EXIT trap
cleanup() {
    rm -f "$WORK"/*.txt "$WORK"/*.R "$WORK"/*.pdf "$WORK"/.Rprofile
    rmdir "$WORK"
}
trap cleanup EXIT
touch "$WORK/.Rprofile"
printf '%s\t%s\n' "$REF" "$QUERY" > "$WORK/pairing.txt"

conpair() {
    singularity exec --cleanenv --bind "$OUTDIR" --pwd "$WORK" "$CONTAINER" "$@" \
        > /dev/null 2>&1 || true
}

conpair python2.7 /usr/bin/conpair/scripts/verify_concordances.py \
    --tumor_pileup="$QUERY_PILEUP" --normal_pileup="$REF_PILEUP" \
    --markers="$MARKERS_TXT" --pairing="$WORK/pairing.txt" \
    --normal_homozygous_markers_only --outpre=hom
conpair python2.7 /usr/bin/conpair/scripts/verify_concordances.py \
    --tumor_pileup="$QUERY_PILEUP" --normal_pileup="$REF_PILEUP" \
    --markers="$MARKERS_TXT" --pairing="$WORK/pairing.txt" --outpre=all
conpair python2.7 /usr/bin/conpair/scripts/estimate_tumor_normal_contaminations.py \
    --tumor_pileup="$QUERY_PILEUP" --normal_pileup="$REF_PILEUP" \
    --markers="$MARKERS_TXT" --pairing="$WORK/pairing.txt" --outpre=contam

for F in hom_concordance.txt all_concordance.txt contam_contamination.txt; do
    [ -s "$WORK/$F" ] || { echo "ERROR: Conpair did not write $F" >&2; exit 1; }
done

CONPAIR_HOM=$(awk -F'\t' 'NR==2{print $2}' "$WORK/hom_concordance.txt")
CONPAIR_ALL=$(awk -F'\t' 'NR==2{print $2}' "$WORK/all_concordance.txt")
CONPAIR_SINGLE=$(awk -F'\t' '$1=="N"{print $3}' "$WORK/contam_contamination.txt")
CONPAIR_PAIRED=$(awk -F'\t' '$1=="T"{print $3}' "$WORK/contam_contamination.txt")

OURS_HOM=$(awk -F'\t' -v r="$REF" -v q="$QUERY" '$1==r && $2==q {print $3}' "$OUTDIR/conpair_concordance.tsv")
OURS_ALL=$(awk -F'\t' -v r="$REF" -v q="$QUERY" '$1==r && $2==q {print $5}' "$OUTDIR/conpair_concordance.tsv")
OURS_SINGLE=$(awk -F'\t' -v r="$REF" '$1==r {print $2}' "$OUTDIR/conpair_contamination.tsv")
OURS_PARTNER=$(awk -F'\t' -v q="$QUERY" '$1==q {print $3}' "$OUTDIR/conpair_contamination.tsv")
if [ "$OURS_PARTNER" == "$REF" ]; then
    OURS_PAIRED=$(awk -F'\t' -v q="$QUERY" '$1==q {print $5}' "$OUTDIR/conpair_contamination.tsv")
else
    OURS_PAIRED="skip"
fi

STATUS=0
compare() {
    local label=$1 conpairValue=$2 oursValue=$3 result
    if [ "$oursValue" == "skip" ]; then
        result="not compared: partner of $QUERY is $OURS_PARTNER"
    elif awk -v a="$conpairValue" -v b="$oursValue" 'BEGIN{exit !(a + 0 == b + 0)}'; then
        result="same"
    else
        result="DIFFERENT"
        STATUS=1
    fi
    printf '  %-34s %10s %14s   %s\n' "$label" "$conpairValue" "$oursValue" "$result"
}

echo
echo "Reference (normal) $REF, query (tumor) $QUERY"
printf '  %-34s %10s %14s\n' "value" "Conpair" "conpairCohort"
compare "concordanceHom (ref -> query)" "$CONPAIR_HOM" "$OURS_HOM"
compare "concordanceAll" "$CONPAIR_ALL" "$OURS_ALL"
compare "contaminationSingle (ref)" "$CONPAIR_SINGLE" "$OURS_SINGLE"
compare "contaminationPaired (query | ref)" "$CONPAIR_PAIRED" "$OURS_PAIRED"
echo
exit $STATUS
