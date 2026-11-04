#!/bin/bash
##
## Output check of a conpair run. slurmTools' checkRun.sh runs it once
## every job of the run has completed; the run is OK only if it exits 0.
## Checks that OUTDIR holds a complete result for exactly the samples in
## OUTDIR/manifest.tsv: a finished pileup per sample, one row per sample
## or per pair in every table, and the report. Lists each problem and
## exits 1 if there is any.
##
## With SINCE (a file; runConpairCohort.sh gives the run's START file) the
## tables and the report must also be newer than SINCE, so results left
## by an earlier run cannot pass. Pileups are exempt: a pileup whose BAM is
## unchanged is reused from an earlier run.
##

set -euo pipefail

SDIR=$(dirname "$(readlink -f "$0")")
VERSION=$(cat "$SDIR/../VERSION")

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo -e "\n   usage: conpairVerify.sh OUTDIR [SINCE]\n" >&2
    exit 1
fi
OUTDIR=$1
SINCE=${2:-}
if [ -n "$SINCE" ] && [ ! -e "$SINCE" ]; then
    echo "SINCE file not found: $SINCE"
    exit 1
fi

PROBLEMS=0
problem() {
    echo "$*"
    PROBLEMS=$((PROBLEMS + 1))
}

MANIFEST=$OUTDIR/manifest.tsv
if [ ! -s "$MANIFEST" ]; then
    echo "manifest missing or empty: $MANIFEST"
    exit 1
fi
SAMPLES=$(tail -n +2 "$MANIFEST" | cut -f1 | sort)
N=$(grep -c . <<< "$SAMPLES" || true)
if [ "$N" -lt 2 ]; then
    echo "manifest has $N samples, at least 2 expected: $MANIFEST"
    exit 1
fi

##
## Pileups. conpairPileup.sh writes the fingerprint only after a pileup
## has succeeded.
##
while read -r S; do
    P=$OUTDIR/pileup/$S.pileup
    if [ ! -s "$P" ]; then
        problem "pileup missing or empty: $P"
    elif ! grep -q '^\[REDUCE RESULT\]' < <(tail -n 1 "$P"); then
        problem "pileup incomplete, no closing [REDUCE RESULT] line: $P"
    fi
    if [ ! -s "$OUTDIR/pileup/$S.fingerprint" ]; then
        problem "fingerprint missing: $OUTDIR/pileup/$S.fingerprint"
    fi
done <<< "$SAMPLES"

## An output file is present, not empty and, with SINCE, newer than it.
## Returns 1 when the file is missing or empty.
output() {
    local file=$OUTDIR/$1
    if [ ! -s "$file" ]; then
        problem "missing or empty: $file"
        return 1
    fi
    if [ -n "$SINCE" ] && [ ! "$file" -nt "$SINCE" ]; then
        problem "not written by this run, older than $SINCE: $file"
    fi
    return 0
}

## Sorted values of the column named COL in a CSV; fails without it.
## FPAT keeps a quoted field with commas in it as one field.
columnValues() {
    awk -v FPAT='([^,]*)|("[^"]*")' -v c="$2" '
        NR == 1 {
            for (i = 1; i <= NF; i++) if ($i == c) k = i
            if (!k) exit 1
            next
        }
        { print $k }
    ' "$1" | sort
}

## A per-sample table: column COL lists every manifest sample once.
perSample() {
    local file=$1 col=$2 values
    output "$file" || return 0
    if ! values=$(columnValues "$OUTDIR/$file" "$col"); then
        problem "$file: no column $col"
        return 0
    fi
    if [ "$values" != "$SAMPLES" ]; then
        problem "$file: column $col does not list each of the $N manifest samples once ($(grep -c . <<< "$values" || true) rows)"
    fi
}

## A pair table: columns A and B hold every pair of distinct manifest
## samples once, as ordered pairs when ORDERED is 1, else unordered.
pairTable() {
    local file=$1 a=$2 b=$3 ordered=$4 expect got
    output "$file" || return 0
    if [ "$ordered" == 1 ]; then
        expect=$((N * (N - 1)))
    else
        expect=$((N * (N - 1) / 2))
    fi
    got=$(awk -v FPAT='([^,]*)|("[^"]*")' -v a="$a" -v b="$b" -v ordered="$ordered" '
        FNR == NR { s[$1] = 1; next }
        FNR == 1 {
            for (i = 1; i <= NF; i++) { if ($i == a) ka = i; if ($i == b) kb = i }
            if (!ka || !kb) { print "no columns " a " and " b; exit }
            next
        }
        {
            x = $ka; y = $kb
            if (!(x in s) || !(y in s) || x == y) { other++; next }
            key = (ordered == 1 || x < y) ? x SUBSEP y : y SUBSEP x
            if (key in seen) { other++; next }
            seen[key] = 1; n++
        }
        END { if (ka && kb) print n + 0 " pairs, " other + 0 " other rows" }
    ' <(echo "$SAMPLES") "$OUTDIR/$file")
    if [ "$got" != "$expect pairs, 0 other rows" ]; then
        problem "$file: $got; expected $expect pairs and no other rows"
    fi
}

perSample conpair_samples.csv sample
pairTable conpair_concordance.csv reference query 1
perSample conpair_contamination.csv sample
perSample conpair_summary.csv sample
pairTable conpair_pairs.csv sampleA sampleB 0

if output conpair_run.csv; then
    RUN_VERSION=$(awk -F, '$1 == "version" { print $2 }' "$OUTDIR/conpair_run.csv")
    RUN_N=$(awk -F, '$1 == "nSamples" { print $2 }' "$OUTDIR/conpair_run.csv")
    [ "$RUN_VERSION" == "$VERSION" ] \
        || problem "conpair_run.csv: version [$RUN_VERSION], this conpair is $VERSION"
    [ "$RUN_N" == "$N" ] \
        || problem "conpair_run.csv: nSamples [$RUN_N], the manifest has $N"
fi

output conpair_report.txt || true
## An .xlsx is a zip archive; an incomplete one has no end-of-archive record.
if output conpair_tables.xlsx; then
    XLSX=$OUTDIR/conpair_tables.xlsx
    [ "$(head -c 2 "$XLSX")" == "PK" ] || problem "not an xlsx (zip) file: $XLSX"
    grep -aq $'PK\x05\x06' < <(tail -c 1024 "$XLSX") || problem "xlsx incomplete, no zip end record: $XLSX"
fi
if output conpair_report.pdf; then
    PDF=$OUTDIR/conpair_report.pdf
    [ "$(head -c 5 "$PDF")" == "%PDF-" ] || problem "not a PDF: $PDF"
    grep -aq '%%EOF' < <(tail -c 1024 "$PDF") || problem "PDF incomplete, no %%EOF at its end: $PDF"
fi

if [ "$PROBLEMS" -gt 0 ]; then
    echo "$PROBLEMS problem(s) in $OUTDIR"
    exit 1
fi
echo "complete: $N samples, $((N * (N - 1) / 2)) pairs, in $OUTDIR"
