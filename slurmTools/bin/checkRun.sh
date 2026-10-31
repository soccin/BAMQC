#!/bin/bash
##
## checkRun.sh [-w] DIR: did every job of a run succeed?
##
## DIR is an output directory (its newest run is checked) or one run
## directory, OUTDIR/SLM/<RUN_ID>. The run is OK when the launcher reached
## slurmFinish, every job in RUN_DIR/jobs ended COMPLETED in sacct, and the
## output check in RUN_DIR/verify, if any, exits 0.
##
## Prints a report ending in STATUS= and REASON= lines. With -w it also
## writes RUN_DIR/RUNSTATUS.txt and RUN_DIR/checkRun.txt, and copies
## RUNSTATUS.txt to OUTDIR when the run is the newest there.
##
## Exit: 0 OK, 1 FAILED, 2 RUNNING, 3 no run found or sacct failed.
##

set -uo pipefail

SDIR=$(dirname "$(readlink -f "$0")")
VERSION=$(cat "$SDIR/../VERSION")

WRITE=""
if [ "${1:-}" == "-w" ]; then
    WRITE=1
    shift
fi
if [ "$#" -ne 1 ]; then
    echo "usage: checkRun.sh [-w] DIR" >&2
    exit 3
fi

## Newest run of an output directory: run ids sort by time.
newestRun() {
    local f latest=""
    for f in "$1"/SLM/[0-9]*_[0-9]*/RUNINFO; do
        [ -f "$f" ] && latest=$(dirname "$(readlink -f "$f")")
    done
    echo "$latest"
}

if [ -f "$1/RUNINFO" ]; then
    RUN_DIR=$(readlink -f "$1")
else
    RUN_DIR=$(newestRun "$1")
fi
if [ -z "$RUN_DIR" ] || [ ! -f "$RUN_DIR/jobs" ]; then
    echo "ERROR: no run in $1" >&2
    exit 3
fi
OUTDIR=$(dirname "$(dirname "$RUN_DIR")")
info() {
    sed -n "s/^$1=//p" "$RUN_DIR/RUNINFO" | tail -n 1
}

##
## One line per task: JobID JobIDRaw JobName State StdOut, '|' separated.
## Only jobs since the run started: SLURM reuses job ids.
##
ACCT=""
IDS=$(paste -sd, "$RUN_DIR/jobs")
if [ -n "$IDS" ]; then
    START=$(info START)
    if ! ACCT=$(sacct -X -n -P -S "${START/ /T}" -j "$IDS" -o JobID,JobIDRaw,JobName,State,StdOut); then
        echo "ERROR: sacct failed" >&2
        exit 3
    fi
fi

## Counts "tasks done active failed", then one line per task that did not
## complete: job, name, state, log. A recorded job sacct does not know
## counts as failed.
TALLY=$(awk -F'|' '
    FILENAME == ARGV[1] { order[++n] = $1; next }
    NF < 5 { next }
    {
        split($1, a, "_"); seen[a[1]] = 1; tasks++
        st = $4; sub(/ .*/, "", st)
        if (st == "COMPLETED") { done++; next }
        if (st ~ /^(PENDING|RUNNING|REQUEUED|SUSPENDED|COMPLETING|CONFIGURING|RESIZING)$/) { active++; next }
        failed++
        lg = $5
        gsub(/%A/, a[1], lg); gsub(/%a/, a[2], lg); gsub(/%j/, $2, lg); gsub(/%x/, $3, lg)
        bad = bad sprintf("%s|%s|%s|%s\n", $1, $3, $4, lg)
    }
    END {
        for (i = 1; i <= n; i++)
            if (!(order[i] in seen)) {
                tasks++; failed++
                bad = bad sprintf("%s|?|NOT_IN_SACCT|\n", order[i])
            }
        printf "%d %d %d %d\n%s", tasks, done, active, failed, bad
    }
' "$RUN_DIR/jobs" <(printf '%s\n' "$ACCT"))
read -r TASKS DONE ACTIVE FAILED <<< "$(head -n 1 <<< "$TALLY")"
BAD=$(tail -n +2 <<< "$TALLY")

VERIFY_OUT=""
if [ -z "$(info SUBMIT_COMPLETE)" ] && [ "$ACTIVE" -eq 0 ]; then
    STATUS=FAILED
    REASON="the launcher stopped before slurmFinish; $TASKS tasks were submitted"
elif [ "$FAILED" -gt 0 ]; then
    STATUS=FAILED
    REASON="$FAILED of $TASKS tasks did not complete; first $(head -n 1 <<< "$BAD" | cut -d'|' -f1-3 | tr '|' ' ')"
elif [ "$ACTIVE" -gt 0 ]; then
    STATUS=RUNNING
    REASON="$ACTIVE of $TASKS tasks pending or running"
elif [ ! -s "$RUN_DIR/verify" ]; then
    STATUS=OK
    REASON="all $TASKS tasks completed; no output check"
else
    mapfile -t VERIFY < "$RUN_DIR/verify"
    if VERIFY_OUT=$("${VERIFY[@]}" 2>&1); then
        STATUS=OK
        REASON="all $TASKS tasks completed; output check passed"
    else
        STATUS=FAILED
        REASON="all $TASKS tasks completed, but the output check failed: $(head -n 1 <<< "$VERIFY_OUT")"
    fi
fi

report() {
    local job name state log
    echo "checkRun.sh $VERSION  $RUN_DIR"
    echo "  $TASKS tasks: $DONE completed, $ACTIVE pending or running, $FAILED not completed"
    if [ -n "$BAD" ]; then
        echo
        echo "  Not completed:"
        while IFS='|' read -r job name state log; do
            printf '    %-14s %-24s %s\n' "$job" "$name" "$state"
            if [ -f "$log" ]; then
                echo "      log: $log"
                grep -aiE 'error|fatal|exception|traceback|killed|oom' "$log" | tail -n 1 | cut -c 1-200 | sed 's/^/      /'
            fi
        done <<< "$BAD"
    fi
    if [ -n "$VERIFY_OUT" ]; then
        echo
        echo "  Output check:"
        head -n 40 <<< "$VERIFY_OUT" | sed 's/^/    | /'
    fi
    echo
    echo "STATUS=$STATUS"
    echo "REASON=$REASON"
}

REPORT=$(report)
echo "$REPORT"

if [ -n "$WRITE" ]; then
    echo "$REPORT" > "$RUN_DIR/checkRun.txt"
    printf 'STATUS=%s\nREASON=%s\nRUN_DIR=%s\nCHECKED=%s\nREPORT=%s\n' \
        "$STATUS" "$REASON" "$RUN_DIR" "$(date '+%F %T')" "$RUN_DIR/checkRun.txt" \
        > "$RUN_DIR/.RUNSTATUS.txt"
    mv -f "$RUN_DIR/.RUNSTATUS.txt" "$RUN_DIR/RUNSTATUS.txt"
    if [ "$(newestRun "$OUTDIR")" == "$RUN_DIR" ]; then
        cp "$RUN_DIR/RUNSTATUS.txt" "$OUTDIR/.RUNSTATUS.txt"
        mv -f "$OUTDIR/.RUNSTATUS.txt" "$OUTDIR/RUNSTATUS.txt"
    fi
fi

case "$STATUS" in
    OK)      exit 0 ;;
    FAILED)  exit 1 ;;
    *)       exit 2 ;;
esac
