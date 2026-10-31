# shellcheck shell=bash
##
## slurmTools.sh: know whether every SLURM job of a run succeeded.
## Sourced by a bash launcher:
##
##   . slurmTools.sh OUTDIR          # start a run; record every sbatch after this
##   sbatch ...                      # the launcher's own calls, unchanged
##   slurmFinish [VERIFY [ARG...]]   # status job: checks the run when all jobs end
##
## The run is recorded in OUTDIR/SLM/<RUN_ID>/, set as RUN_DIR:
##
##   RUNINFO   start time, launcher, SUBMIT_COMPLETE once slurmFinish ran
##   START     empty; its time is the start of the run (for output checks)
##   jobs      the id of every job submitted, one per line
##   verify    the output check, one argument per line
##
## Every job gets --kill-on-invalid-dep=yes, so a job whose dependency
## failed is cancelled instead of pending forever. The status job runs
## bin/checkRun.sh -w, which writes OUTDIR/RUNSTATUS.txt, and SLURM mails
## when it ends: COMPLETED for an OK run, FAILED otherwise.
##
## Settings: SLURMTOOLS_MAIL (default $USER@mskcc.org; empty for no mail),
## SLURMTOOLS_PARTITION (default cpushort,cmobic_short).
##

SLURMTOOLS_DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
SLURMTOOLS_VERSION=$(cat "$SLURMTOOLS_DIR/VERSION")

if [ "$#" -ne 1 ]; then
    echo "ERROR: usage: . slurmTools.sh OUTDIR" >&2
    return 1
fi
RUN_ID=$(date +%Y%m%d_%H%M%S)
RUN_DIR=$(readlink -m "$1")/SLM/$RUN_ID
mkdir -p "$(dirname "$RUN_DIR")"
mkdir "$RUN_DIR" || return 1
touch "$RUN_DIR/START" "$RUN_DIR/jobs"
{
    echo "START=$(date '+%F %T')"
    echo "LAUNCHER=$(readlink -f "$0")"
    echo "SLURMTOOLS_VERSION=$SLURMTOOLS_VERSION"
} > "$RUN_DIR/RUNINFO"
echo "STATUS=SUBMITTING" > "$RUN_DIR/RUNSTATUS.txt"
cp "$RUN_DIR/RUNSTATUS.txt" "$RUN_DIR/../../RUNSTATUS.txt"

## Records the id of every job submitted. Output and exit status are those
## of sbatch.
sbatch() {
    local out id
    out=$(command sbatch --kill-on-invalid-dep=yes "$@") || return
    ## "123", "123;cluster" (--parsable) or "Submitted batch job 123"
    id=${out##* }
    id=${id%%;*}
    if [[ ! "$id" =~ ^[0-9]+$ ]]; then
        echo "ERROR: slurmTools: no job id in sbatch output [$out]" >&2
        return 1
    fi
    echo "$id" >> "$RUN_DIR/jobs"
    echo "$out"
}

## Submits the status job, after every recorded job (afterany). VERIFY is
## the tool's output check: the run is OK only if it exits 0. Sets
## RUN_STATUS_JOBID, and RUN_MAIL (empty for no mail).
slurmFinish() {
    local ids
    local -a opts
    ids=$(paste -sd: "$RUN_DIR/jobs")
    if [ -z "$ids" ]; then
        echo "ERROR: slurmFinish: no jobs were submitted" >&2
        return 1
    fi
    if [ "$#" -gt 0 ]; then
        printf '%s\n' "$@" > "$RUN_DIR/verify"
    fi
    RUN_MAIL=${SLURMTOOLS_MAIL-$USER@mskcc.org}
    opts=(--parsable --job-name="$(basename "$0" .sh)Status"
          --output="$RUN_DIR/status.%j.out" --dependency="afterany:$ids"
          --partition="${SLURMTOOLS_PARTITION:-cpushort,cmobic_short}"
          --cpus-per-task=1 --mem=2G --time=30)
    if [ -n "$RUN_MAIL" ]; then
        opts+=("--mail-type=END,FAIL" --mail-user="$RUN_MAIL")
    fi
    ## The sleep lets accounting (sacct) catch up with the jobs' end.
    # shellcheck disable=SC2034  # read by the launcher
    RUN_STATUS_JOBID=$(command sbatch "${opts[@]}" \
        --wrap="sleep 30; $(printf '%q ' "$SLURMTOOLS_DIR/bin/checkRun.sh" -w "$RUN_DIR")") \
        || return 1
    echo "SUBMIT_COMPLETE=$(date '+%F %T')" >> "$RUN_DIR/RUNINFO"
}
