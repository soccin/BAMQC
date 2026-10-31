#!/bin/bash
##
## runTests.sh ROOT: end-to-end tests with real SLURM jobs. Submits one run
## per scenario into ROOT/<scenario>, waits for them, then checks each
## verdict. Only submits and polls, so it can run on a login node. Takes a
## few minutes, depending on the queue. ROOT must be new and on shared
## storage. No mail is sent.
##

set -uo pipefail

SDIR=$(dirname "$(readlink -f "$0")")
TOOLS=$(dirname "$SDIR")
CHECKRUN=$TOOLS/bin/checkRun.sh

if [ "$#" -ne 1 ] || [ -e "$1" ]; then
    echo "usage: runTests.sh ROOT  (ROOT must not exist)" >&2
    exit 1
fi
ROOT=$(readlink -m "$1")
mkdir -p "$ROOT"
export SLURMTOOLS_MAIL=""
J=("--partition=cpushort,cmobic_short" --time=5 --mem=100M)

## Submits the jobs of scenario $1 into ROOT/$1, in a subshell as a
## launcher would. Prints the job ids to wait for.
launch() (
    local s=$1 out=$ROOT/$1 id
    # shellcheck source-path=SCRIPTDIR source=../slurmTools.sh
    . "$TOOLS/slurmTools.sh" "$out"
    case "$s" in
        ok)
            id=$(sbatch --parsable "${J[@]}" -a 1-3 -o "$out/a.%A_%a.out" --wrap 'exit 0')
            sbatch "${J[@]}" -d "afterok:$id" -o "$out/b.%j.out" --wrap 'sleep 60' > /dev/null
            slurmFinish true ;;
        fail)
            id=$(sbatch --parsable "${J[@]}" -o "$out/a.%j.out" \
                --wrap 'echo "ERROR: step a fails on purpose"; exit 1')
            sbatch "${J[@]}" -d "afterok:$id" -o "$out/b.%j.out" --wrap 'exit 0' > /dev/null
            slurmFinish true ;;
        failtask)
            # shellcheck disable=SC2016  # expanded in the job
            sbatch "${J[@]}" -a 1-3 -o "$out/a.%A_%a.out" \
                --wrap '[ "$SLURM_ARRAY_TASK_ID" -ne 2 ] || { echo "ERROR: task 2 fails on purpose"; exit 1; }' > /dev/null
            slurmFinish ;;
        verifyfail)
            sbatch "${J[@]}" -o "$out/a.%j.out" --wrap 'exit 0' > /dev/null
            slurmFinish sh -c 'echo "an output is missing on purpose"; exit 1' ;;
        nofinish)
            sbatch "${J[@]}" -o "$out/a.%j.out" --wrap 'exit 0' > /dev/null ;;
    esac
    paste -sd, "$RUN_DIR/jobs"
    echo "${RUN_STATUS_JOBID:-}"
)

PASSED=0
FAILED=0
## check DESCRIPTION GOT EXPECTED
check() {
    if [ "$2" == "$3" ]; then
        echo "PASS  $1"
        PASSED=$((PASSED + 1))
    else
        echo "FAIL  $1: got [$2], expected [$3]"
        FAILED=$((FAILED + 1))
    fi
}
status() {
    sed -n 's/^STATUS=//p' "$ROOT/$1/RUNSTATUS.txt"
}
checkRc() {
    "$CHECKRUN" "$ROOT/$1" > /dev/null 2>&1
    echo $?
}

SCENARIOS=(ok fail failtask verifyfail nofinish)
WAIT=""
for S in "${SCENARIOS[@]}"; do
    WAIT="$WAIT,$(launch "$S" | paste -sd,)"
done
WAIT=$(tr -s ',' <<< "${WAIT#,}" | sed 's/,$//')

check "checkRun on a run with jobs queued exits 2" "$(checkRc ok)" 2
check "OUTDIR/RUNSTATUS.txt is SUBMITTING before the status job" "$(status ok)" SUBMITTING

echo "waiting for jobs $WAIT"
for ((I = 0; I < 80; I++)); do
    [ -z "$(squeue -h -o %i -j "$WAIT" 2> /dev/null)" ] && break
    sleep 15
done

declare -A WANT=([ok]=OK [fail]=FAILED [failtask]=FAILED [verifyfail]=FAILED)
for S in "${!WANT[@]}"; do
    check "$S: RUNSTATUS.txt from the status job" "$(status "$S")" "${WANT[$S]}"
done
check "ok: checkRun exits 0" "$(checkRc ok)" 0
check "fail: checkRun exits 1" "$(checkRc fail)" 1
## Lines of the "Not completed" list: job id, name, state.
notCompleted() {
    grep -E '^    [0-9]+(_[0-9]+)? ' "$ROOT/$1/SLM/"*/checkRun.txt | awk '{ print $1, $3 }'
}
check "fail: the failed job and its cancelled dependent reported" \
    "$(notCompleted fail | awk '{ print $2 }' | sort | paste -sd' ')" "CANCELLED FAILED"
check "fail: report shows the error line" \
    "$(grep -c 'step a fails on purpose' "$ROOT/fail/SLM/"*/checkRun.txt)" 1
check "failtask: only task 2 reported" \
    "$(notCompleted failtask | sed -E 's/^[0-9]+_([0-9]+) /task \1 /')" "task 2 FAILED"
check "verifyfail: reason gives the check's message" \
    "$(grep -c 'an output is missing on purpose' "$ROOT/verifyfail/RUNSTATUS.txt")" 1
check "nofinish: checkRun exits 1" "$(checkRc nofinish)" 1

echo
echo "$PASSED passed, $FAILED failed"
[ "$FAILED" -eq 0 ]
