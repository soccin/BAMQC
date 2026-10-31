# HOWTO: add slurmTools to a pipeline

`README.md` says what slurmTools does. This file says how to add it to a
pipeline, new or old. A working example is BAMQC's
`conpair/bin/runConpairCohort.sh`.

slurmTools works with a bash launcher that submits its jobs with `sbatch`.
It does not see jobs submitted from Python or R.

## The pattern

```bash
#!/bin/bash
set -euo pipefail

SDIR=$(dirname "$(readlink -f "$0")")
SLURMTOOLS_DIR=${SLURMTOOLS_DIR:-$(dirname "$(dirname "$SDIR")")/slurmTools}

OUTDIR=$(readlink -m "$2")
# ... check the inputs; stop here on any problem ...

. "$SLURMTOOLS_DIR/slurmTools.sh" "$OUTDIR"    # from here on, every sbatch is recorded

JOB=$(sbatch --parsable --array "1-$N" -o "$RUN_DIR/stepA.%A_%a.out" "$SDIR/stepA.sh" "$OUTDIR")
sbatch --dependency "afterok:$JOB" -o "$RUN_DIR/stepB.%j.out" "$SDIR/stepB.sh" "$OUTDIR"

slurmFinish "$SDIR/verify.sh" "$OUTDIR" "$RUN_DIR/START"

echo "Status job $RUN_STATUS_JOBID. Check with: $SLURMTOOLS_DIR/bin/checkRun.sh $OUTDIR"
```

That is all a launcher needs: one line before the first `sbatch`, one after
the last. Here `SLURMTOOLS_DIR` defaults to a `slurmTools` checkout beside
the pipeline's repository. A pipeline that carries its own copy of
slurmTools (as BAMQC does, in `BAMQC/slurmTools/`) points it there instead.
Copy the whole tree except `.git`, and do not edit the copy: change
slurmTools itself, then copy it again.

## Converting an old pipeline

1. **Source `slurmTools.sh`** after the input checks and before the first
   `sbatch`. Sourcing starts a run, so a launcher that stops on a bad input
   should stop before this line. Without `set -e`, write
   `. "$SLURMTOOLS_DIR/slurmTools.sh" "$OUTDIR" || exit 1`.
2. **Call `slurmFinish`** after the last `sbatch`, with the output check
   and its arguments.
3. **Every `sbatch` must be a plain `sbatch` in the launcher's shell.**
   - Recorded: `sbatch` in the launcher, in its functions and loops, inside
     `$( )` and in pipes.
   - Not recorded: `command sbatch`, `/usr/bin/sbatch`, `xargs sbatch`,
     `parallel sbatch`. Change these to plain `sbatch` in a loop.
   - Not recorded: a helper script run as `./submit.sh`, because it is a new
     process. Source it instead (`. ./submit.sh`), unless it calls `exit`,
     or move its `sbatch` calls into the launcher.
4. **Remove jobs that run only on failure** (`--dependency=afternotok:...`).
   Every job gets `--kill-on-invalid-dep=yes`, so SLURM cancels such a job
   when nothing failed, and the cancel makes the run FAILED. The status job
   already reports failures.
5. **Remove `--mail-type` from the jobs.** The status job sends one mail per
   run.
6. **Make every job exit non-zero when it fails, and 0 when it succeeds.**
   Use `set -euo pipefail` in job scripts. Any job that does not end
   `COMPLETED` makes the run FAILED, including an optional step that is
   allowed to fail; such a step must exit 0.
7. **Logs** (optional): write them to `$RUN_DIR/NAME.%j.out`, or
   `NAME.%A_%a.out` for an array, so each run keeps its own. Use only `%j`,
   `%A`, `%a` and `%x` in log names; `checkRun.sh` fills these in to find
   the log of a failed job.
8. **Write an output check** (next section).
9. **Tell the user:** print the status job and the `checkRun.sh` command in
   the launcher's summary, mention `OUTDIR/RUNSTATUS.txt` in the README, and
   bump the version.

One run per output directory at a time: a second launch is not stopped, and
both would write the same files.

## The output check

A job can exit 0 without doing its work, so each pipeline needs a check of
its outputs. It is any command that exits 0 when the outputs are complete;
`checkRun.sh` runs it once every job has completed.

```bash
#!/bin/bash
## verify.sh OUTDIR SINCE: are the outputs of a run complete?
set -euo pipefail
OUTDIR=$1
SINCE=$2
PROBLEMS=0
problem() { echo "$*"; PROBLEMS=$((PROBLEMS + 1)); }

## Expected outputs follow from the inputs, not from what is there.
while read -r SAMPLE; do
    F=$OUTDIR/stepA/$SAMPLE.txt
    if [ ! -s "$F" ]; then
        problem "missing or empty: $F"
    elif [ ! "$F" -nt "$SINCE" ]; then
        problem "left by an earlier run: $F"
    fi
done < "$OUTDIR/samples.txt"

[ "$PROBLEMS" -eq 0 ] || exit 1
echo "complete: $OUTDIR"
```

- **Check that each expected file is there, written by this run, and
  complete.** "This run" means newer than `$RUN_DIR/START`. Outputs a
  rerun reuses on purpose are exempt; check their content instead.
  "Complete" can mean a closing line, a file footer, or one row per sample.
- **Print one problem per line, then exit 1.** The first line becomes the
  `REASON` in `RUNSTATUS.txt`.
- **Keep it light.** It runs in the status job (1 CPU, 2 GB, 30 min), and on
  the login node whenever someone runs `checkRun.sh`.
- **Use absolute paths.** It may run from any directory.

## Check that it works

Use a small input on shared storage, and turn mail off while testing with
`export SLURMTOOLS_MAIL=`.

1. **A normal run:** `checkRun.sh OUTDIR` exits 0. The task count in the
   report must match the jobs and array tasks you expect; a lower count
   means an `sbatch` was not recorded (step 3).
2. **A failing run:** give one sample a bad input, or `scancel` one task.
   The verdict is FAILED, and the report names the task, its log and its
   error line.
3. **The output check on its own:** on a copy of good outputs, delete one
   file, truncate another and age a third (`touch -d 2000-01-01 FILE`). The
   check must report each one.

## Using the verdict elsewhere

`checkRun.sh` exits 0 OK, 1 FAILED, 2 RUNNING. A pipeline can refuse to
start on input from a run that is not OK:

```bash
"$SLURMTOOLS_DIR/bin/checkRun.sh" "$UPSTREAM" > /dev/null \
    || { echo "ERROR: run in $UPSTREAM is not OK" >&2; exit 1; }
```
