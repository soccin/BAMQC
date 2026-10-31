# slurmTools

Tells you whether every SLURM job of a pipeline run succeeded.

SLURM gives no single answer to "did my run work?": an array task can fail
while the others succeed, the job waiting on it is cancelled, and the logs
look normal. slurmTools records the id of every job a launcher submits,
then a last status job asks SLURM accounting (`sacct`) how each one ended.

## Use

Two lines in a bash launcher:

```bash
. "$SLURMTOOLS_DIR/slurmTools.sh" "$OUTDIR"    # start a run; record every sbatch after this
# ... the launcher's sbatch calls, unchanged ...
slurmFinish "$SDIR/myVerify.sh" "$OUTDIR" "$RUN_DIR/START"   # output check optional
```

Sourcing `slurmTools.sh` defines a bash function `sbatch` that runs the real
`sbatch` and records the job id. `slurmFinish` submits a status job that
starts once every recorded job has ended (`afterany`), runs
`bin/checkRun.sh -w`, writes `OUTDIR/RUNSTATUS.txt`, and mails you. The mail
says `COMPLETED` for an OK run and `FAILED` otherwise.

A run is **OK** when the launcher reached `slurmFinish`, every task of every
job ended `COMPLETED` in `sacct`, and the output check, if given, exits 0.

To add slurmTools to a new or existing pipeline, see [HOWTO.md](HOWTO.md).

## Check a run

```
bin/checkRun.sh [-w] DIR
```

`DIR` is the output directory (its newest run) or a run directory. Exit 0 OK,
1 FAILED, 2 RUNNING, 3 no run found or `sacct` failed. The report lists every
task that did not complete, its log, and the last error-like line in it.

A conpair run given a file that is not a container:

```
checkRun.sh 0.2.0  /scratch/.../conpair/bad2/SLM/20261031_112312
  4 tasks: 0 completed, 0 pending or running, 4 not completed

  Not completed:
    5544200_1      conpairPileup            FAILED
      log: /scratch/.../conpair/bad2/SLM/20261031_112312/conpairPileup.5544200_1.out
      FATAL:   While checking container encryption: could not open image ...
    ...
    5544205        conpairAnalyze           CANCELLED

STATUS=FAILED
REASON=4 of 4 tasks did not complete; first 5544200_1 conpairPileup FAILED
```

## The output check

A job that exits 0 without doing its work ends `COMPLETED`, so give each
tool an output check: a script that exits 0 when the outputs are complete,
and otherwise prints one problem per line and exits 1. Pass it
`$RUN_DIR/START` and require outputs to be newer, so results left by an
earlier run cannot pass. Keep it light: it runs in the status job (1 CPU,
2 GB, 30 min) and on the login node when you run `checkRun.sh`.

## Run directory

```
OUTDIR/RUNSTATUS.txt     verdict of the newest run
OUTDIR/SLM/<RUN_ID>/     RUN_ID is the start time, YYYYmmdd_HHMMSS
  RUNINFO                start, launcher, SUBMIT_COMPLETE
  START                  empty; its time is the start of the run
  jobs                   one job id per line
  verify                 the output check, one argument per line
  checkRun.txt           report behind the verdict
  status.<jobid>.out     the status job's log
```

Job logs stay wherever the launcher puts them; `$RUN_DIR` is a good place.

## Notes

- Only `sbatch` calls made in the launcher's own shell are recorded. A
  separate script it runs, or a launcher in Python or R, does not see the
  function.
- Every job gets `--kill-on-invalid-dep=yes`: a job whose dependency failed
  is cancelled instead of pending forever. A job that waits `afternotok`
  for a failure is therefore cancelled when nothing fails, and the run
  shows FAILED; drop such jobs, the status job reports failures.
- A launcher that stops before `slurmFinish` gets no status job;
  `checkRun.sh` reports it FAILED once its jobs end.
- Settings: `SLURMTOOLS_MAIL` (default `$USER@mskcc.org`; empty for no
  mail), `SLURMTOOLS_PARTITION` (default `cpushort,cmobic_short`).
- Tests: `tests/runTests.sh ROOT` submits small real jobs into a new
  directory ROOT on shared storage and checks the verdicts.
