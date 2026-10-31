# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Cohort-wide QC tools for aligned BAMs, run on a SLURM cluster. Each tool
lives in its own directory (`conpair/` is the only one so far) with its own
`README.md`, `VERSION`, `bin/` and `tests/`. The tool READMEs are the
reference for method, flags and outputs; read them before changing
behavior, and update them in the same change.

Shared pieces at the top level:

- `bin/getGenomeBuildBAM.sh` identifies a build from the md5 of the sorted
  `@SQ` SN/LN lines. The md5 values depend on sorting under
  `en_US.UTF-8`; do not change the locale. New references are added to its
  `case` table.
- `bin/getGenomeFasta.sh` maps a build to a FASTA on IRIS and verifies the
  FASTA's `.dict` hashes back to the same build.
- `R/palettes.R` holds every color used in figures. Do not hard-code colors
  in a tool; use `PAL_QC`, `PAL_CAT`/`PAL_CAT32` via `catColors()`, and
  `NA_GREY`.

No tool takes a genome argument. The build always comes from the BAM header.

Every launcher uses slurmTools, a separate repository found in
`$SLURMTOOLS_DIR`, else `../slurmTools` beside this one. Sourcing its
`slurmTools.sh OUTDIR` records every later `sbatch` call (a bash function
wraps it), and `slurmFinish VERIFY...` adds the status job; see its
README. Call `sbatch` in the launcher's own shell, not from a separate
script, or the job is not recorded. A run is OK only when every task
ended `COMPLETED` in `sacct` and the tool's output check passes, so each
tool needs an output check that fails on missing, partial or stale outputs.

## conpair architecture

The pipeline is split across host and container, and across SLURM jobs:

1. `runConpairCohort.sh` (login node, light only): validates the manifest,
   detects each BAM's build, picks FASTA and markers, writes the normalized
   `OUTDIR/manifest.tsv` (`sample, bam, individual, build, fasta`), then
   submits the jobs through slurmTools. All later steps read that
   normalized manifest, never the user's.
2. `conpairPileup.sh`: SLURM array, one task per manifest row
   (`SLURM_ARRAY_TASK_ID` = row). Runs Conpair's GATK pileup script in the
   container. That script ignores GATK's exit status, so success is checked
   by the trailing `[REDUCE RESULT]` line. A pileup is reused when
   `pileup/<sample>.fingerprint` (BAM size plus md5 of its first and last
   4 MB) still matches the BAM, so moved BAMs keep their pileups.
3. `conpairAnalyze.sh`: depends `afterok` on the array. Runs
   `conpairCohort.py` inside the container, then `conpairReport.R` on the
   host.
4. The slurmTools status job (`afterany` on both) checks the job states,
   then runs `conpairVerify.sh OUTDIR RUN_DIR/START`, and writes
   `OUTDIR/RUNSTATUS.txt`. A change to the outputs (a new file, a changed
   row count) must be reflected in `conpairVerify.sh`.

`conpairMarkers.sh` is sourced (not executed) by all three shell scripts and
is the single place mapping build to Conpair marker files. All BAMs in a run
must resolve to one marker set.

Jobs locate their sibling scripts through `SBATCH_SCRIPT_DIR`, exported at
submit time, because SLURM runs a copy of the script from a spool
directory. Keep that pattern in any new job script.

Constraints that are easy to break:

- `conpairCohort.py` runs under the container's Python 2.7. Use type
  comments (`# type: (...) -> ...`), not annotations, and keep it
  2.7-compatible. It imports Conpair's modules from `/usr/bin/conpair` at
  run time.
- Never copy Conpair code into the repository. Its NYGC license is
  single-user, academic only, and forbids redistribution. Call its
  functions or scripts in the container instead.
- Tools write only into the output directory they are given, never into
  the repository. Outputs contain sample names.
- `conpairReport.R` exits 0 whenever it can read its inputs; it reports,
  it does not gate.
- Thresholds are shared: `conpairReport.R` takes the match and
  min-markers defaults from `conpair_run.tsv` so groups and partners agree
  with what `conpairCohort.py` computed.

## Commands

Run a cohort (from the login node; it only submits):

```
conpair/bin/runConpairCohort.sh [-c FILE] [-p LIST] [-o OUTDIR] MANIFEST
```

`OUTDIR` defaults to `./qc/conpair` and may not be under `/localscratch`
or `/tmp`.

Did the run finish, and did it succeed? Exit 0 OK, 1 FAILED, 2 RUNNING,
3 UNKNOWN; `OUTDIR/RUNSTATUS.txt` holds the status job's verdict:

```
../slurmTools/bin/checkRun.sh OUTDIR
```

Check the outputs of a run directly (exit 0 when complete):

```
conpair/bin/conpairVerify.sh OUTDIR [SINCE]
```

Rerun only the report with other thresholds:

```
Rscript conpair/bin/conpairReport.R DIR [--related 80] [--contam-warn 2] [--contam-fail 5]
```

Test: compare one pair from a finished run against Conpair's stock
scripts. Exits 1 if any value differs. Submit with `sbatch`; do not `srun`
from an interactive allocation.

```
sbatch -o compare.%j.out conpair/tests/compareWithConpair.sh DIR REFERENCE QUERY
```

Lint shell scripts (the scripts carry `shellcheck` directives and should
stay clean):

```
shellcheck -x bin/*.sh conpair/bin/*.sh conpair/tests/*.sh
```

The container defaults to `$CONPAIR_SIF`, else
`/home/soccin/Code/ConPair/conpair_v0.3.3.sif`. `.sif` files are not
committed.

## Conventions

- Conventional commits with the tool as scope:
  `feat(conpair): ...`, `docs(conpair): ...`, `test(conpair): ...`.
- A behavior change bumps the tool's `VERSION`, and the version appears in
  the commit subject (e.g. `feat(conpair): ..., v0.2.0`). The top-level
  `README.md` tool table lists the current version.
- Shell scripts: `set -euo pipefail`, a `##` header block explaining what
  the script does, `die`/`ERROR:` messages to stderr.
- Command lines follow the Unix convention: required inputs are
  positional arguments, options are optional and have defaults. Do not
  copy the required-named-option style of bioinformatics tools.
