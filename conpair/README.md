# conpair: cohort-wide sample identity and contamination

Conpair is normally run on one tumor and its matched normal. This tool runs
it across a whole cohort: every BAM is compared with every other BAM, so it
finds samples that do not match the individual they are labeled as, samples
that match an individual they should not, and contaminated samples.

All genotype and contamination numbers come from Conpair's own code, run
from the Tempo container. Nothing from Conpair is copied into this
repository: its NYGC license is single-user and academic only.

## Layout

```
bin/runConpairCohort.sh    entry point: checks the manifest, submits the jobs
bin/conpairMarkers.sh      genome build -> Conpair marker files (sourced)
bin/conpairPileup.sh       SLURM array task: GATK Pileup at the markers, one BAM
bin/conpairAnalyze.sh      SLURM job: conpairCohort.py, then conpairReport.R
bin/conpairCohort.py       concordance and contamination, in the container
bin/conpairReport.R        groups, flags, tables, PDF
bin/conpairVerify.sh       output check, run by the status job
tests/compareWithConpair.sh  checks one pair against Conpair's stock scripts
VERSION
```

## Usage

```
bin/runConpairCohort.sh [-c FILE] [-p LIST] [-o OUTDIR] MANIFEST
```

Run it on the login node; it only checks inputs and calls `sbatch`.
`OUTDIR` defaults to `./qc/conpair`. It must be on shared storage: an
`OUTDIR` under `/localscratch` or `/tmp` is rejected, because the jobs on
other nodes cannot see it.

### Knowing that a run is done and correct

The jobs are submitted through
[slurmTools](../../slurmTools/README.md): `$SLURMTOOLS_DIR`, by default the
`slurmTools` directory beside the BAMQC repository. Each run is recorded in
`OUTDIR/SLM/<RUN_ID>/`, with its logs. A last status job starts once every
job has ended, whatever their state. It checks each job's SLURM state, then,
if all completed, the outputs (`bin/conpairVerify.sh`):

- every sample in the manifest has a finished pileup (closing
  `[REDUCE RESULT]` line) and its fingerprint;
- every table has one row per sample, or per pair, for exactly the
  manifest's samples: N rows in the per-sample tables, N(N-1) in
  `conpair_concordance.tsv`, N(N-1)/2 in `conpair_pairs.tsv`;
- `conpair_run.tsv` has this version and N samples;
- the text report is not empty and the PDF is complete;
- every table and report was written by this run, not left by an earlier
  one. Pileups are exempt, since unchanged ones are reused.

The verdict, `OK` or `FAILED` with the reason, goes to
`OUTDIR/RUNSTATUS.txt`, and SLURM mails `$USER@mskcc.org` when the status
job ends (`SLURMTOOLS_MAIL` sets another address; empty for no mail). The
mail only says the run is done; read the verdict in `RUNSTATUS.txt`. To check
a run at any time (`runConpairCohort.sh` prints the full path):

```
$SLURMTOOLS_DIR/bin/checkRun.sh OUTDIR
```

Do not start a second run into an `OUTDIR` while one is still running: both
would write the same files.

`MANIFEST` is a TSV or CSV file with a header. A tab in the header line means
TSV, otherwise a comma means CSV; the file name does not matter. Fields may
be in double quotes (`""` inside them is a literal quote), as Excel writes
them. Spaces around fields, a UTF-8 BOM, CRLF line ends and rows with no
values are ignored. Columns, in any order:

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Sample ID, `[A-Za-z0-9._-]`, unique |
| `bam` | yes | Indexed BAM |
| `individual` | no | Expected identity, e.g. cmoPatientId |

Without `individual` the samples are still grouped by genotype, but nothing
can be flagged as a mismatch.

The genome is never given. Each BAM's build is read from its header with
`../bin/getGenomeBuildBAM.sh`, and the build selects the reference FASTA
(`../bin/getGenomeFasta.sh`) and the Conpair markers:

| Builds | Markers |
|---|---|
| b37, b37_dmp, b37_dmp_v2, GRCh37-lite | GRCh37 (contigs `1`, `2`, ...) |
| GRCh38, GRCh38_GDC/GRCh38.d1.vd1, hg38, hg38-igenomes | GRCh38 liftover (contigs `chr1`, `chr2`, ...) |

A BAM whose build is unknown, has no FASTA on IRIS, or has no Conpair
markers (hg19: `chr` names on GRCh37 coordinates) stops the run before
anything is submitted. BAMs on different builds can be run together when
they share the markers; each is piled up against its own FASTA.

The container defaults to `$CONPAIR_SIF`, else
`/home/soccin/Code/ConPair/conpair_v0.3.3.sif`. To pull it again:

```
singularity pull conpair_v0.3.3.sif docker://cmopipeline/conpair:v0.3.3
```

A pileup already in `OUTDIR/pileup` is reused when its BAM's fingerprint is
unchanged, so a rerun after adding samples only piles up the new ones. The
fingerprint is the BAM's size and the md5 of its first and last 4 MB. It
reads 8 MB rather than the whole BAM, and does not include the path, so a
moved BAM keeps its pileup, while a different BAM under the same sample
name is piled up again. The report alone can be rerun with other
thresholds:

```
Rscript bin/conpairReport.R OUTDIR [--related 80] [--contam-warn 2] [--contam-fail 5]
```

## Steps

1. **Pileup** (array job, one task per BAM, about 2 minutes for a 30-100x
   WGS BAM). Conpair's `run_gatk_pileup_for_sample.py` runs GATK 3.8 Pileup
   at the 7,387 Conpair markers. That script ignores GATK's exit status, so
   the task checks for GATK's closing `[REDUCE RESULT]` line itself.
2. **Genotypes and concordance** (seconds). Conpair's
   `genotype_likelihoods_for_markers` is run once per sample. Concordance is
   then computed for every ordered pair with the rule in
   `verify_concordance.py`: a marker counts when both samples have coverage
   of at least 10; the pair agrees there when the called genotypes are
   equal.
3. **Contamination** (about 2 minutes per sample at 30x, longer at higher
   depth, run in parallel on the job's CPUs). Conpair's
   `run_estimate_tumor_normal_contamination` is called directly.
4. **Report** (seconds, host R).
5. **Status** (seconds). The slurmTools status job: job states, then
   `conpairVerify.sh`; see above.

## Method

### Concordance and identity

Two concordance values are reported for every ordered pair (reference,
query):

- `concordanceHom` counts only markers called homozygous in the reference.
  This is Conpair's `--normal_homozygous_markers_only`, which Tempo uses.
  LOH and copy number changes in the query cannot make a homozygous marker
  heterozygous, so they do not lower this value, but it depends on which
  sample is the reference.
- `concordanceAll` counts all covered markers. It is symmetric, but LOH
  lowers it between samples of the same individual.

A pair's `identity` is the larger of its two `concordanceHom` values. LOH
in one sample lowers only the direction in which that sample is the
reference, so the larger value is the one LOH has not affected. Values
based on fewer than 100 markers are not used.

Expected values with these markers (minor allele frequency at least 0.4):
the same individual is near 100%, unrelated people are about 25-40%, and
relatives fall in between.

Samples whose identity is at least the match threshold (95%) are joined
into genotype groups (connected components).

### Contamination

Conpair estimates contamination in two ways, and both are reported:

- `contaminationSingle` is Conpair's normal-sample estimate: the sample's
  own reads against population allele frequencies. It assumes heterozygous
  sites are 50/50. Aneuploidy and LOH break that assumption and inflate the
  estimate, so for tumors and cell lines it overstates contamination.
- `contaminationPaired` is Conpair's tumor-sample estimate. It takes the
  markers homozygous in another sample from the same individual (the
  partner) as the sample's known genotype, so reads of the other allele at
  those markers come from contamination. Copy number changes do not affect
  this.

A sample's partner is the other sample with the highest `concordanceHom`
from partner to sample, provided it is at least the match threshold. The
verdict uses `contaminationPaired` when a partner exists and
`contaminationSingle` otherwise; the `contaminationMethod` column says
which.

When a sample has no partner, its single estimate comes from running
Conpair's function against an empty pileup, which skips the paired part.

The source of contamination is not identified. That is planned for v2.

## Flags and verdicts

Pair status (needs `individual` on both samples):

| Status | Condition | Level |
|---|---|---|
| `EXPECTED_MATCH_FAIL` | same individual, identity below 95% | FAIL |
| `UNEXPECTED_MATCH` | different individuals, identity at least 95% | FAIL |
| `PARTIAL_MATCH` | different individuals, identity 80-95% | WARN |
| `LOW_MARKERS` | fewer than 100 markers, no identity value | none |
| `OK` | none of the above | none |

Sample verdict is the worst of:

- its pair flags;
- contamination: at least 2% is WARN, at least 5% is FAIL (Tempo's
  cutoffs);
- fewer than 100 covered markers: WARN.

`INCOMPLETE` means no contamination estimate could be made. It is never
reported as `PASS`.

`PARTIAL_MATCH` is deliberately neutral. Identity between 80% and 95% can
come from a contaminated or mixed sample, low coverage, or close relatives.

`asymmetry` is the difference between the two directed concordances. It
only matters when identity is high: between unrelated samples a few points
of asymmetry is normal, because the two samples have different numbers of
homozygous markers. A large asymmetry in a matching pair means one of two
things. Either one sample has LOH: its own direction is
lower. Or one sample is a mixture: the direction starting from the mixture
is high, and the other direction is low. The contamination estimates tell
these apart.

## Outputs

Written into `OUTDIR`:

| File | Content |
|---|---|
| `manifest.tsv` | Normalized manifest used by all jobs, with each BAM's detected `build` and reference `fasta` |
| `pileup/<sample>.pileup` | GATK Pileup at the markers |
| `pileup/<sample>.fingerprint` | Fingerprint of the BAM the pileup came from: size and md5 |
| `conpair_samples.tsv` | Markers covered, % heterozygous, median coverage |
| `conpair_concordance.tsv` | Every ordered pair: both concordance values and marker counts |
| `conpair_contamination.tsv` | Single and paired estimate, partner |
| `conpair_run.tsv` | Version, parameters, run time |
| `conpair_summary.tsv` | One row per sample: verdict, group, contamination, flags |
| `conpair_pairs.tsv` | One row per pair: identity, both directions, status |
| `conpair_report.txt` | Text report |
| `conpair_report.pdf` | Concordance heatmaps (homozygous-marker and all-marker) and contamination plot |
| `RUNSTATUS.txt` | Verdict on the newest run: `OK`, `FAILED`, `RUNNING`, ... and why |
| `SLM/<RUN_ID>/` | One run: its record (`RUNINFO`, `jobs.tsv`), SLURM logs, `RUNSTATUS.txt` and the report behind it, `checkRun.txt` |
| `tmp/` | Java temporary space for the pileups; empty after a successful run |

All of these carry sample names. They stay in `OUTDIR`, outside the
repository.

## Limits

- Human only: the builds in the table under Usage. The container also has
  GRCm38 mouse exome markers, which this tool does not use.
- Markers are autosomal, so sex is not checked.
- The container's Conpair predates the upstream fix for very high coverage
  (above about 1000x at a marker, likelihoods underflow). WGS and exome
  depths are unaffected; GATK downsamples to 1000 reads.
- With many samples from one individual, every sample gets one partner,
  and the partner's single estimate is recomputed for each sample that
  chose it.

## Validation

| Cohort | n | Result |
|---|---|---|
| Proj_18079 cell lines, WGS 29-103x, GRCh37 | 12 | 6 genotype groups of 2. Five match their patient labels; one pair labeled as two different patients is genetically identical (99.93%), and both samples FAIL as `UNEXPECTED_MATCH`. No paired contamination estimate above 0.4%. |

Details from that run (2026-10-30):

- Same-individual `concordanceHom` 99.3-99.95%; unrelated pairs 26-42%.
- Run time: pileups about 2 minutes per BAM, in parallel; analysis 5.7
  minutes on 8 CPUs, almost all of it contamination (105-206 s per
  sample, longest for the 103x sample).
- The 103x sample has a single contamination estimate of 21.3% and a
  paired estimate of 0.36%. The paired estimate rules out DNA from another
  person, so the single estimate is the overstatement described under
  Contamination.
- `tests/compareWithConpair.sh` reproduces Conpair's stock scripts exactly
  for a same-individual pair (all four values) and an unrelated pair
  (concordance and single estimate). For the unrelated pair Conpair's
  paired estimate is 100%, which is why a paired estimate is only taken
  from a same-individual partner.

To check a pair from any run (it compares the paired estimate only when
REFERENCE is QUERY's partner):

```
sbatch -o compare.%j.out tests/compareWithConpair.sh OUTDIR REFERENCE QUERY
```
