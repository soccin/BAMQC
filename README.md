# BAMQC

Tools for quality control of aligned BAM files, applied to a whole cohort
rather than one sample or one pair at a time.

## Tools

| Tool | Question it answers | Status |
|---|---|---|
| [`conpair/`](conpair/README.md) | Which samples come from the same individual, and how contaminated is each sample? | 0.2.0 |

## Layout

```
bin/getGenomeBuildBAM.sh  genome build of a BAM, from its @SQ header lines
bin/getGenomeFasta.sh     reference FASTA on IRIS for a build
R/palettes.R              color palettes shared by every tool's figures
conpair/                  cohort-wide Conpair: sample identity and contamination
```

Each tool lives in its own directory with its own `README.md`, `VERSION`,
`bin/` and `tests/`. Tools write only into the output directory they are
given; nothing is written into the repository.

## Genome build

No tool asks for the genome. The build of each BAM is read from its header
by `bin/getGenomeBuildBAM.sh`: the md5 of its sorted `@SQ` SN/LN lines,
looked up in a table of known references. `bin/getGenomeFasta.sh` then
gives the reference FASTA for that build and checks that the FASTA's
`.dict` hashes to the same build. A BAM on a reference that is not in the
table is an error; add its md5 (printed in the error) to the table.

## Requirements

- SLURM (`sbatch`), apptainer/singularity, samtools, R 4.5 with
  tidyverse, fs, glue, igraph, ComplexHeatmap and circlize.
- Tool-specific containers are listed in each tool's README. Containers are
  not stored in the repository.
