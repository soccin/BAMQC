# BAMQC

Tools for quality control of aligned BAM files, applied to a whole cohort
rather than one sample or one pair at a time.

## Tools

| Tool | Question it answers | Status |
|---|---|---|
| [`conpair/`](conpair/README.md) | Which samples come from the same individual, and how contaminated is each sample? | 0.1.0 |

## Layout

```
R/palettes.R     color palettes shared by every tool's figures
conpair/         cohort-wide Conpair: sample identity and contamination
```

Each tool lives in its own directory with its own `README.md`, `VERSION`,
`bin/` and `tests/`. Tools write only into the output directory they are
given; nothing is written into the repository.

## Requirements

- SLURM (`sbatch`), apptainer/singularity, R 4.5 with tidyverse, fs, glue,
  igraph, ComplexHeatmap and circlize.
- Tool-specific containers are listed in each tool's README. Containers are
  not stored in the repository.
