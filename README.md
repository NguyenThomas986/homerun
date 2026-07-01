# homerun

<div align="center">

[![Issues](https://img.shields.io/github/issues/NguyenThomas986/homerun.svg)](https://github.com/NguyenThomas986/homerun/issues)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

</div>

This repository is the SLURM-oriented wrapper around the Python package in [csrnaseq](csrnaseq). It is meant to help you launch a complete csRNA-seq analysis workflow from a project directory on a cluster.


## What this repo contains

- [submit_array.sh](submit_array.sh): submits the full three-phase workflow
- [prepare.sbatch](prepare.sbatch): prepares the project and runs the initial setup
- [align_array.sbatch](align_array.sbatch): runs trim + align for one sample per array task
- [collect.sbatch](collect.sbatch): runs downstream analyses after the array finishes
- [csrnaseq](csrnaseq): the Python package that actually performs the pipeline steps

## Typical workflow

1. Create a project directory and place your FASTQ files in it.
2. Run the submission script with the required paths and cluster settings.
3. Let the prepare, align-array, and collect stages run in sequence.
4. Inspect the outputs in the project tree, especially the report and TSS directories.

## Quick start (w/o copying files and having fastqs in the project dir)

```bash
path/to/homerun/submit_array.sh \
  --project /PATH/TO/PROJECT \
  --partition kamiak \
  --conda-env CONDA_ENV_NAME \
  --genome-index /path/to/STARIndex \
  --genome hg38
```

## Quick start (with copying files via copy-src)

```bash
path/to/homerun/submit_array.sh \
  --project /PATH/TO/PROJECT \
  --partition kamiak \
  --conda-env CONDA_ENV_NAME \
  --genome-index /path/to/STARIndex \
  --genome hg38 \
  --copy-src /PATH/TO/*.FASTQ.GZ
```

The controller script will:

- stage raw FASTQs into the project
- submit a prepare job
- submit a sample-array job for trim + align
- submit a collect job that runs tagdirs, bedGraphs, TSS, QC, stability, and report

## All flags

Required:

| Flag | Description |
|---|---|
| `--project` | Path to your project directory |
| `--partition` | SLURM partition name |
| `--conda-env` | Conda environment with pipeline dependencies |
| `--genome-index` | Path to STAR `--genomeDir` or HISAT2 `-x` index prefix |
| `--genome` | HOMER genome name (e.g. `hg38`) |

Optional:

| Flag | Default | Description |
|---|---|---|
| `--conda-module` | `anaconda3` | Cluster module that provides Conda |
| `--aligner` | `star` | Aligner to use: `star` or `hisat2` |
| `--throttle` | `16` | Max array tasks running at once |
| `--email` | | Email address for SLURM notifications |
| `--copy-src` | | Glob path to FASTQs to copy into `RawData/` if none are present |
| `--starindex-url` | | URL to auto-download a STAR index tarball if none is found |
| `--ntag-threshold` | `7` | Minimum tags to call a TSS cluster |
| `--trim-min` | `20` | Discard reads shorter than this after trimming |
| `--trim-max` | `58` | Discard reads longer than this after trimming |
| `--trim-adapter` | `AGATCGGAAGAGCACACGTCT` | 3′ adapter sequence |
| `--skip-chr` | `chrEBV` | Chromosome to exclude from bedGraphs |

To pass extra flags directly to the Python pipeline (e.g. STAR tuning):

```bash
path/to/homerun/submit_array.sh \
  --project /PATH/TO/PROJECT \
  --partition kamiak \
  --conda-env CONDA_ENV_NAME \
  --genome-index /path/to/STARIndex \
  --genome hg38 \
  -- --star-filter-multimap 5000
```

To see all the flags do:

```bash
path/to/homerun/submit_array.sh \
    --h (or --help)
```

## Requirements

The conda environment passed via `--conda-env` must provide:

| Tool | Used for | Notes |
|---|---|---|
| [HOMER](http://homer.ucsd.edu/homer/) | Tag directories, bedGraphs, TSS calling/annotation | Must include `homerTools` on `PATH`; genome data for `--genome` (e.g. `hg38`) must be installed via `perl configureHomer.pl -install hg38` |
| Python ≥3.9 | Runs the `csrnaseq` package | See [csrnaseq/README.md](csrnaseq/README.md) for Python package dependencies |
| [STAR](https://github.com/alexdobin/STAR) or [HISAT2](http://daehwankimlab.github.io/hisat2/) | Read alignment | Select via `--aligner`; a matching pre-built genome index is required (`--genome-index`) |
| SLURM | Job scheduling | Only required if using `submit_array.sh`. Tested on partitions with `sbatch`/`squeue`; `--partition` must be a valid partition on your cluster |

**Not on a SLURM cluster?** The same HOMER/Python/STAR-or-HISAT2 requirements apply — SLURM is only used by `submit_array.sh` to schedule the jobs. The underlying `csrnaseq` Python package can be run directly (`python -m csrnaseq ...`) on any machine that meets the requirements above, without SLURM at all. See [csrnaseq/README.md](csrnaseq/README.md) for running it standalone.

Quick check that required tools are available in your environment:

```bash
conda activate CONDA_ENV_NAME
for t in homerTools STAR hisat2; do
    command -v "$t" >/dev/null && echo "OK: $t" || echo "MISSING: $t"
done
```

Note: only one of `STAR`/`hisat2` is required, depending on `--aligner`.

## Project layout

The submission flow expects a project directory with these standard folders:

- `RawData/`: input FASTQ files
- `Trimmed/`, `Aligned/`: intermediate files
- `TagDirs/`, `bedGraphs/`, `TSS/`: analysis outputs
- `QC/`, `Reports/`: QC plots and HTML reports
- `logs_slurm/`: SLURM logs

## Notes

- The repo is built around the `csrnaseq` CLI, so most of the actual processing logic lives in the package under [csrnaseq](csrnaseq).
- For package-specific installation details, module descriptions, and CLI internals, see [csrnaseq/README.md](csrnaseq/README.md).
