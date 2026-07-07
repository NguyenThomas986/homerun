# homerun

<div align="center">

[![Open Issues](https://img.shields.io/github/issues/NguyenThomas986/homerun.svg)](https://github.com/NguyenThomas986/homerun/issues)
[![Closed Issues](https://img.shields.io/github/issues-closed/NguyenThomas986/homerun.svg)](https://github.com/NguyenThomas986/homerun/issues?q=is%3Aissue+is%3Aclosed)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

</div>

This repository is the SLURM-oriented wrapper around the Python package in [csrnaseq](csrnaseq). It is meant to help you launch a complete csRNA-seq analysis workflow from a project directory on a cluster.

## What this repo contains

- [submit_array.sh](submit_array.sh): submits the full three-phase workflow
- [prepare.sbatch](prepare.sbatch): prepares the project and runs the initial setup
- [align_array.sbatch](align_array.sbatch): runs trim + align for one sample per array task
- [tagdir_array.sbatch](tagdir_array.sbatch): builds leaf-level TagDirs in parallel, one array task per leaf
- [collect.sbatch](collect.sbatch): combines TagDirs and runs downstream analyses after the array finishes
- [csrnaseq](csrnaseq): the Python package that actually performs the pipeline steps

## Project structure

Projects are organized as a nested `Species/Sample/Leaf` directory layout. Rather than writing outputs to shared top-level directories, each sample generates and keeps its own intermediate and final outputs alongside it. Legacy flat output directories (`TagDirs`, `bedGraphs`, `TSS`, `QC`, `Reports`) are no longer used.

Species and sample names are parsed straight out of the FASTQ filenames. For example, `homo_sapiens_K562_csRNA-r1_DB422_S1_R1_001.fastq.gz` and `homo_sapiens_K562_csRNA-r2_DB423_S2_R1_001.fastq.gz` become species `homo_sapiens`, sample `K562`, with leaves `csRNA_r1` and `csRNA_r2`. Each distinct assay/replicate combination (`csRNA_r1`, `sRNA_r1`, `RNA_r1`, ...) becomes its own leaf directory, and paired-end reads (files sharing a leaf but differing only in `_R1`/`_R2`) are staged into that same leaf.

This produces a tree like:

    homerun_test/
    └── homo_sapiens/
        ├── HepG2/
        │   ├── csRNA_r1/
        │   │   ├── RawData/
        │   │   ├── Trimmed/
        │   │   ├── Aligned/
        │   │   └── TagDirs/         # leaf-level TagDir
        │   ├── csRNA_r2/
        │   │   └── ...
        │   ├── csRNA-combo/          # combined TagDir across csRNA_r1 + csRNA_r2
        │   ├── sRNA_r1/
        │   │   └── ...
        │   ├── sRNA_r2/
        │   │   └── ...
        │   ├── sRNA-combo/           # combined TagDir across sRNA_r1 + sRNA_r2
        │   ├── bedGraphs/
        │   ├── TSS/
        │   ├── QC/
        │   └── qc_report.html
        └── K562/
            ├── RNA_r1/                # paired-end (R1 + R2)
            │   ├── RawData/
            │   ├── Trimmed/
            │   ├── Aligned/
            │   └── TagDirs/
            ├── csRNA_r1/
            │   └── ...
            ├── csRNA_r2/
            │   └── ...
            ├── csRNA-combo/           # combined TagDir across csRNA_r1 + csRNA_r2
            ├── sRNA_r1/
            │   └── ...
            ├── sRNA_r2/
            │   └── ...
            ├── sRNA_combo/            # combined TagDir across sRNA_r1 + sRNA_r2
            ├── bedGraphs/
            ├── TSS/
            ├── QC/
            └── qc_report.html

- **Trimmed reads** are written into each leaf's local `Trimmed/` directory.
- **Alignments** are written into each leaf's local `Aligned/` directory.
- **Leaf-level TagDirs** are built inside each replicate leaf (e.g. `K562/csRNA_r1/TagDirs/`).
- **Combined TagDirs** are built as their own sibling directory per assay, named `<assay>_combo/` (e.g. `K562/csRNA_combo/`, `K562/sRNA_combo/`), merging that assay's replicate leaves.
- **bedGraphs, TSS results, and QC outputs** are built within each sample's directory (e.g. `homo_sapiens/K562/`), downstream of the combo TagDirs.
- **Reporting** is simplified to a single `qc_report.html` per sample.

## Typical workflow

1. Create a project directory and place your FASTQ files in it — they'll be parsed and staged into the nested `Species/Sample/Leaf` layout automatically (e.g. `homo_sapiens/K562/csRNA_r1/RawData/`).
2. Run the submission script with the required paths and cluster settings.
3. Let the prepare, align-array, tagdir-array, and collect stages run in sequence.
4. Inspect the outputs alongside each sample, especially the `qc_report.html` and `TSS/` directory (e.g. `homo_sapiens/HepG2/qc_report.html`).

### Stability handling

The stability step is automatically skipped for samples without matching total RNA data. For example, if a sample has `csRNA_r1`/`csRNA_r2` and `sRNA_r1`/`sRNA_r2` leaves but no `RNA_r*` leaf, stability analysis is skipped for that sample. When stability analysis isn't performed, stability-specific columns are omitted from the report.

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
| `--copy-src` | | Glob path to FASTQs to copy into the nested sample layout if none are present |
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

## Notes

- The repo is built around the `csrnaseq` CLI, so most of the actual processing logic lives in the package under [csrnaseq](csrnaseq).
- Outputs live per-sample (nested under `Species/Sample/`, e.g. `homo_sapiens/K562/`) rather than in shared top-level directories, and each sample produces a single `qc_report.html` with concise, focused figures.
- Combined TagDirs are per-assay `<assay>_combo/` directories (e.g. `csRNA_combo`, `sRNA_combo`), not a single catch-all `TagDirs/` folder — this keeps replicate merges scoped to their own assay.
- For package-specific installation details, module descriptions, and CLI internals, see [csrnaseq/README.md](csrnaseq/README.md).
