# homerun

<div align="center">

[![Open Issues](https://img.shields.io/github/issues/NguyenThomas986/homerun.svg)](https://github.com/NguyenThomas986/homerun/issues)
[![Closed Issues](https://img.shields.io/github/issues-closed/NguyenThomas986/homerun.svg)](https://github.com/NguyenThomas986/homerun/issues?q=is%3Aissue+is%3Aclosed)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

</div>

This repository is the SLURM-oriented wrapper around the Python package in [csrnaseq](csrnaseq). It is meant to help you launch a complete csRNA-seq analysis workflow from a project directory on a cluster.

## What this repo contains

- [submit_array.sh](submit_array.sh): submits the full 7-phase workflow
- [prepare.sbatch](prepare.sbatch): prepares the project and runs the initial setup
- [align_array.sbatch](align_array.sbatch): runs trim + align for one leaf run per array task
- [tagdir_array.sbatch](tagdir_array.sbatch): builds leaf-level TagDirs in parallel, one array task per leaf
- [tagdirs_combo_array.sbatch](tagdirs_combo_array.sbatch): merges replicates into combo TagDirs in parallel, one array task per Species/Sample
- [tss_array.sbatch](tss_array.sbatch): calls TSRs with `findcsRNATSS.pl` in parallel, one array task per Species/Sample
- [bedgraphs_array.sbatch](bedgraphs_array.sbatch): builds bedGraphs in parallel, one array task per Species/Sample
- [collect.sbatch](collect.sbatch): runs `ritrie`, QC, stability, and reporting after every array phase finishes
- [csrnaseq](csrnaseq): the Python package that actually performs the pipeline steps

## Project structure

Projects are organized as a nested `Species/Sample/Leaf` directory layout. Rather than writing outputs to shared top-level directories, each leaf run generates and keeps its own intermediate outputs alongside it, and each sample keeps its own combined/downstream outputs one level up.

Species and sample names are parsed straight out of the FASTQ filenames. For example, `homo_sapiens_K562_csRNA-r1_DB422_S1_R1_001.fastq.gz` and `homo_sapiens_K562_csRNA-r2_DB423_S2_R1_001.fastq.gz` become species `homo_sapiens`, sample `K562`, with leaves `csRNA_r1` and `csRNA_r2`. Each distinct assay/replicate combination becomes its own leaf directory, and paired-end reads (files sharing a leaf but differing only in `_R1`/`_R2`) are staged into that same leaf.

This produces a tree like:

    homerun_test/
    ├── config.txt                        # auto-generated summary of every run's config + samples
    └── homo_sapiens/
        └── K562/
            ├── RNA_r1/                    # paired-end (R1 + R2)
            │   ├── RawData/
            │   ├── Trimmed/
            │   ├── Aligned/
            │   └── TagDir/
            ├── csRNA_r1/
            │   ├── RawData/
            │   ├── Trimmed/
            │   ├── Aligned/
            │   └── TagDir/
            ├── csRNA_r2/
            │   └── ...
            ├── csRNA-combo/               # combined TagDir/bedGraph across csRNA_r1 + csRNA_r2
            │   ├── TagDir/
            │   └── bedGraph/
            ├── sRNA_r1/
            │   └── ...
            ├── sRNA_r2/
            │   └── ...
            ├── sRNA-combo/
            │   ├── TagDir/
            │   └── bedGraph/
            ├── TSS/                       # called TSRs for this sample
            ├── QC/                        # QC plots, ritrie output, qc_report.html
            └── QC/qc_report.html

- **Trimmed reads** are written into each leaf's local `Trimmed/` directory.
- **Alignments** are written into each leaf's local `Aligned/` directory.
- **Leaf-level TagDirs** are built inside each replicate leaf (e.g. `K562/csRNA_r1/TagDir/`).
- **Combined TagDirs/bedGraphs** are built as their own sibling directory per assay, named `<assay>-combo/` (e.g. `K562/csRNA-combo/`, `K562/sRNA-combo/`), merging that assay's replicate leaves.
- **TSS, QC, RITRIE, and reporting** are built within each sample's directory (e.g. `homo_sapiens/K562/`), downstream of the combo TagDirs.
- **Reporting** is simplified to a single `qc_report.html` per sample.
- **`config.txt`** at the project root records every config value used for the run (genome, aligner, thresholds, `--gtf`, etc.) and every sample/RawData file discovered — see [Run configuration (config.txt)](#run-configuration-configtxt) below.

## Typical workflow

1. Create a project directory and place your FASTQ files in it — they'll be parsed and staged into the nested `Species/Sample/Leaf` layout automatically (e.g. `homo_sapiens/K562/csRNA_r1/RawData/`).
2. Run the submission script with the required paths and cluster settings.
3. Let the prepare, align-array, tagdir/tagdirs-combo-array, tss/bedgraphs-array, and collect stages run in sequence (see dependency graph below).
4. Inspect the outputs alongside each sample, especially `QC/qc_report.html` and the `TSS/` directory (e.g. `homo_sapiens/K562/QC/qc_report.html`), and check `config.txt` at the project root to confirm what the run was actually configured with.

### Job-array dependency graph

```
prepare
  └─afterok─> align_array[0..N-1]
                ├─afterok─> tagdir_array[0..N-1]          (leaf TagDirs)
                └─afterok─> tagdirs_combo_array[0..S-1]   (combo TagDirs)
                              ├─afterok─> tss_array[0..S-1]
                              └─(+ tagdir_array)─afterok─> bedgraphs_array[0..S-1]
                                                              └─afterok─(+tss_array)─> collect
```

- **N** = number of leaf runs (R1 files in RawData) — `align_array`/`tagdir_array` indexing, one task per leaf run.
- **S** = number of Species/Sample groups — `tagdirs_combo_array`/`bedgraphs_array`/`tss_array` indexing, one task per Species/Sample (usually smaller than N, since a sample typically has several leaf runs).

`tagdirs_combo_array` only needs `align_array` to finish (not `tagdir_array`), so it runs in parallel with the leaf TagDir build rather than after it — combo TagDirs are built straight from the aligned SAM files, not from the leaf TagDirs. `tss_array` needs the combo TagDirs. `bedgraphs_array` needs **both** `tagdir_array` and `tagdirs_combo_array`, since it builds bedGraphs for every TagDir (leaf and combo alike) under a sample. `collect` then runs `ritrie`, QC, stability, and reporting over everything the array phases built.

### TagDir generation

- **Leaf TagDirs** — built in parallel via `tagdir_array.sbatch`, one array task per leaf (e.g. one task each for `K562/csRNA_r1`, `K562/csRNA_r2`, `K562/sRNA_r1`, `K562/sRNA_r2`, `K562/RNA_r1`).
- **Combined TagDirs** — built in parallel via `tagdirs_combo_array.sbatch`, one array task per Species/Sample, merging that sample's replicate leaves into a `<assay>-combo/` directory (e.g. `K562/csRNA-combo/`, `K562/sRNA-combo/`).

### RITRIE (RIT/RIE QC metric)

`ritrie` computes RIT/RIE — reads at called TSRs vs. reads at exons (excluding TSR overlap) — a RIN-free proxy for RNA degradation. Higher values generally indicate stronger enrichment of transcription-initiation signal relative to exonic RNA; lower values may indicate more degraded RNA. Thresholds haven't yet been formally established.

This step requires a genome annotation GTF, passed via `--gtf` (or `CSRNA_GTF`). If it's unset, `ritrie` is simply skipped — that's not an error. If it *is* set but points at a missing/unreadable file, `prepare` fails immediately with a clear error rather than letting the pipeline run for hours and fail on `ritrie` at the very end.

### Stability handling

The stability step is automatically skipped for samples without matching total RNA data. For example, if a sample has `csRNA_r1`/`csRNA_r2` and `sRNA_r1`/`sRNA_r2` leaves but no `RNA_r*` leaf, stability analysis is skipped for that sample. When stability analysis isn't performed, stability-specific columns are omitted from the report.

### Run configuration (config.txt)

Every `prepare` run writes (or refreshes) `<project>/config.txt` — a plain-text summary of the run, meant to make a project directory self-documenting without digging through `config.env`, CLI history, or logs:

```
# HomeRun csRNA-seq pipeline — run configuration
# Generated: 2026-07-14T12:01:49

[Config]
project = /path/to/homerun_test
aligner = star
genome_index = /path/to/STARIndex
genome = hg38
gtf = /path/to/gencode.v44.annotation.gtf
threads = 8
...

[Samples]
count = 2
homo_sapiens/HepG2
homo_sapiens/K562

[RawData]
count = 3
homo_sapiens/HepG2/sRNA_r1/RawData/homo_sapiens_HepG2_sRNA-r1_DB477_S7_R1_001.fastq.gz
homo_sapiens/K562/csRNA_r1/RawData/homo_sapiens_K562_csRNA-r1_DB422_S1_R1_001.fastq.gz
homo_sapiens/K562/csRNA_r2/RawData/homo_sapiens_K562_csRNA-r2_DB423_S2_R1_001.fastq.gz
```

`[Config]` lists every setting the run actually used — whatever was passed via CLI flag or `CSRNA_*` env var, or its built-in default if neither was set. `[Samples]` and `[RawData]` list every Species/Sample discovered and every staged FASTQ, so you can confirm what the pipeline saw before it spends hours processing it. It's safe to re-run `prepare` on a project that already has staged data or an old `config.txt` sitting in it — the file is simply overwritten with current information each time, and any leftovers from an older/different layout are ignored.

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

### Quick start (with RITRIE enabled)

```bash
path/to/homerun/submit_array.sh \
  --project /PATH/TO/PROJECT \
  --partition kamiak \
  --conda-env CONDA_ENV_NAME \
  --genome-index /path/to/STARIndex \
  --genome hg38 \
  -- --gtf /path/to/genome.gtf
```

The controller will:

- parse and stage raw FASTQs into the nested `Species/Sample/Leaf` project layout (e.g. `homo_sapiens_K562_csRNA-r1_DB422_S1_R1_001.fastq.gz` → `homo_sapiens/K562/csRNA_r1/RawData/`)
- write/refresh `config.txt` at the project root
- submit a prepare job (which also validates `--gtf` up front, if given)
- submit a sample-array job for trim + align
- submit a tagdir-array job (leaf TagDirs) and a tagdirs-combo-array job (combo TagDirs) in parallel
- submit a tss-array job and a bedgraphs-array job, each depending on the combo TagDirs (bedgraphs also needs the leaf TagDirs)
- submit a collect job that runs `ritrie`, QC, stability, and reporting

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
| `--gtf` | | Genome annotation GTF; enables the `ritrie` step. Skipped (not an error) if unset |

To pass extra flags directly to the Python pipeline (e.g. STAR tuning, or `--gtf`):

```bash
path/to/homerun/submit_array.sh \
  --project /PATH/TO/PROJECT \
  --partition kamiak \
  --conda-env CONDA_ENV_NAME \
  --genome-index /path/to/STARIndex \
  --genome hg38 \
  -- --star-filter-multimap 5000 --gtf /path/to/genome.gtf
```

To see all the flags do:

```bash
path/to/homerun/submit_array.sh \
    --h (or --help)
```

Or call the Python package directly for the same output:

```
$ python -m csrnaseq --help
usage: csrnaseq [-h] [--project PROJECT] [--log-path LOG_PATH]
                 [--steps {trim,align,tagdirs,tagdirs-combo,bedgraphs,tss,ritrie,qc,stability,report} [{trim,align,tagdirs,tagdirs-combo,bedgraphs,tss,ritrie,qc,stability,report} ...]]
                 [--sample-index SAMPLE_INDEX] [--group-index GROUP_INDEX] [--skip-prepare] [--only-prepare] [--count-samples] [--count-groups]
                 [--stage-raw] [--aligner {star,hisat2}] [--genome-index GENOME_INDEX] [--genome GENOME] [--gtf GTF] [--copy-src COPY_SRC]
                 [--threads THREADS] [--trim-adapter TRIM_ADAPTER] [--trim-min TRIM_MIN] [--trim-max TRIM_MAX] [--ntag-threshold NTAG_THRESHOLD]
                 [--skip-chr SKIP_CHR] [--star-filter-multimap STAR_FILTER_MULTIMAP] [--star-multimap-out STAR_MULTIMAP_OUT]
                 [--star-multimap-order {Random,Old_2.4}] [--hisat2-strandness {F,R,FR,RF}]

    __  __                                    
   / / / /___  ____ ___  ___  _______  ______ 
  / /_/ / __ \/ __ `__ \/ _ \/ ___/ / / / __ \
 / __  / /_/ / / / / / /  __/ /  / /_/ / / / /
/_/ /_/\____/_/ /_/ /_/\___/_/   \__,_/_/ /_/ 

RNA-seq analysis pipeline for HPC clusters.
Version: 1.0.0

optional arguments:
  -h, --help            show this help message and exit
  --project PROJECT     Project root (default: $CSRNA_PROJECT or CWD).
  --log-path LOG_PATH   Pipeline log file path (overrides CSRNA_LOG; else timestamped file under <project>/logs/).
  --steps {trim,align,tagdirs,tagdirs-combo,bedgraphs,tss,ritrie,qc,stability,report} [{trim,align,tagdirs,tagdirs-combo,bedgraphs,tss,ritrie,qc,stability,report} ...]
                        Run only these steps (still executed in canonical order).
  --sample-index SAMPLE_INDEX
                        0-based index into RawData R1 files. Restricts trim/align/tagdirs to one leaf run. Used by SLURM_ARRAY_TASK_ID.
  --group-index GROUP_INDEX
                        0-based index into list_samples(cfg). Restricts tagdirs-combo/bedgraphs/tss to one Species/Sample group. Used by
                        SLURM_ARRAY_TASK_ID.
  --skip-prepare        Skip folder creation/raw copy/STARIndex setup.
  --only-prepare        Run prepare and exit.
  --count-samples       Print number of leaf runs (R1 files in RawData) and exit.
  --count-groups        Print number of Species/Sample groups and exit.
  --stage-raw           Move loose *_R1*/*_R2* FASTQs into RawData/ and exit.

config overrides (override config.env when given):
  --aligner {star,hisat2}
                        Aligner (overrides CSRNA_ALIGNER).
  --genome-index GENOME_INDEX
                        STAR genomeDir or HISAT2 prefix (overrides CSRNA_GENOME_INDEX).
  --genome GENOME       HOMER genome (overrides CSRNA_GENOME).
  --gtf GTF             GTF annotation file for RIT/RIE metric (overrides CSRNA_GTF).
  --copy-src COPY_SRC   FASTQ copy source (overrides CSRNA_COPY_SRC).
  --threads THREADS     Threads (overrides CSRNA_THREADS / SLURM_CPUS_PER_TASK).
  --trim-adapter TRIM_ADAPTER
                        3' adapter sequence.
  --trim-min TRIM_MIN   Minimum read length after trimming.
  --trim-max TRIM_MAX   Maximum read length after trimming.
  --ntag-threshold NTAG_THRESHOLD
                        Minimum tags for TSS calling.
  --skip-chr SKIP_CHR   Chromosome excluded from bedGraphs.

alignment overrides:
  --star-filter-multimap STAR_FILTER_MULTIMAP
                        STAR outFilterMultimapNmax.
  --star-multimap-out STAR_MULTIMAP_OUT
                        STAR outSAMmultNmax.
  --star-multimap-order {Random,Old_2.4}
                        STAR outMultimapperOrder.
  --hisat2-strandness {F,R,FR,RF}
                        HISAT2 RNA strandness.
```

## Requirements

The conda environment passed via `--conda-env` must provide:

| Tool | Used for | Notes |
|---|---|---|
| [HOMER](http://homer.ucsd.edu/homer/) | Tag directories, bedGraphs, TSS calling/annotation, RITRIE | Must include `homerTools` on `PATH`; genome data for `--genome` (e.g. `hg38`) must be installed via `perl configureHomer.pl -install hg38` |
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
- Combined TagDirs/bedGraphs are per-assay `<assay>-combo/` directories (e.g. `csRNA-combo`, `sRNA-combo`), not a single catch-all folder — this keeps replicate merges scoped to their own assay.
- `tagdirs-combo`, `bedgraphs`, and `tss` each run as their own parallel SLURM array (indexed by `--group-index`, one task per Species/Sample) rather than looping serially inside `collect` — see the [job-array dependency graph](#job-array-dependency-graph) above.
- `config.txt` at the project root is regenerated at the end of every `prepare` run — check it first if a run behaves unexpectedly, to confirm the pipeline actually saw the config/samples you intended.
- For package-specific installation details, module descriptions, and CLI internals, see [csrnaseq/README.md](csrnaseq/README.md).
