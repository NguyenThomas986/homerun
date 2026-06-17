# csRNA-seq Pipeline

A SLURM-based pipeline for processing csRNA-seq data — trimming, alignment, tag directory construction, TSS calling, QC, and HTML report generation.

---

## Requirements

- HOMER (`homerTools`, `makeTagDirectory`, `makeUCSCfile`, `findcsRNATSS.pl`)
- STAR or HISAT2
- Python ≥ 3.10 with the `csrnaseq` package installed (see [Package README](csrnaseq/README.md))
- Conda environment with all dependencies (default: `miniComputer`)
- SLURM cluster

---

## Project layout

```
<project>/
├── config.env              # all pipeline settings (edit this before running)
├── submit_array.sh         # per-sample trim + align (SLURM array)
├── collect.sh              # downstream steps: tagdirs → TSS → QC → report
├── RawData/                # input FASTQ files (*_R1*.fastq.gz)
├── Trimmed/                # trimmed reads
├── Aligned/                # SAM files from STAR / HISAT2
├── TagDirs/                # HOMER tag directories (*-combo)
├── bedGraphs/              # strand-specific bedGraph files
├── TSS/                    # findcsRNATSS output (*.tss.txt, *.alltss.txt, *.stats.txt)
├── QC/                     # QC plots (PNG) and mapping stats
├── Reports/                # generated HTML reports
└── logs_slurm/             # SLURM stdout / stderr logs
```

---

## Quick start

### 1. Configure

Edit `config.env` — at minimum set these four variables:

```bash
export CSRNA_ALIGNER="star"           # "star" or "hisat2"
export CSRNA_GENOME_INDEX="/path/to/STARIndex"
export CSRNA_GENOME="hg38"
export CSRNA_PROJECT="/path/to/your/project"
```

See [Configuration reference](#configuration-reference) below for all options.

### 2. Run

```bash
# Step 1 — trim + align (one SLURM job per sample, runs in parallel)
sbatch submit_array.sh

# Step 2 — everything downstream (runs once after the array finishes)
sbatch collect.sh
```

`collect.sh` runs these steps in order: `tagdirs → bedgraphs → tss → qc → stability → report`

To run a single step manually:

```bash
csrnaseq --steps report
csrnaseq --steps qc stability report
```

---

## Input file naming

The pipeline classifies each sample by its filename. Files **must** contain one of these tags:

| Tag | Library type | Notes |
|---|---|---|
| `_csRNA` | csRNA-seq | single-end, trimmed with homerTools |
| `_sRNA` | spike-in control | single-end, trimmed with homerTools |
| `_totalRNA` or `_RNA` | total RNA-seq | paired-end, trimmed with skewer |

Example: `K562_csRNA_rep1_R1.fastq.gz`, `K562_sRNA_rep1_R1.fastq.gz`

---

## Configuration reference

All variables are set in `config.env` and loaded automatically by both SLURM scripts.

### Aligner & genome

| Variable | Default | Description |
|---|---|---|
| `CSRNA_ALIGNER` | `star` | Aligner to use: `star` or `hisat2` |
| `CSRNA_GENOME_INDEX` | *(required)* | STAR `--genomeDir` directory or HISAT2 `-x` index prefix |
| `CSRNA_GENOME` | *(required)* | HOMER `-genome` value (e.g. `hg38`) |
| `CSRNA_STARINDEX_URL` | *(optional)* | URL to auto-download a STAR index tarball if none is found |

### Project

| Variable | Default | Description |
|---|---|---|
| `CSRNA_PROJECT` | *(required)* | Absolute path to the project root directory |
| `CSRNA_COPY_SRC` | *(empty)* | Glob pattern for raw FASTQs to copy into `RawData/` before running |
| `CSRNA_ARRAY_THROTTLE` | `16` | Max number of array tasks running simultaneously |

### Trimming

| Variable | Default | Description |
|---|---|---|
| `CSRNA_TRIM_ADAPTER` | `AGATCGGAAGAGCACACGTCT` | 3′ adapter sequence |
| `CSRNA_TRIM_MINLEN` | `20` | Discard reads shorter than this after trimming |
| `CSRNA_TRIM_MAXLEN` | `58` | Discard reads longer than this after trimming |

### TSS calling

| Variable | Default | Description |
|---|---|---|
| `CSRNA_NTAG_THRESHOLD` | `7` | Minimum tags to call a TSS cluster |
| `CSRNA_SKIP_CHR` | `chrEBV` | Chromosome to exclude from bedGraphs |

### Cluster

| Variable | Default | Description |
|---|---|---|
| `CSRNA_CONDA_MODULE` | `anaconda3` | Cluster module providing Conda |
| `CSRNA_CONDA_ENV` | `miniComputer` | Conda environment with pipeline dependencies |
| `CSRNA_PARTITION` | `kamiak` | SLURM partition for job submission |
| `CSRNA_EMAIL` | *(empty)* | Email address for SLURM job notifications |

---

## Output files

| Path | Contents |
|---|---|
| `TSS/<sample>.tss.txt` | Filtered TSS clusters with annotation |
| `TSS/<sample>.alltss.txt` | All candidate TSS before score filtering |
| `TSS/<sample>.stats.txt` | Run statistics from `findcsRNATSS.pl` |
| `QC/*.png` | All QC plots |
| `Reports/tss_clusters.html` | Interactive TSS cluster tables |
| `Reports/alltss.html` | All pre-filter TSS candidates |
| `Reports/stats.html` | Run statistics per sample |
| `Reports/qc_report.html` | All QC plots + data files in one page |

---

## Re-running steps

All steps are idempotent — outputs are skipped if they already exist. To re-run a step, delete its output first:

```bash
# Re-run QC and report only
rm -rf QC/* Reports/*
csrnaseq --steps qc stability report

# Re-run TSS calling for one sample
rm TSS/mysample.*
csrnaseq --steps tss
```

---

## Package documentation

See [csrnaseq/README.md](csrnaseq/README.md) for installation instructions and a description of every module in the Python package.
