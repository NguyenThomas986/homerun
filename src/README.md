# Homerun — CSV-Driven Pipeline

**Supplementary Software for HOMER · Duttke Lab**

---

## Overview

Homerun automates the full RNA/csRNA-seq analysis pipeline:

```
Trim → STAR alignment → Tag directory → QC → DESeq2
```

The pipeline is **entirely driven by a CSV file** that acts as the single source of truth for all sample state. Every step is tracked per-sample, and the pipeline can be safely interrupted and resumed at any point.

---

## Quick Start

### 1. Generate a CSV from your FASTQ directory

```bash
python Homerun.py --generate-csv /path/to/fastqs \
                  --genome Homo_sapiens \
                  --out homerun_pipeline.csv
```

This scans the directory for `*.fastq.gz` files and creates one row per sample with all steps set to `PENDING`.

### 2. Run the pipeline

```bash
python Homerun.py --csv homerun_pipeline.csv --working-path /path/to/project
```

### 3. On SLURM

```bash
sbatch run_homerun.sh
```

---

## CSV Format

The pipeline CSV must contain these columns (in any order):

| Column   | Description |
|----------|-------------|
| `sample` | Unique sample identifier |
| `fastq`  | Absolute path to FASTQ file or directory |
| `genome` | Genome name (e.g. `Homo_sapiens`) |
| `trim`   | Output path or step status |
| `star`   | SAM/BAM output path or step status |
| `tagdir` | Tag directory path or step status |
| `qc`     | QC output directory or step status |
| `deseq2` | DESeq2 output path or step status |
| `status` | Global sample state |
| `notes`  | SLURM job IDs, errors, comments |

### Step column values

| Value | Meaning |
|-------|---------|
| `PENDING` or empty | Not yet run — will be executed |
| `DONE` | Completed — will be skipped |
| `FAILED` | Failed — skipped unless `--retry-failed` |
| Any file/dir path | Treated as DONE (output already exists) |

### Global status values

| Value | Meaning |
|-------|---------|
| `NOT_STARTED` | No steps have run yet |
| `RUNNING` | At least one step has completed |
| `FAILED` | A step failed; downstream steps were skipped |
| `DONE` | All steps complete |

---

## Resume Behavior

The CSV is updated **after every step**. If the pipeline is killed mid-run:

- Completed steps retain their output paths
- Incomplete steps remain `PENDING`
- Simply re-run the same command to pick up where it left off

---

## CLI Reference

```
python Homerun.py --csv pipeline.csv [options]

Required:
  --csv / -c PATH         Pipeline CSV file

Execution control:
  --dry-run               Print what would run without executing
  --retry-failed          Reset FAILED steps to PENDING and re-run
  --steps STEP [STEP...]  Only run these steps: trim star tagdir qc deseq2
  --samples NAME [...]    Only run these specific samples

Resources:
  --cpus N                CPUs to use (default: $SLURM_CPUS_PER_TASK or 1)
  --working-path PATH     Project root directory (default: current directory)

CSV generation:
  --generate-csv FASTQ_DIR   Auto-generate CSV from FASTQ directory
  --genome NAME              Genome name for generated CSV (default: Homo_sapiens)
  --out PATH                 Output path for generated CSV (default: homerun_pipeline.csv)
```

---

## Example: Resume after QC failure

```csv
sample,fastq,genome,trim,star,tagdir,qc,deseq2,status,notes
THP1_rep1,/data/THP1_R1.fastq.gz,Homo_sapiens,/out/THP1.trimmed,/out/THP1.sam,/tagDirs/THP1,FAILED,PENDING,FAILED,QC plot error
```

Re-run with `--retry-failed` to attempt QC again:

```bash
python Homerun.py --csv pipeline.csv --retry-failed --steps qc
```

---

## File Structure

The pipeline expects and creates this directory layout under `--working-path`:

```
<working-path>/
├── data/
│   └── <genome>/
│       ├── fastq/
│       │   ├── csRNA/       ← csRNA FASTQs + trimmed outputs
│       │   ├── sRNA/        ← sRNA FASTQs + trimmed outputs
│       │   └── totalRNA/    ← totalRNA FASTQs + trimmed outputs
│       └── tagDirs/         ← HOMER tag directories
├── genomes/
│   └── <genome>/
│       ├── *.fa             ← genome FASTA
│       └── STARIndex/       ← STAR genome index
├── files/
│   ├── QC/                  ← QC plots (.png / .svg)
│   ├── mappingStats/        ← STAR mapping stats per sample
│   ├── TSR/                 ← peak calling outputs
│   └── homerun_run_summary.tsv  ← end-of-run report
└── analysis/
    └── DESeq2/              ← differential expression outputs
```

---

## FASTQ Naming Convention

For `--generate-csv` to correctly parse sample names, FASTQ files should follow:

```
<sample>_csRNA_R1.fastq.gz
<sample>_sRNA_R1.fastq.gz
<sample>_totalRNA_R1.fastq.gz   (or _RNA_)
```

Paired-end R2 files are auto-detected from the R1 filename.

---

## Dependencies

| Tool | Version |
|------|---------|
| Python | ≥ 3.8 |
| STAR | ≥ 2.7 |
| HOMER | ≥ 4.11 |
| pandas | ≥ 1.3 |
| seaborn | ≥ 0.12 |
| matplotlib | ≥ 3.5 |
| numpy | ≥ 1.21 |

Optional: HISAT2, skewer (for totalRNA paired-end trimming)

---

## Module Summary

| File | Role |
|------|------|
| `Homerun.py` | CLI entry point and per-sample dispatcher |
| `HomerunCSV.py` | CSV read/write with file locking |
| `HomerunState.py` | Decides which steps need to run |
| `HomerunTrim.py` | Trimming (homerTools / skewer) |
| `HomerunSTAR.py` | STAR alignment |
| `HomerunTagdir.py` | HOMER tag directory creation |
| `HomerunQC.py` | QC plots (all original functions) |
| `HomerunDeseq2.py` | DESeq2 differential expression |
| `HomerunReport.py` | Console summary + TSV report |
| `HomerunGenCSV.py` | Auto-generate CSV from FASTQ directory |
| `run_homerun.sh` | SLURM submission script |
