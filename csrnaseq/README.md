# csrnaseq package reference

This package implements the csRNA-seq pipeline as a Python CLI. It handles preparation, trimming, alignment, tag directory creation, TSS calling, QC, stability analysis, and report generation.

## Installation

```bash
pip install -e .
```

The console entry point is defined in [pyproject.toml](../pyproject.toml) and points to `csrnaseq.pipeline:main`.

## Command-line usage

```bash
python -m csrnaseq --help
# or
csrnaseq --help
```

Useful flags include:

- `--project`: choose the project root
- `--steps`: run only a subset of the pipeline stages
- `--sample-index`: restrict trim + align to one sample
- `--only-prepare`: create folders and exit
- `--skip-prepare`: skip setup on later runs
- `--count-samples` / `--list-samples`: inspect input files

## Configuration

The package loads settings from environment variables prefixed with `CSRNA_` and from CLI flags. A minimal setup typically includes:

```bash
export CSRNA_PROJECT=/path/to/project
export CSRNA_ALIGNER=star
export CSRNA_GENOME_INDEX=/path/to/STARIndex
export CSRNA_GENOME=hg38
```

The main configuration logic lives in [config.py](config.py).

## Module overview

- [pipeline.py](pipeline.py): orchestrates the workflow in canonical order
- [prepare.py](prepare.py): creates output directories and stages inputs
- [trim.py](trim.py): trims adapters and filters reads
- [mapping.py](mapping.py): runs the aligner
- [tagdirs.py](tagdirs.py): builds HOMER tag directories
- [bedgraphs.py](bedgraphs.py): writes strand-specific bedGraph files
- [tss.py](tss.py): calls TSS clusters
- [qc.py](qc.py): generates QC plots and summaries
- [stability.py](stability.py): evaluates stability-related output
- [report.py](report.py): assembles HTML reports
- [utils.py](utils.py): shared helpers and logging

## Typical execution

```bash
csrnaseq \
  --project /path/to/project \
  --genome-index /path/to/STARIndex \
  --genome hg38 \
  --steps trim align tagdirs bedgraphs tss qc stability report
```

This package is intended to be the engine behind the repo-level SLURM helpers in the parent directory.
