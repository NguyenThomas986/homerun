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
- `--check-rerun`: preflight exit code for wrapper scripts (see below);
  prints nothing on success
- `--force` / `--overwrite`: allow rerunning on a project that already has
  pipeline outputs (see below)

### Rerunning on an existing project

Every step already skips its own finished work (via a `done()`-style check
before each command), so re-running the pipeline on a project that has some
samples finished and others new or in-progress is always safe and does
exactly what you'd want — it fills in whatever's missing and leaves
completed samples untouched. **No flag is needed for that.**

HOMERun only refuses to start in one specific case: every sample currently
staged in the project already has a populated `QC/` (i.e. the pipeline has
already run to completion for everything present) and no new or
partially-processed sample was found. That means the invocation would touch
nothing new — almost always a sign you pointed `--project` at the wrong
place, or forgot you'd already run this project. In that case it exits
immediately with a message listing the existing output it found, before
`--only-prepare` or any pipeline step runs. Pass `--force` (or `--overwrite`)
to proceed anyway (e.g. to force a full redo). This check doesn't apply to
`--count-samples`, `--count-groups`, or `--stage-raw`, which never modify or
delete existing outputs and are used as a cheap pre-flight before SLURM
array submission.

If you're driving the pipeline through `submit_array.sh` (or any wrapper
that submits each step as its own `sbatch` job), the in-process check above
runs too late to help — by the time a submitted job reaches
`--only-prepare`/`run_pipeline`, every other job is already queued.
`submit_array.sh` instead calls `--check-rerun` as a lightweight preflight,
right after it reports the sample/group counts and before it submits
anything: it does the same "is there anything new to do" check as above and
exits 1 with the same message if not, which (under `set -euo pipefail`)
stops the script before any `sbatch` call. Pass `--force`/`--overwrite` after
`--` to `submit_array.sh` to bypass it the same way.

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

## QC: sample-level vs. per-replicate

`qc.py` generates two kinds of plots into each sample's `QC/` folder, and
`qc_report.html` presents them in two separate sections:

- **Sample-Level QC** — built from the merged **combo** TagDir(s) (one per
  assay: `csRNA-combo`, `sRNA-combo`, `totalRNA-combo`). This describes the
  final, merged sample and is unchanged from before: read length
  distribution, nucleotide frequency, autocorrelation, tag directory stats,
  median tags/position, the A-plot and tags-vs-fraction-of-positions
  overlays, plus everything derived from TSS calling (threshold
  optimization, TSR summary/annotation, TSS nucleotide frequency,
  stability/location, and the distal-vs-proximal pie chart) — these are all
  inherently per-sample, since `findcsRNATSS.pl` only ever runs once per
  sample, on the combo TagDirs, not per replicate.
- **Per-Replicate QC** — built from each **individual** leaf TagDir
  (`csRNA_r1`, `csRNA_r2`, `sRNA_r1`, ...), so a problem specific to one
  replicate (a bad library, an alignment issue, an unusual length or
  nucleotide profile) is visible even after replicates are merged into the
  combo. Only metrics that are actually meaningful per-replicate are
  duplicated here: read length distribution, nucleotide frequency,
  autocorrelation, and tag directory stats. Rendered as a compact grid (one
  small panel per replicate) rather than one crowded overlay, so it scales
  to samples with many replicates.

## Typical execution

```bash
csrnaseq \
  --project /path/to/project \
  --genome-index /path/to/STARIndex \
  --genome hg38 \
  --steps trim align tagdirs bedgraphs tss qc stability report
```

This package is intended to be the engine behind the repo-level SLURM helpers in the parent directory.