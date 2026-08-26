# CLI Reference

!!! warning "Alpha"
    HOMERun is in alpha testing. Nothing here is finalized, and some flags and
    behaviors are still rough or subject to change.

Full usage for the `homerun` command. This page mirrors `homerun --help`.

On start, `homerun` prints:

```text
    __  __
   / / / /___  ____ ___  ___  _______  ______
  / /_/ / __ \/ __ `__ \/ _ \/ ___/ / / / __ \
 / __  / /_/ / / / / / /  __/ /  / /_/ / / / /
/_/ /_/\____/_/ /_/ /_/\___/_/   \__,_/_/ /_/

RNA-seq analysis pipeline for HPC clusters.
Version: 1.0.0
```
<!-- TODO: banner still says "RNA-seq"; source string in build_parser/_BANNER should read "csRNA-seq" — update there and this block together. Version is hardcoded here; it will drift when __version__ changes. -->

## Usage

```text
homerun [-h] [--project PROJECT] [--log-path LOG_PATH]
        [--steps {trim,align,tagdirs,tagdirs-combo,bedgraphs,tss,ritrie,qc,stability,report} ...]
        [--sample-index SAMPLE_INDEX] [--group-index GROUP_INDEX]
        [--skip-prepare] [--force] [--only-prepare]
        [--count-samples] [--count-groups] [--check-rerun] [--stage-raw]
        [--aligner {star,hisat2}] [--genome-index GENOME_INDEX] [--genome GENOME]
        [--gtf GTF] [--copy-src COPY_SRC] [--threads THREADS]
        [--trim-adapter TRIM_ADAPTER] [--trim-min TRIM_MIN] [--trim-max TRIM_MAX]
        [--ntag-threshold NTAG_THRESHOLD] [--skip-chr SKIP_CHR]
        [--cleanup-intermediates]
        [--star-filter-multimap STAR_FILTER_MULTIMAP]
        [--star-multimap-out STAR_MULTIMAP_OUT]
        [--star-multimap-order {Random,Old_2.4}]
        [--hisat2-strandness {F,R,FR,RF}]
```

## Options

| Flag | Description |
| --- | --- |
| `-h`, `--help` | Show the help message and exit. |
| `--project PROJECT` | Project root (default: `$CSRNA_PROJECT` or CWD). |
| `--log-path LOG_PATH` | Pipeline log file path (overrides `CSRNA_LOG`; else a timestamped file under `<project>/logs/`). |
| `--steps STEP [STEP ...]` | Run only these steps (still executed in canonical order). Choices: `trim`, `align`, `tagdirs`, `tagdirs-combo`, `bedgraphs`, `tss`, `ritrie`, `qc`, `stability`, `report`. |
| `--sample-index N` | 0-based index into RawData R1 files. Restricts `trim`/`align`/`tagdirs` to one leaf run. Used by `SLURM_ARRAY_TASK_ID`. |
| `--group-index N` | 0-based index into the sample list. Restricts `tagdirs-combo`/`bedgraphs`/`tss` to one Species/Sample group. Used by `SLURM_ARRAY_TASK_ID`. |
| `--skip-prepare` | Skip folder creation / raw copy / STARIndex setup. |
| `--force`, `--overwrite` | Allow rerunning on a project that already has pipeline outputs (`Trimmed/`, `Aligned/`, `TagDirs/`, `bedGraphs/`, `TSS/`, `QC/`). Without this flag, HOMERun refuses to start (and touches nothing) if any of those exist, to avoid overwriting a previous run. No effect on `--count-samples`/`--count-groups`/`--stage-raw`. |
| `--only-prepare` | Run prepare and exit. |
| `--count-samples` | Print the number of leaf runs (R1 files in RawData) and exit. |
| `--count-groups` | Print the number of Species/Sample groups and exit. |
| `--check-rerun` | Preflight for wrapper scripts: exit 1 (with a message) if every staged sample already has a populated `QC/` and `--force` was not given — i.e. the run would touch nothing new. Exit 0 otherwise. Submits/runs nothing itself; call it *before* submitting SLURM jobs. |
| `--stage-raw` | Move loose `*_R1*` / `*_R2*` FASTQs into `RawData/` and exit. |

## Config overrides

These override `config.env` values when given.
<!-- TODO: config.env is being deprecated in favor of flags — reword once the precedence story is finalized. -->

| Flag | Description |
| --- | --- |
| `--aligner {star,hisat2}` | Aligner (overrides `CSRNA_ALIGNER`). |
| `--genome-index GENOME_INDEX` | STAR genomeDir or HISAT2 prefix (overrides `CSRNA_GENOME_INDEX`). |
| `--genome GENOME` | HOMER genome (overrides `CSRNA_GENOME`). |
| `--gtf GTF` | GTF annotation file for the RIT/RIE metric (overrides `CSRNA_GTF`). |
| `--copy-src COPY_SRC` | FASTQ copy source (overrides `CSRNA_COPY_SRC`). |
| `--threads THREADS` | Threads (overrides `CSRNA_THREADS` / `SLURM_CPUS_PER_TASK`). |
| `--trim-adapter TRIM_ADAPTER` | 3′ adapter sequence. |
| `--trim-min TRIM_MIN` | Minimum read length after trimming. |
| `--trim-max TRIM_MAX` | Maximum read length after trimming. |
| `--ntag-threshold NTAG_THRESHOLD` | Minimum tags for TSS calling. |
| `--skip-chr SKIP_CHR` | Chromosome excluded from bedGraphs. |
| `--cleanup-intermediates` | Delete each sample's `Trimmed/` and `Aligned/` directories once the `qc` step has generated its tables from them (overrides `CSRNA_CLEANUP_INTERMEDIATES`). Off by default. |

## Alignment overrides

| Flag | Description |
| --- | --- |
| `--star-filter-multimap N` | STAR `outFilterMultimapNmax`. |
| `--star-multimap-out N` | STAR `outSAMmultNmax`. |
| `--star-multimap-order {Random,Old_2.4}` | STAR `outMultimapperOrder`. |
| `--hisat2-strandness {F,R,FR,RF}` | HISAT2 RNA strandness. |
