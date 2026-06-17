# csrnaseq — Python package

The `csrnaseq` package contains all pipeline logic for the csRNA-seq workflow. It is designed to be installed once and then driven entirely through `config.env` and SLURM scripts — no code editing required for routine runs.

---

## Installation

From the project root (where `pyproject.toml` or `setup.py` lives):

```bash
pip install -e .
```

Or inside your Conda environment:

```bash
conda activate miniComputer
pip install -e .
```

After installation the pipeline is available as a command:

```bash
csrnaseq --help
```

---

## Running

```bash
# Full pipeline
csrnaseq

# Specific steps only (executed in canonical order regardless of input order)
csrnaseq --steps tagdirs bedgraphs tss qc stability report

# Single sample (SLURM array usage)
csrnaseq --steps trim align --sample-index 3

# Count samples in RawData (used by submit_array.sh)
csrnaseq --count-samples

# Prepare directories only
csrnaseq --only-prepare
```

All settings are read from environment variables set in `config.env`. Override any variable on the command line with `--project`:

```bash
csrnaseq --project /data/lab/myproject --steps report
```

---

## Step order

```
trim → align → tagdirs → bedgraphs → tss → qc → stability → report
```

Steps fail fast — if any step raises an exception the pipeline stops immediately so downstream steps never run on incomplete inputs.

`trim` and `align` accept `--sample-index` for SLURM array parallelism. All other steps run once across all samples.

---

## Module reference

### `__main__.py`
Entry point for `python -m csrnaseq`. Delegates to `pipeline.main()`.

---

### `pipeline.py`
**Orchestrator.** Parses CLI arguments, loads config, and runs the selected steps in canonical order.

Key functions:
- `build_parser()` — constructs the `argparse` parser
- `run_pipeline(cfg, steps, skip_prepare, sample_index)` — executes steps in order, fail-fast
- `main(argv)` — CLI entry point; returns exit code

Step registry (`STEP_ORDER` / `STEP_FUNCS`) maps step names to their `run_*` functions. Adding a new step means adding one entry to each dict.

---

### `config.py`
**Central configuration.** All settings live here as a `Config` dataclass. Every field has a default; values are overridden by `CSRNA_*` environment variables at runtime.

Key items:
- `Config` dataclass — all pipeline settings as typed fields
- `load_config(project)` — reads env vars, returns a populated `Config`
- Derived path properties: `rawdata`, `trimmed`, `aligned`, `tagdirs`, `bedgraphs`, `tss`, `qc`, `logs_dir`, `reports`, `starindex`

No genome is assumed — `genome_index` and `genome` must be set explicitly or the pipeline refuses to run.

---

### `prepare.py`
**Preparation.** Creates output directories, optionally copies raw FASTQs, and downloads/extracts the STAR index if needed.

Functions:
- `setup_dirs(cfg)` — creates all directories from `cfg.output_dirs()`
- `copy_raw(cfg)` — copies files matching `cfg.copy_src` into `RawData/`
- `ensure_starindex(cfg)` — downloads and extracts the STAR index from `cfg.starindex_url` only when `cfg.aligner == "star"` and the index is missing. No-ops silently for HISAT2.
- `prepare(cfg)` — calls the above three in order

---

### `trim.py`
**Step 1 — Trimming.** csRNA/sRNA samples are trimmed single-end with `homerTools trim`; totalRNA samples are trimmed paired-end with `skewer`.

Functions:
- `trim_one(cfg, r1)` — trims one sample; detects library type from filename
- `run_trim(cfg, sample_index)` — trims all (or one) samples; safe for SLURM array use

Output goes to `Trimmed/`. Already-trimmed files are skipped.

---

### `mapping.py`
**Step 2 — Alignment.** Aligns trimmed reads with STAR (default) or HISAT2. Both tools write uncompressed SAM files named `<prefix>.Aligned.out.sam` into `Aligned/`.

Functions:
- `_star_cmd(cfg, reads_in, out_prefix)` — builds the STAR command string
- `_hisat2_cmd(cfg, reads_flag, out_sam)` — builds the HISAT2 command string; writes mapping stats to `QC/<sample>_mappingstats.txt`
- `map_one(cfg, r1)` — aligns one sample; handles SE (csRNA/sRNA) and PE (totalRNA)
- `run_mapping(cfg, sample_index)` — aligns all (or one) samples

---

### `tagdirs.py`
**Step 3 — HOMER tag directories.** Merges all replicates for each sample into a single `<sample>-combo` tag directory using `makeTagDirectory`.

Functions:
- `run_tagdirs(cfg)` — finds all `*[_-]r1*.Aligned.out.sam` files, groups by sample prefix, and builds one combo tag directory per sample

Library type is detected from the sample name and the appropriate `makeTagDirectory` flags are applied (`-omitSN` for csRNA/sRNA, `-read2` for totalRNA).

---

### `bedgraphs.py`
**Step 4 — Genome-browser tracks.** Generates strand-specific uncompressed bedGraph files from each `*-combo` tag directory using `makeUCSCfile`.

Functions:
- `run_bedgraphs(cfg)` — iterates all `*-combo` dirs and writes `.posStrand.bedGraph` and `.negStrand.bedGraph` into `bedGraphs/`

Style is `tss` for csRNA/sRNA and `rnaseq` for totalRNA.

---

### `tss.py`
**Step 5 — TSS calling.** Runs `findcsRNATSS.pl` for each csRNA-combo tag directory, using the matched sRNA-combo as the input control and the totalRNA-combo (if present) as the `-rna` reference.

Functions:
- `run_tss(cfg)` — finds all `*_csRNA-combo` dirs, pairs them with their sRNA and optional totalRNA counterparts, and calls `findcsRNATSS.pl`

Outputs per sample in `TSS/`: `<sample>.tss.txt`, `<sample>.alltss.txt`, `<sample>.stats.txt`, `<sample>.inputDistribution.txt`, `<sample>.freq.tsv`.

Total RNA is optional — TSS are still called without it; the `Stable/Unstable` column will be absent.

---

### `qc.py`
**Step 6 — Quality control plots.** Generates a comprehensive set of PNG plots in `QC/` from tag directory internals and TSS output files.

| Output file | What it shows |
|---|---|
| `tagdir_stats.png` | Key stats table from each tag directory's `tagInfo.txt` |
| `read_length_distribution.png` | Read length distribution across all combo tag dirs |
| `nucleotide_frequency.png` | Per-nucleotide frequency from `tagFreqUniq.txt` |
| `autocorrelation.png` | Strand autocorrelation from `tagAutocorrelation.txt` |
| `median_tags_per_position.png` | Median tags per position for csRNA, sRNA, totalRNA (should be ≈ 1) |
| `csRNA_tagsPer_Vs_FracofPos.png` | Log-log tags vs fraction of positions — csRNA |
| `sRNA_tagsPer_Vs_FracofPos.png` | Log-log tags vs fraction of positions — sRNA |
| `totalRNA_tagsPer_Vs_FracofPos.png` | Log-log tags vs fraction of positions — totalRNA |
| `combined_tagsPer_Vs_FracofPos.png` | csRNA + sRNA overlaid |
| `csRNA_Aplot.png` | A-frequency near TSS — csRNA |
| `sRNA_Aplot.png` | A-frequency near TSS — sRNA |
| `totalRNA_Aplot.png` | A-frequency near TSS — totalRNA |
| `combined_Aplot.png` | csRNA + sRNA overlaid |
| `threshold_optimization.png` | TSS score threshold optimization from `*.inputDistribution.txt` |
| `tss_nucleotide_frequency.png` | Nucleotide frequency at primary TSS from `*.freq.tsv` |
| `tsr_summary.png` | Summary table parsed from `*.stats.txt` |
| `tsr_annotation.png` | Stacked bar of TSR annotation categories |

All plots use `bbox_inches="tight"` so titles are never clipped.

---

### `stability.py`
**Step 7 — TSR characterization.** Classifies TSRs as stable/unstable (requires totalRNA) and proximal/distal from the columns in `*.tss.txt` files. Outputs summary PNGs to `QC/`.

| Output file | What it shows |
|---|---|
| `stability_by_location_stacked_bar.png` | Stable vs unstable, split by proximal/distal (when totalRNA present) |
| `location_stacked_bar.png` | Proximal vs distal per sample (when no totalRNA) |
| `tsr_pie.png` | Pooled pie chart of stable/unstable or distal/proximal |

Column detection is flexible — the module searches by column name first, then by value patterns, then falls back to a numeric totalRNA ratio threshold.

---

### `report.py`
**Step 8 — HTML reports.** Writes four self-contained HTML files into `Reports/`. All images are base64-embedded so the files open anywhere without a server.

| Output file | Source |
|---|---|
| `tss_clusters.html` | `TSS/*.tss.txt` — filtered TSS clusters, scrollable tables with colour-coded annotation |
| `alltss.html` | `TSS/*.alltss.txt` — all pre-filter TSS candidates |
| `stats.html` | `TSS/*.stats.txt` — run statistics rendered as formatted text |
| `qc_report.html` | `QC/*.png` (embedded) + `QC/*.txt/.tsv/.csv`, images ordered by pipeline stage |

---

### `utils.py`
**Shared utilities** used by every other module.

| Function | Purpose |
|---|---|
| `setup_logging(cfg)` | Configures logging to stdout and a timestamped log file in `logs/` |
| `run(cmd, label, check)` | Runs a shell command, logs output, raises on failure |
| `done(path)` | Returns `True` if a file is non-empty or a directory exists (skip logic) |
| `seq_type(name)` | Classifies a filename as `"csRNA"`, `"sRNA"`, `"totalRNA"`, or `None` |
| `list_r1(cfg)` | Returns sorted list of `*_R1*.fastq[.gz]` files in `RawData/` |
| `check_tools(required, optional)` | Logs tool availability; returns list of missing required tools |
