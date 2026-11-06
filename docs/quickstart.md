# Quick Start

## Prepare your project

Copy your FASTQ files into your project folder.

Then run HomeRun using one of the options below.

## Run on SLURM (recommended)

If your system uses SLURM, submit the pipeline as SLURM job arrays:

```bash
path/to/submit_array.sh \
  --project /path/to/project \
  --partition partition_name \
  --conda-env conda_env_name \
  --genome-index /path/to/STARindex \
  --genome hg38
```

<!-- TODO: confirm exact script name and flag names against `submit_array.sh --help` -->

## Run without SLURM

If you do not have access to a SLURM cluster, run HOMERun directly:

```bash
homerun --project /path/to/project --genome hg38
```

Or run it as a Python module:

```bash
python -m homerun --project /path/to/project --genome hg38
```

In these examples, `hg38` is the genome identifier (homo_sapiens). Replace it with the genome appropriate for your dataset.

<!-- TODO: confirm module name (homerun vs csrnaseq) and the required flags -->

## Use a CSV sample sheet

Filename parsing remains the default. If your filenames do not contain the
species, sample, condition, assay, and replicate identity, generate a starter
CSV after placing the FASTQs in the project:

```bash
homerun --project /path/to/project --init-csv
```

HomeRun creates `/path/to/project/samples.csv` with the exact discovered FASTQ
names. It leaves the biological fields blank rather than guessing them. Edit
the file so each FASTQ has its required values:

```csv
FASTQ,Species,Sample,Condition,Assay,Replicate
homo_sapiens_K562_RNA-r1_R1.fastq,homo_sapiens,K562,,RNA,r1
```

Then run with that file:

```bash
homerun \
  --project /path/to/project \
  --csv /path/to/project/samples.csv \
  --genome-index /path/to/STARindex \
  --genome hg38
```

The same option supports `.xlsx` files. When using `submit_array.sh`, forward
it after `--`:

```bash
path/to/submit_array.sh \
  --project /path/to/project \
  --partition PARTITION_NAME \
  --conda-env CONDA_ENV_NAME \
  --genome-index /path/to/STARindex \
  --genome hg38 \
  -- --csv /path/to/project/samples.csv
```

See [Input files](input-files.md) for all columns and validation rules.

## Something went wrong?

See [Report an Issue](issue.md) for instructions on reporting bugs or unexpected output.