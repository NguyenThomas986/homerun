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
  --genome mm10
```

<!-- TODO: confirm exact script name and flag names against `submit_array.sh --help` -->

## Run without SLURM

If you do not have access to a SLURM cluster, run HomeRun directly:

```bash
homerun --project /path/to/project --genome mm10
```

Or run it as a Python module:

```bash
python -m homerun --project /path/to/project --genome mm10
```

<!-- TODO: confirm module name (homerun vs csrnaseq) and the required flags -->

## Something went wrong?

See [Report an Issue](issue.md) for instructions on reporting bugs or unexpected output.