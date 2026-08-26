# Quick Start

## Run on SLURM (recommended)

Submit the pipeline as SLURM job arrays:

```bash
path/to/submit_array.sh \
  --project /path/to/project \
  --partition (parition_Name) \
  --conda-env (conda_ENV_Name) \
  --genome-index /path/to/STARindex \
  --genome mm10
```
<!-- TODO: confirm exact script name and flag names against `submit_array.sh --help` -->

## Run directly

For a single node without job arrays:

```bash
homerun --project /path/to/project --genome mm10
```

or as a module:

```bash
python -m homerun --project /path/to/project --genome mm10
```
<!-- TODO: confirm module name (homerun vs csrnaseq) and the required flags -->

## Something went wrong?

See [Report an Issue](issue.md) for how to report bugs and unexpected output.
