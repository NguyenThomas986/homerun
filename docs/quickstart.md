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

See [Report an Issue](issue.md) for how to report bugs and unexpected output.
