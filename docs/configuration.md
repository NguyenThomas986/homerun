# Configuration

Configuration values use this precedence order:

1. an explicit command-line flag;
2. a `CSRNA_*` environment variable;
3. the built-in default.

`--project` falls back to `CSRNA_PROJECT`, then the current directory. Threads
fall back to `CSRNA_THREADS`, then `SLURM_CPUS_PER_TASK`, then `20`.

## Required analysis values

The SLURM controller requires these values as flags:

| Flag | Environment variable | Meaning |
| --- | --- | --- |
| `--project` | `CSRNA_PROJECT` | Project root |
| `--genome-index` | `CSRNA_GENOME_INDEX` | STAR genome directory or HISAT2 index prefix |
| `--genome` | `CSRNA_GENOME` | HOMER genome name or genome FASTA path |
| `--aligner` | `CSRNA_ALIGNER` | `star` (default) or `hisat2` |

`submit_array.sh` additionally requires `--partition` and `--conda-env`; those
control SLURM and are not Python pipeline settings.

## Common optional values

| Flag | Environment variable | Default |
| --- | --- | --- |
| `--gtf` | `CSRNA_GTF` | unset; RIT/RIE is skipped |
| `--metadata` | `CSRNA_METADATA` | unset; parse identities from FASTQ filenames |
| `--threads` | `CSRNA_THREADS` | `SLURM_CPUS_PER_TASK` or `20` |
| `--trim-min` | `CSRNA_TRIM_MINLEN` | `20` |
| `--trim-max` | `CSRNA_TRIM_MAXLEN` | `58` |
| `--ntag-threshold` | `CSRNA_NTAG_THRESHOLD` | `7` |
| `--skip-chr` | `CSRNA_SKIP_CHR` | `chrEBV` |
| `--cleanup-intermediates` | `CSRNA_CLEANUP_INTERMEDIATES` | off |

See the [CLI reference](cli.md) for every flag and its meaning.

## Generated run record

Preparation writes `<project>/config.txt`. It records the effective settings,
the samples HOMERun discovered, and the staged FASTQ paths. Check this file
before troubleshooting later stages; it shows what the pipeline actually used.
When metadata mode is active, the resolved metadata path appears as
`metadata = /path/to/samples.xlsx`, and the sample list uses metadata-defined
species/sample identities.

!!! note "GTF is optional"
    `--gtf` enables the RIT/RIE metric. Omitting it is valid and causes that
    step to be skipped. Supplying an unreadable path is an error caught during
    preparation.
