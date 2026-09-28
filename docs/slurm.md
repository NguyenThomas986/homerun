# Running with SLURM

`submit_array.sh` coordinates preparation and five worker stages. It activates
the requested conda environment and submits jobs from the repository checkout,
so an editable package installation is not required on each compute node.

## Submit a run

```bash
path/to/homerun/submit_array.sh \
  --project /path/to/project \
  --partition kamiak \
  --conda-env homerun \
  --genome-index /path/to/STARIndex \
  --genome hg38
```

Use `--` to pass additional Python pipeline flags:

```bash
path/to/homerun/submit_array.sh \
  --project /path/to/project \
  --partition kamiak \
  --conda-env homerun \
  --genome-index /path/to/STARIndex \
  --genome hg38 \
  -- --gtf /path/to/annotation.gtf --cleanup-intermediates
```

## Dependency graph

```text
prepare (waits before array sizing)
  └─ align_array[one task per R1]
       └─ tagdir_array[one task per R1]
            ├─ tss_array[one task per species/sample]
            └─ bedgraphs_array[one task per species/sample]
                 └─ collect (after both arrays succeed)
```

The tag-directory array waits for the complete alignment array because one
canonical task may combine SAMs produced by several lane or replicate tasks.
It uses an `afterany` scheduling dependency so the validation worker always
runs; that worker fails if any expected SAM is missing, preventing partial
condition combos from reaching TSS calling.

## Concurrency

- `--throttle` limits align, tag-directory, and bedGraph arrays (default `16`).
- `--tss-throttle` limits sample-level TSS array tasks (default `4`).

Logs are written to `<project>/logs_slurm/`. The final submission summary prints
every SLURM job ID.

## Safe reruns

Before submitting jobs, the controller stages loose FASTQs and runs the
read-only `--check-rerun` preflight. A completed project is rejected unless
`--force` (or `--overwrite`) is passed after `--`.
