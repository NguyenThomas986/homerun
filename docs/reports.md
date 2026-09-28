# Reports and Outputs

HOMERun keeps analysis categories directly under each species directory.
Sample identity remains in filenames and TagDir names so multiple samples can
share the same category directories safely.

```text
project/
├── config.txt
├── logs/
├── logs_slurm/
└── homo_sapiens/
    ├── RawData/
    ├── Trimmed/
    ├── Aligned/
    ├── TagDirs/
    ├── bedGraphs/
    ├── TSS/
    ├── RITRIE/
    └── QC/
        └── K562/
```

## Naming

Individual TagDirs and bedGraph directories use
`<sample>_<condition>_<assay>_<replicate>`. Replicates that differ only by the
final `rN` marker are also merged into a condition-preserving
`<sample>_<condition>_<assay>-combo` directory.

For example, `K562_D2_Post_csRNA_r1` and `K562_D2_Post_csRNA_r2` produce
`K562_D2_Post_csRNA-combo`; they are not merged with a `D2_Pre` condition.

TSS calls are written for matched individual replicates and their matched
condition combos. csRNA is paired with the corresponding sRNA control and,
when available, the corresponding total-RNA reference.

## What to keep

Keep `config.txt`, TagDirs, bedGraphs, TSS results, and QC output with an
archived analysis. If `--cleanup-intermediates` is enabled, `Trimmed/` and
`Aligned/` data may be removed after QC has extracted its tables.
