# Quality Control

HOMERun writes detailed per-sample QC under `Species/QC/<sample>/` and
cross-condition plots directly under `Species/QC/`.

## Per-sample output

Available inputs determine which plots are produced. Outputs can include:

- trimming and alignment summaries;
- read-length, nucleotide-frequency, and autocorrelation plots;
- combined and per-replicate tag-directory statistics;
- TSR summary and annotation plots;
- promoter-proximal versus distal TSS counts;
- RIT/RIE metrics when a GTF is configured.

Missing optional inputs cause the relevant plot to be skipped with a log
message; they do not invalidate unrelated QC.

## Cross-condition nucleotide divergence

For each species, assay, and nucleotide, HOMERun compares condition-preserving
combo TagDirs against the global mean. It writes:

```text
Species/QC/csRNA_A_DivergentPlot.png
Species/QC/csRNA_A_DivergentPlot.svg
Species/QC/csRNA_A_Divergent_Data.tsv
```

Equivalent `C`, `G`, and `T` files are produced, along with sRNA files when
sRNA combos exist. Rows retain the sample and full condition-combo label, and
day labels such as `D2` and `D10` are sorted naturally.

## Workbook

The QC step also writes a species-level Excel workbook containing the
underlying summary tables, including detailed TSR statistics. Use the workbook
for filtering and downstream review; use the PNG/SVG files for quick visual
inspection and reports.
