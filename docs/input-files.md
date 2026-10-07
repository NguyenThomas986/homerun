# Input Files

HOMERun discovers samples from FASTQ filenames. It does not use a sample sheet.

## Filename pattern

Each filename must contain, in order:

1. a two-token species name;
2. an optional sample name;
3. optional condition tokens;
4. an assay token: `csRNA`, `sRNA`, `totalRNA`, or `RNA`;
5. a lowercase replicate marker such as `r1`, `r2`, or `rep1`;
6. optional lane, accession, and read metadata.

Hyphens and underscores are treated as separators. The Illumina read markers
`R1` and `R2` are uppercase and are not interpreted as replicate markers.

Examples:

| Filename | Species | Sample | Replicate identity |
| --- | --- | --- | --- |
| `homo_sapiens_K562_csRNA-r1_DB422_S1_R1_001.fastq.gz` | `homo_sapiens` | `K562` | `csRNA_r1` |
| `homo_sapiens_K562_D2_Post_csRNA_r2_R1.fastq.gz` | `homo_sapiens` | `K562` | `D2_Post_csRNA_r2` |
| `Apis_mellifera_csRNA_r1_R1.fastq.gz` | `apis_mellifera` | `apis_mellifera` | `csRNA_r1` |

When no distinct sample token exists, HOMERun reuses the species name as the
sample name.

## Placement

FASTQs may start in either location:

- directly in the project root; `--stage-raw` and `submit_array.sh` move them
  into the correct species directory;
- in `Species/RawData/`, ready for discovery.

The `--copy-src` option accepts a quoted shell glob and copies matching FASTQs
during preparation. Quote the glob so the controller, rather than your current
shell, receives it.

```bash
submit_array.sh \
  ... \
  --copy-src '/data/run42/*.fastq.gz'
```

Paired total-RNA mates are matched by parsed species, sample, condition, assay,
and replicate identity. Their download accessions do not need to match.

!!! warning "Validate names before a long run"
    A filename without a recognized assay or lowercase replicate marker is
    rejected. Run `homerun --project PROJECT --stage-raw` first and inspect the
    resulting `Species/RawData/` directories.
