# Input Files

HOMERun supports two additive input modes. If `--metadata` is omitted, the
existing filename parser is used exactly as before. If `--metadata` is given,
the CSV/XLSX values define the sample identity and filenames are used only to
match files and distinguish `R1` from `R2`.

## Mode 1: filename-based parsing

This is the default and remains fully backward compatible. Species, sample,
condition, assay, and replicate identity are parsed from each FASTQ filename.

### Filename pattern

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

Species must be present as the first two filename tokens in this mode. For
example, `K562_D2_csRNA_r1_R1.fastq.gz` would not identify
`homo_sapiens`; use metadata mode when filenames do not carry species.

## Mode 2: CSV or Excel metadata

Pass a `.csv` or `.xlsx` file with `--metadata`:

```bash
homerun \
  --project /data/my-project \
  --metadata /data/my-project/samples.xlsx \
  --genome-index /indexes/hg38 \
  --genome hg38
```

To create a starter CSV automatically, place the FASTQs in the project root
or existing `Species/RawData/` directories and run:

```bash
homerun --project /data/my-project --init-metadata
```

This creates `/data/my-project/samples.csv`, fills its `FASTQ` column with the
exact discovered filenames, and leaves the biological fields blank for the
user to complete. It never guesses biological metadata and refuses to
overwrite an existing `samples.csv`.

`CSRNA_METADATA` provides the same setting. An explicit `--metadata` value
overrides the environment variable. The default is unset.

When using the SLURM controller, pass the Python option after its `--`
separator so it is forwarded to every job:

```bash
submit_array.sh \
  --project /data/my-project \
  --partition compute \
  --conda-env homerun \
  --genome-index /indexes/hg38 \
  --genome hg38 \
  -- --metadata /data/my-project/samples.xlsx
```

Required columns are `FASTQ`, `Species`, `Assay`, and `Replicate`. `Sample`
and `Condition` are optional. Column names are case-insensitive and tolerate
spaces, underscores, or hyphens.

CSV example:

```csv
FASTQ,Species,Sample,Condition,Assay,Replicate
K562_D2_csRNA_r1_R1.fastq.gz,homo_sapiens,K562,D2,csRNA,r1
K562_D2_sRNA_r1_R1.fastq.gz,homo_sapiens,K562,D2,sRNA,r1
```

For a conventional total-RNA filename with no condition, leave the
`Condition` cell empty:

```csv
FASTQ,Species,Sample,Condition,Assay,Replicate
homo_sapiens_K562_RNA-r1_R1.fastq,homo_sapiens,K562,,RNA,r1
```

The equivalent Excel sheet looks like this:

| FASTQ | Species | Sample | Condition | Assay | Replicate |
| --- | --- | --- | --- | --- | --- |
| `K562_D2_csRNA_r1_R1.fastq.gz` | `homo_sapiens` | `K562` | `D2` | `csRNA` | `r1` |
| `K562_D2_sRNA_r1_R1.fastq.gz` | `homo_sapiens` | `K562` | `D2` | `sRNA` | `r1` |

For the first row, HOMERun produces
`("homo_sapiens", "K562", "D2_csRNA_r1")`. A blank `Condition` produces
`csRNA_r1`. A blank or omitted `Sample` reuses the normalized species name,
so a row for `apis_mellifera` produces sample `apis_mellifera`.

Metadata rules:

- `FASTQ` is the exact basename of a FASTQ HomeRun can discover. Do not put a
  directory path in the cell.
- Every discovered FASTQ must have a row, and every metadata row must match a
  discovered file. For paired-end data, include both `R1` and `R2` rows with
  the same biological fields.
- `Species` is required and normalized to lowercase.
- `Assay` accepts `csRNA`, `sRNA`, `totalRNA`, or `RNA` (case-insensitive).
- `Replicate` accepts lowercase or uppercase forms such as `r1`, `r2`,
  `rep1`, and `rep2`, and is normalized to lowercase.
- Duplicate FASTQ rows are rejected. Duplicate output identities are also
  rejected, except for the intentional `R1`/`R2` pair of one library.

When metadata is supplied, HomeRun does not fall back to filename parsing for
an unlisted FASTQ. Validation stops early with a message identifying missing
rows, missing files, invalid fields, or output-name collisions.

## FASTQ placement

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

Paired total-RNA mates are matched by species, sample, condition, assay, and
replicate identity, whether that identity came from the filename or metadata.
Their download accessions do not need to match.

!!! warning "Validate names before a long run"
    In filename mode, a filename without a recognized assay or lowercase
    replicate marker is rejected. In metadata mode, the manifest is validated
    strictly before staging. Run `homerun --project PROJECT --stage-raw`
    (adding `--metadata PATH` when applicable) and inspect the resulting
    `Species/RawData/` directories.
