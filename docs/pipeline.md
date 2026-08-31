# Pipeline

HOMERun analyzes csRNA-seq data as a series of steps, each handled by an established bioinformatics tool and submitted as SLURM jobs.

## 1. Adapter trimming

Sequencing reads carry adapter sequence that must be removed, and csRNA/sRNA reads are size-selected (roughly 20–58 nt).

- **Tool:** HOMER `homerTools trim` (csRNA, sRNA); skewer for paired-end totalRNA
<!-- TODO: confirm trimming tool(s) and the exact size range -->

## 2. Alignment to the genome

Trimmed reads are mapped to a reference genome to determine where each read originated, including strand.

- **Tool:** STAR (default), HISAT2 (optional)
- **Output:** `Aligned/<sample>.Aligned.out.sam`

## 3. Tag directory creation

Alignments are converted into HOMER tag directories — position-sorted, depth-normalized read counts plus precomputed QC statistics.

- **Tool:** HOMER `makeTagDirectory`
- Built per replicate (leaf) and merged per assay (combo)

## 4. BedGraph generation

Genome-browser tracks are generated from the tag directories for visualizing signal along the genome.

- **Tool:** HOMER `makeUCSCfile`

## 5. TSS / TSR calling

Capped-RNA peaks are called to identify active transcription start regions, using sRNA as an input control and totalRNA as a stability reference.

- **Tool:** HOMER `findcsRNATSS.pl`
- **Output:** `Species/TSS/<sample>.tss.txt`

## 6. Quality control & reporting

Per-sample QC is written under `Species/QC/<sample>/`. Compact cross-sample A/C/G/T nucleotide-divergence heatmaps (PNG/SVG plus TSV matrices) are written directly under `Species/QC/` for csRNA and sRNA.
<!-- TODO: confirm report format (HTML / PDF) -->

## 7. Stability analysis

Called TSRs are classified as stable or unstable transcripts using the totalRNA reference (HOMER `Log2Ratio vs. stableRNA`).
<!-- TODO: the stability scatter plot was removed from QC output — confirm whether the stability step itself remains in the pipeline, and adjust or delete this section accordingly. -->
