"""Step 2 — Mapping (alignment). STAR (default) or HISAT2, per R1 file.

The aligner is chosen by cfg.aligner ("star" | "hisat2"). Both tools write a
standard, UNCOMPRESSED SAM named "<prefix>.Aligned.out.sam" into each
sample's own Aligned/ (next to its Trimmed/ and RawData/), so the tagdir
step and single-command re-runs work the same regardless of aligner.

cfg.genome_index must be set explicitly (no genome is assumed):
  - STAR   → used as --genomeDir (a directory)
  - HISAT2 → used as -x (an index prefix, e.g. /path/idx/genome)
"""
from __future__ import annotations
import os
from pathlib import Path
from .utils import run, log, seq_type, done, list_r1, leaf_dir


def _check_index(cfg) -> None:
    """Fail fast, with a clear message, on a missing OR wrong-shaped index —
    rather than letting STAR/HISAT2 die deep inside a SLURM array task and
    leave nothing behind but a log. In particular this catches the easy
    mistake of pointing --genome-index at a STAR --genomeDir while running
    with aligner=hisat2 (or vice versa): the two tools' index formats are
    NOT interchangeable (STAR wants a directory; HISAT2 wants a file-prefix
    with matching *.ht2 files next to it), and HISAT2 in particular fails
    near-instantly on a directory path with no SAM/alignment output at all,
    which looks identical to "alignment silently did nothing".
    """
    if cfg.aligner not in ("star", "hisat2"):
        raise ValueError(f"aligner must be 'star' or 'hisat2', got '{cfg.aligner}'")
    if not cfg.genome_index:
        raise ValueError(
            "genome_index is not set. Set CSRNA_GENOME_INDEX to your "
            + ("STAR --genomeDir directory" if cfg.aligner == "star"
               else "HISAT2 -x index prefix") + " (no genome is assumed)."
        )

    idx = Path(cfg.genome_index)

    if cfg.aligner == "star":
        # STAR --genomeDir must be an existing directory; SAindex is the
        # standard marker file every `STAR --runMode genomeGenerate` writes.
        if not idx.is_dir():
            raise ValueError(
                f"aligner is 'star' but genome_index '{cfg.genome_index}' is not a "
                f"directory. STAR needs a pre-built --genomeDir (from `STAR "
                f"--runMode genomeGenerate`), not a raw genome FASTA and not a "
                f"HISAT2 index prefix."
            )
        if not (idx / "SAindex").exists():
            raise ValueError(
                f"aligner is 'star' but '{cfg.genome_index}' doesn't look like a "
                f"STAR genomeDir (no SAindex file inside it). Point "
                f"--genome-index/CSRNA_GENOME_INDEX at a directory built with "
                f"`STAR --runMode genomeGenerate`."
            )
    else:  # hisat2
        # HISAT2 -x is a PREFIX, not a directory: the real index files are
        # <prefix>.1.ht2 ... <prefix>.8.ht2 (or .ht2l for large genomes).
        if idx.is_dir():
            raise ValueError(
                f"aligner is 'hisat2' but genome_index '{cfg.genome_index}' is a "
                f"directory. HISAT2 needs an index PREFIX (e.g. "
                f"/path/to/index/genome), built with `hisat2-build <genome.fa> "
                f"<prefix>`, pointing at files like <prefix>.1.ht2 — not a STAR "
                f"--genomeDir."
            )
        ht2_files = list(idx.parent.glob(f"{idx.name}.*.ht2")) + \
                    list(idx.parent.glob(f"{idx.name}.*.ht2l"))
        if not ht2_files:
            raise ValueError(
                f"aligner is 'hisat2' but no '{idx.name}.*.ht2[l]' files were found "
                f"next to '{cfg.genome_index}'. Build one with `hisat2-build "
                f"<genome.fa> {cfg.genome_index}`, then point "
                f"--genome-index/CSRNA_GENOME_INDEX at that same prefix."
            )


def _star_cmd(cfg, reads_in: str, out_prefix: str) -> str:
    # STAR default --outSAMtype is uncompressed SAM → <out_prefix>.Aligned.out.sam
    return (f"STAR --genomeDir {cfg.genome_index} "
            f"--readFilesIn {reads_in} "
            f"--outFileNamePrefix {out_prefix}. "
            f"--runThreadN {cfg.threads} "
            f"--outSAMstrandField intronMotif "
            f"--outMultimapperOrder Random "
            f"--outSAMmultNmax 1 "
            f"--outFilterMultimapNmax 10000 "
            f"--limitOutSAMoneReadBytes 10000000")


def _hisat2_cmd(cfg, reads_flag: str, out_sam: str, aligned_dir) -> str:
    # derive sample name from the output sam filename; mapping stats land
    # in this run's own Aligned/ dir (a future commit that nests qc.py per
    # Species/Sample can pool these there instead, same as it already
    # pools TSS/stability output).
    sample = os.path.basename(out_sam).replace(".Aligned.out.sam", "")
    stats = aligned_dir / f"{sample}_mappingstats.txt"
    return (f"hisat2 -p {cfg.threads} --rna-strandness F --dta "
            f"-x {cfg.genome_index} "
            f"{reads_flag} -S {out_sam} 2> {stats}")

def map_one(cfg, r1) -> None:
    st = seq_type(r1.name)
    trimmed_dir = leaf_dir(r1) / "Trimmed"
    aligned_dir = leaf_dir(r1) / "Aligned"
    aligned_dir.mkdir(parents=True, exist_ok=True)
    if st in ("csRNA", "sRNA"):                              # single-end
        trimmed = trimmed_dir / f"{r1.name}.trimmed"
        if not trimmed.exists():
            log.warning("mapping: trimmed file missing for %s (%s)", r1.name, trimmed.name); return
        prefix = r1.name.split("_R1")[0]
        sam = aligned_dir / f"{prefix}.Aligned.out.sam"
        if done(sam):
            log.info("  skip (done): %s", sam.name); return
        if cfg.aligner == "hisat2":
            run(_hisat2_cmd(cfg, f"-U {trimmed}", str(sam), aligned_dir), label=f"HISAT2 SE {r1.name}")
        else:
            run(_star_cmd(cfg, str(trimmed), str(aligned_dir / prefix)), label=f"STAR SE {r1.name}")
    elif st == "totalRNA":                                   # paired-end
        base = r1.name.split("_R1")[0]
        p1 = trimmed_dir / f"{base}-trimmed-pair1.fastq"
        p2 = trimmed_dir / f"{base}-trimmed-pair2.fastq"
        if not p1.exists():
            log.warning("mapping: trimmed pair missing for %s (%s)", base, p1.name); return
        sam = aligned_dir / f"{base}.Aligned.out.sam"
        if done(sam):
            log.info("  skip (done): %s", sam.name); return
        if cfg.aligner == "hisat2":
            reads = f"-1 {p1} -2 {p2}" if p2.exists() else f"-U {p1}"
            run(_hisat2_cmd(cfg, reads, str(sam), aligned_dir), label=f"HISAT2 PE {base}")
        else:
            reads = f"{p1} {p2}" if p2.exists() else f"{p1}"
            run(_star_cmd(cfg, reads, str(aligned_dir / base)), label=f"STAR PE {base}")
    else:
        log.warning("mapping: skipping untyped file %s", r1.name)


def run_mapping(cfg, sample_index=None) -> None:
    _check_index(cfg)
    log.info("mapping with %s (index: %s)", cfg.aligner.upper(), cfg.genome_index)
    r1s = list_r1(cfg)
    if not r1s:
        log.info("mapping: no *_R1*.fastq[.gz] under nested RawData/ dirs in %s", cfg.project); return
    if sample_index is not None:
        if not (0 <= sample_index < len(r1s)):
            raise IndexError(f"sample_index {sample_index} out of range (0-{len(r1s)-1})")
        r1s = [r1s[sample_index]]
    for r1 in r1s:
        map_one(cfg, r1)