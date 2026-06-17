"""Step 2 — Mapping (alignment). STAR (default) or HISAT2, per R1 file.

The aligner is chosen by cfg.aligner ("star" | "hisat2"). Both tools write a
standard, UNCOMPRESSED SAM named "<prefix>.Aligned.out.sam" into Aligned/, so
the tagdir step and single-command re-runs work the same regardless of aligner.

cfg.genome_index must be set explicitly (no genome is assumed):
  - STAR   → used as --genomeDir (a directory)
  - HISAT2 → used as -x (an index prefix, e.g. /path/idx/genome)
"""
from __future__ import annotations
import os
from .utils import run, log, seq_type, done, list_r1


def _check_index(cfg) -> None:
    if cfg.aligner not in ("star", "hisat2"):
        raise ValueError(f"aligner must be 'star' or 'hisat2', got '{cfg.aligner}'")
    if not cfg.genome_index:
        raise ValueError(
            "genome_index is not set. Set CSRNA_GENOME_INDEX to your "
            + ("STAR --genomeDir directory" if cfg.aligner == "star"
               else "HISAT2 -x index prefix") + " (no genome is assumed)."
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


def _hisat2_cmd(cfg, reads_flag: str, out_sam: str) -> str:
    # derive sample name from the output sam filename
    sample = os.path.basename(out_sam).replace(".Aligned.out.sam", "")
    stats = cfg.qc / f"{sample}_mappingstats.txt"
    return (f"hisat2 -p {cfg.threads} --rna-strandness F --dta "
            f"-x {cfg.genome_index} "
            f"{reads_flag} -S {out_sam} 2> {stats}")

def map_one(cfg, r1) -> None:
    st = seq_type(r1.name)
    if st in ("csRNA", "sRNA"):                              # single-end
        trimmed = cfg.trimmed / f"{r1.name}.trimmed"
        if not trimmed.exists():
            log.warning("mapping: trimmed file missing for %s (%s)", r1.name, trimmed.name); return
        prefix = r1.name.split("_R1")[0]
        sam = cfg.aligned / f"{prefix}.Aligned.out.sam"
        if done(sam):
            log.info("  skip (done): %s", sam.name); return
        if cfg.aligner == "hisat2":
            run(_hisat2_cmd(cfg, f"-U {trimmed}", str(sam)), label=f"HISAT2 SE {r1.name}")
        else:
            run(_star_cmd(cfg, str(trimmed), str(cfg.aligned / prefix)), label=f"STAR SE {r1.name}")
    elif st == "totalRNA":                                   # paired-end
        base = r1.name.split("_R1")[0]
        p1 = cfg.trimmed / f"{base}-trimmed-pair1.fastq"
        p2 = cfg.trimmed / f"{base}-trimmed-pair2.fastq"
        if not p1.exists():
            log.warning("mapping: trimmed pair missing for %s (%s)", base, p1.name); return
        sam = cfg.aligned / f"{base}.Aligned.out.sam"
        if done(sam):
            log.info("  skip (done): %s", sam.name); return
        if cfg.aligner == "hisat2":
            reads = f"-1 {p1} -2 {p2}" if p2.exists() else f"-U {p1}"
            run(_hisat2_cmd(cfg, reads, str(sam)), label=f"HISAT2 PE {base}")
        else:
            reads = f"{p1} {p2}" if p2.exists() else f"{p1}"
            run(_star_cmd(cfg, reads, str(cfg.aligned / base)), label=f"STAR PE {base}")
    else:
        log.warning("mapping: skipping untyped file %s", r1.name)


def run_mapping(cfg, sample_index=None) -> None:
    _check_index(cfg)
    log.info("mapping with %s (index: %s)", cfg.aligner.upper(), cfg.genome_index)
    r1s = list_r1(cfg)
    if not r1s:
        log.info("mapping: no *_R1*.fastq[.gz] in %s", cfg.rawdata); return
    if sample_index is not None:
        if not (0 <= sample_index < len(r1s)):
            raise IndexError(f"sample_index {sample_index} out of range (0-{len(r1s)-1})")
        r1s = [r1s[sample_index]]
    for r1 in r1s:
        map_one(cfg, r1)
