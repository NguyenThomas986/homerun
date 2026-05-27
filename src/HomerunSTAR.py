"""
HomerunSTAR.py
──────────────
STAR alignment step for the CSV-driven Homerun pipeline.

Reads the trimmed file(s) from row["trim"], runs STAR, moves the
output SAM to the species fastq directory, and returns CSV updates.
"""

from __future__ import annotations

import glob
import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, Any

log = logging.getLogger(__name__)


# ── Public entry point ────────────────────────────────────────────────────────
def run_star(row: Dict[str, Any], working_path: str, cpus: int = 1) -> Dict[str, Any]:
    """
    Run STAR alignment for the sample described by *row*.

    Returns dict with at minimum:
        {"star": <sam_output_path>, "status": "RUNNING"}
    """
    sample    = row["sample"]
    genome    = row["genome"]
    trim_val  = row.get("trim", "")
    seq_type  = _infer_seq_type(sample)

    star_index = Path(working_path) / "genomes" / genome / "STARIndex"
    fastq_dir  = Path(working_path) / "data" / genome / "fastq" / seq_type
    sam_dir    = Path(working_path) / "data" / genome / "fastq"
    mapping_dir = Path(working_path) / "files" / "mappingStats"
    mapping_dir.mkdir(parents=True, exist_ok=True)

    # Find trimmed files
    trimmed_files = _find_trimmed(trim_val, str(fastq_dir), seq_type)
    if not trimmed_files:
        raise FileNotFoundError(
            f"No trimmed FASTQ files found for {sample} "
            f"(trim='{trim_val}', dir='{fastq_dir}')"
        )

    # Determine R1 files to drive alignment (one job per R1)
    r1_files = [f for f in trimmed_files if "_R1" in f or len(trimmed_files) == 1]
    if not r1_files:
        r1_files = trimmed_files  # single-end: use all

    # Run STAR per R1 (sequential — STAR is RAM-hungry)
    sam_outputs = []
    for r1 in r1_files:
        sam_out = _run_star_one(
            r1_file=r1,
            trimmed_files=trimmed_files,
            star_index=str(star_index),
            sam_dir=str(sam_dir),
            mapping_dir=str(mapping_dir),
            sample=sample,
            cpus=cpus,
        )
        sam_outputs.append(sam_out)

    result_path = sam_outputs[0] if len(sam_outputs) == 1 else str(sam_dir)
    log.info("[%s] STAR alignment complete → %s", sample, result_path)

    return {"star": result_path, "status": "RUNNING"}


# ── Core STAR call ────────────────────────────────────────────────────────────
def _run_star_one(
    r1_file: str,
    trimmed_files: list,
    star_index: str,
    sam_dir: str,
    mapping_dir: str,
    sample: str,
    cpus: int,
) -> str:
    stem    = Path(r1_file).stem.split(".fastq")[0]
    out_pfx = str(Path(sam_dir) / stem) + "."
    r2_file = r1_file.replace("_R1", "_R2")

    if Path(r2_file).exists():
        reads_arg = f"--readFilesIn {r1_file} {r2_file}"
    else:
        reads_arg = f"--readFilesIn {r1_file}"

    # Use --readFilesCommand zcat for compressed inputs
    read_cmd = "--readFilesCommand zcat" if r1_file.endswith(".gz") else ""

    star_cmd = (
        f"STAR --genomeDir {star_index} "
        f"--runThreadN {cpus} "
        f"{reads_arg} "
        f"{read_cmd} "
        f"--outFileNamePrefix {out_pfx} "
        f"--genomeLoad NoSharedMemory "
        f"--outSAMstrandField intronMotif "
        f"--outMultimapperOrder Random "
        f"--outSAMmultNmax 1 "
        f"--outFilterMultimapNmax 10000 "
        f"--limitOutSAMoneReadBytes 10000000"
    )
    log.debug("STAR: %s", star_cmd)

    result = subprocess.run(star_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"STAR failed for {r1_file}:\n{result.stderr}")

    # Move mapping stats
    log_final = out_pfx + "Log.final.out"
    if Path(log_final).exists():
        dest = Path(mapping_dir) / f"{stem}.mappingstats.txt"
        Path(log_final).rename(dest)
        log.debug("Mapping stats → %s", dest)

    # Expected SAM output
    sam_path = out_pfx + "Aligned.out.sam"
    if not Path(sam_path).exists():
        raise FileNotFoundError(
            f"STAR finished but SAM not found at expected path: {sam_path}"
        )
    return sam_path


# ── Helpers ───────────────────────────────────────────────────────────────────
def _infer_seq_type(sample: str) -> str:
    if "csRNA" in sample:
        return "csRNA"
    if "_sRNA" in sample:
        return "sRNA"
    return "totalRNA"


def _find_trimmed(trim_val: str, fastq_dir: str, seq_type: str) -> list:
    """
    Locate trimmed FASTQ files.
    trim_val is either a file path, a directory, or a status string like DONE.
    """
    # If trim_val is an actual path (file or dir), use it
    if os.path.isfile(trim_val):
        return [trim_val]
    if os.path.isdir(trim_val):
        src = trim_val
    else:
        src = fastq_dir

    if seq_type in ("csRNA", "sRNA"):
        return sorted(glob.glob(os.path.join(src, "*.trimmed")))
    return sorted(
        glob.glob(os.path.join(src, "*-trimmed*.fastq"))
        + glob.glob(os.path.join(src, "*-trimmed*.fastq.gz"))
    )
