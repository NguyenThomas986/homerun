"""
HomerunTrim.py
──────────────
Trim step for the CSV-driven Homerun pipeline.

Wraps `homerTools trim` (csRNA / sRNA) and `skewer` (totalRNA / paired-end).
Returns a dict of CSV column updates on success; raises on failure.
"""

from __future__ import annotations

import glob
import logging
import multiprocessing
import os
import subprocess
from pathlib import Path
from typing import Dict, Any

log = logging.getLogger(__name__)


# ── Public entry point ────────────────────────────────────────────────────────
def run_trim(row: Dict[str, Any], working_path: str, cpus: int = 1) -> Dict[str, Any]:
    """
    Trim the FASTQ(s) described by *row*.

    Returns a dict with at minimum:
        {"trim": <output_path_or_DONE>, "status": "RUNNING"}
    """
    sample    = row["sample"]
    fastq_src = row["fastq"]        # file or directory
    genome    = row["genome"]
    seq_type  = _infer_seq_type(sample)

    output_dir = Path(working_path) / "data" / genome / "fastq" / seq_type
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect FASTQ files
    fastqs = _collect_fastqs(fastq_src)
    if not fastqs:
        raise FileNotFoundError(f"No FASTQ files found at: {fastq_src}")

    log.info("[%s] trimming %d FASTQ(s) in %s mode (type=%s)",
             sample, len(fastqs), "homer" if seq_type in ("csRNA", "sRNA") else "skewer",
             seq_type)

    if seq_type in ("csRNA", "sRNA"):
        _homer_trim(fastqs, cpus=cpus)
        # homerTools trim writes .trimmed files alongside the input — look there
        src_dir = fastq_src if os.path.isdir(fastq_src) else os.path.dirname(fastq_src)
        trimmed = _collect_trimmed(src_dir, seq_type)
    else:
        _skewer_trim(fastqs, output_prefix=str(output_dir / sample), cpus=cpus)
        trimmed = _collect_trimmed(str(output_dir), seq_type)

    if not trimmed:
        raise RuntimeError(f"Trimming produced no output files for sample {sample}")

    log.info("[%s] trim complete — %d trimmed file(s).", sample, len(trimmed))
    return {
        "trim":   trimmed[0] if len(trimmed) == 1 else str(output_dir),
        "status": "RUNNING",
    }


# ── Homer trim (csRNA / sRNA) ─────────────────────────────────────────────────
def _homer_trim(fastqs: list, cpus: int) -> None:
    if cpus > 1:
        pool = multiprocessing.Pool(min(cpus, len(fastqs)))
        pool.map(_homer_trim_one, fastqs)
        pool.close()
        pool.join()
    else:
        for fq in fastqs:
            _homer_trim_one(fq)


def _homer_trim_one(fastq_path: str) -> None:
    cmd = f"homerTools trim -mis 2 -minMatchLength 4 -min 20 {fastq_path}"
    log.debug("TRIM: %s", cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"homerTools trim failed for {fastq_path}:\n{result.stderr}"
        )


# ── Skewer trim (totalRNA / paired-end) ───────────────────────────────────────
def _skewer_trim(fastqs: list, output_prefix: str, cpus: int) -> None:
    r1_files = [f for f in fastqs if "_R1" in f]
    r2_files = [f for f in fastqs if "_R2" in f]

    cmds = []
    for r1 in r1_files:
        prefix = output_prefix + "_" + Path(r1).stem.split("_R1")[0]
        r2 = r1.replace("_R1", "_R2")
        if r2 in r2_files:
            cmds.append(f"skewer -m mp {r1} {r2} -t {cpus} -o {prefix}")
        else:
            log.warning("No R2 found for %s — trimming single-end.", r1)
            cmds.append(f"skewer -m any {r1} -t {cpus} -o {prefix}")

    for cmd in cmds:
        log.debug("TRIM: %s", cmd)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"skewer failed:\n{result.stderr}")


# ── Helpers ───────────────────────────────────────────────────────────────────
def _infer_seq_type(sample: str) -> str:
    """Infer csRNA / sRNA / totalRNA from sample name."""
    if "csRNA" in sample:
        return "csRNA"
    if "_sRNA" in sample:
        return "sRNA"
    return "totalRNA"


def _collect_fastqs(src: str) -> list:
    """Return list of FASTQ paths from a file or directory."""
    p = Path(src)
    if p.is_file():
        return [str(p)]
    if p.is_dir():
        found = (
            glob.glob(str(p / "*.fastq.gz"))
            + glob.glob(str(p / "*.fastq"))
            + glob.glob(str(p / "*.fq.gz"))
        )
        return sorted(found)
    return []


def _collect_trimmed(directory: str, seq_type: str) -> list:
    """Collect trimmed output files from a directory."""
    if seq_type in ("csRNA", "sRNA"):
        return sorted(glob.glob(os.path.join(directory, "*.trimmed")))
    return sorted(
        glob.glob(os.path.join(directory, "*-trimmed*.fastq"))
        + glob.glob(os.path.join(directory, "*-trimmed*.fastq.gz"))
    )
