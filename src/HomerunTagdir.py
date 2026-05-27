"""
HomerunTagdir.py
────────────────
Tag directory creation step for the CSV-driven Homerun pipeline.

Reads the SAM file from row["star"], builds sampleInfo.txt for this
sample, runs batchMakeTagDirectory.pl, and returns CSV updates.
"""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, Any

log = logging.getLogger(__name__)


# ── Public entry point ────────────────────────────────────────────────────────
def run_tagdir(row: Dict[str, Any], working_path: str, cpus: int = 1) -> Dict[str, Any]:
    """
    Create a HOMER tag directory for the sample described by *row*.

    Returns dict with at minimum:
        {"tagdir": <tagdir_path>, "status": "RUNNING"}
    """
    sample    = row["sample"]
    genome    = row["genome"]
    star_val  = row.get("star", "")

    species_path = Path(working_path) / "data" / genome
    tagdir_path  = species_path / "tagDirs"
    fastq_path   = species_path / "fastq"
    genome_fa    = str(Path(working_path) / "genomes" / genome / "*.fa")

    tagdir_path.mkdir(parents=True, exist_ok=True)

    # Locate SAM file(s)
    sam_file = _find_sam(star_val, str(fastq_path), sample)
    if not sam_file:
        raise FileNotFoundError(
            f"No SAM file found for sample {sample} "
            f"(star='{star_val}', fastq_dir='{fastq_path}')"
        )

    # Derive tag directory name from SAM filename
    sam_stem = Path(sam_file).stem          # e.g. "THP1_rep1_csRNA-r1"
    tagdir_name = _derive_tagdir_name(sam_stem)
    tagdir_full = str(tagdir_path / tagdir_name)

    # Write a minimal sampleInfo.txt for this one sample
    sample_info = fastq_path / f"sampleInfo_{sample}.txt"
    with open(sample_info, "w") as f:
        f.write(f"{tagdir_full}\t{sam_file}\n")

    # Build tag directory
    cmd = (
        f"batchMakeTagDirectory.pl {sample_info} "
        f"-cpu {cpus} "
        f"-genome {genome_fa} "
        f"-omitSN -checkGC -fragLength 150 -r"
    )
    log.info("[%s] building tag directory → %s", sample, tagdir_full)
    log.debug("TAGDIR: %s", cmd)

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                            cwd=str(tagdir_path))
    if result.returncode != 0:
        raise RuntimeError(
            f"batchMakeTagDirectory.pl failed for {sample}:\n{result.stderr}"
        )

    # Verify output
    if not (Path(tagdir_full) / "tagInfo.txt").exists():
        raise RuntimeError(
            f"Tag directory created but tagInfo.txt missing at {tagdir_full}"
        )

    log.info("[%s] tag directory complete → %s", sample, tagdir_full)
    return {"tagdir": tagdir_full, "status": "RUNNING"}


# ── Helpers ───────────────────────────────────────────────────────────────────
def _find_sam(star_val: str, fastq_dir: str, sample: str) -> str | None:
    """Resolve the SAM file path from the star column value."""
    import glob

    if os.path.isfile(star_val) and star_val.endswith(".sam"):
        return star_val

    # Search fastq_dir for a SAM matching the sample name
    candidates = glob.glob(os.path.join(fastq_dir, f"*{sample}*.sam"))
    if not candidates:
        candidates = glob.glob(os.path.join(fastq_dir, "*.sam"))
    return candidates[0] if candidates else None


def _derive_tagdir_name(sam_stem: str) -> str:
    """
    Convert a SAM filename stem to the expected tagdir naming convention.

    e.g. "THP1_rep1_csRNA-r1.Aligned.out" → "THP1_rep1_csRNA-r1"
         "THP1_rep1_csRNA-r1"             → "THP1_rep1_csRNA-r1"
    """
    # Strip STAR output suffixes
    name = sam_stem
    for suffix in (".Aligned.out", ".Aligned", ".sorted"):
        name = name.replace(suffix, "")
    return name
