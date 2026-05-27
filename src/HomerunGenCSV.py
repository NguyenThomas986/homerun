"""
HomerunGenCSV.py
────────────────
Auto-generate a pipeline CSV from a directory of FASTQ files.

Naming convention expected:
    <sample>_<type>_R1.fastq.gz    (type: csRNA, sRNA, totalRNA / RNA)

Each unique sample name (everything before _csRNA / _sRNA / _RNA / _totalRNA)
becomes one row in the CSV, with all step columns set to PENDING.

Usage (CLI):
    python Homerun.py --generate-csv /path/to/fastqs --genome Homo_sapiens --out pipeline.csv

Usage (API):
    from HomerunGenCSV import generate_csv
    generate_csv("/data/fastqs", genome="Homo_sapiens", output_path="pipeline.csv")
"""

from __future__ import annotations

import glob
import logging
import os
import re
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

STEP_COLUMNS     = ["trim", "star", "tagdir", "qc", "deseq2"]
REQUIRED_COLUMNS = ["sample", "fastq", "genome"] + STEP_COLUMNS + ["status", "notes"]

# Regex to extract sample name from typical Homerun FASTQ filenames
_SAMPLE_RE = re.compile(
    r"^(?P<sample>.+?)(?:_csRNA|_sRNA|_totalRNA|_RNA).*?(?:_R[12])?\.fastq(?:\.gz)?$",
    re.IGNORECASE,
)


def generate_csv(fastq_dir: str, genome: str, output_path: str) -> pd.DataFrame:
    """
    Scan *fastq_dir* for FASTQ files and write a pipeline CSV to *output_path*.

    Returns the generated DataFrame.
    """
    fastq_dir = Path(fastq_dir).resolve()
    if not fastq_dir.is_dir():
        raise NotADirectoryError(f"FASTQ directory not found: {fastq_dir}")

    # Gather all FASTQ files (recursively one level deep)
    fastq_files = sorted(
        glob.glob(str(fastq_dir / "*.fastq.gz"))
        + glob.glob(str(fastq_dir / "*.fastq"))
        + glob.glob(str(fastq_dir / "**" / "*.fastq.gz"))
        + glob.glob(str(fastq_dir / "**" / "*.fastq"))
    )
    fastq_files = list(dict.fromkeys(fastq_files))  # deduplicate, preserve order

    if not fastq_files:
        log.warning("No FASTQ files found in: %s", fastq_dir)

    # Map sample name → first R1 FASTQ file found
    samples: dict[str, str] = {}
    unmatched: list[str] = []

    for fpath in fastq_files:
        fname = Path(fpath).name
        m = _SAMPLE_RE.match(fname)
        if m:
            sample = m.group("sample")
            # Use R1 as the canonical fastq path; fall back to whatever we find
            if "_R2" in fname:
                continue   # skip R2; R1 is the driver
            if sample not in samples:
                samples[sample] = fpath
        else:
            # Fallback: use the full stem as the sample name
            sample = fname.split(".fastq")[0]
            if "_R2" in sample:
                continue
            if sample not in samples:
                samples[sample] = fpath
                unmatched.append(fname)

    if unmatched:
        log.warning(
            "Could not parse sample names from %d file(s) — used full stem: %s",
            len(unmatched), ", ".join(unmatched[:5]),
        )

    rows = []
    for sample, fastq_path in sorted(samples.items()):
        rows.append({
            "sample":  sample,
            "fastq":   fastq_path,
            "genome":  genome,
            "trim":    "PENDING",
            "star":    "PENDING",
            "tagdir":  "PENDING",
            "qc":      "PENDING",
            "deseq2":  "PENDING",
            "status":  "NOT_STARTED",
            "notes":   "",
        })

    df = pd.DataFrame(rows, columns=REQUIRED_COLUMNS) if rows else \
         pd.DataFrame(columns=REQUIRED_COLUMNS)

    df.to_csv(output_path, index=False)
    log.info("Generated CSV with %d sample(s) → %s", len(df), output_path)
    return df
