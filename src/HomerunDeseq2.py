"""
HomerunDeseq2.py
────────────────
Downstream analysis / DESeq2 step for the CSV-driven Homerun pipeline.

Calls getDiffExpression.pl (HOMER) on the merged TSR counts file.
Returns CSV column updates on success; raises on failure.
"""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, Any

import pandas as pd

log = logging.getLogger(__name__)


# ── Public entry point ────────────────────────────────────────────────────────
def run_deseq2(row: Dict[str, Any], working_path: str, cpus: int = 1) -> Dict[str, Any]:
    """
    Run DESeq2 differential expression analysis for the sample set.

    Expects:
        working_path/analysis/peakCalling/keyFiles/allTSSmerged_anoRaw.txt

    Returns dict with at minimum:
        {"deseq2": <output_path>, "status": "RUNNING"}
    """
    sample = row["sample"]
    genome = row["genome"]

    analysis_dir = Path(working_path) / "analysis" / "DESeq2"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    raw_tsr_path = (
        Path(working_path)
        / "analysis"
        / "peakCalling"
        / "keyFiles"
        / "allTSSmerged_anoRaw.txt"
    )

    if not raw_tsr_path.exists():
        raise FileNotFoundError(
            f"DESeq2 input not found: {raw_tsr_path}\n"
            "Run the 'stat' step first to generate merged TSR counts."
        )

    # Read raw TSR file and derive sample conditions
    raw_tsrs = pd.read_csv(raw_tsr_path, sep="\t", dtype=str)
    raw_tsrs.columns = [
        c.split("_")[1] if "Total" in c else c
        for c in raw_tsrs.columns
    ]

    # Extract condition columns (columns 19 onward, 0-indexed)
    condition_cols = list(raw_tsrs.columns[19:29]) if len(raw_tsrs.columns) > 19 else []
    if not condition_cols:
        log.warning("[%s] DESeq2: no condition columns found beyond index 19.", sample)
        condition_cols = list(raw_tsrs.columns[1:])   # fallback: use all but ID

    conditions_str = " ".join(condition_cols)
    output_prefix  = str(analysis_dir / "ALL_conditions")
    output_tsv     = output_prefix + "_DeSeq2.tsv"

    cmd = (
        f"getDiffExpression.pl {raw_tsr_path} "
        f"{conditions_str} "
        f"-export ALL_conditions -AvsA "
        f"> {output_tsv}"
    )
    log.info("[%s] running DESeq2 → %s", sample, output_tsv)
    log.debug("DESEQ2: %s", cmd)

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                            cwd=str(analysis_dir))
    if result.returncode != 0:
        raise RuntimeError(f"getDiffExpression.pl failed:\n{result.stderr}")

    if not Path(output_tsv).exists():
        raise RuntimeError(f"DESeq2 finished but output not found: {output_tsv}")

    log.info("[%s] DESeq2 complete → %s", sample, output_tsv)
    return {"deseq2": output_tsv, "status": "RUNNING"}
