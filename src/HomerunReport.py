"""
HomerunReport.py
────────────────
Summary and reporting utilities for the CSV-driven Homerun pipeline.

print_summary()  – human-readable console table at start/end of a run
write_report()   – write a TSV summary report to working_path/files/
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

STEP_COLUMNS = ["trim", "star", "tagdir", "qc", "deseq2"]


# ── Console summary ───────────────────────────────────────────────────────────
def print_summary(df: pd.DataFrame, label: str = "") -> None:
    """Print a compact pipeline status table to the logger."""
    total = len(df)
    if total == 0:
        log.info("[%s] CSV is empty — nothing to process.", label)
        return

    lines = [f"\n{'─' * 60}", f"  {label}  ({total} sample(s))", f"{'─' * 60}"]

    # Per-step counts
    for step in STEP_COLUMNS:
        if step not in df.columns:
            continue
        done    = (df[step].str.upper() == "DONE").sum() + \
                  (~df[step].isin(["", "PENDING", "FAILED", "DONE", None])).sum()
        pending = df[step].isin(["", "PENDING"]).sum()
        failed  = (df[step].str.upper() == "FAILED").sum()
        lines.append(
            f"  {step:<8}  done={done:<4} pending={pending:<4} failed={failed}"
        )

    # Global status counts
    lines.append(f"{'─' * 60}")
    for status in ["NOT_STARTED", "RUNNING", "DONE", "FAILED"]:
        n = (df.get("status", pd.Series()) == status).sum()
        if n:
            lines.append(f"  status={status:<12} {n} sample(s)")

    lines.append(f"{'─' * 60}\n")
    log.info("\n".join(lines))


# ── File report ───────────────────────────────────────────────────────────────
def write_report(df: pd.DataFrame, output_dir: str) -> None:
    """
    Write a tab-separated summary report.

    Columns: sample, status, per-step status, notes
    """
    report_dir = Path(output_dir) / "files"
    report_dir.mkdir(parents=True, exist_ok=True)
    report_path = report_dir / "homerun_run_summary.tsv"

    cols = ["sample"] + STEP_COLUMNS + ["status", "notes"]
    out  = df[[c for c in cols if c in df.columns]].copy()

    # Simplify step columns to DONE / PENDING / FAILED for readability
    for step in STEP_COLUMNS:
        if step not in out.columns:
            continue
        out[step] = out[step].apply(_simplify_status)

    out.to_csv(report_path, sep="\t", index=False)
    log.info("Run summary written to: %s", report_path)


def _simplify_status(val) -> str:
    v = str(val).strip().upper()
    if v in ("DONE", "PENDING", "FAILED", ""):
        return v if v else "PENDING"
    return "DONE"   # any non-empty, non-keyword value is an output path → DONE
