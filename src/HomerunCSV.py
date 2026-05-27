"""
HomerunCSV.py
─────────────
Thread- and process-safe CSV read/write layer.
The CSV is the sole source of truth for all pipeline state.

Public API
----------
PipelineCSV(path)           – load (or create) the CSV
csv.dataframe               – current pandas DataFrame
csv.update_row(sample, row) – persist one updated row to disk
csv.reload()                – re-read from disk (e.g. after SLURM steps)

CSV Schema
----------
sample, fastq, genome, trim, star, tagdir, qc, deseq2, status, notes
"""

from __future__ import annotations

import fcntl
import logging
import os
import time
from pathlib import Path
from typing import Dict, Any

import pandas as pd

log = logging.getLogger(__name__)

# ── Required columns (in display order) ──────────────────────────────────────
REQUIRED_COLUMNS = [
    "sample", "fastq", "genome",
    "trim", "star", "tagdir", "qc", "deseq2",
    "status", "notes",
]

STEP_COLUMNS   = ["trim", "star", "tagdir", "qc", "deseq2"]
STATUS_OPTIONS = {"NOT_STARTED", "RUNNING", "FAILED", "DONE"}
PENDING        = "PENDING"


class PipelineCSV:
    """Manages the pipeline CSV with safe concurrent writes."""

    def __init__(self, path: str | Path):
        self.path = Path(path).resolve()
        self._df: pd.DataFrame | None = None
        self.reload()

    # ── Load ──────────────────────────────────────────────────────────────────
    def reload(self) -> None:
        """Read CSV from disk (or create skeleton if file is new)."""
        if not self.path.exists():
            log.warning("CSV not found — creating empty skeleton at %s", self.path)
            self._df = _empty_dataframe()
            self._write(self._df)
        else:
            raw = pd.read_csv(self.path, dtype=str).fillna("")
            self._df = _normalise(raw)
            log.debug("Loaded %d rows from %s", len(self._df), self.path)

    @property
    def dataframe(self) -> pd.DataFrame:
        return self._df.copy()

    # ── Update ────────────────────────────────────────────────────────────────
    def update_row(self, sample: str, updated: Dict[str, Any]) -> None:
        """
        Persist one updated sample row to disk.

        Steps:
        1. Acquire an exclusive file lock (fcntl — works on Linux/SLURM).
        2. Re-read CSV to pick up any concurrent changes.
        3. Merge the updated dict into the matching row.
        4. Write back.
        5. Release lock.
        """
        retries = 5
        for attempt in range(retries):
            try:
                lock_path = str(self.path) + ".lock"
                with open(lock_path, "w") as lf:
                    try:
                        fcntl.flock(lf, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    except BlockingIOError:
                        log.debug("CSV locked, retry %d/%d…", attempt + 1, retries)
                        time.sleep(0.5 * (attempt + 1))
                        continue

                    # Re-read from disk (another process may have written)
                    current = pd.read_csv(self.path, dtype=str).fillna("")
                    current = _normalise(current)

                    mask = current["sample"] == sample
                    if not mask.any():
                        log.warning("Sample '%s' not found in CSV — appending.", sample)
                        new_row = _row_defaults()
                        new_row.update({k: str(v) for k, v in updated.items()})
                        new_row["sample"] = sample
                        current = pd.concat(
                            [current, pd.DataFrame([new_row])], ignore_index=True
                        )
                    else:
                        for key, val in updated.items():
                            if key in current.columns:
                                current.loc[mask, key] = str(val) if val is not None else ""

                    self._write(current)
                    self._df = current
                    log.debug("Row updated: %s", sample)
                    return  # success

            except Exception as exc:
                log.error("Failed to update CSV (attempt %d): %s", attempt + 1, exc)
                time.sleep(1)

        log.error("Gave up updating CSV row for sample '%s' after %d attempts.", sample, retries)

    # ── Internal write ────────────────────────────────────────────────────────
    def _write(self, df: pd.DataFrame) -> None:
        """Write DataFrame to CSV atomically via a temp file."""
        tmp = str(self.path) + ".tmp"
        df[REQUIRED_COLUMNS].to_csv(tmp, index=False)
        os.replace(tmp, self.path)   # atomic on POSIX


# ── Helpers ───────────────────────────────────────────────────────────────────
def _empty_dataframe() -> pd.DataFrame:
    return pd.DataFrame(columns=REQUIRED_COLUMNS)


def _row_defaults() -> Dict[str, str]:
    return {col: "" for col in REQUIRED_COLUMNS}


def _normalise(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure all required columns exist; add missing ones as empty strings."""
    for col in REQUIRED_COLUMNS:
        if col not in df.columns:
            log.warning("Column '%s' missing from CSV — adding empty column.", col)
            df[col] = ""
    # Keep only known columns + any extras the user added, but put required first
    extra = [c for c in df.columns if c not in REQUIRED_COLUMNS]
    return df[REQUIRED_COLUMNS + extra].copy()
