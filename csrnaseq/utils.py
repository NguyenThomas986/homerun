"""Shared utilities: logging setup, command runner, helpers."""
from __future__ import annotations

import datetime
import logging
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger("csrnaseq")


def setup_logging(cfg) -> None:
    """Log to stdout (→ SLURM .out) and to a logfile under the project."""
    cfg.logs_dir.mkdir(parents=True, exist_ok=True)
    if not cfg.log_path:
        cfg.log_path = str(cfg.logs_dir / f"csrna_{datetime.datetime.now():%Y%m%d_%H%M%S}.log")

    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s  %(levelname)-7s  %(message)s", "%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(); sh.setFormatter(fmt); log.addHandler(sh)
    fh = logging.FileHandler(cfg.log_path); fh.setFormatter(fmt); log.addHandler(fh)
    log.info("Logging to %s", cfg.log_path)


def run(cmd: str, label: str = "", check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command, log output, and raise on failure (fail-fast)."""
    log.info("\u25b6 %s", label or cmd)
    job = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if job.stdout and job.stdout.strip():
        log.info(job.stdout.rstrip())
    if job.returncode != 0:
        log.error("FAILED (%d): %s", job.returncode, cmd)
        if job.stderr and job.stderr.strip():
            log.error(job.stderr.rstrip())
        if check:
            raise RuntimeError(f"Command failed ({job.returncode}): {cmd}")
    elif job.stderr and job.stderr.strip():
        log.debug(job.stderr.rstrip())
    return job


def done(path) -> bool:
    """True if an output already exists (dir present, or file non-empty)."""
    p = Path(path)
    return p.is_dir() or (p.is_file() and p.stat().st_size > 0)


def seq_type(name: str) -> str | None:
    """Classify a filename/sample by library tag. Handles _RNA and _totalRNA."""
    if "_csRNA" in name:
        return "csRNA"
    if "_sRNA" in name:
        return "sRNA"
    if "_totalRNA" in name or "_RNA" in name:
        return "totalRNA"
    return None


def list_r1(cfg):
    """Sorted list of R1 FASTQs in RawData — the unit of array parallelism.
    The array task index maps 1:1 to this (deterministic) ordering."""
    raw = cfg.rawdata
    return sorted(p for p in raw.glob("*_R1*")
                  if p.name.endswith(".fastq") or p.name.endswith(".fastq.gz"))


def check_tools(required=(), optional=()) -> list:
    """Log tool availability; return list of missing required tools."""
    log.info("Tool availability:")
    missing = []
    for t in required:
        w = shutil.which(t)
        log.info("  %-22s %s", t, w or "NOT FOUND")
        if not w:
            missing.append(t)
    for t in optional:
        w = shutil.which(t)
        log.info("  %-22s %s", t, (w or "not found") + " (optional)")
    return missing
