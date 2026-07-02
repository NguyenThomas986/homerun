"""Shared utilities: logging setup, command runner, helpers."""
from __future__ import annotations

import datetime
import logging
import re
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


def run(cmd: str, label: str = "", check: bool = True, cwd=None) -> subprocess.CompletedProcess:
    """Run a shell command, log output, and raise on failure (fail-fast).

    cwd: working directory for the command. Use this for tools (like
    findcsRNATSS.pl) that write temp files to the current directory, so
    those temp files land somewhere the user actually has write access to
    (e.g. the project's TSS/ dir) instead of wherever the job happened to
    start from.
    """
    log.info("▶ %s", label or cmd)
    job = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
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


def parse_sample_name(filename: str) -> tuple[str, str, str]:
    """Parse 'homo_sapiens_K562_csRNA_r1...' -> ('homo_sapiens', 'K562', 'csRNA_r1').

    Splits on both '-' and '_' (so '-r2' and '_r2' are equivalent), then
    truncates at the first replicate marker (r1, r2, rep1...). Everything
    after the marker (lane/index/sequencer metadata) is discarded. Raises
    ValueError if no replicate marker is found.
    """
    stem = filename.split(".")[0]
    tokens = re.split(r"[-_]", stem)
    # Only lowercase 'r1'/'rep2' etc. count as a replicate marker. Illumina
    # read-tags (R1/R2) are conventionally uppercase and must NOT match here,
    # or a filename with no real replicate marker before its R1/R2 read tag
    # would silently misparse (see test_illumina_r1_tag_is_not_a_replicate).
    rep_idx = next(
        (i for i, t in enumerate(tokens) if re.fullmatch(r"r(ep)?\d+", t)),
        None,
    )
    if rep_idx is None:
        raise ValueError(f"parse_sample_name: no replicate marker (r1, rep2...) in '{filename}'")
    relevant = tokens[: rep_idx + 1]
    if len(relevant) < 3:
        raise ValueError(f"parse_sample_name: not enough tokens before replicate marker in '{filename}'")
    species = f"{relevant[0]}_{relevant[1]}".lower()
    sample = relevant[2]
    leaf_name = "_".join(relevant[3:]) or relevant[-1]
    return species, sample, leaf_name


def seq_type(name: str) -> str | None:
    """Classify a filename/sample by library tag. Handles _RNA and _totalRNA."""
    if "_csRNA" in name:
        return "csRNA"
    if "_sRNA" in name:
        return "sRNA"
    if "_totalRNA" in name or "_RNA" in name:
        return "totalRNA"
    return None


def leaf_dir(r1: Path) -> Path:
    """Given an R1 fastq Path at .../Species/Sample/<leaf>/RawData/<file>,
    return the <leaf> directory itself (.../Species/Sample/<leaf>/).

    Used by trim/mapping/tagdirs so each sample's Trimmed/Aligned/TagDir/
    bedGraph outputs land next to its own RawData/, without re-parsing the
    filename again — the directory structure already encodes it.
    """
    return r1.parent.parent


def list_r1(cfg):
    """Sorted list of R1 FASTQs under every nested Species/Sample/<leaf>/RawData/ —
    the unit of array parallelism. The array task index maps 1:1 to this
    (deterministic) ordering."""
    return sorted(
        p for p in cfg.project.glob("*/*/*/RawData/*_R1*")
        if p.name.endswith(".fastq") or p.name.endswith(".fastq.gz")
    )


def assay_of_leaf(leaf_name: str) -> str | None:
    """Classify a leaf directory name like 'csRNA_r1' or 'csRNAseq_r2' into
    csRNA / sRNA / totalRNA. Unlike seq_type(), this does NOT require a
    leading underscore before the assay token, since leaf names (produced by
    parse_sample_name) already have the species/sample prefix stripped off —
    e.g. 'csRNA_r1' rather than '..._csRNA_r1'.
    """
    base = re.sub(r"_r(ep)?\d+$", "", leaf_name)
    low = base.lower()
    if low.startswith("csrna"):
        return "csRNA"
    if low.startswith("srna"):
        return "sRNA"
    if low.startswith("totalrna") or low.startswith("rna"):
        return "totalRNA"
    return None


def replicate_of_leaf(leaf_name: str) -> str | None:
    """Extract the replicate marker ('r1', 'r2', ...) from a leaf dir name."""
    m = re.search(r"_(r(ep)?\d+)$", leaf_name)
    return m.group(1) if m else None


def iter_leaf_dirs(cfg):
    """Yield (species, sample, leaf_dir) for every nested
    Species/Sample/<leaf>/ directory under the project that has a RawData/
    subdir — the same set of runs list_r1() draws its R1 files from, but
    exposed as directories rather than individual FASTQs (used by steps that
    key off the leaf dir itself, e.g. tagdirs/bedgraphs)."""
    for rawdata in sorted(cfg.project.glob("*/*/*/RawData")):
        if not rawdata.is_dir():
            continue
        leaf = rawdata.parent
        sample = leaf.parent
        species = sample.parent
        yield species.name, sample.name, leaf


def iter_samples(cfg):
    """Yield (species, sample) pairs for every nested Species/Sample/ that
    has at least one leaf run dir — the unit of iteration for sample-level
    outputs (QC/, TSS/)."""
    seen = set()
    for species, sample, _leaf in iter_leaf_dirs(cfg):
        key = (species, sample)
        if key not in seen:
            seen.add(key)
            yield key


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
        log.info("  %-22s %s", t, w or "not found (optional)")
    return missing