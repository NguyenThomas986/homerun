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

def read_homer_table(path):
    """Read any HOMER TSV output (peak files, mergePeaks output, annotatePeaks.pl
    output, tag directory stat files...) into a DataFrame.

    Some HOMER tools (e.g. findcsRNATSS.pl's .tss.txt, mergePeaks,
    annotatePeaks.pl) write just ONE '#'-prefixed header line. Others (e.g.
    findPeaks) write a whole metadata preamble first — parameters, genome,
    tag directory, thresholds — as several more '#' lines, with the real
    tab-separated column header further down. Reading naively with
    pd.read_csv(sep='\t') treats whichever line comes first as the header;
    for a multi-line preamble that's metadata text, not real columns, so it
    correctly detects 1 field there and then fails ('Expected 1 fields...')
    the moment it hits a real data row with many tab-separated fields.

    This scans every leading '#'-prefixed line and treats the LAST one as
    the actual column header (stripping its '#'), skipping everything above
    it — correct for both single-header-line files (where that's just the
    first and only line) and multi-line-preamble files alike.
    """
    import pandas as pd
    header_line, header_idx = None, -1
    with open(path) as fh:
        for i, line in enumerate(fh):
            if line.startswith("#"):
                header_line, header_idx = line, i
            else:
                break
    if header_line is None:
        # No '#'-prefixed line at all — fall back to plain read (unexpected,
        # but don't crash on a file that happens to have no HOMER header).
        return pd.read_csv(path, sep="\t")
    columns = header_line.lstrip("#").rstrip("\n").split("\t")
    return pd.read_csv(path, sep="\t", skiprows=header_idx + 1, names=columns)

def done(path) -> bool:
    """True if an output already exists (dir present, or file non-empty)."""
    p = Path(path)
    return p.is_dir() or (p.is_file() and p.stat().st_size > 0)


def _find_assay(tokens: list[str]) -> tuple[int | None, str | None]:
    """Search a list of tokens for an assay-type token, returning
    (index, assay_name) for the first match, or (None, None) if none found.

    Shared by parse_sample_name() and assay_of_leaf() so both agree on what
    counts as an assay token — checked by *position-independent* search
    (any token, not just the first/last), since a condition like 'p53KO'
    may sit before the assay token in a real filename (e.g.
    'p53KO_csRNA_r1'). Specific forms (csRNA/sRNA/totalRNA) are checked
    before the generic 'RNA' fallback so e.g. a 'totalRNA' token isn't
    mistaken for anything else.
    """
    for i, t in enumerate(tokens):
        low = t.lower()
        if low.startswith("csrna"):
            return i, "csRNA"
        if low.startswith("srna"):
            return i, "sRNA"
        if low.startswith("totalrna"):
            return i, "totalRNA"
    for i, t in enumerate(tokens):
        if t.lower().startswith("rna"):
            return i, "totalRNA"
    return None, None


def parse_sample_name(filename: str) -> tuple[str, str, str]:
    """Parse 'homo_sapiens_K562_csRNA_r1...' -> ('homo_sapiens', 'K562', 'csRNA_r1').
    Also handles a condition token before the assay, e.g.
    'homo_sapiens_K562_p53KO_csRNA_r1...' -> ('homo_sapiens', 'K562', 'p53KO_csRNA_r1').

    Splits on both '-' and '_' (so '-r2' and '_r2' are equivalent), then
    truncates at the first replicate marker (r1, r2, rep1...). Everything
    after the marker (lane/index/sequencer metadata) is discarded. Raises
    ValueError if no replicate marker is found, or if no assay-type token
    (csRNA/sRNA/totalRNA/RNA) is found between the sample and the replicate
    marker — the assay token's *position* isn't assumed (a condition like
    'p53KO' may come before it), only that one exists somewhere in there,
    via the same search _find_assay() uses for assay_of_leaf().
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
    _, assay = _find_assay(relevant[3:])
    if assay is None:
        raise ValueError(
            f"parse_sample_name: no assay type (csRNA/sRNA/totalRNA) found in '{filename}'"
        )
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
    """Classify a leaf directory name like 'csRNA_r1', 'csRNAseq_r2', or
    'p53KO_csRNA_r1' into csRNA / sRNA / totalRNA. Unlike seq_type(), this
    does NOT require a leading underscore before the assay token, since leaf
    names (produced by parse_sample_name) already have the species/sample
    prefix stripped off — e.g. 'csRNA_r1' rather than '..._csRNA_r1'.

    The assay token is searched for anywhere in the (underscore-split) leaf
    name via _find_assay(), not assumed to be the first token — a condition
    like 'p53KO' may legitimately sit in front of it (e.g. 'p53KO_csRNA_r1'),
    and a startswith()-on-the-whole-string check would silently miss those
    and return None, causing real samples to be skipped as unclassifiable.
    """
    base = re.sub(r"_r(ep)?\d+$", "", leaf_name)
    _, assay = _find_assay(base.split("_"))
    return assay


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


def list_samples(cfg):
    """Sorted, materialized (species, sample) pairs — the unit of array
    parallelism for GROUP_STEPS (tagdirs-combo, bedgraphs, tss). A SLURM
    array task's --group-index N maps 1:1 onto list_samples(cfg)[N], so this
    must be deterministic across calls; sorted() guarantees that regardless
    of filesystem iteration order."""
    return sorted(set(iter_samples(cfg)))


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