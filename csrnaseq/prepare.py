"""Preparation: create folders, stage loose FASTQs, copy raw FASTQs, ensure STARIndex exists.

FASTQs always land under the flat Species/Sample/RawData/ layout (shared by
every assay and replicate of that sample — see Config.rawdata_dir /
utils.parse_sample_name) — there is no per-assay subfolder and no
project-level fallback.
"""
from __future__ import annotations
import dataclasses
import glob
import shutil
from datetime import datetime
from pathlib import Path
from .utils import run, log, parse_sample_name, list_samples

def setup_dirs(cfg) -> None:
    # Only the project-wide logs/ dir is created up front; per-sample
    # RawData/Trimmed/Aligned/TagDir/bedGraph/QC/TSS dirs are all created on
    # demand as each sample's files are discovered (their paths depend on
    # parsing the filename, which we don't know until we see it).
    d = cfg.logs_dir
    existed = d.is_dir()
    d.mkdir(parents=True, exist_ok=True)
    log.info("  %s  %s", "exists " if existed else "CREATED", d)

def _stage_one(cfg, src: Path) -> None:
    """Parse src's filename and move/copy it into the sample's shared RawData/
    dir (Species/Sample/RawData/ — shared across every replicate of every
    assay in that sample, not one folder per assay or per replicate; the
    filename itself still uniquely identifies both downstream)."""
    species, sample, _leaf = parse_sample_name(src.name)
    dst_dir = cfg.rawdata_dir(species, sample)
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / src.name
    if dst.exists():
        log.warning("stage: %s already exists at %s — leaving source in place.",
                    src.name, dst_dir)
        return dst
    shutil.move(str(src), str(dst)) if src.exists() else None
    log.info("stage: moved %s -> %s/", src.name, dst_dir)
    return dst

def copy_raw(cfg) -> None:
    """Copy FASTQs matched by cfg.copy_src directly into their nested RawData/ dirs.

    Each matched file is parsed individually (species/sample/leaf), unlike a
    flat `cp -r glob dest/`, since the destination now depends on the
    filename itself.
    """
    if not cfg.copy_src:
        log.info("copy_src empty — skipping raw copy.")
        return
    matches = sorted(Path(p) for p in glob.glob(cfg.copy_src) if Path(p).is_file())
    if not matches:
        log.warning("copy_src '%s' matched no files.", cfg.copy_src)
        return
    for src in matches:
        try:
            species, sample, _leaf = parse_sample_name(src.name)
        except ValueError as exc:
            log.warning("copy_raw: skipping %s (%s)", src.name, exc)
            continue
        dst_dir = cfg.rawdata_dir(species, sample)
        dst_dir.mkdir(parents=True, exist_ok=True)
        dst = dst_dir / src.name
        if dst.exists():
            log.warning("copy_raw: %s already exists at %s — skipping.", src.name, dst_dir)
            continue
        run(f"cp {src} {dst}", label=f"copy raw {src.name}")

def stage_loose_fastqs(cfg) -> None:
    """Move loose *_R1*/*_R2* FASTQs sitting in the project ROOT into their
    Species/Sample/RawData/ dir.

    Non-recursive (only the project root is scanned, never subdirs), so files
    already staged are untouched. If a same-named file already exists at the
    destination, the loose copy is LEFT IN PLACE (never clobbered) and a
    warning is logged. Filenames that don't parse (no replicate marker) are
    skipped with a warning rather than crashing the whole prepare step.
    Safe to call repeatedly — a no-op once everything is staged.
    """
    loose = sorted(
        p for p in cfg.project.glob("*")
        if p.is_file()
        and ("_R1" in p.name or "_R2" in p.name)
        and (p.name.endswith(".fastq") or p.name.endswith(".fastq.gz"))
    )
    if not loose:
        log.info("stage: no loose FASTQs in project root — nothing to move.")
        return
    for src in loose:
        try:
            _stage_one(cfg, src)
        except ValueError as exc:
            log.warning("stage: skipping %s (%s)", src.name, exc)

def ensure_starindex(cfg) -> None:
    if cfg.aligner != "star":
        log.info("ensure_starindex: aligner is '%s' — skipping.", cfg.aligner)
        return
    si = cfg.starindex
    if si.is_dir() and any(si.iterdir()):
        log.info("STARIndex present: %s", si)
        return
    if not cfg.starindex_url:
        raise ValueError(
            "STARIndex not found and CSRNA_STARINDEX_URL is not set. "
            "Either set CSRNA_GENOME_INDEX to an existing index, or set "
            "CSRNA_STARINDEX_URL to download it automatically."
        )
    tarball = cfg.project / "GSE287021_STARIndex_hg38.tar.gz"
    run(f"wget -O {tarball} '{cfg.starindex_url}'", label="download STARIndex")
    run(f"tar -xvzf {tarball} -C {cfg.project}", label="extract STARIndex")
    run(f"rm -f {tarball}", label="cleanup tarball")
    log.info("STARIndex extracted to %s", si)

def validate_gtf(cfg) -> None:
    """If --gtf/CSRNA_GTF is set, confirm it actually points at a real,
    readable file NOW (in the 'prepare' job — the first phase), rather than
    discovering a typo'd/missing path 3 jobs later when 'ritrie' (the very
    last collect step) finally tries to read it. If --gtf is unset, ritrie is
    simply skipped later — that's fine, not an error — so this only raises
    when a value WAS given but doesn't check out."""
    if not cfg.gtf:
        return
    p = Path(cfg.gtf)
    if not p.is_file():
        raise ValueError(
            f"--gtf/CSRNA_GTF is set to '{cfg.gtf}' but that file does not exist "
            f"(or isn't visible from this node). Double-check the exact path and "
            f"extension — e.g. with: ls -la {cfg.gtf}"
        )
    log.info("GTF found: %s", p)

def write_config_summary(cfg) -> None:
    """Write <project>/config.txt: every Config field (whatever was passed in,
    or its default if not) plus every sample discovered so far. Meant to make
    a project directory self-documenting — what genome/aligner/thresholds/gtf
    it was run with, and which samples it saw — without digging through
    config.env, CLI history, or logs. Safe to call repeatedly (overwrites);
    called at the end of every prepare() so it stays current as new samples
    get staged across separate submit_array.sh runs."""
    lines = [
        "# HomeRun csRNA-seq pipeline — run configuration",
        f"# Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "[Config]",
    ]
    for f in dataclasses.fields(cfg):
        lines.append(f"{f.name} = {getattr(cfg, f.name)}")
 
    lines += ["", "[Samples]"]
    samples = list_samples(cfg)
    if not samples:
        lines.append("(none discovered yet)")
    else:
        lines.append(f"count = {len(samples)}")
        for species, sample in samples:
            lines.append(f"{species}/{sample}")
 
    lines += ["", "[RawData]"]
    raw_files = sorted(p for p in cfg.project.glob("*/*/RawData/*") if p.is_file())
    if not raw_files:
        lines.append("(none staged yet)")
    else:
        lines.append(f"count = {len(raw_files)}")
        for p in raw_files:
            lines.append(str(p.relative_to(cfg.project)))
 
    out = cfg.project / "config.txt"
    out.write_text("\n".join(lines) + "\n")
    log.info("Wrote run configuration (%d sample(s), %d raw file(s)) to %s",
              len(samples), len(raw_files), out)

def prepare(cfg) -> None:
    log.info("=== PREPARE: folders / stage loose / raw copy / STARIndex ===")
    validate_gtf(cfg)
    setup_dirs(cfg)
    stage_loose_fastqs(cfg)
    copy_raw(cfg)
    if cfg.aligner == "star":
        ensure_starindex(cfg)
    else:
        log.info("Aligner is %s — skipping STARIndex.", cfg.aligner)
    write_config_summary(cfg)