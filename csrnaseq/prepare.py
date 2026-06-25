"""Preparation: create folders, stage loose FASTQs, copy raw FASTQs, ensure STARIndex exists."""
from __future__ import annotations
import shutil
from .utils import run, log

def setup_dirs(cfg) -> None:
    for d in cfg.output_dirs():
        existed = d.is_dir()
        d.mkdir(parents=True, exist_ok=True)
        log.info("  %s  %s", "exists " if existed else "CREATED", d)

def copy_raw(cfg) -> None:
    if cfg.copy_src:
        run(f"cp -r {cfg.copy_src} {cfg.rawdata}/", label="copy raw FASTQs")
    else:
        log.info("copy_src empty — skipping raw copy.")

def stage_loose_fastqs(cfg) -> None:
    """Move loose *_R1*/*_R2* FASTQs sitting in the project ROOT into RawData/.

    Non-recursive (only the project root is scanned, never subdirs), so files
    already in RawData/ are untouched. If a same-named file already exists in
    RawData/, the loose copy is LEFT IN PLACE (never clobbered) and a warning is
    logged. Safe to call repeatedly — a no-op once everything is staged.
    """
    cfg.rawdata.mkdir(parents=True, exist_ok=True)
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
        dst = cfg.rawdata / src.name
        if dst.exists():
            log.warning("stage: %s already exists in RawData/ — leaving loose file in place.",
                        src.name)
            continue
        shutil.move(str(src), str(dst))
        log.info("stage: moved %s -> RawData/", src.name)

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

def prepare(cfg) -> None:
    log.info("=== PREPARE: folders / stage loose / raw copy / STARIndex ===")
    setup_dirs(cfg)
    stage_loose_fastqs(cfg)
    copy_raw(cfg)
    if cfg.aligner == "star":
        ensure_starindex(cfg)
    else:
        log.info("Aligner is %s — skipping STARIndex.", cfg.aligner)