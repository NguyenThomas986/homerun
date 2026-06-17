"""Preparation: create folders, copy raw FASTQs, ensure STARIndex exists."""
from __future__ import annotations
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
    log.info("=== PREPARE: folders / raw copy / STARIndex ===")
    setup_dirs(cfg)
    copy_raw(cfg)
    if cfg.aligner == "star":
        ensure_starindex(cfg)
    else:
        log.info("Aligner is %s — skipping STARIndex.", cfg.aligner)