"""Step 1 — Trim. csRNA/sRNA (SE) via homerTools; totalRNA (PE) via skewer.

Operates per R1 file so it can run inside a SLURM array (one task per sample).
homerTools writes outputs next to the input (i.e. into the sample's own
nested RawData/), so each call moves ONLY its own sample's *.trimmed/*.lengths
into that same sample's Trimmed/ — safe under concurrency.
"""
from __future__ import annotations

import shutil

from .utils import run, log, seq_type, done, list_r1, leaf_dir


def trim_one(cfg, r1) -> None:
    st = seq_type(r1.name)
    trimmed_dir = leaf_dir(r1) / "Trimmed"
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    if st in ("csRNA", "sRNA"):                              # single-end
        out = trimmed_dir / f"{r1.name}.trimmed"
        if done(out):
            log.info("  skip (done): %s", r1.name); return
        run(f"homerTools trim -3 {cfg.trim_adapter} -mis {cfg.trim_mis} "
            f"-minMatchLength {cfg.trim_minmatch} -min {cfg.trim_min} "
            f"-max {cfg.trim_max} {r1}", label=f"trim SE {r1.name}")
        for suffix in (".trimmed", ".lengths"):              # move only THIS sample's outputs
            src = r1.parent / f"{r1.name}{suffix}"           # homerTools writes next to r1 (RawData/)
            if src.exists():
                shutil.move(str(src), str(trimmed_dir / src.name))
    elif st == "totalRNA":                                   # paired-end
        out_prefix = trimmed_dir / r1.name.split("_R1")[0]
        if done(f"{out_prefix}-trimmed-pair1.fastq"):
            log.info("  skip (done): %s", r1.name); return
        r2 = r1.parent / r1.name.replace("_R1", "_R2")        # r2 lives alongside r1 in the same RawData/
        run(f"skewer -m pe {r1} {r2} -t {cfg.threads} -o {out_prefix}",
            label=f"trim PE {r1.name}")
    else:
        log.warning("trim: skipping untyped file %s", r1.name)


def run_trim(cfg, sample_index=None) -> None:
    r1s = list_r1(cfg)
    if not r1s:
        log.info("trim: no *_R1*.fastq[.gz] under nested RawData/ dirs in %s", cfg.project); return
    if sample_index is not None:
        if not (0 <= sample_index < len(r1s)):
            raise IndexError(f"sample_index {sample_index} out of range (0-{len(r1s)-1})")
        r1s = [r1s[sample_index]]
    for r1 in r1s:
        trim_one(cfg, r1)