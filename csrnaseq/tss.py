"""Step 5 — Find TSRs with findcsRNATSS.pl.

For each _csRNA-combo, uses the matched _sRNA-combo as input control and, if present,
the totalRNA combo (_totalRNA-combo or _RNA-combo) as the -rna reference. The -rna
reference is optional: without total RNA the TSRs are still called, just without
stable/unstable assignment.
"""
from __future__ import annotations

from .utils import run, log, done


def run_tss(cfg) -> None:
    cs = sorted(d for d in cfg.tagdirs.glob("*_csRNA-combo") if d.is_dir())
    if not cs:
        log.info("TSS: no *_csRNA-combo dirs in %s", cfg.tagdirs)
        return

    for d in cs:
        prefix   = d.name.split("_csRNA")[0]
        sRNA_dir = cfg.tagdirs / f"{prefix}_sRNA-combo"
        out      = cfg.tss / prefix
        if done(f"{out}.tss.txt"):
            log.info("  skip (done): %s.tss.txt", prefix); continue

        rna_dir = None
        for cand in (f"{prefix}_totalRNA-combo", f"{prefix}_RNA-combo"):
            p = cfg.tagdirs / cand
            if p.is_dir():
                rna_dir = p; break

        cmd = (f"findcsRNATSS.pl {d} -o {out} -genome {cfg.genome} "
               f"-ntagThreshold {cfg.ntag_threshold} -i {sRNA_dir}")
        if rna_dir:
            cmd += f" -rna {rna_dir}"
            log.info("  using total-RNA reference: %s", rna_dir.name)
        else:
            log.info("  no total-RNA reference for %s (stability will be unavailable)", prefix)
        run(cmd, label=f"findcsRNATSS {prefix}")
