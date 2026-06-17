"""Step 4 — Genome-browser bedGraphs (strand-specific) per -combo tag dir.

Writes plain UNCOMPRESSED .bedGraph files with makeUCSCfile -o (stable names,
no piping/gzip) so the output structure stays the same across runs and a single
re-run command keeps working.
"""
from __future__ import annotations
from .utils import run, log, seq_type, done


def run_bedgraphs(cfg) -> None:
    combos = sorted(d for d in cfg.tagdirs.glob("*-combo") if d.is_dir())
    if not combos:
        log.info("bedGraph: no *-combo tag dirs in %s", cfg.tagdirs)
        return
    skip = f"-skipChr {cfg.skip_chr} " if cfg.skip_chr else ""
    for d in combos:
        prefix = d.name.split("-combo")[0]
        style  = "rnaseq" if seq_type(d.name) == "totalRNA" else "tss"
        pos = cfg.bedgraphs / f"{prefix}.posStrand.bedGraph"
        neg = cfg.bedgraphs / f"{prefix}.negStrand.bedGraph"
        if not done(pos):
            run(f"makeUCSCfile {d} -style {style} -strand + {skip}-o {pos}",
                label=f"bedGraph + {prefix}")
        if not done(neg):
            run(f"makeUCSCfile {d} -style {style} -strand - -neg {skip}-o {neg}",
                label=f"bedGraph - {prefix}")
