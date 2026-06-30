"""Step 4 — Genome-browser bedGraphs (strand-specific) per tag dir.

Generates bedGraphs for both:
  • <sample>-combo      — merged replicates (existing behavior)
  • <sample>_r<N>        — individual replicate tag dirs

Writes plain UNCOMPRESSED .bedGraph files with makeUCSCfile -o (stable names,
no piping/gzip) so the output structure stays the same across runs and a single
re-run command keeps working.
"""
from __future__ import annotations
import re
from .utils import run, log, seq_type, done


def _tagdir_prefix(name: str) -> str:
    """Strip '-combo' or '_r<N>' suffix to get the bare sample name for seq_type()."""
    if name.endswith("-combo"):
        return name[: -len("-combo")]
    m = re.match(r"^(.*)_r\d+$", name)
    if m:
        return m.group(1)
    return name


def run_bedgraphs(cfg) -> None:
    combo_dirs = sorted(d for d in cfg.tagdirs.glob("*-combo") if d.is_dir())
    rep_dirs   = sorted(d for d in cfg.tagdirs.glob("*_r[0-9]*") if d.is_dir())
    all_dirs   = combo_dirs + rep_dirs
    if not all_dirs:
        log.info("bedGraph: no tag dirs in %s", cfg.tagdirs)
        return
    skip = f"-skipChr {cfg.skip_chr} " if cfg.skip_chr else ""
    for d in all_dirs:
        sample_prefix = _tagdir_prefix(d.name)
        style  = "rnaseq" if seq_type(sample_prefix) == "totalRNA" else "tss"
        pos = cfg.bedgraphs / f"{d.name}.posStrand.bedGraph"
        neg = cfg.bedgraphs / f"{d.name}.negStrand.bedGraph"
        if not done(pos):
            run(f"makeUCSCfile {d} -style {style} -strand + {skip}-o {pos}",
                label=f"bedGraph + {d.name}")
        if not done(neg):
            run(f"makeUCSCfile {d} -style {style} -strand - -neg {skip}-o {neg}",
                label=f"bedGraph - {d.name}")
