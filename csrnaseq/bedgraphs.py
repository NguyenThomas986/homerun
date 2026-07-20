"""Step 4 — Genome-browser bedGraphs (strand-specific), nested next to each
TagDir under Species/Sample/<assay>/TagDirs/ instead of a flat
project-level bedGraphs/.

Generates a bedGraph folder next to every TagDir built by tagdirs.py:
  • Species/Sample/<assay>/TagDirs/<assay>-combo -> Species/Sample/<assay>/bedGraphs/<assay>-combo/
  • Species/Sample/<assay>/TagDirs/<leaf_name>   -> Species/Sample/<assay>/bedGraphs/<leaf_name>/

Writes plain UNCOMPRESSED .bedGraph files with makeUCSCfile -o (stable names,
no piping/gzip) so a re-run command keeps working.
"""
from __future__ import annotations
from .utils import run, log, done


def run_bedgraphs(cfg, group=None) -> None:
    """Array-capable via --group-index (group=(species, sample) restricts to
    just that one Species/Sample's TagDirs — both its leaf and combo TagDirs,
    since this needs whichever of each already exist), or all Species/Sample
    at once when group=None."""
    # Species/Sample/<assay>/TagDirs/<leaf_or_combo>/ — the assay is the
    # directory's own name two levels up, no classification needed (unlike
    # the old flat layout, where a folder like '<assay_rep>' or
    # '<assay>-combo' had to be parsed to recover the assay).
    tagdirs = sorted(p for p in cfg.project.glob("*/*/*/TagDirs/*") if p.is_dir())
    if group is not None:
        sp, sa = group
        tagdirs = [td for td in tagdirs
                   if td.parent.parent.name == sa and td.parent.parent.parent.name == sp]
    if not tagdirs:
        log.info("bedGraph: no TagDirs/* under %s", cfg.project)
        return

    skip = f"-skipChr {cfg.skip_chr} " if cfg.skip_chr else ""
    for td in tagdirs:
        assay_dir = td.parent.parent                     # Species/Sample/<assay>/
        sample = assay_dir.parent.name
        species = assay_dir.parent.parent.name
        assay = assay_dir.name
        species_sample_run = f"{species}/{sample}/{assay}/{td.name}"

        bedgraph_dir = assay_dir / "bedGraphs" / td.name
        bedgraph_dir.mkdir(parents=True, exist_ok=True)
        style = "rnaseq" if assay == "totalRNA" else "tss"
        pos = bedgraph_dir / "posStrand.bedGraph"
        neg = bedgraph_dir / "negStrand.bedGraph"

        if not done(pos):
            run(f"makeUCSCfile {td} -style {style} -strand + {skip}-o {pos}",
                label=f"bedGraph + {species_sample_run}")
        if not done(neg):
            run(f"makeUCSCfile {td} -style {style} -strand - -neg {skip}-o {neg}",
                label=f"bedGraph - {species_sample_run}")