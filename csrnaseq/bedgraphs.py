"""Step 4 — Genome-browser bedGraphs (strand-specific), nested next to each
TagDir under Species/Sample/ instead of a flat project-level bedGraphs/.

Generates a bedGraph/ folder next to every TagDir built by tagdirs.py:
  • Species/Sample/<assay>-combo/TagDir  -> Species/Sample/<assay>-combo/bedGraph/
  • Species/Sample/<assay_rep>/TagDir    -> Species/Sample/<assay_rep>/bedGraph/

Writes plain UNCOMPRESSED .bedGraph files with makeUCSCfile -o (stable names,
no piping/gzip) so a re-run command keeps working.
"""
from __future__ import annotations
from .utils import run, log, done, assay_of_leaf


def _assay_for(run_dir_name: str) -> str | None:
    """Classify the folder a TagDir sits in ('<assay>-combo' or '<assay_rep>')."""
    if run_dir_name.endswith("-combo"):
        return run_dir_name[: -len("-combo")]
    return assay_of_leaf(run_dir_name)


def run_bedgraphs(cfg, group=None) -> None:
    """Array-capable via --group-index (group=(species, sample) restricts to
    just that one Species/Sample's TagDirs — both its leaf and combo TagDirs,
    since this needs whichever of each already exist), or all Species/Sample
    at once when group=None."""
    tagdirs = sorted(cfg.project.glob("*/*/*/TagDir"))
    if group is not None:
        sp, sa = group
        tagdirs = [td for td in tagdirs
                   if td.parent.parent.name == sa and td.parent.parent.parent.name == sp]
    if not tagdirs:
        log.info("bedGraph: no nested TagDir/ dirs under %s", cfg.project)
        return

    skip = f"-skipChr {cfg.skip_chr} " if cfg.skip_chr else ""
    for td in tagdirs:
        run_dir = td.parent                      # Species/Sample/<assay>-combo or <assay_rep>
        species_sample_run = f"{run_dir.parent.parent.name}/{run_dir.parent.name}/{run_dir.name}"
        assay = _assay_for(run_dir.name)
        if not assay:
            log.warning("bedGraph: could not classify assay for %s", species_sample_run)
            continue

        bedgraph_dir = run_dir / "bedGraph"
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