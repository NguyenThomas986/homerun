"""Step 4 — Genome-browser bedGraphs (strand-specific), written into each
sample's own flat Species/Sample/bedGraphs/ next to Species/Sample/TagDirs/.

Generates a bedGraph folder next to every TagDir built by tagdirs.py:
  • Species/Sample/TagDirs/<assay>-combo -> Species/Sample/bedGraphs/<assay>-combo/
  • Species/Sample/TagDirs/<leaf_name>   -> Species/Sample/bedGraphs/<leaf_name>/

Writes plain UNCOMPRESSED .bedGraph files with makeUCSCfile -o (stable names,
no piping/gzip) so a re-run command keeps working.
"""
from __future__ import annotations
from .utils import run, log, done, assay_of_leaf


def _assay_of_tagdir(name: str) -> str | None:
    """Recover the assay from a TagDir's own name, since there's no longer a
    per-assay parent folder to read it off of. Combo dirs are literally
    '<assay>-combo'; leaf dirs (e.g. 'csRNA_r1') are classified the same way
    assay_of_leaf() already classifies any other leaf name."""
    if name.endswith("-combo"):
        return name[: -len("-combo")]
    return assay_of_leaf(name)


def run_bedgraphs(cfg, group=None) -> None:
    """Array-capable via --group-index (group=(species, sample) restricts to
    just that one Species/Sample's TagDirs — both its leaf and combo TagDirs,
    since this needs whichever of each already exist), or all Species/Sample
    at once when group=None."""
    # Species/Sample/TagDirs/<leaf_or_combo>/ — every assay's TagDirs sit
    # together directly under one sample now, so the assay is recovered from
    # the TagDir's own name (via _assay_of_tagdir), not from a per-assay
    # parent folder.
    all_tagdirs = sorted(p for p in cfg.project.glob("*/*/TagDirs/*") if p.is_dir())
    if not all_tagdirs:
        log.info("bedGraph: no TagDirs/* under %s", cfg.project)
        return

    tagdirs = all_tagdirs
    if group is not None:
        sp, sa = group
        # td = Species/Sample/TagDirs/<leaf_or_combo> so td.parent.parent is
        # the SAMPLE dir, and td.parent.parent.parent is the species dir.
        tagdirs = [td for td in all_tagdirs
                   if td.parent.parent.name == sa
                   and td.parent.parent.parent.name == sp]
        if not tagdirs:
            log.info("bedGraph: %d TagDirs exist under %s, but none matched "
                     "group %s/%s", len(all_tagdirs), cfg.project, sp, sa)
            return

    skip = f"-skipChr {cfg.skip_chr} " if cfg.skip_chr else ""
    for td in tagdirs:
        sample_dir = td.parent.parent                    # Species/Sample/
        sample = sample_dir.name
        species = sample_dir.parent.name
        assay = _assay_of_tagdir(td.name)
        if not assay:
            log.warning("bedGraph: could not classify assay for TagDir %s — skipping.", td)
            continue
        species_sample_run = f"{species}/{sample}/{td.name}"

        bedgraph_dir = sample_dir / "bedGraphs" / td.name
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