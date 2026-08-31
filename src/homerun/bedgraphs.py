"""Step 4 — Genome-browser bedGraphs (strand-specific), written into
Species/bedGraphs/ next to Species/TagDirs/.

Generates a bedGraph folder next to every TagDir built by tagdirs.py:
  • Species/TagDirs/<sample>_<assay>-combo -> Species/bedGraphs/<same-name>/
  • Species/TagDirs/<sample>_<leaf_name>   -> Species/bedGraphs/<same-name>/

Writes plain UNCOMPRESSED .bedGraph files with makeUCSCfile -o (stable names,
no piping/gzip) so a re-run command keeps working.
"""
from __future__ import annotations
from .utils import run, log, done, assay_of_leaf, list_samples


def _assay_of_tagdir(name: str) -> str | None:
    """Recover the assay from a TagDir's own name, since there's no longer a
    per-assay parent folder to read it off of. TagDir names now carry a
    <sample>_ prefix (e.g. 'IMR90_csRNA-combo', 'IMR90_csRNA_r1' — see
    Config.leaf_tagdir/combo_tagdir), so a plain '-combo' strip alone would
    leave the sample prefix stuck to the assay ('IMR90_csRNA' instead of
    'csRNA'). Strip '-combo' first if present, then run the same
    position-independent token search assay_of_leaf() already uses for
    leaf names — it already correctly ignores non-assay tokens (species,
    sample, condition), so it handles the sample-prefixed combo case too."""
    if name.endswith("-combo"):
        name = name[: -len("-combo")]
    return assay_of_leaf(name)


def run_bedgraphs(cfg, group=None) -> None:
    """Array-capable via --group-index (group=(species, sample) restricts to
    just that one Species/Sample's TagDirs — both its leaf and combo TagDirs,
    since this needs whichever of each already exist), or all Species/Sample
    at once when group=None."""
    # Species/TagDirs/<sample>_<leaf_or_combo>/ — every sample and assay sits
    # together under its species, so the assay is recovered from
    # the TagDir's own name (via _assay_of_tagdir), not from a per-assay
    # parent folder.
    all_tagdirs = sorted(p for p in cfg.project.glob("*/TagDirs/*") if p.is_dir())
    if not all_tagdirs:
        log.info("bedGraph: no TagDirs/* under %s", cfg.project)
        return

    tagdirs = all_tagdirs
    if group is not None:
        sp, sa = group
        tagdirs = [td for td in all_tagdirs
                   if td.parent.parent.name == sp
                   and td.name.startswith(f"{sa}_")]
        if not tagdirs:
            log.info("bedGraph: %d TagDirs exist under %s, but none matched "
                     "group %s/%s", len(all_tagdirs), cfg.project, sp, sa)
            return

    skip = f"-skipChr {cfg.skip_chr} " if cfg.skip_chr else ""
    for td in tagdirs:
        species_dir = td.parent.parent                   # Species/
        species = species_dir.name
        candidates = [sample for sp, sample in list_samples(cfg)
                      if sp == species and td.name.startswith(f"{sample}_")]
        if not candidates:
            log.warning("bedGraph: could not identify sample for TagDir %s — skipping.", td)
            continue
        sample = max(candidates, key=len)
        assay = _assay_of_tagdir(td.name)
        if not assay:
            log.warning("bedGraph: could not classify assay for TagDir %s — skipping.", td)
            continue
        species_sample_run = f"{species}/{sample}/{td.name}"

        bedgraph_dir = species_dir / "bedGraphs" / td.name
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
