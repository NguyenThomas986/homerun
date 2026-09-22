"""Step 4 — Genome-browser bedGraphs (strand-specific), written into
Species/bedGraphs/ next to Species/TagDirs/.

Generates a bedGraph folder next to every TagDir built by tagdirs.py:

  • Species/TagDirs/<sample>_<assay>-combo -> Species/bedGraphs/<same-name>/
  • Species/TagDirs/<sample>_<leaf_name>   -> Species/bedGraphs/<same-name>/

Each output keeps the full TagDir identity in the filename, for example:

  Testicle_D2_Post_BDAD_NO_RA_csRNA-combo.posStrand.bedGraph.gz
  Testicle_D2_Post_BDAD_NO_RA_csRNA-combo.negStrand.bedGraph.gz

makeUCSCfile adds gzip compression itself, so -o is given the uncompressed
base filename while the completion check looks for the resulting .gz file.
"""

from __future__ import annotations

from .utils import run, log, done, assay_of_leaf, list_samples


def _assay_of_tagdir(name: str) -> str | None:
    """Recover the assay from a TagDir's own name.

    TagDir names carry a <sample>_ prefix and may also contain condition
    information, for example:

        Testicle_D2_Post_BDAD_NO_RA_csRNA-combo
        Testicle_D2_Post_BDAD_NO_RA_csRNA_r1

    Strip '-combo' first when present, then use assay_of_leaf() to identify
    csRNA, sRNA, or totalRNA from the remaining tokens.
    """
    if name.endswith("-combo"):
        name = name[:-len("-combo")]

    return assay_of_leaf(name)


def run_bedgraphs(cfg, group=None) -> None:
    """Create strand-specific bedGraphs for every existing TagDir.

    When group=(species, sample) is supplied, only TagDirs belonging to that
    Species/Sample are processed. Otherwise all existing TagDirs are used.
    """

    # Species/TagDirs/<sample>_<leaf_or_combo>/
    all_tagdirs = sorted(
        p
        for p in cfg.project.glob("*/TagDirs/*")
        if p.is_dir()
    )

    if not all_tagdirs:
        log.info(
            "bedGraph: no TagDirs/* under %s",
            cfg.project,
        )
        return

    tagdirs = all_tagdirs

    if group is not None:
        sp, sa = group

        tagdirs = [
            td
            for td in all_tagdirs
            if td.parent.parent.name == sp
            and td.name.startswith(f"{sa}_")
        ]

        if not tagdirs:
            log.info(
                "bedGraph: %d TagDirs exist under %s, but none matched "
                "group %s/%s",
                len(all_tagdirs),
                cfg.project,
                sp,
                sa,
            )
            return

    skip = (
        f"-skipChr {cfg.skip_chr} "
        if cfg.skip_chr
        else ""
    )

    for td in tagdirs:
        species_dir = td.parent.parent
        species = species_dir.name

        candidates = [
            sample
            for sp, sample in list_samples(cfg)
            if sp == species
            and td.name.startswith(f"{sample}_")
        ]

        if not candidates:
            log.warning(
                "bedGraph: could not identify sample for TagDir %s — skipping.",
                td,
            )
            continue

        # Prefer the longest matching sample name in case one sample name is
        # a prefix of another.
        sample = max(
            candidates,
            key=len,
        )

        assay = _assay_of_tagdir(td.name)

        if not assay:
            log.warning(
                "bedGraph: could not classify assay for TagDir %s — skipping.",
                td,
            )
            continue

        species_sample_run = (
            f"{species}/{sample}/{td.name}"
        )

        bedgraph_dir = (
            species_dir
            / "bedGraphs"
            / td.name
        )

        bedgraph_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        style = (
            "rnaseq"
            if assay == "totalRNA"
            else "tss"
        )

        # Important:
        # makeUCSCfile appends .gz itself.
        #
        # Give -o the .bedGraph base path:
        #
        #   ...posStrand.bedGraph
        #
        # HOMER then creates:
        #
        #   ...posStrand.bedGraph.gz
        #
        # This avoids the old:
        #
        #   negStrand.gz.gz
        #
        pos_base = (
            bedgraph_dir
            / f"{td.name}.posStrand.bedGraph"
        )

        neg_base = (
            bedgraph_dir
            / f"{td.name}.negStrand.bedGraph"
        )

        pos_output = (
            bedgraph_dir
            / f"{td.name}.posStrand.bedGraph.gz"
        )

        neg_output = (
            bedgraph_dir
            / f"{td.name}.negStrand.bedGraph.gz"
        )

        if done(pos_output):
            log.info(
                "  skip (done): %s",
                pos_output,
            )
        else:
            run(
                (
                    f"makeUCSCfile {td} "
                    f"-style {style} "
                    f"-strand + "
                    f"{skip}"
                    f"-o {pos_base}"
                ),
                label=(
                    f"bedGraph + "
                    f"{species_sample_run}"
                ),
            )

        if done(neg_output):
            log.info(
                "  skip (done): %s",
                neg_output,
            )
        else:
            run(
                (
                    f"makeUCSCfile {td} "
                    f"-style {style} "
                    f"-strand - "
                    f"-neg "
                    f"{skip}"
                    f"-o {neg_base}"
                ),
                label=(
                    f"bedGraph - "
                    f"{species_sample_run}"
                ),
            )