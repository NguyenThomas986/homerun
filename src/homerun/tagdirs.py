"""Step 3 — build one HOMER tag directory per biological replicate.

FASTQ filenames are reduced to a replicate identity by ``parse_sample_name``.
Every aligned SAM with the same ``(species, sample, leaf_name)`` is supplied
to one ``makeTagDirectory`` call. Thus lane/accession suffixes are combined,
while ``r1`` and ``r2`` remain separate. Assay-wide uber combos are not built.

The SLURM array remains indexed by R1 FASTQ for compatibility with trim and
align. Only the first sorted R1 in a replicate group performs the build; the
other lane tasks log a skip.
"""
from __future__ import annotations

from collections import defaultdict

from .utils import assay_of_leaf, done, iter_leaf_dirs, list_r1, log, parse_sample_name, run


def _make_tagdir(input_sams: list, tagdir, assay: str, label: str, cfg) -> None:
    if done(tagdir):
        log.info("  skip (done): %s", tagdir)
        return

    sams = " ".join(str(path) for path in sorted(set(input_sams)))
    if assay in ("csRNA", "sRNA"):
        cmd = (f"makeTagDirectory {tagdir} {sams} "
               f"-genome {cfg.genome} -checkGC -fragLength 150 -omitSN")
    elif assay == "totalRNA":
        cmd = (f"makeTagDirectory {tagdir} {sams} "
               f"-genome {cfg.genome} -checkGC -fragLength 150 -read2")
    else:
        log.warning("tagdir: unrecognized assay '%s' for %s", assay, tagdir)
        return
    run(cmd, label=label)


def _sam_for_r1(cfg, species, sample, r1):
    """Return the SAM produced by mapping.py for one R1 FASTQ."""
    prefix = r1.name.split("_R1")[0]
    return cfg.aligned_dir(species, sample) / f"{prefix}.Aligned.out.sam"


def _replicate_groups(cfg):
    """Map one parsed replicate identity to all R1 files sharing it."""
    groups = defaultdict(list)
    for species, sample, leaf_name, r1 in iter_leaf_dirs(cfg):
        groups[(species, sample, leaf_name)].append(r1)
    return {key: sorted(paths) for key, paths in groups.items()}


def _build_replicate_tagdir(cfg, key, r1s) -> None:
    species, sample, leaf_name = key
    assay = assay_of_leaf(leaf_name)
    if not assay:
        log.warning("tagdir: could not classify assay for %s/%s/%s",
                    species, sample, leaf_name)
        return

    expected = [_sam_for_r1(cfg, species, sample, r1) for r1 in r1s]
    missing = [sam for sam in expected if not sam.exists()]
    if missing:
        log.warning(
            "tagdir: %s/%s/%s is missing %d of %d aligned SAM file(s); "
            "not building a partial TagDir: %s",
            species, sample, leaf_name, len(missing), len(expected),
            ", ".join(path.name for path in missing),
        )
        return

    _make_tagdir(
        expected,
        cfg.leaf_tagdir(species, sample, leaf_name),
        assay,
        f"tagdir {species}/{sample}/{leaf_name} ({len(expected)} file(s))",
        cfg,
    )


def run_leaf_tagdirs(cfg, sample_index=None) -> None:
    """Build replicate TagDirs, merging files with an identical r-X identity."""
    r1s = list_r1(cfg)
    if not r1s:
        log.info("tagdir: no *_R1*.fastq[.gz] under nested RawData/ dirs in %s",
                 cfg.project)
        return

    groups = _replicate_groups(cfg)
    if sample_index is not None:
        if not (0 <= sample_index < len(r1s)):
            raise IndexError(f"sample_index {sample_index} out of range (0-{len(r1s)-1})")
        selected = r1s[sample_index]
        key = parse_sample_name(selected.name)
        members = groups[key]
        if selected != members[0]:
            log.info(
                "tagdir: %s belongs to %s/%s/%s; canonical task is %s — skip",
                selected.name, *key, members[0].name,
            )
            return
        groups = {key: members}

    for key in sorted(groups):
        _build_replicate_tagdir(cfg, key, groups[key])


def run_combo_tagdirs(cfg, group=None) -> None:
    """Compatibility no-op: assay-wide uber-combo TagDirs were removed."""
    log.info(
        "tagdirs-combo: skipped; only files sharing the same replicate "
        "identity are merged (r1 with r1, r2 with r2)."
    )
