"""Step 3 — HOMER tag directories.

Builds two kinds of tag directories per sample, nested under
Species/Sample/ instead of a flat project-level TagDirs/:
  • Species/Sample/<assay>-combo/TagDir   — all replicates of that assay merged
  • Species/Sample/<assay_rep>/TagDir     — one tag dir per individual replicate

Both are built from the same aligned SAM files (Species/Sample/<assay_rep>/Aligned/).
"""
from __future__ import annotations
from collections import defaultdict
from .utils import run, log, done, iter_leaf_dirs, assay_of_leaf


def _make_tagdir(cmd_input_sams: str, tagdir, assay: str, label: str, cfg) -> None:
    if done(tagdir):
        log.info("  skip (done): %s", tagdir)
        return
    if assay in ("csRNA", "sRNA"):
        cmd = (f"makeTagDirectory {tagdir} {cmd_input_sams} "
               f"-genome {cfg.genome} -checkGC -fragLength 150 -omitSN")
    elif assay == "totalRNA":
        cmd = (f"makeTagDirectory {tagdir} {cmd_input_sams} "
               f"-genome {cfg.genome} -checkGC -fragLength 150 -read2")
    else:
        log.warning("tagdir: unrecognized assay '%s' for %s", assay, tagdir)
        return
    run(cmd, label=label)


def run_tagdirs(cfg) -> None:
    combo_groups: dict[tuple[str, str, str], list] = defaultdict(list)
    any_found = False

    for species, sample, leaf_dir in iter_leaf_dirs(cfg):
        sams = sorted((leaf_dir / "Aligned").glob("*.Aligned.out.sam"))
        if not sams:
            continue
        any_found = True

        leaf_name = leaf_dir.name
        assay = assay_of_leaf(leaf_name)
        if not assay:
            log.warning("tagdir: could not classify assay for %s/%s/%s",
                        species, sample, leaf_name)
            continue

        # ── Per-replicate: one tag dir per individual leaf run ──────────────
        sams_str = " ".join(str(s) for s in sams)
        leaf_td = cfg.leaf_tagdir(species, sample, leaf_name)
        _make_tagdir(sams_str, leaf_td, assay,
                     f"tagdir {species}/{sample}/{leaf_name}", cfg)

        combo_groups[(species, sample, assay)].extend(sams)

    if not any_found:
        log.info("tagdir: no *.Aligned.out.sam under nested Aligned/ dirs in %s", cfg.project)
        return

    # ── Combo: all replicates of an assay merged, at the sample level ───────
    for (species, sample, assay), sams in combo_groups.items():
        sams_str = " ".join(str(s) for s in sorted(set(sams)))
        combo_td = cfg.combo_tagdir(species, sample, assay)
        _make_tagdir(sams_str, combo_td, assay,
                     f"tagdir {species}/{sample}/{assay}-combo", cfg)