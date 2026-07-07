"""Step 3 — HOMER tag directories, split into two phases so leaf builds can
run in parallel across a SLURM array:

  run_leaf_tagdirs(cfg, sample_index=...)  — ARRAY-CAPABLE, one leaf TagDir
                                              per call, same --sample-index
                                              indexing as trim/align (one
                                              task per R1 file).
  run_combo_tagdirs(cfg)                   — runs once (in collect), merges
                                              every replicate's raw SAM files
                                              per assay into a combo TagDir.

Both build under Species/Sample/ instead of a flat project-level TagDirs/:
  • Species/Sample/<assay>-combo/TagDir   — all replicates of that assay merged
  • Species/Sample/<assay_rep>/TagDir     — one tag dir per individual replicate

Both are built from the same aligned SAM files (Species/Sample/<assay_rep>/Aligned/).
The combo build does NOT depend on the leaf TagDirs existing — it merges the
raw SAM files directly — so run_combo_tagdirs only needs the align phase to
be done, not run_leaf_tagdirs.
"""
from __future__ import annotations
from collections import defaultdict
from .utils import run, log, done, iter_leaf_dirs, assay_of_leaf, list_r1, leaf_dir


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


def _leaf_tagdir_for_r1(cfg, r1) -> None:
    run_dir = leaf_dir(r1)                     # .../Species/Sample/<assay_rep>/
    leaf_name = run_dir.name
    sample = run_dir.parent.name
    species = run_dir.parent.parent.name

    sams = sorted((run_dir / "Aligned").glob("*.Aligned.out.sam"))
    if not sams:
        log.warning("tagdir: no aligned SAMs yet for %s/%s/%s (run 'align' first)",
                    species, sample, leaf_name)
        return

    assay = assay_of_leaf(leaf_name)
    if not assay:
        log.warning("tagdir: could not classify assay for %s/%s/%s",
                    species, sample, leaf_name)
        return

    sams_str = " ".join(str(s) for s in sams)
    leaf_td = cfg.leaf_tagdir(species, sample, leaf_name)
    _make_tagdir(sams_str, leaf_td, assay,
                 f"tagdir {species}/{sample}/{leaf_name}", cfg)


def run_leaf_tagdirs(cfg, sample_index=None) -> None:
    """Build individual-replicate leaf TagDirs. Array-capable via
    --sample-index, indexed the same way as trim/align (one R1 file = one
    leaf run), so a tagdir_array.sbatch task maps 1:1 onto the align_array
    task that produced its SAM file."""
    r1s = list_r1(cfg)
    if not r1s:
        log.info("tagdir: no *_R1*.fastq[.gz] under nested RawData/ dirs in %s", cfg.project)
        return
    if sample_index is not None:
        if not (0 <= sample_index < len(r1s)):
            raise IndexError(f"sample_index {sample_index} out of range (0-{len(r1s)-1})")
        r1s = [r1s[sample_index]]
    for r1 in r1s:
        _leaf_tagdir_for_r1(cfg, r1)


def run_combo_tagdirs(cfg) -> None:
    """Merge every replicate's raw SAM files per assay into one combo TagDir
    per Species/Sample. Runs once (not array-capable) — typically in collect,
    after every leaf run's align phase has finished. Independent of whether
    run_leaf_tagdirs has run: combo TagDirs are built straight from the SAM
    files, not from the leaf TagDirs."""
    combo_groups: dict[tuple[str, str, str], list] = defaultdict(list)
    any_found = False

    for species, sample, ld in iter_leaf_dirs(cfg):
        sams = sorted((ld / "Aligned").glob("*.Aligned.out.sam"))
        if not sams:
            continue
        any_found = True

        assay = assay_of_leaf(ld.name)
        if not assay:
            log.warning("tagdir-combo: could not classify assay for %s/%s/%s",
                        species, sample, ld.name)
            continue

        combo_groups[(species, sample, assay)].extend(sams)

    if not any_found:
        log.info("tagdir-combo: no *.Aligned.out.sam under nested Aligned/ dirs in %s", cfg.project)
        return

    for (species, sample, assay), sams in combo_groups.items():
        sams_str = " ".join(str(s) for s in sorted(set(sams)))
        combo_td = cfg.combo_tagdir(species, sample, assay)
        _make_tagdir(sams_str, combo_td, assay,
                     f"tagdir {species}/{sample}/{assay}-combo", cfg)