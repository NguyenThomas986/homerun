"""Step 3 — HOMER tag directories, split into two phases so leaf builds can
run in parallel across a SLURM array:

  run_leaf_tagdirs(cfg, sample_index=...)  — ARRAY-CAPABLE, one leaf TagDir
                                              per call, same --sample-index
                                              indexing as trim/align (one
                                              task per R1 file).
  run_combo_tagdirs(cfg, group=...)        — ARRAY-CAPABLE via --group-index,
                                              one Species/Sample per call
                                              (group=(species, sample)), or
                                              all of them at once when
                                              group=None. Merges every
                                              replicate's raw SAM files per
                                              assay into a combo TagDir.

Built under Species/Sample/TagDirs/ (all assays of a sample together):
  • Species/Sample/TagDirs/<assay>-combo  — all replicates of that assay merged
  • Species/Sample/TagDirs/<leaf_name>     — one tag dir per individual replicate

Both are built from the same aligned SAM files, which live in the shared
Species/Sample/Aligned/ (one folder per sample, not per assay or per
replicate — see Config.aligned_dir). Since replicates (and assays) no longer
get their own directory, a replicate's SAM is located by filename prefix
rather than by listing an Aligned/ folder that belongs to just that one
replicate.

The combo build does NOT depend on the leaf TagDirs existing — it merges the
raw SAM files directly — so run_combo_tagdirs only needs the align phase to
be done, not run_leaf_tagdirs.
"""
from __future__ import annotations
from collections import defaultdict
from .utils import run, log, done, iter_leaf_dirs, assay_of_leaf, list_r1, parse_sample_name


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


def _sam_for_r1(cfg, species, sample, r1):
    """The aligned SAM this replicate's R1 produced, in the shared
    Species/Sample/Aligned/ folder — same naming mapping.py already
    uses (<prefix>.Aligned.out.sam where prefix = r1.name up to '_R1')."""
    prefix = r1.name.split("_R1")[0]
    return cfg.aligned_dir(species, sample) / f"{prefix}.Aligned.out.sam"


def _leaf_tagdir_for_r1(cfg, r1) -> None:
    species, sample, leaf_name = parse_sample_name(r1.name)
    assay = assay_of_leaf(leaf_name)
    if not assay:
        log.warning("tagdir: could not classify assay for %s", r1.name)
        return

    sam = _sam_for_r1(cfg, species, sample, r1)
    if not sam.exists():
        log.warning("tagdir: no aligned SAM yet for %s/%s/%s (run 'align' first)",
                    species, sample, leaf_name)
        return

    leaf_td = cfg.leaf_tagdir(species, sample, leaf_name)
    _make_tagdir(str(sam), leaf_td, assay,
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


def run_combo_tagdirs(cfg, group=None) -> None:
    """Merge every replicate's raw SAM files per assay into one combo TagDir
    per Species/Sample. Array-capable via --group-index (group=(species,
    sample) restricts to just that one Species/Sample), or runs once for
    every Species/Sample when group=None. Only needs the align phase to be
    done, not run_leaf_tagdirs — combo TagDirs are built straight from the
    SAM files, not from the leaf TagDirs."""
    combo_groups: dict[tuple[str, str, str], list] = defaultdict(list)
    any_found = False

    expected_leaves: dict[tuple[str, str, str], int] = defaultdict(int)

    for species, sample, leaf_name, r1 in iter_leaf_dirs(cfg):
        if group is not None and (species, sample) != group:
            continue

        assay = assay_of_leaf(leaf_name)
        if not assay:
            log.warning("tagdir-combo: could not classify assay for %s/%s/%s",
                        species, sample, leaf_name)
            continue
        expected_leaves[(species, sample, assay)] += 1

        sam = _sam_for_r1(cfg, species, sample, r1)
        if not sam.exists():
            log.warning("tagdir-combo: no aligned SAM for %s/%s/%s — this replicate "
                        "will be MISSING from the combo TagDir (check align logs).",
                        species, sample, leaf_name)
            continue
        any_found = True

        combo_groups[(species, sample, assay)].append(sam)

    for key, n_expected in expected_leaves.items():
        n_found = len(set(combo_groups.get(key, [])))
        if 0 < n_found < n_expected:
            species, sample, assay = key
            log.warning("tagdir-combo: %s/%s/%s-combo built from only %d of %d "
                        "expected replicate(s) — one or more alignments failed/missing.",
                        species, sample, assay, n_found, n_expected)

    if not any_found:
        log.info("tagdir-combo: no aligned SAM files found under %s", cfg.project)
        return

    for (species, sample, assay), sams in combo_groups.items():
        sams_str = " ".join(str(s) for s in sorted(set(sams)))
        combo_td = cfg.combo_tagdir(species, sample, assay)
        _make_tagdir(sams_str, combo_td, assay,
                     f"tagdir {species}/{sample}/{assay}-combo", cfg)