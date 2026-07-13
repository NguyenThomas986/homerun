"""Pipeline orchestrator.

Runs the csRNA-seq steps in strict dependency order:

    trim → align → tagdirs → tagdirs-combo → bedgraphs → tss → qc → stability → report

Steps run sequentially and FAIL FAST: if any step raises, the pipeline stops
immediately with a non-zero exit, so a downstream step never runs on missing
inputs. Even when a subset is requested with --steps, they execute in this
canonical order.

For SLURM job arrays, --sample-index restricts trim/align/tagdirs to a single
leaf run (the Nth R1 file in RawData) — tagdirs builds that one leaf's
TagDir, letting the (slow) makeTagDirectory step run in parallel across an
array the same way trim/align do.

tagdirs-combo, bedgraphs, and tss are each keyed by Species/Sample rather
than by individual leaf/R1 file, so they honor a separate --group-index
instead: the Nth entry of utils.list_samples(cfg) (a sorted, deduplicated
list of every Species/Sample pair). This lets each of those three steps also
run as its own parallel SLURM array, one task per Species/Sample, rather
than looping over every sample serially inside one job. qc, stability, and
report remain whole-project, run-once steps.
"""
from __future__ import annotations

import argparse
import sys

from .config import load_config
from .utils import setup_logging, log, check_tools, list_r1, list_samples
from . import prepare, trim, mapping, tagdirs, bedgraphs, tss, qc, stability, report

STEP_ORDER = ["trim", "align", "tagdirs", "tagdirs-combo", "bedgraphs", "tss", "qc", "stability", "report"]
PER_SAMPLE = {"trim", "align", "tagdirs"}              # steps that honor --sample-index (one R1 file each)
GROUP_STEPS = {"tagdirs-combo", "bedgraphs", "tss"}    # steps that honor --group-index (one Species/Sample each)

STEP_FUNCS = {
    "trim":           trim.run_trim,
    "align":          mapping.run_mapping,
    "tagdirs":        tagdirs.run_leaf_tagdirs,     # array-capable: one leaf TagDir per call
    "tagdirs-combo":  tagdirs.run_combo_tagdirs,    # runs once: merges replicates per assay
    "bedgraphs":      bedgraphs.run_bedgraphs,
    "tss":            tss.run_tss,
    "qc":             qc.run_qc,
    "stability":      stability.run_stability,
    "report":         report.run_report,
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="csrnaseq",
        description="csRNA-seq pipeline (trim → STAR → tagdirs → bedGraphs → TSS → QC → stability)",
    )
    p.add_argument("--project", help="Project root (default: $CSRNA_PROJECT or CWD).")
    p.add_argument("--log-path", default=None,
                   help="Pipeline log file path (overrides CSRNA_LOG; else a timestamped "
                        "file under <project>/logs/).")
    p.add_argument("--steps", nargs="+", choices=STEP_ORDER,
                   help="Run only these steps (still executed in canonical order).")
    p.add_argument("--sample-index", type=int, default=None,
                   help="0-based index into RawData R1 files; restrict trim/align/tagdirs to "
                        "that one leaf run (used by SLURM array tasks via $SLURM_ARRAY_TASK_ID).")
    p.add_argument("--group-index", type=int, default=None,
                   help="0-based index into utils.list_samples(cfg) (sorted Species/Sample "
                        "pairs); restrict tagdirs-combo/bedgraphs/tss to that one Species/Sample "
                        "(used by SLURM array tasks via $SLURM_ARRAY_TASK_ID).")
    p.add_argument("--skip-prepare", action="store_true",
                   help="Skip folder creation / raw copy / STARIndex setup.")
    p.add_argument("--only-prepare", action="store_true",
                   help="Run prepare (folders/copy/STARIndex) and exit.")
    p.add_argument("--count-samples", action="store_true",
                   help="Print the number of leaf runs (R1 files in RawData) and exit.")
    p.add_argument("--count-groups", action="store_true",
                   help="Print the number of Species/Sample groups and exit.")
    p.add_argument("--stage-raw", action="store_true",
                   help="Move loose *_R1*/*_R2* FASTQs from the project root into RawData/ and exit.")

    # ── Config overrides (each beats its CSRNA_* env var; unset → env/default) ──
    g = p.add_argument_group("config overrides (override config.env when given)")
    g.add_argument("--aligner", choices=["star", "hisat2"], default=None,
                   help="Aligner (overrides CSRNA_ALIGNER).")
    g.add_argument("--genome-index", default=None,
                   help="STAR --genomeDir dir or HISAT2 -x prefix (overrides CSRNA_GENOME_INDEX).")
    g.add_argument("--genome", default=None,
                   help="HOMER -genome for tagdirs/TSS (overrides CSRNA_GENOME).")
    g.add_argument("--copy-src", default=None,
                   help="Glob of FASTQs to copy into RawData/ (overrides CSRNA_COPY_SRC).")
    g.add_argument("--threads", type=int, default=None,
                   help="Thread count (overrides CSRNA_THREADS / SLURM_CPUS_PER_TASK).")
    g.add_argument("--trim-adapter", default=None,
                   help="3' adapter sequence (overrides CSRNA_TRIM_ADAPTER).")
    g.add_argument("--trim-min", default=None,
                   help="Min read length after trimming (overrides CSRNA_TRIM_MINLEN).")
    g.add_argument("--trim-max", default=None,
                   help="Max read length after trimming (overrides CSRNA_TRIM_MAXLEN).")
    g.add_argument("--ntag-threshold", default=None,
                   help="Min tags to call a TSS cluster (overrides CSRNA_NTAG_THRESHOLD).")
    g.add_argument("--skip-chr", default=None,
                   help="Chromosome to exclude from bedGraphs (overrides CSRNA_SKIP_CHR).")

    # ── Alignment overrides (csRNA-tuned defaults stay unless set) ──────────────
    a = p.add_argument_group("alignment overrides (override config.env when given)")
    a.add_argument("--star-filter-multimap", default=None,
                   help="STAR --outFilterMultimapNmax: max loci a read may map to "
                        "(default 10000; overrides CSRNA_STAR_FILTER_MULTIMAP).")
    a.add_argument("--star-multimap-out", default=None,
                   help="STAR --outSAMmultNmax: alignments written per multimapping read "
                        "(default 1; overrides CSRNA_STAR_MULTIMAP_OUT).")
    a.add_argument("--star-multimap-order", choices=["Random", "Old_2.4"], default=None,
                   help="STAR --outMultimapperOrder (default Random; "
                        "overrides CSRNA_STAR_MULTIMAP_ORDER).")
    a.add_argument("--hisat2-strandness", choices=["F", "R", "FR", "RF"], default=None,
                   help="HISAT2 --rna-strandness (default F; overrides CSRNA_HISAT2_STRANDNESS).")
    return p


def run_pipeline(cfg, steps=None, skip_prepare=False, sample_index=None, group_index=None) -> None:
    if not skip_prepare:
        prepare.prepare(cfg)

    missing = check_tools(
        required=["homerTools", "STAR", "makeTagDirectory", "makeUCSCfile", "findcsRNATSS.pl"],
        optional=["skewer"],
    )
    if missing:
        log.warning("Missing required tools: %s (steps using them will fail).", ", ".join(missing))

    group = None
    if group_index is not None:
        groups = list_samples(cfg)
        if not (0 <= group_index < len(groups)):
            raise IndexError(f"group_index {group_index} out of range (0-{len(groups)-1})")
        group = groups[group_index]
        log.info("Group %d/%d: %s/%s", group_index, len(groups) - 1, *group)

    selected = [s for s in STEP_ORDER if (not steps or s in steps)]
    log.info("Running steps in order: %s", ", ".join(selected))
    for step in selected:
        log.info("=== STEP: %s ===", step)
        if step in PER_SAMPLE:
            STEP_FUNCS[step](cfg, sample_index=sample_index)
        elif step in GROUP_STEPS:
            STEP_FUNCS[step](cfg, group=group)
        else:
            STEP_FUNCS[step](cfg)      # raises on failure → fail-fast, order preserved
    log.info("Pipeline complete.")


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    cfg = load_config(args)

    # --count-samples / --count-groups: clean stdout (no logging) for the array controller
    if args.count_samples:
        print(len(list_r1(cfg)))
        return 0
    if args.count_groups:
        print(len(list_samples(cfg)))
        return 0

    setup_logging(cfg)
    log.info("Project: %s | genome=%s | threads=%d", cfg.project, cfg.genome, cfg.threads)

    # --stage-raw: move loose root FASTQs into RawData/ and exit (used by the
    # controller BEFORE --count-samples so the count sees the moved files).
    if args.stage_raw:
        prepare.stage_loose_fastqs(cfg)
        return 0

    if args.sample_index is not None:
        log.info("Sample index: %d", args.sample_index)
    if args.group_index is not None:
        log.info("Group index: %d", args.group_index)

    try:
        if args.only_prepare:
            prepare.prepare(cfg)
            log.info("Prepare complete.")
            return 0
        run_pipeline(cfg, steps=args.steps, skip_prepare=args.skip_prepare,
                     sample_index=args.sample_index, group_index=args.group_index)
    except Exception as exc:                 # fail-fast: surface and exit non-zero
        log.error("Pipeline aborted: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())