"""Pipeline orchestrator.

Runs the csRNA-seq steps in strict dependency order:

    trim → align → tagdirs → bedgraphs → tss → qc → stability

Steps run sequentially and FAIL FAST: if any step raises, the pipeline stops
immediately with a non-zero exit, so a downstream step never runs on missing
inputs. Even when a subset is requested with --steps, they execute in this
canonical order.

For SLURM job arrays, --sample-index restricts trim/align to a single sample
(the Nth R1 file in RawData); the downstream collect steps run once afterward.
"""
from __future__ import annotations

import argparse
import sys

from .config import load_config
from .utils import setup_logging, log, check_tools, list_r1
from . import prepare, trim, mapping, tagdirs, bedgraphs, tss, qc, stability, report

STEP_ORDER = ["trim", "align", "tagdirs", "bedgraphs", "tss", "qc", "stability", "report"]
PER_SAMPLE = {"trim", "align"}          # steps that honor --sample-index

STEP_FUNCS = {
    "trim":       trim.run_trim,
    "align":      mapping.run_mapping,
    "tagdirs":    tagdirs.run_tagdirs,
    "bedgraphs":  bedgraphs.run_bedgraphs,
    "tss":        tss.run_tss,
    "qc":         qc.run_qc,
    "stability":  stability.run_stability,
    "report":     report.run_report,
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
                   help="0-based index into RawData R1 files; restrict trim/align to that one "
                        "sample (used by SLURM array tasks via $SLURM_ARRAY_TASK_ID).")
    p.add_argument("--skip-prepare", action="store_true",
                   help="Skip folder creation / raw copy / STARIndex setup.")
    p.add_argument("--only-prepare", action="store_true",
                   help="Run prepare (folders/copy/STARIndex) and exit.")
    p.add_argument("--count-samples", action="store_true",
                   help="Print the number of samples (R1 files in RawData) and exit.")
    p.add_argument("--list-samples", action="store_true",
                   help="Print each sample index and name and exit.")
    return p


def run_pipeline(cfg, steps=None, skip_prepare=False, sample_index=None) -> None:
    if not skip_prepare:
        prepare.prepare(cfg)

    missing = check_tools(
        required=["homerTools", "STAR", "makeTagDirectory", "makeUCSCfile", "findcsRNATSS.pl"],
        optional=["skewer"],
    )
    if missing:
        log.warning("Missing required tools: %s (steps using them will fail).", ", ".join(missing))

    selected = [s for s in STEP_ORDER if (not steps or s in steps)]
    log.info("Running steps in order: %s", ", ".join(selected))
    for step in selected:
        log.info("=== STEP: %s ===", step)
        if step in PER_SAMPLE:
            STEP_FUNCS[step](cfg, sample_index=sample_index)
        else:
            STEP_FUNCS[step](cfg)      # raises on failure → fail-fast, order preserved
    log.info("Pipeline complete.")


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    cfg = load_config(args)

    # --count-samples: clean stdout (no logging) for the array controller
    if args.count_samples:
        print(len(list_r1(cfg)))
        return 0

    # --list-samples: print index: samplename for the array controller
    if args.list_samples:
        for i, r1 in enumerate(list_r1(cfg)):
            name = r1.name.split("_R1")[0]
            print(f"{i}: {name}")
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

    try:
        if args.only_prepare:
            prepare.prepare(cfg)
            log.info("Prepare complete.")
            return 0
        run_pipeline(cfg, steps=args.steps, skip_prepare=args.skip_prepare,
                     sample_index=args.sample_index)
    except Exception as exc:                 # fail-fast: surface and exit non-zero
        log.error("Pipeline aborted: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
