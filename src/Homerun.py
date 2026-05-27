#!/usr/bin/env python3
"""
Homerun — CSV-driven bioinformatics pipeline
Duttke Lab

Drives trimming → STAR alignment → tag directory → QC → DESeq2
from a single CSV file that serves as the sole source of truth.

Usage:
    python Homerun.py --csv /path/to/homerun_pipeline.csv
    python Homerun.py --csv pipeline.csv --dry-run
    python Homerun.py --generate-csv /path/to/fastqs --genome Homo_sapiens --out pipeline.csv
"""

import os
import sys
import time
import argparse
import logging
import traceback

# ── Bootstrap: install missing packages before anything else ──────────────────
def _bootstrap():
    required = ["pandas", "numpy"]
    for pkg in required:
        try:
            __import__(pkg)
        except ModuleNotFoundError:
            os.system(f"pip install {pkg} --quiet")

_bootstrap()

import pandas as pd  # noqa: E402  (after bootstrap)

from HomerunCSV     import PipelineCSV           # CSV read/write/lock
from HomerunState   import resolve_steps         # decide what to run
from HomerunTrim    import run_trim              # trimming
from HomerunSTAR    import run_star              # STAR alignment
from HomerunTagdir  import run_tagdir            # tag directory
from HomerunQC      import run_qc               # quality control
from HomerunDeseq2  import run_deseq2           # DESeq2 / downstream
from HomerunReport  import print_summary, write_report  # reporting
from HomerunGenCSV  import generate_csv         # auto-generate CSV helper

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)
log = logging.getLogger("Homerun")

# ── CLI ───────────────────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="Homerun",
        description="CSV-driven RNA/csRNA-seq pipeline (trim → STAR → tagdir → QC → DESeq2)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Primary mode
    p.add_argument("-c", "--csv",
                   help="Path to the pipeline CSV (single source of truth).")

    # Optional flags
    p.add_argument("--dry-run", action="store_true",
                   help="Print what would run without executing anything.")
    p.add_argument("--retry-failed", action="store_true",
                   help="Reset FAILED steps to PENDING before running.")
    p.add_argument("--steps", nargs="+",
                   choices=["trim", "star", "tagdir", "qc", "deseq2"],
                   help="Limit execution to specific steps only.")
    p.add_argument("--samples", nargs="+",
                   help="Limit execution to specific sample names only.")

    # CSV generation helper
    gen = p.add_argument_group("CSV generation (--generate-csv)")
    gen.add_argument("--generate-csv", metavar="FASTQ_DIR",
                     help="Auto-generate a pipeline CSV from a FASTQ directory.")
    gen.add_argument("--genome", default="Homo_sapiens",
                     help="Genome name to use when generating CSV (default: Homo_sapiens).")
    gen.add_argument("--out", default="homerun_pipeline.csv",
                     help="Output path for the generated CSV (default: homerun_pipeline.csv).")

    # SLURM / resource
    p.add_argument("--cpus", type=int,
                   default=int(os.getenv("SLURM_CPUS_PER_TASK", 1)),
                   help="CPUs to use (default: SLURM_CPUS_PER_TASK or 1).")
    p.add_argument("--working-path", default=os.getcwd(),
                   help="Root working directory for outputs.")

    return p


# ── Per-sample pipeline ───────────────────────────────────────────────────────
STEP_RUNNERS = {
    "trim":   run_trim,
    "star":   run_star,
    "tagdir": run_tagdir,
    "qc":     run_qc,
    "deseq2": run_deseq2,
}

STEP_ORDER = ["trim", "star", "tagdir", "qc", "deseq2"]


def process_sample(row: dict, csv: "PipelineCSV", args) -> dict:
    """
    Execute all incomplete steps for one sample.
    Returns the (possibly updated) row dict.
    """
    sample = row["sample"]
    steps_to_run = resolve_steps(row,
                                 allowed_steps=args.steps,
                                 retry_failed=args.retry_failed)

    if not steps_to_run:
        log.info("[%s] all steps complete — skipping.", sample)
        return row

    log.info("[%s] steps to run: %s", sample, ", ".join(steps_to_run))

    for step in STEP_ORDER:
        if step not in steps_to_run:
            continue

        log.info("[%s] ▶ starting step: %s", sample, step)

        if args.dry_run:
            log.info("[%s] (dry-run) would execute: %s", sample, step)
            continue

        try:
            runner = STEP_RUNNERS[step]
            result = runner(row, working_path=args.working_path, cpus=args.cpus)
            # result is a dict with keys matching CSV columns + optional "notes"
            row.update(result)
            row[step] = row.get(step, "DONE")   # runner should set this; belt+braces
            row["status"] = "RUNNING"
            log.info("[%s] ✔ step %s complete.", sample, step)

        except Exception as exc:
            tb = traceback.format_exc()
            row[step] = "FAILED"
            row["notes"] = f"{step} error: {exc}"
            row["status"] = "FAILED"
            log.error("[%s] ✘ step %s FAILED: %s\n%s", sample, step, exc, tb)
            break   # do not proceed to later steps for this sample

        finally:
            # Always persist state after every step attempt
            csv.update_row(sample, row)

    # Mark DONE only when every step is done or skipped
    all_done = all(
        row.get(s) not in (None, "", "PENDING", "FAILED")
        for s in STEP_ORDER
    )
    if all_done and row["status"] != "FAILED":
        row["status"] = "DONE"
        csv.update_row(sample, row)
        log.info("[%s] ✔✔ all steps complete.", sample)

    return row


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = build_parser()
    args = parser.parse_args()

    # ── CSV generation mode ───────────────────────────────────────────────────
    if args.generate_csv:
        log.info("Generating CSV from: %s", args.generate_csv)
        generate_csv(
            fastq_dir=args.generate_csv,
            genome=args.genome,
            output_path=args.out,
        )
        log.info("CSV written to: %s", args.out)
        return

    # ── Normal pipeline mode ──────────────────────────────────────────────────
    if not args.csv:
        parser.error("Provide --csv <path> to run the pipeline, "
                     "or --generate-csv <fastq_dir> to create a new CSV.")

    if not os.path.isfile(args.csv):
        log.error("CSV not found: %s", args.csv)
        sys.exit(1)

    timer_start = time.perf_counter()
    csv = PipelineCSV(args.csv)

    # Optionally filter to specific samples
    samples_to_run = args.samples  # None → all

    print_summary(csv.dataframe, label="Pipeline start")

    rows = csv.dataframe.to_dict(orient="records")
    processed = 0

    for row in rows:
        sample = row.get("sample", "")
        if samples_to_run and sample not in samples_to_run:
            log.debug("Skipping sample (not in --samples list): %s", sample)
            continue
        if row.get("status") == "DONE" and not args.retry_failed:
            log.info("[%s] status=DONE — skipping.", sample)
            continue

        process_sample(row, csv, args)
        processed += 1

    elapsed = time.perf_counter() - timer_start
    log.info("Pipeline finished. %d sample(s) processed in %.1fs.", processed, elapsed)

    print_summary(csv.dataframe, label="Pipeline end")

    if not args.dry_run:
        write_report(csv.dataframe, output_dir=args.working_path)


if __name__ == "__main__":
    main()
