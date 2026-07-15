"""Pipeline orchestrator.

Runs the csRNA-seq steps in strict dependency order:

    trim → align → tagdirs → tagdirs-combo → bedgraphs → tss → ritrie → qc → stability → report

Steps run sequentially and FAIL FAST: if any step raises, the pipeline stops
immediately with a non-zero exit, so a downstream step never runs on missing
inputs.

SLURM array support:
    --sample-index
        Restricts trim/align/tagdirs to one leaf run (Nth R1 FASTQ).

    --group-index
        Restricts tagdirs-combo/bedgraphs/tss to one Species/Sample group.

QC/report steps remain project-level single jobs.
"""

from __future__ import annotations

import argparse
import sys

from .config import load_config
from .utils import (
    setup_logging,
    log,
    check_tools,
    list_r1,
    list_samples,
)

from . import (
    __version__,
    prepare,
    trim,
    mapping,
    tagdirs,
    bedgraphs,
    tss,
    ritrie,
    qc,
    stability,
    report,
)


_BANNER = r"""
    __  __                                    
   / / / /___  ____ ___  ___  _______  ______ 
  / /_/ / __ \/ __ `__ \/ _ \/ ___/ / / / __ \
 / __  / /_/ / / / / / /  __/ /  / /_/ / / / /
/_/ /_/\____/_/ /_/ /_/\___/_/   \__,_/_/ /_/ 
"""


STEP_ORDER = [
    "trim",
    "align",
    "tagdirs",
    "tagdirs-combo",
    "bedgraphs",
    "tss",
    "ritrie",
    "qc",
    "stability",
    "report",
]


# SLURM array steps: one RawData FASTQ per task
PER_SAMPLE = {
    "trim",
    "align",
    "tagdirs",
}


# SLURM array steps: one Species/Sample per task
GROUP_STEPS = {
    "tagdirs-combo",
    "bedgraphs",
    "tss",
}


STEP_FUNCS = {
    "trim": trim.run_trim,
    "align": mapping.run_mapping,
    "tagdirs": tagdirs.run_leaf_tagdirs,
    "tagdirs-combo": tagdirs.run_combo_tagdirs,
    "bedgraphs": bedgraphs.run_bedgraphs,
    "tss": tss.run_tss,
    "ritrie": ritrie.run_ritrie,
    "qc": qc.run_qc,
    "stability": stability.run_stability,
    "report": report.run_report,
}



def build_parser() -> argparse.ArgumentParser:

    p = argparse.ArgumentParser(
        prog="csrnaseq",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            _BANNER
            + "\nRNA-seq analysis pipeline for HPC clusters."
            + f"\nVersion: {__version__}"
        ),
    )

    p.add_argument(
        "--project",
        help="Project root (default: $CSRNA_PROJECT or CWD).",
    )

    p.add_argument(
        "--log-path",
        default=None,
        help=(
            "Pipeline log file path (overrides CSRNA_LOG; "
            "else timestamped file under <project>/logs/)."
        ),
    )


    p.add_argument(
        "--steps",
        nargs="+",
        choices=STEP_ORDER,
        help=(
            "Run only these steps "
            "(still executed in canonical order)."
        ),
    )


    p.add_argument(
        "--sample-index",
        type=int,
        default=None,
        help=(
            "0-based index into RawData R1 files. "
            "Restricts trim/align/tagdirs to one leaf run. "
            "Used by SLURM_ARRAY_TASK_ID."
        ),
    )


    p.add_argument(
        "--group-index",
        type=int,
        default=None,
        help=(
            "0-based index into list_samples(cfg). "
            "Restricts tagdirs-combo/bedgraphs/tss "
            "to one Species/Sample group. "
            "Used by SLURM_ARRAY_TASK_ID."
        ),
    )


    p.add_argument(
        "--skip-prepare",
        action="store_true",
        help="Skip folder creation/raw copy/STARIndex setup.",
    )


    p.add_argument(
        "--only-prepare",
        action="store_true",
        help="Run prepare and exit.",
    )


    p.add_argument(
        "--count-samples",
        action="store_true",
        help=(
            "Print number of leaf runs "
            "(R1 files in RawData) and exit."
        ),
    )


    p.add_argument(
        "--count-groups",
        action="store_true",
        help=(
            "Print number of Species/Sample groups "
            "and exit."
        ),
    )


    p.add_argument(
        "--stage-raw",
        action="store_true",
        help=(
            "Move loose *_R1*/*_R2* FASTQs "
            "into RawData/ and exit."
        ),
    )


    # ----------------------------------------------------------
    # Config overrides
    # ----------------------------------------------------------

    g = p.add_argument_group(
        "config overrides (override config.env when given)"
    )


    g.add_argument(
        "--aligner",
        choices=["star", "hisat2"],
        default=None,
        help="Aligner (overrides CSRNA_ALIGNER).",
    )


    g.add_argument(
        "--genome-index",
        default=None,
        help=(
            "STAR genomeDir or HISAT2 prefix "
            "(overrides CSRNA_GENOME_INDEX)."
        ),
    )


    g.add_argument(
        "--genome",
        default=None,
        help="HOMER genome (overrides CSRNA_GENOME).",
    )


    g.add_argument(
        "--gtf",
        default=None,
        help=(
            "GTF annotation file for RIT/RIE metric "
            "(overrides CSRNA_GTF)."
        ),
    )


    g.add_argument(
        "--copy-src",
        default=None,
        help="FASTQ copy source (overrides CSRNA_COPY_SRC).",
    )


    g.add_argument(
        "--threads",
        type=int,
        default=None,
        help=(
            "Threads "
            "(overrides CSRNA_THREADS / SLURM_CPUS_PER_TASK)."
        ),
    )


    g.add_argument(
        "--trim-adapter",
        default=None,
        help="3' adapter sequence.",
    )


    g.add_argument(
        "--trim-min",
        default=None,
        help="Minimum read length after trimming.",
    )


    g.add_argument(
        "--trim-max",
        default=None,
        help="Maximum read length after trimming.",
    )


    g.add_argument(
        "--ntag-threshold",
        default=None,
        help="Minimum tags for TSS calling.",
    )


    g.add_argument(
        "--skip-chr",
        default=None,
        help="Chromosome excluded from bedGraphs.",
    )



    # ----------------------------------------------------------
    # Alignment overrides
    # ----------------------------------------------------------

    a = p.add_argument_group(
        "alignment overrides"
    )


    a.add_argument(
        "--star-filter-multimap",
        default=None,
        help="STAR outFilterMultimapNmax.",
    )


    a.add_argument(
        "--star-multimap-out",
        default=None,
        help="STAR outSAMmultNmax.",
    )


    a.add_argument(
        "--star-multimap-order",
        choices=["Random", "Old_2.4"],
        default=None,
        help="STAR outMultimapperOrder.",
    )


    a.add_argument(
        "--hisat2-strandness",
        choices=["F", "R", "FR", "RF"],
        default=None,
        help="HISAT2 RNA strandness.",
    )


    return p

def run_pipeline(
    cfg,
    steps=None,
    skip_prepare=False,
    sample_index=None,
    group_index=None,
) -> None:

    if not skip_prepare:
        prepare.prepare(cfg)


    missing = check_tools(
        required=[
            "homerTools",
            "STAR",
            "makeTagDirectory",
            "makeUCSCfile",
            "findcsRNATSS.pl",
        ],
        optional=[
            "skewer",
            "findPeaks",
            "mergePeaks",
            "annotatePeaks.pl",
            "parseGTF.pl",
        ],
    )


    if missing:
        log.warning(
            "Missing required tools: %s "
            "(steps using them may fail).",
            ", ".join(missing),
        )


    group = None

    if group_index is not None:

        groups = list_samples(cfg)

        if not (0 <= group_index < len(groups)):
            raise IndexError(
                f"group_index {group_index} out of range "
                f"(0-{len(groups)-1})"
            )


        group = groups[group_index]

        log.info(
            "Group %d/%d: %s/%s",
            group_index,
            len(groups)-1,
            *group,
        )


    selected = [
        step
        for step in STEP_ORDER
        if not steps or step in steps
    ]


    log.info(
        "Running steps in order: %s",
        ", ".join(selected),
    )


    for step in selected:

        log.info(
            "=== STEP: %s ===",
            step,
        )


        if step in PER_SAMPLE:

            STEP_FUNCS[step](
                cfg,
                sample_index=sample_index,
            )


        elif step in GROUP_STEPS:

            STEP_FUNCS[step](
                cfg,
                group=group,
            )


        else:

            # qc, ritrie, stability, report
            # run once on the full project
            STEP_FUNCS[step](cfg)


    log.info(
        "Pipeline complete."
    )



def print_banner():

    print(_BANNER)
    print(
        "RNA-seq analysis pipeline for HPC clusters."
    )
    print(
        f"Version: {__version__}"
    )



def main(argv=None) -> int:

    args = build_parser().parse_args(argv)

    cfg = load_config(args)



    # ----------------------------------------------------------
    # Array controller helpers
    # Keep stdout clean for SLURM submit scripts
    # ----------------------------------------------------------

    if args.count_samples:

        print(
            len(list_r1(cfg))
        )

        return 0



    if args.count_groups:

        print(
            len(list_samples(cfg))
        )

        return 0



    print_banner()



    setup_logging(cfg)


    log.info(
        "Project: %s | genome=%s | threads=%d",
        cfg.project,
        cfg.genome,
        cfg.threads,
    )



    # ----------------------------------------------------------
    # Stage raw FASTQs then exit
    # ----------------------------------------------------------

    if args.stage_raw:

        # Regenerate config.txt here too, not just from the full prepare().
        # --stage-raw is the narrow, fast path (move loose FASTQs and exit)
        # that submit_array.sh's local pre-flight step uses before counting
        # samples/groups and submitting the SLURM array — it deliberately
        # skips copy_raw()/ensure_starindex()/validate_gtf(), but skipping
        # write_config_summary() too meant config.txt silently went stale
        # (still showing whatever the last full `prepare()` run wrote) even
        # though the array jobs themselves read the filesystem directly via
        # list_r1()/list_samples() and always saw the newly staged samples
        # correctly. Cheap to call, so just always keep config.txt honest.
        prepare.stage_loose_fastqs(cfg)

        prepare.write_config_summary(cfg)

        return 0



    if args.sample_index is not None:

        log.info(
            "Sample index: %d",
            args.sample_index,
        )


    if args.group_index is not None:

        log.info(
            "Group index: %d",
            args.group_index,
        )



    try:

        if args.only_prepare:

            prepare.prepare(cfg)

            log.info(
                "Prepare complete."
            )

            return 0



        run_pipeline(
            cfg,
            steps=args.steps,
            skip_prepare=args.skip_prepare,
            sample_index=args.sample_index,
            group_index=args.group_index,
        )


    except Exception as exc:

        log.error(
            "Pipeline aborted: %s",
            exc,
        )

        return 1



    return 0



if __name__ == "__main__":

    sys.exit(main())