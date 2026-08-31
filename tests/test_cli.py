"""CLI tests: argument parsing, --help, exit codes, and the
--count-samples/--count-groups/--check-rerun "array controller" fast paths.

None of these invoke run_pipeline (no HOMER/STAR/etc. needed) — they only
exercise build_parser()/main()'s pre-pipeline branches, which is everything
a SLURM wrapper script actually depends on.
"""
from __future__ import annotations

import pytest

from homerun.pipeline import build_parser, main, STEP_ORDER


# ── argparse construction ────────────────────────────────────────────────────

def test_build_parser_returns_parser():
    p = build_parser()
    assert p.prog == "homerun"


def test_help_exits_zero(capsys):
    p = build_parser()
    with pytest.raises(SystemExit) as exc_info:
        p.parse_args(["--help"])
    assert exc_info.value.code == 0
    out = capsys.readouterr().out
    assert "--project" in out
    assert "--steps" in out


def test_no_args_parses_with_defaults():
    p = build_parser()
    args = p.parse_args([])
    assert args.project is None
    assert args.steps is None
    assert args.sample_index is None
    assert args.group_index is None
    assert args.force is False


def test_steps_accepts_valid_step_names():
    p = build_parser()
    args = p.parse_args(["--steps", "trim", "align"])
    assert args.steps == ["trim", "align"]


def test_steps_rejects_invalid_step_name():
    p = build_parser()
    with pytest.raises(SystemExit):
        p.parse_args(["--steps", "not-a-real-step"])


def test_sample_index_must_be_int():
    p = build_parser()
    args = p.parse_args(["--sample-index", "3"])
    assert args.sample_index == 3
    with pytest.raises(SystemExit):
        p.parse_args(["--sample-index", "not-an-int"])


def test_group_index_must_be_int():
    p = build_parser()
    args = p.parse_args(["--group-index", "2"])
    assert args.group_index == 2


def test_aligner_choices_are_restricted():
    p = build_parser()
    args = p.parse_args(["--aligner", "hisat2"])
    assert args.aligner == "hisat2"
    with pytest.raises(SystemExit):
        p.parse_args(["--aligner", "bwa"])


def test_force_and_overwrite_are_aliases():
    p = build_parser()
    assert p.parse_args(["--force"]).force is True
    assert p.parse_args(["--overwrite"]).force is True
    assert p.parse_args([]).force is False


def test_star_multimap_order_choices():
    p = build_parser()
    args = p.parse_args(["--star-multimap-order", "Random"])
    assert args.star_multimap_order == "Random"
    with pytest.raises(SystemExit):
        p.parse_args(["--star-multimap-order", "bogus"])


def test_hisat2_strandness_choices():
    p = build_parser()
    args = p.parse_args(["--hisat2-strandness", "FR"])
    assert args.hisat2_strandness == "FR"
    with pytest.raises(SystemExit):
        p.parse_args(["--hisat2-strandness", "bogus"])


# ── main() fast paths (no pipeline execution) ────────────────────────────────

def test_count_samples_on_empty_project_prints_zero(project_dir, capsys):
    rc = main(["--project", str(project_dir), "--count-samples"])
    assert rc == 0
    out = capsys.readouterr().out.strip()
    assert out == "0"


def test_count_samples_counts_staged_r1_files(project_dir, capsys, make_fastq):
    rawdata = project_dir / "homo_sapiens" / "RawData"
    rawdata.mkdir(parents=True)
    make_fastq(rawdata / "homo_sapiens_K562_csRNA_r1_R1.fastq.gz")
    make_fastq(rawdata / "homo_sapiens_K562_csRNA_r2_R1.fastq.gz")

    rc = main(["--project", str(project_dir), "--count-samples"])
    assert rc == 0
    assert capsys.readouterr().out.strip() == "2"


def test_count_groups_on_empty_project_prints_zero(project_dir, capsys):
    rc = main(["--project", str(project_dir), "--count-groups"])
    assert rc == 0
    assert capsys.readouterr().out.strip() == "0"


def test_count_groups_counts_distinct_species_sample_pairs(project_dir, capsys, make_fastq):
    rawdata = project_dir / "homo_sapiens" / "RawData"
    rawdata.mkdir(parents=True)
    # two replicates of the SAME sample -> still one group
    make_fastq(rawdata / "homo_sapiens_K562_csRNA_r1_R1.fastq.gz")
    make_fastq(rawdata / "homo_sapiens_K562_csRNA_r2_R1.fastq.gz")

    rc = main(["--project", str(project_dir), "--count-groups"])
    assert rc == 0
    assert capsys.readouterr().out.strip() == "1"


def test_check_rerun_on_empty_project_exits_zero(project_dir):
    # Nothing staged, nothing produced -> not a "would touch nothing new"
    # situation (there's nothing existing to conflict with), so this must
    # NOT be blocked.
    rc = main(["--project", str(project_dir), "--check-rerun"])
    assert rc == 0


def test_check_rerun_blocks_when_everything_already_finished(project_dir, make_fastq):
    rawdata = project_dir / "homo_sapiens" / "RawData"
    rawdata.mkdir(parents=True)
    make_fastq(rawdata / "homo_sapiens_K562_csRNA_r1_R1.fastq.gz")
    qc_dir = project_dir / "homo_sapiens" / "QC" / "K562"
    qc_dir.mkdir(parents=True)
    (qc_dir / "qc_report.html").write_text("<html></html>")
    # Some other output dir must also exist for find_existing_outputs() to
    # register anything at all (QC/ itself is one of OUTPUT_DIR_NAMES).

    rc = main(["--project", str(project_dir), "--check-rerun"])
    assert rc == 1


def test_check_rerun_force_bypasses_the_block(project_dir, make_fastq):
    rawdata = project_dir / "homo_sapiens" / "RawData"
    rawdata.mkdir(parents=True)
    make_fastq(rawdata / "homo_sapiens_K562_csRNA_r1_R1.fastq.gz")
    qc_dir = project_dir / "homo_sapiens" / "QC" / "K562"
    qc_dir.mkdir(parents=True)
    (qc_dir / "qc_report.html").write_text("<html></html>")

    rc = main(["--project", str(project_dir), "--check-rerun", "--force"])
    assert rc == 0


def test_stage_raw_moves_loose_fastqs_and_exits(project_dir, make_fastq):
    make_fastq(project_dir / "homo_sapiens_K562_csRNA_r1_R1.fastq.gz")
    rc = main(["--project", str(project_dir), "--stage-raw"])
    assert rc == 0
    dest = project_dir / "homo_sapiens" / "RawData" / "homo_sapiens_K562_csRNA_r1_R1.fastq.gz"
    assert dest.is_file()
    assert (project_dir / "config.txt").is_file()


def test_sample_index_out_of_range_returns_nonzero(project_dir, make_fastq):
    # run_trim (and friends) only bounds-check --sample-index AFTER
    # confirming list_r1(cfg) is non-empty (an empty list returns early with
    # a log message, never reaching the bounds check) — so this needs at
    # least one real staged R1 file to actually exercise the IndexError path.
    rawdata = project_dir / "homo_sapiens" / "RawData"
    rawdata.mkdir(parents=True)
    make_fastq(rawdata / "homo_sapiens_K562_csRNA_r1_R1.fastq.gz")

    rc = main([
        "--project", str(project_dir),
        "--skip-prepare",
        "--steps", "trim",
        "--sample-index", "99",
    ])
    assert rc == 1
