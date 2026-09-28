"""Metadata-backed discovery, staging, grouping, and TSS integration tests."""
from __future__ import annotations

import csv

import pytest

from homerun import prepare, tagdirs, tss
from homerun.utils import assay_for_fastq, iter_leaf_dirs, parse_sample_name


HEADERS = ["FASTQ", "Species", "Sample", "Condition", "Assay", "Replicate"]


def _manifest(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADERS)
        writer.writerows(rows)


def test_metadata_stages_filename_without_species(make_cfg, make_fastq, project_dir):
    name = "K562_D2_csRNA_r1_R1.fastq.gz"
    manifest = project_dir / "samples.csv"
    _manifest(manifest, [[name, "homo_sapiens", "K562", "D2", "csRNA", "r1"]])
    make_fastq(project_dir / name)
    cfg = make_cfg(metadata=str(manifest))

    prepare.stage_loose_fastqs(cfg)

    staged = project_dir / "homo_sapiens" / "RawData" / name
    assert staged.is_file()
    assert list(iter_leaf_dirs(cfg)) == [
        ("homo_sapiens", "K562", "D2_csRNA_r1", staged)
    ]
    assert assay_for_fastq(cfg, staged) == "csRNA"
    # Without a manifest, the unchanged legacy parser interprets the first
    # two filename tokens as species; it cannot recover the intended species.
    assert parse_sample_name(name) != ("homo_sapiens", "K562", "D2_csRNA_r1")


def test_discovered_fastq_missing_metadata_row_raises(make_cfg, make_fastq, project_dir):
    listed = "listed_R1.fastq.gz"
    unlisted = "unlisted_R1.fastq.gz"
    manifest = project_dir / "samples.csv"
    _manifest(manifest, [[listed, "homo_sapiens", "K562", "", "csRNA", "r1"]])
    make_fastq(project_dir / listed)
    make_fastq(project_dir / unlisted)
    cfg = make_cfg(metadata=str(manifest))

    with pytest.raises(ValueError, match="missing from metadata.*unlisted_R1"):
        prepare.stage_loose_fastqs(cfg)


def test_metadata_fastq_that_cannot_be_found_raises(make_cfg, project_dir):
    manifest = project_dir / "samples.csv"
    _manifest(manifest, [["absent_R1.fastq.gz", "homo_sapiens", "K562", "", "csRNA", "r1"]])
    cfg = make_cfg(metadata=str(manifest))

    with pytest.raises(ValueError, match="do(es)? not match any discovered FASTQ.*absent_R1"):
        prepare.stage_loose_fastqs(cfg)


def test_config_summary_records_metadata_identity(make_cfg, make_fastq, project_dir):
    name = "opaque_R1.fastq.gz"
    manifest = project_dir / "samples.csv"
    _manifest(manifest, [[name, "apis_mellifera", "", "", "csRNA", "r1"]])
    make_fastq(project_dir / name)
    cfg = make_cfg(metadata=str(manifest))
    prepare.stage_loose_fastqs(cfg)
    prepare.write_config_summary(cfg)

    summary = (project_dir / "config.txt").read_text()
    assert f"metadata = {manifest}" in summary
    assert "apis_mellifera/apis_mellifera" in summary


def test_metadata_conditions_remain_separate_combo_groups(make_cfg, make_fastq, project_dir):
    rows = []
    for condition in ("D2_Pre", "D2_Post"):
        for replicate in ("r1", "r2"):
            name = f"{condition}_{replicate}_R1.fastq.gz"
            rows.append([name, "homo_sapiens", "K562", condition, "csRNA", replicate])
            make_fastq(project_dir / name)
    manifest = project_dir / "samples.csv"
    _manifest(manifest, rows)
    cfg = make_cfg(metadata=str(manifest))
    prepare.stage_loose_fastqs(cfg)

    groups = tagdirs._combo_groups(cfg)
    assert set(groups) == {
        ("homo_sapiens", "K562", "D2_Pre_csRNA"),
        ("homo_sapiens", "K562", "D2_Post_csRNA"),
    }
    assert all(len(members) == 2 for members in groups.values())


def test_tss_matches_metadata_condition_and_replicate(
    make_cfg, make_fastq, project_dir, monkeypatch
):
    rows = []
    for assay in ("csRNA", "sRNA"):
        name = f"input_{assay}_R1.fastq.gz"
        rows.append([name, "homo_sapiens", "K562", "D2_Post", assay, "r1"])
        make_fastq(project_dir / name)
    manifest = project_dir / "samples.csv"
    _manifest(manifest, rows)
    cfg = make_cfg(metadata=str(manifest), threads=1)
    prepare.stage_loose_fastqs(cfg)

    for species, sample, leaf, _r1 in iter_leaf_dirs(cfg):
        cfg.leaf_tagdir(species, sample, leaf).mkdir(parents=True)

    commands = []
    monkeypatch.setattr(tss, "run", lambda command, **kwargs: commands.append(command))
    tss.run_tss(cfg)

    assert len(commands) == 1
    assert "K562_D2_Post_csRNA_r1" in commands[0]
    assert "K562_D2_Post_sRNA_r1" in commands[0]
