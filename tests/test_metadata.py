"""Unit tests for CSV/XLSX sample metadata loading and validation."""
from __future__ import annotations

import csv

import pytest

from homerun.metadata import load_sample_metadata


HEADERS = ["FASTQ", "Species", "Sample", "Condition", "Assay", "Replicate"]


def _write_csv(path, rows, headers=HEADERS):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)


def test_load_csv_metadata_builds_identity(tmp_path):
    path = tmp_path / "samples.csv"
    _write_csv(path, [["K562_D2_csRNA_r1_R1.fastq.gz", "HOMO_SAPIENS", "K562", "D2", "csRNA", "r1"]])

    record = load_sample_metadata(path)["K562_D2_csRNA_r1_R1.fastq.gz"]
    assert record == {
        "species": "homo_sapiens",
        "sample": "K562",
        "condition": "D2",
        "assay": "csRNA",
        "replicate": "r1",
        "leaf_name": "D2_csRNA_r1",
    }


def test_load_xlsx_metadata_builds_same_identity(tmp_path):
    from openpyxl import Workbook

    path = tmp_path / "samples.xlsx"
    workbook = Workbook()
    sheet = workbook.active
    sheet.append(HEADERS)
    sheet.append(["K562_D2_csRNA_r1_R1.fastq.gz", "homo_sapiens", "K562", "D2", "csRNA", "r1"])
    workbook.save(path)
    workbook.close()

    record = load_sample_metadata(path)["K562_D2_csRNA_r1_R1.fastq.gz"]
    assert (record["species"], record["sample"], record["leaf_name"]) == (
        "homo_sapiens", "K562", "D2_csRNA_r1"
    )


def test_blank_sample_falls_back_to_species(tmp_path):
    path = tmp_path / "samples.csv"
    _write_csv(path, [["bee_csRNA_r1_R1.fastq.gz", "apis_mellifera", "", "", "csRNA", "r1"]])
    record = load_sample_metadata(path)["bee_csRNA_r1_R1.fastq.gz"]
    assert record["sample"] == "apis_mellifera"
    assert record["leaf_name"] == "csRNA_r1"


def test_blank_condition_omits_condition_prefix(tmp_path):
    path = tmp_path / "samples.csv"
    _write_csv(path, [["K562_csRNA_r1_R1.fastq.gz", "homo_sapiens", "K562", "", "csRNA", "r1"]])
    assert load_sample_metadata(path)["K562_csRNA_r1_R1.fastq.gz"]["leaf_name"] == "csRNA_r1"


def test_duplicate_fastq_rows_raise(tmp_path):
    path = tmp_path / "samples.csv"
    row = ["K562_csRNA_r1_R1.fastq.gz", "homo_sapiens", "K562", "", "csRNA", "r1"]
    _write_csv(path, [row, row])
    with pytest.raises(ValueError, match="Duplicate FASTQ metadata entry"):
        load_sample_metadata(path)


def test_invalid_assay_raises(tmp_path):
    path = tmp_path / "samples.csv"
    _write_csv(path, [["K562_chip_r1_R1.fastq.gz", "homo_sapiens", "K562", "", "ChIP", "r1"]])
    with pytest.raises(ValueError, match="Invalid assay 'ChIP'"):
        load_sample_metadata(path)


def test_missing_required_columns_raise(tmp_path):
    path = tmp_path / "samples.csv"
    _write_csv(path, [["x.fastq.gz", "homo_sapiens"]], headers=["FASTQ", "Species"])
    with pytest.raises(ValueError, match="missing required column.*Assay.*Replicate"):
        load_sample_metadata(path)


def test_duplicate_identity_raises_but_r1_r2_pair_is_allowed(tmp_path):
    path = tmp_path / "collision.csv"
    shared = ["homo_sapiens", "K562", "", "totalRNA", "r1"]
    _write_csv(path, [
        ["library_lane1_R1.fastq.gz", *shared],
        ["library_lane2_R1.fastq.gz", *shared],
    ])
    with pytest.raises(ValueError, match="Duplicate metadata identity would collide"):
        load_sample_metadata(path)

    pair_path = tmp_path / "pair.csv"
    _write_csv(pair_path, [
        ["library_R1.fastq.gz", *shared],
        ["library_R2.fastq.gz", *shared],
    ])
    assert len(load_sample_metadata(pair_path)) == 2
