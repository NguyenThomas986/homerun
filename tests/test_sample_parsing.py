"""Sample/input tests for v1.0.0's filename-based identity system —
parse_sample_name, assay_of_leaf, seq_type, and R2 mate-pairing.

Explicitly NOT covered here (per v1.1.0 scope): any Excel/CSV/TSV manifest
input. v1.0.0 identifies every sample purely from its FASTQ filename.
"""
from __future__ import annotations

import pytest

from homerun.utils import (
    parse_sample_name,
    assay_of_leaf,
    seq_type,
    replicate_of_leaf,
    find_r2_for_r1,
    _find_assay,
)


# ── parse_sample_name: species/sample/leaf extraction ────────────────────────

def test_parse_standard_csrna_filename():
    species, sample, leaf = parse_sample_name("homo_sapiens_K562_csRNA_r1_R1.fastq.gz")
    assert species == "homo_sapiens"
    assert sample == "K562"
    assert leaf == "csRNA_r1"


def test_parse_srna_filename():
    species, sample, leaf = parse_sample_name("homo_sapiens_K562_sRNA_r2_R1.fastq.gz")
    assert species == "homo_sapiens"
    assert sample == "K562"
    assert leaf == "sRNA_r2"


def test_parse_totalrna_filename():
    species, sample, leaf = parse_sample_name("homo_sapiens_K562_totalRNA_r1_R1.fastq.gz")
    assert leaf == "totalRNA_r1"


def test_parse_generic_rna_token_maps_to_totalrna():
    # bare "RNA" (not csRNA/sRNA/totalRNA) falls back to totalRNA per
    # _find_assay's second pass.
    species, sample, leaf = parse_sample_name("homo_sapiens_IMR90_RNA-r1_ENCFF000HAZ_R1.fastq.gz")
    assert species == "homo_sapiens"
    assert sample == "IMR90"
    assert leaf == "RNA_r1"
    assert assay_of_leaf(leaf) == "totalRNA"


def test_parse_condition_token_before_assay():
    species, sample, leaf = parse_sample_name("homo_sapiens_K562_p53KO_csRNA_r1_R1.fastq.gz")
    assert species == "homo_sapiens"
    assert sample == "K562"
    assert leaf == "p53KO_csRNA_r1"
    assert assay_of_leaf(leaf) == "csRNA"


def test_parse_no_distinct_sample_token_reuses_species():
    # e.g. Apis_mellifera_csRNA_r1... -> species reused as sample so
    # same-species replicates across assays still land in one Species/Sample/.
    species, sample, leaf = parse_sample_name("Apis_mellifera_csRNA_r1_R1.fastq.gz")
    assert species == "apis_mellifera"
    assert sample == "apis_mellifera"
    assert leaf == "csRNA_r1"


def test_parse_hyphenated_replicate_marker_equivalent_to_underscore():
    a = parse_sample_name("homo_sapiens_K562_csRNA-r1_R1.fastq.gz")
    b = parse_sample_name("homo_sapiens_K562_csRNA_r1_R1.fastq.gz")
    assert a == b


def test_parse_repN_style_marker():
    species, sample, leaf = parse_sample_name("homo_sapiens_K562_csRNA_rep3_R1.fastq.gz")
    assert leaf == "csRNA_rep3"


def test_parse_species_is_lowercased():
    species, _sample, _leaf = parse_sample_name("HOMO_SAPIENS_K562_csRNA_r1_R1.fastq.gz")
    assert species == "homo_sapiens"


def test_parse_everything_after_replicate_marker_is_discarded():
    # Lane/index/accession metadata after the replicate marker must not
    # affect the parsed identity — this is exactly what makes
    # find_r2_for_r1's identity-matching approach work.
    a = parse_sample_name("homo_sapiens_K562_csRNA_r1_L001_R1_001.fastq.gz")
    b = parse_sample_name("homo_sapiens_K562_csRNA_r1_R1.fastq.gz")
    assert a == b


# ── error cases ───────────────────────────────────────────────────────────────

def test_parse_raises_without_replicate_marker():
    with pytest.raises(ValueError, match="no replicate marker"):
        parse_sample_name("homo_sapiens_K562_csRNA_R1.fastq.gz")


def test_parse_raises_without_assay_token():
    with pytest.raises(ValueError, match="no assay type"):
        parse_sample_name("homo_sapiens_K562_r1_R1.fastq.gz")


def test_parse_raises_with_too_few_tokens_before_marker():
    with pytest.raises(ValueError):
        parse_sample_name("r1_R1.fastq.gz")


def test_illumina_r1_tag_is_not_a_replicate_marker():
    """Regression test: only a LOWERCASE 'r1'/'rep2'-style token counts as a
    replicate marker. Illumina's own read-tag 'R1' is conventionally
    uppercase and must not be mistaken for one, or a filename with no real
    replicate marker before its R1 read tag would silently misparse."""
    with pytest.raises(ValueError):
        # No lowercase replicate marker anywhere in this filename — must
        # raise, not silently treat the uppercase "R1" tag as a match.
        parse_sample_name("homo_sapiens_K562_csRNA_R1.fastq.gz")


# ── assay_of_leaf / seq_type / replicate_of_leaf ─────────────────────────────

@pytest.mark.parametrize("leaf,expected", [
    ("csRNA_r1", "csRNA"),
    ("sRNA_r2", "sRNA"),
    ("totalRNA_r1", "totalRNA"),
    ("p53KO_csRNA_r1", "csRNA"),
    ("csRNAseq_r2", "csRNA"),
    ("not_a_real_assay_r1", None),
])
def test_assay_of_leaf(leaf, expected):
    assert assay_of_leaf(leaf) == expected


@pytest.mark.parametrize("name,expected", [
    ("sample_csRNA_r1_R1.fastq.gz", "csRNA"),
    ("sample_sRNA_r1_R1.fastq.gz", "sRNA"),
    ("sample_totalRNA_r1_R1.fastq.gz", "totalRNA"),
    ("sample_RNA_r1_R1.fastq.gz", "totalRNA"),
    ("sample_unknown_r1_R1.fastq.gz", None),
])
def test_seq_type(name, expected):
    assert seq_type(name) == expected


def test_replicate_of_leaf():
    assert replicate_of_leaf("csRNA_r1") == "r1"
    assert replicate_of_leaf("p53KO_csRNA_rep12") == "rep12"
    assert replicate_of_leaf("csRNA") is None


def test_find_assay_prefers_specific_form_over_generic_rna():
    idx, assay = _find_assay(["totalRNA", "r1"])
    assert assay == "totalRNA"
    assert idx == 0


# ── R2 mate-pairing (identity-based, not filename substitution) ─────────────

def test_find_r2_matches_same_basename_mates(paired_totalrna_sample):
    # sanity check on the simple/typical case even though the fixture uses
    # mismatched accessions on purpose for the ENCODE regression below
    _cfg, r1, r2 = paired_totalrna_sample
    found = find_r2_for_r1(r1)
    assert found == r2


def test_find_r2_handles_encode_style_mismatched_accessions(make_cfg, make_fastq):
    """Regression test: r1.name.replace("_R1", "_R2") assumed both mates
    share a basename apart from the R1/R2 tag. ENCODE downloads break that
    assumption — each mate has its own independent accession number. This
    verifies find_r2_for_r1 correctly matches by parsed identity instead."""
    cfg = make_cfg()
    rawdata = cfg.rawdata_dir("homo_sapiens", "IMR90")
    rawdata.mkdir(parents=True)
    r1 = rawdata / "homo_sapiens_IMR90_RNA-r1_ENCFF000HAZ_R1.fastq.gz"
    r2 = rawdata / "homo_sapiens_IMR90_RNA-r1_ENCFF000HBG_R2.fastq.gz"
    make_fastq(r1)
    make_fastq(r2)

    found = find_r2_for_r1(r1)
    assert found == r2

    # the naive substitution would have looked for a file that never exists
    naive_guess = rawdata / r1.name.replace("_R1", "_R2")
    assert naive_guess != r2
    assert not naive_guess.exists()


def test_find_r2_does_not_cross_match_between_replicates(make_cfg, make_fastq):
    """Two replicates sharing one flat RawData/ dir must each find their OWN
    mate, never the other replicate's."""
    cfg = make_cfg()
    rawdata = cfg.rawdata_dir("homo_sapiens", "IMR90")
    rawdata.mkdir(parents=True)
    files = {
        "r1_r1": rawdata / "homo_sapiens_IMR90_RNA-r1_ENCFF000HAZ_R1.fastq.gz",
        "r1_r2": rawdata / "homo_sapiens_IMR90_RNA-r1_ENCFF000HBG_R2.fastq.gz",
        "r2_r1": rawdata / "homo_sapiens_IMR90_RNA-r2_ENCFF000HBA_R1.fastq.gz",
        "r2_r2": rawdata / "homo_sapiens_IMR90_RNA-r2_ENCFF000HBI_R2.fastq.gz",
    }
    for p in files.values():
        make_fastq(p)

    assert find_r2_for_r1(files["r1_r1"]) == files["r1_r2"]
    assert find_r2_for_r1(files["r2_r1"]) == files["r2_r2"]


def test_find_r2_returns_none_when_no_mate_exists(make_cfg, make_fastq):
    cfg = make_cfg()
    rawdata = cfg.rawdata_dir("homo_sapiens", "IMR90")
    rawdata.mkdir(parents=True)
    r1 = rawdata / "homo_sapiens_IMR90_totalRNA_r1_R1.fastq.gz"
    make_fastq(r1)
    assert find_r2_for_r1(r1) is None


def test_find_r2_returns_none_on_ambiguous_multiple_matches(make_cfg, make_fastq):
    cfg = make_cfg()
    rawdata = cfg.rawdata_dir("homo_sapiens", "IMR90")
    rawdata.mkdir(parents=True)
    r1 = rawdata / "homo_sapiens_IMR90_totalRNA_r1_R1.fastq.gz"
    # two candidate R2 files that would BOTH parse to the same identity as r1
    dup1 = rawdata / "homo_sapiens_IMR90_totalRNA_r1_R2.fastq.gz"
    dup2 = rawdata / "homo_sapiens_IMR90_totalRNA_r1_lane2_R2.fastq"
    for p in (r1, dup1, dup2):
        make_fastq(p)
    assert find_r2_for_r1(r1) is None
