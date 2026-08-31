"""Tests for compact species-level QC summaries."""
from __future__ import annotations

import pandas as pd
import pytest

from homerun.qc import _sample_tss_files, qc_nucleotide_divergence_heatmaps


def _write_tag_freq(path, adjustment: float) -> None:
    path.mkdir(parents=True)
    pd.DataFrame({
        "Offset": [-1, 0, 1],
        "A": [0.20 + adjustment, 0.25 + adjustment, 0.30 + adjustment],
        "C": [0.30 - adjustment, 0.25 - adjustment, 0.20 - adjustment],
        "G": [0.25, 0.30 + adjustment, 0.25],
        "T": [0.25, 0.20, 0.25 - adjustment],
    }).to_csv(path / "tagFreqUniq.txt", sep="\t", index=False)


def test_nucleotide_divergence_heatmaps_cover_all_bases_and_samples(cfg):
    samples = [("homo_sapiens", "K562"), ("homo_sapiens", "IMR90")]
    _write_tag_freq(cfg.combo_tagdir("homo_sapiens", "K562", "csRNA"), 0.02)
    _write_tag_freq(cfg.combo_tagdir("homo_sapiens", "IMR90", "csRNA"), -0.02)

    qc_nucleotide_divergence_heatmaps(cfg, samples)

    qc_root = cfg.species_qc("homo_sapiens")
    assert (qc_root / "csRNA_nucleotide_divergence_heatmap.png").is_file()
    assert (qc_root / "csRNA_nucleotide_divergence_heatmap.svg").is_file()
    for nucleotide in "ACGT":
        data_path = qc_root / f"csRNA_{nucleotide}_nucleotide_divergence.tsv"
        assert data_path.is_file()
        data = pd.read_csv(data_path, sep="\t", index_col=0)
        assert set(data.index) == {"K562", "IMR90"}
        assert data.to_numpy().mean() == pytest.approx(0.0, abs=1e-12)


def test_nucleotide_divergence_heatmaps_skip_missing_assay(cfg):
    qc_nucleotide_divergence_heatmaps(cfg, [("homo_sapiens", "K562")])
    assert not cfg.species_qc("homo_sapiens").exists()


def test_shared_tss_directory_does_not_prefix_match_other_sample(cfg):
    tss_dir = cfg.sample_tss("homo_sapiens", "K5")
    tss_dir.mkdir(parents=True)
    (tss_dir / "K5.tss.txt").write_text("K5")
    (tss_dir / "K562.tss.txt").write_text("K562")

    found = _sample_tss_files(cfg, "homo_sapiens", "K5", ".tss.txt")
    assert [path.name for path in found] == ["K5.tss.txt"]
