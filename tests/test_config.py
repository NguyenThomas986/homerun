"""Configuration tests: Config defaults, CLI-flag > env-var > default
precedence, path-builder correctness, and validation of the documented
"no genome is assumed" behavior (validated at the point of use by
mapping._check_index / prepare.ensure_starindex, not by Config itself —
Config never raises on construction, since threads/samples/etc. all still
need a valid Config to run --count-samples style preflight checks even
before a genome is configured).
"""
from __future__ import annotations

import argparse

import pytest

from homerun.config import Config, load_config, _pick


# ── defaults ──────────────────────────────────────────────────────────────────

def test_config_requires_only_project(tmp_path):
    cfg = Config(project=tmp_path)
    assert cfg.aligner == "star"
    assert cfg.genome_index == ""
    assert cfg.genome == ""
    assert cfg.threads == 20
    assert cfg.force is False


def test_load_config_with_no_args_uses_cwd_or_env(monkeypatch, tmp_path):
    monkeypatch.delenv("CSRNA_PROJECT", raising=False)
    monkeypatch.chdir(tmp_path)
    cfg = load_config(None)
    assert cfg.project == tmp_path.resolve()


def test_load_config_project_flag_overrides_cwd(tmp_path, monkeypatch):
    other = tmp_path / "elsewhere"
    other.mkdir()
    monkeypatch.delenv("CSRNA_PROJECT", raising=False)
    args = argparse.Namespace(project=str(other))
    cfg = load_config(args)
    assert cfg.project == other.resolve()


# ── env var precedence ───────────────────────────────────────────────────────

def test_env_var_used_when_no_cli_flag(tmp_path, monkeypatch):
    monkeypatch.setenv("CSRNA_GENOME", "mm10")
    args = argparse.Namespace(project=str(tmp_path), genome=None)
    cfg = load_config(args)
    assert cfg.genome == "mm10"


def test_cli_flag_overrides_env_var(tmp_path, monkeypatch):
    monkeypatch.setenv("CSRNA_GENOME", "mm10")
    args = argparse.Namespace(project=str(tmp_path), genome="hg38")
    cfg = load_config(args)
    assert cfg.genome == "hg38"


def test_default_used_when_neither_cli_nor_env_set(tmp_path, monkeypatch):
    monkeypatch.delenv("CSRNA_ALIGNER", raising=False)
    args = argparse.Namespace(project=str(tmp_path), aligner=None)
    cfg = load_config(args)
    assert cfg.aligner == "star"


def test_pick_precedence_helper_directly():
    args = argparse.Namespace(foo="from_cli")
    assert _pick(args, "foo", "SOME_ENV_VAR_THAT_DOES_NOT_EXIST", "default") == "from_cli"
    args2 = argparse.Namespace(foo=None)
    assert _pick(args2, "foo", "SOME_ENV_VAR_THAT_DOES_NOT_EXIST", "default") == "default"
    assert _pick(None, "foo", "SOME_ENV_VAR_THAT_DOES_NOT_EXIST", "default") == "default"


def test_threads_falls_back_to_slurm_cpus_per_task(tmp_path, monkeypatch):
    monkeypatch.delenv("CSRNA_THREADS", raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    args = argparse.Namespace(project=str(tmp_path), threads=None)
    cfg = load_config(args)
    assert cfg.threads == 8


def test_aligner_is_lowercased(tmp_path):
    args = argparse.Namespace(project=str(tmp_path), aligner="STAR")
    cfg = load_config(args)
    assert cfg.aligner == "star"


def test_force_flag_wired_through_load_config(tmp_path):
    """Regression test: Config previously had no `force` field at all, so
    prepare.prepare()'s `getattr(cfg, "force", False)` check always
    evaluated False and --force never actually wiped outputs (see
    prepare.wipe_outputs and pipeline.main). Fixed by adding Config.force
    and wiring it through load_config()."""
    args = argparse.Namespace(project=str(tmp_path), force=True)
    cfg = load_config(args)
    assert cfg.force is True

    args_no_force = argparse.Namespace(project=str(tmp_path), force=False)
    cfg2 = load_config(args_no_force)
    assert cfg2.force is False


# ── path builders ─────────────────────────────────────────────────────────────

def test_sample_dir_layout(tmp_path):
    cfg = Config(project=tmp_path)
    assert cfg.species_dir("homo_sapiens") == tmp_path / "homo_sapiens"
    assert cfg.sample_dir("homo_sapiens", "K562") == tmp_path / "homo_sapiens"


def test_flat_species_category_dirs(tmp_path):
    cfg = Config(project=tmp_path)
    base = tmp_path / "homo_sapiens"
    assert cfg.rawdata_dir("homo_sapiens", "K562") == base / "RawData"
    assert cfg.trimmed_dir("homo_sapiens", "K562") == base / "Trimmed"
    assert cfg.aligned_dir("homo_sapiens", "K562") == base / "Aligned"
    assert cfg.species_qc("homo_sapiens") == base / "QC"
    assert cfg.sample_qc("homo_sapiens", "K562") == base / "QC" / "K562"
    assert cfg.sample_tss("homo_sapiens", "K562") == base / "TSS"


def test_tagdir_names_are_sample_prefixed(tmp_path):
    cfg = Config(project=tmp_path)
    leaf = cfg.leaf_tagdir("homo_sapiens", "IMR90", "csRNA_r1")
    combo = cfg.combo_tagdir("homo_sapiens", "IMR90", "csRNA")
    assert leaf.name == "IMR90_csRNA_r1"
    assert combo.name == "IMR90_csRNA-combo"
    assert leaf.parent.name == "TagDirs"
    assert combo.parent.name == "TagDirs"


def test_bedgraph_names_match_tagdir_naming(tmp_path):
    cfg = Config(project=tmp_path)
    leaf = cfg.leaf_bedgraph("homo_sapiens", "IMR90", "csRNA_r1")
    combo = cfg.combo_bedgraph("homo_sapiens", "IMR90", "csRNA")
    assert leaf.name == "IMR90_csRNA_r1"
    assert combo.name == "IMR90_csRNA-combo"


def test_ritrie_paths(tmp_path):
    cfg = Config(project=tmp_path)
    species_exons = cfg.species_ritrie_gtf_exons("homo_sapiens")
    assert species_exons == tmp_path / "homo_sapiens" / "RITRIE" / "parsed_gtf_exons.tsv"
    leaf = cfg.leaf_ritrie("homo_sapiens", "IMR90", "csRNA_r1")
    assert leaf == tmp_path / "homo_sapiens" / "RITRIE" / "IMR90_csRNA_r1"


def test_logs_dir_and_starindex_properties(tmp_path):
    cfg = Config(project=tmp_path, genome_index=str(tmp_path / "idx"))
    assert cfg.logs_dir == tmp_path / "logs"
    assert cfg.starindex == tmp_path / "idx"


def test_relative_project_path_is_resolved_to_absolute(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "sub").mkdir()
    args = argparse.Namespace(project="sub")
    cfg = load_config(args)
    assert cfg.project.is_absolute()
    assert cfg.project == (tmp_path / "sub").resolve()
