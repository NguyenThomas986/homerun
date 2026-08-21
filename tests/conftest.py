"""Shared fixtures for the homerun test suite.

Design notes:
- Every fixture that touches the filesystem uses pytest's tmp_path, never a
  real/absolute path — tests must be runnable on any machine, in any repo
  checkout location, with no network and no HOMER/STAR/HISAT2/skewer
  installed.
- `cfg`/`make_cfg` build a real homerun.config.Config pointed at a tmp_path
  project root. No subprocess, filesystem-outside-tmp_path, or SLURM access
  is ever needed to construct one.
- Fixtures that need external tool output (FASTQ, SAM, HOMER tag directory
  files) build minimal, syntactically-valid fakes by hand rather than
  shelling out to real tools.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Make sure "homerun" resolves to the package under test even if the test
# runner's cwd/sys.path doesn't already include it (pip install -e . also
# achieves this, but conftest doing it too means `pytest` "just works" from
# a fresh checkout without requiring the editable install first).
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from homerun.config import Config  # noqa: E402


@pytest.fixture
def project_dir(tmp_path) -> Path:
    """Empty project root directory."""
    d = tmp_path / "project"
    d.mkdir()
    return d


@pytest.fixture
def make_cfg(project_dir):
    """Factory for a Config pointed at project_dir, with sane, harmless
    defaults (genome/genome_index set to dummy strings so code paths that
    just read cfg.genome/cfg.genome_index don't need a real index) and any
    field overridable via kwargs."""

    def _make(**overrides) -> Config:
        kwargs = dict(
            project=project_dir,
            genome="hg38",
            genome_index=str(project_dir / "fake_index"),
            aligner="star",
        )
        kwargs.update(overrides)
        return Config(**kwargs)

    return _make


@pytest.fixture
def cfg(make_cfg) -> Config:
    """A default Config for tests that don't care about specific overrides."""
    return make_cfg()


def _write_fastq(path: Path, n_reads: int = 4) -> None:
    """Minimal but syntactically valid FASTQ content — enough for any code
    path that just checks a file's existence/non-emptiness (utils.done())
    or parses its filename; nothing in this suite runs a real aligner
    against these reads."""
    lines = []
    for i in range(n_reads):
        lines += [f"@read{i}", "ACGT", "+", "IIII"]
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def make_fastq():
    """Factory: make_fastq(path, n_reads=4) writes a minimal valid FASTQ."""
    return _write_fastq


@pytest.fixture
def staged_sample(make_cfg, make_fastq):
    """A Config plus one already-staged csRNA replicate
    (Species/Sample/RawData/<file>_R1.fastq.gz), for tests that need at
    least one real, discoverable sample without exercising prepare.py's
    staging logic itself."""
    cfg = make_cfg()
    r1_name = "homo_sapiens_K562_csRNA_r1_R1.fastq.gz"
    rawdata = cfg.rawdata_dir("homo_sapiens", "K562")
    rawdata.mkdir(parents=True)
    r1 = rawdata / r1_name
    make_fastq(r1)
    return cfg, r1


@pytest.fixture
def paired_totalrna_sample(make_cfg, make_fastq):
    """A Config plus one staged totalRNA replicate with BOTH mates present,
    using ENCODE-style mismatched accessions on purpose — this is the exact
    shape of input that once broke naive '_R1'->'_R2' filename substitution
    (see utils.find_r2_for_r1's docstring and tests/test_utils.py's
    regression test for it)."""
    cfg = make_cfg()
    rawdata = cfg.rawdata_dir("homo_sapiens", "IMR90")
    rawdata.mkdir(parents=True)
    r1 = rawdata / "homo_sapiens_IMR90_totalRNA_r1_ENCFF000HAZ_R1.fastq.gz"
    r2 = rawdata / "homo_sapiens_IMR90_totalRNA_r1_ENCFF000HBG_R2.fastq.gz"
    make_fastq(r1)
    make_fastq(r2)
    return cfg, r1, r2


@pytest.fixture
def non_interactive_matplotlib():
    """Force the non-interactive 'Agg' backend before any test imports a
    plotting module, so QC/plotting tests never try to open a GUI window
    (would hang/fail in a headless CI environment)."""
    import matplotlib
    matplotlib.use("Agg", force=True)
    yield
    # Nothing to tear down — backend selection is process-global and Agg is
    # a safe default for the rest of the test session either way.