"""Static safeguards for the SLURM controller and worker wrappers."""
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
WORKERS = (
    "prepare.sbatch",
    "align_array.sbatch",
    "tagdir_array.sbatch",
    "tagdirs_combo_array.sbatch",
    "bedgraphs_array.sbatch",
    "tss_array.sbatch",
    "collect.sbatch",
)


def _text(name: str) -> str:
    return (ROOT / name).read_text(encoding="utf-8")


def test_slurm_scripts_import_src_layout():
    for name in ("submit_array.sh", *WORKERS):
        assert '/src:${PYTHONPATH:-}' in _text(name)


def test_controller_sizes_arrays_after_prepare_and_fails_closed():
    script = _text("submit_array.sh")
    assert script.index("sbatch --wait --parsable") < script.index("--count-samples")
    assert "--dependency=afterany" not in script
    assert '--array=0-$((S-1))%"${TSS_THROTTLE}"' in script


def test_collect_never_reruns_prepare():
    script = _text("collect.sbatch")
    assert "--skip-prepare" in script


def test_all_worker_files_are_tracked_by_gitignore_policy():
    ignore = _text(".gitignore")
    assert "*.sbatch" not in ignore
    assert "*.sh" not in ignore
