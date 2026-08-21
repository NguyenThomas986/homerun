"""Package/import tests: homerun and its important public modules import
cleanly, and the package version is accessible."""
from __future__ import annotations

import importlib

import pytest


def test_homerun_imports():
    import homerun
    assert homerun is not None


def test_version_is_accessible():
    import homerun
    assert hasattr(homerun, "__version__")
    assert isinstance(homerun.__version__, str)
    assert homerun.__version__ == "1.0.0"


@pytest.mark.parametrize("module_name", [
    "homerun.pipeline",
    "homerun.config",
    "homerun.utils",
    "homerun.prepare",
    "homerun.trim",
    "homerun.mapping",
    "homerun.tagdirs",
    "homerun.bedgraphs",
    "homerun.tss",
    "homerun.ritrie",
    "homerun.qc",
    "homerun.stability",
    "homerun.report",
])
def test_public_modules_import(module_name):
    mod = importlib.import_module(module_name)
    assert mod is not None


def test_pipeline_exposes_main():
    from homerun.pipeline import main
    assert callable(main)


def test_pipeline_exposes_step_registry():
    from homerun.pipeline import STEP_ORDER, STEP_FUNCS
    assert isinstance(STEP_ORDER, list) and STEP_ORDER
    assert isinstance(STEP_FUNCS, dict) and STEP_FUNCS


def test_config_exposes_load_config():
    from homerun.config import Config, load_config
    assert callable(load_config)
    assert Config is not None
