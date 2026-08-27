# Testing

!!! warning "Alpha"
    HOMERun is in alpha testing. Test layout and coverage are still evolving.

The test suite is **pure Python**: it needs no HOMER, STAR, HISAT2, or skewer, no network, and no SLURM. Every filesystem fixture uses a temporary directory, so the tests run anywhere from a fresh checkout. Fixtures that need external-tool output (FASTQ, SAM, HOMER tag directories) build minimal syntactically-valid fakes by hand instead of shelling out to real tools.

## Install test dependencies

```bash
pip install -e ".[test]"
```
<!-- TODO: confirm the extras name is `test` (it matches [project.optional-dependencies] in pyproject.toml). -->

## Run

```bash
pytest
```

Options are preconfigured in `pyproject.toml` (`testpaths = ["tests"]`, `-ra --strict-markers`), so a bare `pytest` from the repo root runs everything.

Run a subset:

```bash
pytest tests/test_cli.py        # one file
pytest -k precedence            # tests matching a keyword
```

Slow tests are marked; select or skip them:

```bash
pytest -m "not slow"            # skip slow tests
pytest -m slow                  # run only slow tests
```

Because `--strict-markers` is on, any marker not registered in `pyproject.toml` is an error — register new markers there before using them.

## What's covered

- **CLI parsing and `main()` fast paths** — argparse construction, `--help`, exit codes, and the `--count-samples` / `--count-groups` / `--check-rerun` / `--stage-raw` controller paths that a SLURM wrapper depends on, none of which run the pipeline.
- **Configuration** — flag > env > default precedence, the `_pick()` helper, `SLURM_CPUS_PER_TASK` fallback, and path-builder correctness.
- **Filename identity** — `parse_sample_name`, `assay_of_leaf` / `seq_type`, replicate markers, and identity-based R2 mate-pairing (including the ENCODE mismatched-accession regression).
- **Package imports** — every public module imports cleanly and `__version__` is accessible.
<!-- TODO: confirm test file names (e.g. the filename-identity tests and the import tests). -->

## Fixtures (`tests/conftest.py`)

- `project_dir` — an empty temporary project root.
- `make_cfg` / `cfg` — a real `Config` pointed at a temporary project, any field overridable.
- `make_fastq` — writes a minimal valid FASTQ.
- `staged_sample` — a `Config` plus one already-staged csRNA replicate.
- `paired_totalrna_sample` — a staged totalRNA replicate with both mates present.
- `non_interactive_matplotlib` — forces matplotlib's headless `Agg` backend so plotting tests never open a GUI window.

## Writing new tests

- Use `tmp_path` and the provided fixtures; never touch absolute paths or the network.
- Build fake tool outputs by hand (see `make_fastq`) rather than calling real tools.
- Register any new marker in `pyproject.toml`, or `--strict-markers` will fail the run.
