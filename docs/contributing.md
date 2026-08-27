# Contributing

!!! warning "Alpha"
    HOMERun is in alpha. Interfaces and layout may change, so check with the
    maintainer before building anything large on top of the current API.

Thanks for helping improve HOMERun.

## Development setup

```bash
git clone https://github.com/NguyenThomas986/homerun.git
cd homerun
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -e ".[test]"
```

The editable install puts the `homerun` command on your `PATH` and lets the test suite import the package. (The test suite also adds the repo root to `sys.path` itself, so `pytest` works even before an editable install.)

## Before opening a pull request

- Run the full suite and make sure it's green:

  ```bash
  pytest
  ```

- Add tests for new behavior. The suite runs without any external tools, so new logic should be unit-tested against the `tmp_path` fixtures rather than a live cluster.

- **Keep output directory and file names stable.** Single-command reruns and SLURM array indexing depend on the existing layout (`RawData/`, `Trimmed/`, `Aligned/`, `TagDirs/`, `bedGraphs/`, `TSS/`, `QC/`) and on tag-directory / bedGraph naming (`<sample>_<assay>_<replicate>`, `<sample>_<assay>-combo`). Renaming any of these breaks reruns, so avoid it unless that's the explicit point of the change.
<!-- TODO: confirm the exact set of names to treat as frozen. -->

## Reporting issues

See [Report an Issue](issue.md).

## Code style

<!-- TODO: state the formatter/linter if one is used (e.g. black, ruff) and any conventions to follow, so contributors match the existing style. -->

## Scope notes

- HOMERun currently identifies every sample purely from its FASTQ filename. Manifest input (CSV / TSV / Excel) is out of scope for now.
<!-- TODO: confirm current version scope and roadmap wording. -->
