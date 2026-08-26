# Installation

HOMERun can be installed three ways. **Conda is recommended**, because it also installs the external bioinformatics tools HOMERun depends on (HOMER, STAR, HISAT2, samtools, skewer). The other two methods install only the Python package — you provide the tools yourself.

## Recommended: Conda (Bioconda)

```bash
conda install -c bioconda homerun
```
<!-- TODO: Bioconda package not published yet — this command works once the recipe is merged. Confirm final package name if "homerun" is already taken on Bioconda/conda-forge. -->

This pulls in every required tool automatically. Verify the install:

```bash
homerun --help
```

## PyPI (pip)

```bash
pip install homerun
```
<!-- TODO: confirm PyPI distribution name; "homerun" may be taken. Import name can stay `homerun` even if the published name differs. -->

pip installs the Python package only. The external tools must already be on your `PATH`. The easiest way is a dedicated conda environment:

```bash
conda create -n homerun-tools -c bioconda homer star hisat2 samtools skewer
```
<!-- TODO: confirm exact tool list and any minimum versions -->

## From source (git clone)

```bash
git clone https://github.com/NguyenThomas986/homerun.git
cd homerun
pip install -e .
```
<!-- TODO: on Kamiak the pipeline is run via PYTHONPATH + `python -m` rather than `pip install -e .`. Confirm which method to document here, and the exact module name (homerun vs csrnaseq). -->

## Requirements

- Python >= 3.9
- A SLURM-based HPC cluster
- External tools: HOMER, STAR, HISAT2, samtools, skewer
<!-- TODO: confirm minimum tool versions -->
