"""Central configuration for the csRNA-seq pipeline.

All settings come from defaults here, overridable by CSRNA_* environment
variables (so the SLURM script can drive everything without editing code).

REQUIRED, no defaults (the pipeline refuses to guess a genome):
  CSRNA_ALIGNER        "star" (default) or "hisat2"
  CSRNA_GENOME_INDEX   path to the aligner index for the CHOSEN tool:
                         - STAR   → the --genomeDir DIRECTORY
                         - HISAT2 → the -x index PREFIX (e.g. /path/idx/genome)
  CSRNA_GENOME         HOMER -genome for tag dirs / TSS (genome name like 'hg38'
                       if installed in HOMER, or a path to a genome FASTA)
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


def _env(name: str, default):
    v = os.getenv(name)
    return v if v not in (None, "") else default


@dataclass
class Config:
    project: Path

    # ── Aligner (REQUIRED to be set explicitly; no genome is assumed) ─────────
    aligner: str = "star"             # "star" or "hisat2"
    genome_index: str = ""            # STAR --genomeDir OR hisat2 -x prefix (REQUIRED)
    genome: str = ""                  # HOMER -genome for tagdirs/tss (REQUIRED)

    threads: int = 20

    # Trimming (homerTools, single-end csRNA/sRNA)
    trim_adapter: str = "AGATCGGAAGAGCACACGTCT"
    trim_mis: str = "2"
    trim_minmatch: str = "4"
    trim_min: str = "20"
    trim_max: str = "58"

    # TSS calling / browser tracks
    ntag_threshold: str = "7"
    skip_chr: str = "chrEBV"          # "" to skip nothing

    # Data sources
    copy_src: str = ""                # optional cp glob into RawData/

    # Stability/location detection (findcsRNATSS column names; override via env if needed)
    stability_col: str = "Stable/Unstable"
    rna_col: str = ""
    rna_stable_threshold: float = 0.0
    distal_col: str = "Promoter Proximal/Distal"

    log_path: str = ""
    sample_filter: str = ""           # restrict trim/mapping to one sample (job arrays)

    # STARIndex auto-download (only used when aligner == "star" and genome_index is missing)
    starindex_url: str = ""           # set via CSRNA_STARINDEX_URL in config.env

    # ── Derived directories ───────────────────────────────────────────────────
    @property
    def rawdata(self) -> Path:   return self.project / "RawData"
    @property
    def trimmed(self) -> Path:   return self.project / "Trimmed"
    @property
    def aligned(self) -> Path:   return self.project / "Aligned"
    @property
    def tagdirs(self) -> Path:   return self.project / "TagDirs"
    @property
    def bedgraphs(self) -> Path: return self.project / "bedGraphs"
    @property
    def tss(self) -> Path:       return self.project / "TSS"
    @property
    def qc(self) -> Path:        return self.project / "QC"
    @property
    def logs_dir(self) -> Path:  return self.project / "logs"
    @property
    def reports(self) -> Path:   return self.project / "Reports"
    @property
    def starindex(self) -> Path: return self.project / "STARIndex"

    def output_dirs(self):
        return [self.rawdata, self.trimmed, self.aligned, self.tagdirs,
                self.bedgraphs, self.tss, self.qc, self.reports]


def load_config(project: str | None = None) -> Config:
    project_path = Path(_env("CSRNA_PROJECT", project or os.getcwd())).resolve()
    return Config(
        project=project_path,
        aligner=_env("CSRNA_ALIGNER", "star").lower(),
        genome_index=_env("CSRNA_GENOME_INDEX", ""),
        genome=_env("CSRNA_GENOME", ""),
        threads=int(_env("CSRNA_THREADS", os.getenv("SLURM_CPUS_PER_TASK", "20"))),
        trim_adapter=_env("CSRNA_TRIM_ADAPTER", "AGATCGGAAGAGCACACGTCT"),
        trim_min=_env("CSRNA_TRIM_MINLEN", "20"),
        trim_max=_env("CSRNA_TRIM_MAXLEN", "58"),
        ntag_threshold=_env("CSRNA_NTAG_THRESHOLD", "7"),
        skip_chr=_env("CSRNA_SKIP_CHR", "chrEBV"),
        copy_src=_env("CSRNA_COPY_SRC", ""),
        stability_col=_env("CSRNA_STABILITY_COL", "Stable/Unstable"),
        rna_col=_env("CSRNA_RNA_COL", ""),
        rna_stable_threshold=float(_env("CSRNA_RNA_STABLE_THRESHOLD", "0") or 0),
        distal_col=_env("CSRNA_DISTAL_COL", "Promoter Proximal/Distal"),
        log_path=_env("CSRNA_LOG", ""),
        starindex_url=_env("CSRNA_STARINDEX_URL", ""),
    )
