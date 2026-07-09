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

OPTIONAL — enables extra features when set:
  CSRNA_GTF            path to a GTF annotation file. Only needed for the
                       'ritrie' step (RIT/RIE QC metric); without it, ritrie
                       is skipped with a log message rather than failing.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


def _env(name: str, default):
    v = os.getenv(name)
    return v if v not in (None, "") else default


def _pick(args, attr, env_name, default):
    """Precedence: explicit CLI flag (if given) > CSRNA_* env var > built-in default.

    A flag left unset by argparse is None, so it falls through to the env var,
    keeping config.env behavior identical when no flag is passed.
    """
    cli = getattr(args, attr, None) if args is not None else None
    if cli is not None:
        return cli
    return _env(env_name, default)


@dataclass
class Config:
    project: Path

    # ── Aligner (REQUIRED to be set explicitly; no genome is assumed) ─────────
    aligner: str = "star"             # "star" or "hisat2"
    genome_index: str = ""            # STAR --genomeDir OR hisat2 -x prefix (REQUIRED)
    genome: str = ""                  # HOMER -genome for tagdirs/tss (REQUIRED)
    gtf: str = ""                     # GTF annotation, only needed for the ritrie step (OPTIONAL)

    threads: int = 20

    # Alignment — STAR (csRNA-tuned defaults; override via flag/env if needed)
    star_filter_multimap: str = "10000"   # --outFilterMultimapNmax
    star_multimap_out: str = "1"          # --outSAMmultNmax
    star_multimap_order: str = "Random"   # --outMultimapperOrder
    # Alignment — HISAT2
    hisat2_strandness: str = "F"          # --rna-strandness

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
    # NOTE: RawData/Trimmed/Aligned/TagDir/bedGraph/QC/TSS are all nested under
    # Species/Sample/ (see below) — there is no flat project-level equivalent.
    @property
    def logs_dir(self) -> Path:  return self.project / "logs"
    @property
    def starindex(self) -> Path: return Path(self.genome_index)

    # ── Species/Sample nested-layout helpers (always used, not opt-in) ────────
    def sample_dir(self, species: str, sample: str) -> Path:
        """Species/Sample/ — where per-sample QC and TSS live."""
        return self.project / species / sample

    def run_dir(self, species: str, sample: str, leaf_name: str) -> Path:
        """Species/Sample/<assay_rep>/ — where per-run RawData/Trimmed/Aligned/
        TagDir/bedGraph live."""
        return self.sample_dir(species, sample) / leaf_name

    # ── Nested TagDir/bedGraph/QC/TSS paths (replace the old flat dirs) ───────
    def leaf_tagdir(self, species: str, sample: str, leaf_name: str) -> Path:
        """Species/Sample/<assay_rep>/TagDir — one per individual replicate."""
        return self.run_dir(species, sample, leaf_name) / "TagDir"

    def leaf_bedgraph(self, species: str, sample: str, leaf_name: str) -> Path:
        """Species/Sample/<assay_rep>/bedGraph — one per individual replicate."""
        return self.run_dir(species, sample, leaf_name) / "bedGraph"

    def combo_dir(self, species: str, sample: str, assay: str) -> Path:
        """Species/Sample/<assay>-combo/ — merged-replicate run for one assay."""
        return self.sample_dir(species, sample) / f"{assay}-combo"

    def combo_tagdir(self, species: str, sample: str, assay: str) -> Path:
        return self.combo_dir(species, sample, assay) / "TagDir"

    def combo_bedgraph(self, species: str, sample: str, assay: str) -> Path:
        return self.combo_dir(species, sample, assay) / "bedGraph"

    def sample_qc(self, species: str, sample: str) -> Path:
        """Species/Sample/QC — one QC dir per sample, covering all assays."""
        return self.sample_dir(species, sample) / "QC"

    def sample_tss(self, species: str, sample: str) -> Path:
        """Species/Sample/TSS — one TSS dir per sample."""
        return self.sample_dir(species, sample) / "TSS"

    # ── RIT/RIE (Reads in TSR / Reads in Exon) QC metric ──────────────────────
    def species_ritrie_gtf_exons(self, species: str) -> Path:
        """Species/RITRIE/parsed_gtf_exons.tsv — parsed once per species (built
        from that species' --gtf) and shared across every sample of that
        species, rather than mixed across species at the flat project root."""
        return self.project / species / "RITRIE" / "parsed_gtf_exons.tsv"

    def leaf_ritrie(self, species: str, sample: str, leaf_name: str) -> Path:
        """Species/Sample/<assay_rep>/RITRIE — working dir for one csRNA
        replicate's RIT/RIE intermediates (iTSS peaks, merges, annotations)."""
        return self.run_dir(species, sample, leaf_name) / "RITRIE"


def load_config(args=None) -> Config:
    # --project flag (if given) overrides CSRNA_PROJECT; else env; else CWD.
    cli_project = getattr(args, "project", None) if args is not None else None
    if cli_project:
        project_path = Path(cli_project).resolve()
    else:
        project_path = Path(_env("CSRNA_PROJECT", os.getcwd())).resolve()
    return Config(
        project=project_path,
        # ── Core (flag > env > default) ──────────────────────────────────────
        aligner=_pick(args, "aligner", "CSRNA_ALIGNER", "star").lower(),
        genome_index=_pick(args, "genome_index", "CSRNA_GENOME_INDEX", ""),
        genome=_pick(args, "genome", "CSRNA_GENOME", ""),
        gtf=_pick(args, "gtf", "CSRNA_GTF", ""),
        copy_src=_pick(args, "copy_src", "CSRNA_COPY_SRC", ""),
        # ── Alignment (flag > env > default; csRNA-tuned defaults) ───────────
        star_filter_multimap=_pick(args, "star_filter_multimap",
                                   "CSRNA_STAR_FILTER_MULTIMAP", "10000"),
        star_multimap_out=_pick(args, "star_multimap_out",
                                "CSRNA_STAR_MULTIMAP_OUT", "1"),
        star_multimap_order=_pick(args, "star_multimap_order",
                                  "CSRNA_STAR_MULTIMAP_ORDER", "Random"),
        hisat2_strandness=_pick(args, "hisat2_strandness",
                                "CSRNA_HISAT2_STRANDNESS", "F"),
        # ── Trimming / TSS (flag > env > default) ────────────────────────────
        threads=int(_pick(args, "threads", "CSRNA_THREADS",
                          os.getenv("SLURM_CPUS_PER_TASK", "20"))),
        trim_adapter=_pick(args, "trim_adapter", "CSRNA_TRIM_ADAPTER",
                          "AGATCGGAAGAGCACACGTCT"),
        trim_min=_pick(args, "trim_min", "CSRNA_TRIM_MINLEN", "20"),
        trim_max=_pick(args, "trim_max", "CSRNA_TRIM_MAXLEN", "58"),
        ntag_threshold=_pick(args, "ntag_threshold", "CSRNA_NTAG_THRESHOLD", "7"),
        skip_chr=_pick(args, "skip_chr", "CSRNA_SKIP_CHR", "chrEBV"),
        # ── Env-only (no flags) ──────────────────────────────────────────────
        stability_col=_env("CSRNA_STABILITY_COL", "Stable/Unstable"),
        rna_col=_env("CSRNA_RNA_COL", ""),
        rna_stable_threshold=float(_env("CSRNA_RNA_STABLE_THRESHOLD", "0") or 0),
        distal_col=_env("CSRNA_DISTAL_COL", "Promoter Proximal/Distal"),
        log_path=_pick(args, "log_path", "CSRNA_LOG", ""),
        starindex_url=_env("CSRNA_STARINDEX_URL", ""),
    )