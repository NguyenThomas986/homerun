"""RIT/RIE (Reads in TSR / Reads in Exon) — a QC metric that estimates csRNA
sample quality (RNA degradation) without relying on RIN scores.

How it works: for each individual csRNA replicate, find its own TSS peaks
("iTSS"), then classify each iTSS peak as falling inside the sample's called
TSR regions, or inside an exon (but NOT inside a TSR — the two are made
mutually exclusive), then pull each iTSS peak's real read strength straight
from that replicate's own tag directory. Summing those reads on each side and
dividing (TSR reads / exon-only reads) gives one RIT/RIE ratio per replicate.
A low ratio suggests a higher proportion of degraded RNA relative to intact
transcription-start signal, without needing a RIN score.

Requires --gtf (CSRNA_GTF): without it, this step is skipped entirely with a
log message rather than failing the pipeline, the same way 'stability' skips
without total RNA.

Per-replicate intermediates and the one-row result land in each replicate's
own Species/Sample/<assay_rep>/RITRIE/. Per-sample results are aggregated
into Species/Sample/QC/ritrie_summary.tsv (+ ritrie.png), so they show up in
that sample's qc_report.html automatically.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt   # noqa: E402
import pandas as pd                     # noqa: E402

from .utils import run, log, done, iter_leaf_dirs, iter_samples, assay_of_leaf, read_homer_table  # noqa: E402

_PEAK_COLS = ["TSS_ID", "chr", "start", "end", "strand"]


# ── HOMER command wrappers ────────────────────────────────────────────────────

def _ensure_gtf_exons(cfg) -> bool:
    """Build the project-wide parsed-GTF-exons file once (cached across every
    sample). Returns False (logging why) if it can't be built."""
    if not cfg.gtf:
        log.info("ritrie: CSRNA_GTF / --gtf not set — skipping ritrie for the whole project.")
        return False
    out = cfg.ritrie_gtf_exons
    if done(out):
        return True
    out.parent.mkdir(parents=True, exist_ok=True)
    job = run(f"parseGTF.pl {cfg.gtf} exons > {out}", label="ritrie: parseGTF exons", check=False)
    if job.returncode != 0 or not done(out):
        log.warning("ritrie: parseGTF.pl failed or produced no output — skipping ritrie.")
        return False
    return True


def _build_itss(cfg, tagdir, out_dir, label: str):
    """findPeaks -style tss on one replicate's own tag directory."""
    out = out_dir / "itss.txt"
    if done(out):
        return out
    out_dir.mkdir(parents=True, exist_ok=True)
    job = run(f"findPeaks {tagdir} -style tss -o {out}", label=f"ritrie: findPeaks {label}", check=False)
    if job.returncode != 0 or not done(out):
        log.warning("ritrie: findPeaks failed for %s — skipping this replicate.", label)
        return None
    return out


def _overlapping_ids(cfg, reference_file, itss_file, out_path, label: str) -> set[str] | None:
    """mergePeaks(reference, itss) -strand, then return the set of iTSS peak
    IDs that landed in a merged region reference also contributed to (i.e.
    iTSS peaks overlapping the reference file). mergePeaks always writes one
    trailing ID column per input file; since itss_file is always passed
    second, its contributed IDs are always the LAST column."""
    if not done(out_path):
        job = run(f"mergePeaks {reference_file} {itss_file} -strand > {out_path}",
                  label=f"ritrie: mergePeaks {label}", check=False)
        if job.returncode != 0 or not done(out_path):
            log.warning("ritrie: mergePeaks failed for %s", label)
            return None
    merged = read_homer_table(out_path)
    if merged.shape[1] < 2:
        log.warning("ritrie: unexpected mergePeaks output for %s (too few columns)", label)
        return None
    ref_col, itss_col = merged.columns[-2], merged.columns[-1]
    overlaps = merged.dropna(subset=[ref_col, itss_col])
    ids = (overlaps[itss_col].astype(str).str.split(",")
           .explode().str.strip())
    return set(ids[ids != ""])


def _annotate_raw(cfg, peak_file, tagdir, out_path, label: str):
    """annotatePeaks.pl <peaks> <genome> -d <tagdir> -raw — pulls each peak's
    real read count straight from that replicate's own tag directory."""
    if not done(out_path):
        job = run(f"annotatePeaks.pl {peak_file} {cfg.genome} -strand + -fragLength 1 "
                  f"-d {tagdir} -raw > {out_path}",
                  label=f"ritrie: annotatePeaks {label}", check=False)
        if job.returncode != 0 or not done(out_path):
            log.warning("ritrie: annotatePeaks.pl failed for %s", label)
            return None
    df = read_homer_table(out_path)
    if df.empty:
        return df
    return df


def _write_peak_subset(itss_df, ids: set[str], out_path) -> int:
    subset = itss_df[itss_df["TSS_ID"].astype(str).isin(ids)]
    subset[_PEAK_COLS].to_csv(out_path, sep="\t", index=False, header=False)
    return len(subset)


# ── Per-replicate RIT/RIE ─────────────────────────────────────────────────────

def _ritrie_for_leaf(cfg, species, sample, leaf_dir, tsr_file) -> dict | None:
    leaf_name = leaf_dir.name
    label = f"{species}/{sample}/{leaf_name}"
    tagdir = cfg.leaf_tagdir(species, sample, leaf_name)
    if not tagdir.is_dir():
        log.info("ritrie: no leaf TagDir yet for %s — run 'tagdirs' first.", label)
        return None

    ritrie_dir = cfg.leaf_ritrie(species, sample, leaf_name)
    out_row = ritrie_dir / "ritrie.tsv"
    if done(out_row):
        return pd.read_csv(out_row, sep="\t").iloc[0].to_dict()

    # Step 1: build this replicate's own TSS peaks, then find which overlap the
    # sample's called TSR regions.
    itss_file = _build_itss(cfg, tagdir, ritrie_dir, label)
    if itss_file is None:
        return None
    itss_df = read_homer_table(itss_file)
    itss_df = itss_df.rename(columns={itss_df.columns[0]: "TSS_ID"})
    if itss_df.empty:
        log.info("ritrie: no TSS peaks found for %s — skipping.", label)
        return None

    tsr_ids = _overlapping_ids(cfg, tsr_file, itss_file, ritrie_dir / "itss_in_TSR.merged.tsv",
                               f"{label} vs TSR")
    if tsr_ids is None:
        return None

    # Step 2: find which iTSS peaks overlap an exon, then exclude anything
    # that's already counted on the TSR side, so the two sets are disjoint.
    exon_ids = _overlapping_ids(cfg, cfg.ritrie_gtf_exons, itss_file, ritrie_dir / "itss_in_exon.merged.tsv",
                                f"{label} vs exons")
    if exon_ids is None:
        return None
    exon_only_ids = exon_ids - tsr_ids

    tsr_peaks = ritrie_dir / "itss_in_TSR.peaks.tsv"
    exon_peaks = ritrie_dir / "itss_in_exon_only.peaks.tsv"
    _write_peak_subset(itss_df, tsr_ids, tsr_peaks)
    _write_peak_subset(itss_df, exon_only_ids, exon_peaks)

    # Step 3: pull each peak's real read strength from this replicate's own
    # tag directory, and sum each side.
    tsr_anno = _annotate_raw(cfg, tsr_peaks, tagdir, ritrie_dir / "itss_in_TSR.anno.tsv", f"{label} TSR")
    exon_anno = _annotate_raw(cfg, exon_peaks, tagdir, ritrie_dir / "itss_in_exon_only.anno.tsv", f"{label} exon")
    if tsr_anno is None or exon_anno is None:
        return None

    tss_in_tsr, tss_in_exon = len(tsr_anno), len(exon_anno)
    reads_in_tsr = float(tsr_anno.iloc[:, -1].sum()) if not tsr_anno.empty else 0.0
    reads_in_exon = float(exon_anno.iloc[:, -1].sum()) if not exon_anno.empty else 0.0

    if reads_in_exon <= 0:
        log.warning("ritrie: %s has zero exon-side reads — RIT/RIE undefined.", label)
        rit_rie = None
    else:
        rit_rie = reads_in_tsr / reads_in_exon
    tit_tie = (tss_in_tsr / tss_in_exon) if tss_in_exon > 0 else None

    row = {
        "Species": species, "Sample": sample, "Library": leaf_name,
        "TSS_in_TSR": tss_in_tsr, "TSS_in_Exon": tss_in_exon,
        "Reads_in_TSR": reads_in_tsr, "Reads_in_Exon": reads_in_exon,
        "RIT_RIE": rit_rie, "TIT_TIE": tit_tie,
    }
    pd.DataFrame([row]).to_csv(out_row, sep="\t", index=False)
    log.info("ritrie: %s → RIT/RIE = %s", label, f"{rit_rie:.3f}" if rit_rie is not None else "NA")
    return row


# ── Per-sample aggregation ────────────────────────────────────────────────────

def _plot_ritrie(cfg, species, sample, df, qc_dir) -> None:
    plotted = df.dropna(subset=["RIT_RIE"])
    if plotted.empty:
        return
    plt.figure(figsize=(max(6, 1.2 * len(plotted)), 5))
    plt.bar(plotted["Library"], plotted["RIT_RIE"], color="steelblue", edgecolor=".4")
    plt.axhline(1, color="gray", ls="--", lw=1, label="RIT/RIE = 1")
    plt.ylabel("RIT / RIE"); plt.xlabel("")
    plt.title(f"{species}/{sample} — RIT/RIE per csRNA replicate")
    plt.xticks(rotation=30, ha="right"); plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(qc_dir / "ritrie.png", dpi=150, bbox_inches="tight")
    plt.close()
    log.info("ritrie: ritrie.png (%s/%s)", species, sample)


def _run_ritrie_sample(cfg, species, sample) -> None:
    tsr_file = cfg.sample_tss(species, sample) / f"{sample}.tss.txt"
    if not tsr_file.exists():
        log.info("ritrie: no %s.tss.txt for %s/%s yet — run 'tss' first.", sample, species, sample)
        return

    leaf_runs = [ld for sp, sa, ld in iter_leaf_dirs(cfg)
                 if sp == species and sa == sample and assay_of_leaf(ld.name) == "csRNA"]
    if not leaf_runs:
        log.info("ritrie: no csRNA leaf runs for %s/%s", species, sample)
        return

    rows = []
    for leaf_dir in leaf_runs:
        row = _ritrie_for_leaf(cfg, species, sample, leaf_dir, tsr_file)
        if row is not None:
            rows.append(row)

    if not rows:
        log.info("ritrie: no results for %s/%s", species, sample)
        return

    df = pd.DataFrame(rows)
    qc_dir = cfg.sample_qc(species, sample)
    qc_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(qc_dir / "ritrie_summary.tsv", sep="\t", index=False)

    valid = df["RIT_RIE"].dropna()
    if len(valid):
        log.info("ritrie: %s/%s sample RIT/RIE (mean across %d replicate(s)) = %.3f",
                 species, sample, len(valid), valid.mean())
    _plot_ritrie(cfg, species, sample, df, qc_dir)


def run_ritrie(cfg) -> None:
    if not _ensure_gtf_exons(cfg):
        return
    samples = list(iter_samples(cfg))
    if not samples:
        log.info("ritrie: no Species/Sample dirs found under %s", cfg.project)
        return
    for species, sample in samples:
        _run_ritrie_sample(cfg, species, sample)
