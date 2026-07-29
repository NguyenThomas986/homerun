"""QC plots from -combo tag directories, one set per Species/Sample.

Each sample gets its own QC/ dir (Species/Sample/QC/) combining csRNA, sRNA,
and totalRNA combo TagDirs and that sample's TSS output — instead of one
flat project-level QC/ mixing every sample together.
"""
from __future__ import annotations

import math
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402
import numpy as np                     # noqa: E402
import pandas as pd                    # noqa: E402
import seaborn as sns                  # noqa: E402

from .utils import log, iter_samples, iter_leaf_dirs, assay_of_leaf  # noqa: E402
from .stability import _read_homer_tss, _location, DISTAL_C, PROX_C  # noqa: E402

_ASSAYS = ("csRNA", "sRNA", "totalRNA")
_ASSAY_COLORS = {"csRNA": "#2c7fb8", "sRNA": "#de7c00", "totalRNA": "#636363"}


def _replicate_grid(n: int, ncols: int = 4):
    """Fixed-grid figure/axes for one small subplot per replicate (used by
    the per-replicate read-length/nucleotide-frequency/autocorrelation
    plots), instead of a single heatmap — this keeps each replicate's full
    curve visible and scales to any number of replicates by adding more
    rows, wrapping at `ncols` columns per row. Unused axes (when n isn't a
    multiple of ncols) are hidden rather than left blank/empty-looking."""
    ncols = min(ncols, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 2.8 * nrows), squeeze=False)
    flat = axes.flat
    for ax in list(flat)[n:]:
        ax.axis("off")
    return fig, axes


# ── Helpers ───────────────────────────────────────────────────────────────────

def _combo(cfg, species, sample, assay):
    """Combo TagDir for one assay of one sample, or None if it doesn't exist."""
    d = cfg.combo_tagdir(species, sample, assay)
    return d if d.is_dir() else None


def _all_combos(cfg, species, sample):
    return [d for a in _ASSAYS if (d := _combo(cfg, species, sample, a))]


def _label(d) -> str:
    """Readable name for a TagDir, e.g. 'csRNA-combo' or 'csRNA_r1'. The tag
    directory itself is named that (Species/Sample/TagDirs/<name>), not
    nested one level further under a generic 'TagDir' folder."""
    return d.name


def _leaf_tagdirs_for_sample(cfg, species, sample):
    """(leaf_name, tagdir_path) for every individual (non-combo) replicate
    TagDir under this sample, across all assays — sorted csRNA replicates
    first, then sRNA, then totalRNA, so a multi-assay grid/table groups
    naturally instead of interleaving by whatever order the filesystem
    happens to return."""
    seen = set()
    out = []
    for sp, sa, leaf_name, _r1 in iter_leaf_dirs(cfg):
        if sp != species or sa != sample or leaf_name in seen:
            continue
        seen.add(leaf_name)
        td = cfg.leaf_tagdir(species, sample, leaf_name)
        if td.is_dir():
            out.append((leaf_name, td))

    order = {"csRNA": 0, "sRNA": 1, "totalRNA": 2}

    def _sort_key(item):
        leaf_name, _td = item
        return (order.get(assay_of_leaf(leaf_name), 3), leaf_name)

    return sorted(out, key=_sort_key)



# ── Existing plots ────────────────────────────────────────────────────────────

def qc_read_length(cfg, species, sample, qc_dir) -> None:
    combos = [d for d in _all_combos(cfg, species, sample)
              if (d / "tagLengthDistribution.txt").exists()]
    if not combos:
        log.info("QC read-length: no tagLengthDistribution.txt yet."); return
    plt.figure(figsize=(12, 6))
    for d in combos:
        df = pd.read_csv(d / "tagLengthDistribution.txt", sep="\t", index_col=0)
        sns.lineplot(x=df.index, y=df.iloc[:, 0], label=_label(d))
    plt.xlabel("Read length (nt)"); plt.ylabel("Fraction of reads")
    plt.title("Read length distribution"); plt.tight_layout()
    plt.savefig(qc_dir / "read_length_distribution.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: read_length_distribution.png")


def qc_nucleotide_freq(cfg, species, sample, qc_dir) -> None:
    have = {a: d for a in _ASSAYS
            if (d := _combo(cfg, species, sample, a)) and (d / "tagFreqUniq.txt").exists()}
    if not have:
        log.info("QC nt-freq: no tagFreqUniq.txt yet."); return
    fig, axes = plt.subplots(len(have), 1, figsize=(12, 4 * len(have)),
                             sharex=True, squeeze=False)
    for ax, (label, d) in zip(axes[:, 0], have.items()):
        df = pd.read_csv(d / "tagFreqUniq.txt", index_col=0, sep="\t")
        sns.lineplot(df.iloc[:, 0:4], ax=ax)
        ax.set_title(label); ax.set_ylabel("Nucleotide frequency")
    plt.xlabel("Distance from 5' end of read"); plt.tight_layout()
    plt.savefig(qc_dir / "nucleotide_frequency.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: nucleotide_frequency.png")


def qc_autocorrelation(cfg, species, sample, qc_dir) -> None:
    have = {a: d for a in _ASSAYS
            if (d := _combo(cfg, species, sample, a)) and (d / "tagAutocorrelation.txt").exists()}
    if not have:
        log.info("QC autocorrelation: no tagAutocorrelation.txt yet."); return
    fig, axes = plt.subplots(len(have), 1, figsize=(12, 4 * len(have)),
                             sharex=True, squeeze=False)
    for ax, (label, d) in zip(axes[:, 0], have.items()):
        df = pd.read_csv(d / "tagAutocorrelation.txt", sep="\t", index_col=0, header=0,
                         names=["Same Strand", "Opposite Strand"]).iloc[1400:2601]
        sns.lineplot(data=df, ax=ax)
        ax.set_title(label); ax.set_ylabel("Read counts rel. to 5' end")
    plt.xlabel("Distance from 5' end of read"); plt.tight_layout()
    plt.savefig(qc_dir / "autocorrelation.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: autocorrelation.png")


# ── New plots (ported from old pipeline) ─────────────────────────────────────

def _median_tags_bar_merged(cfg, species, sample, qc_dir) -> None:
    """Single side-by-side bar chart of median tags per position for csRNA, sRNA, totalRNA."""
    panels = []
    for assay in _ASSAYS:
        d = _combo(cfg, species, sample, assay)
        if not d or not (d / "tagCountDistribution.txt").exists():
            continue
        df = pd.read_csv(d / "tagCountDistribution.txt", sep="\t")
        median_val = str(list(df.columns)).split("=")[1].split(",")[0].strip()
        row = pd.DataFrame([{"Library": _label(d), "Median": float(median_val), "Type": assay}])
        panels.append((assay, row))

    if not panels:
        log.info("QC median-tags: no tagCountDistribution.txt found."); return

    n = len(panels)
    fig, axes = plt.subplots(n, 1, figsize=(8, max(3, 4 * n)), squeeze=False)
    for ax, (label, mdf) in zip(axes[:, 0], panels):
        sns.barplot(data=mdf, y="Library", x="Median",
                    capsize=0.4, linewidth=1.5, edgecolor=".4", color="steelblue", ax=ax)
        ax.axvline(1,   color="green",  ls="--", lw=1.2, label="ideal (1)")
        ax.axvline(1.2, color="red",    ls="-",  lw=1.2, label="warn (1.2)")
        ax.set_xlabel("Median tags per position (should be =1)")
        ax.set_title(f"{label}", fontsize=11, fontweight="bold")
        ax.legend(fontsize=8)

    fig.suptitle("Median Tags Per Position", fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(qc_dir / "median_tags_per_position.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: median_tags_per_position.png")


def qc_threshold_optimization(cfg, species, sample, qc_dir) -> None:
    """Threshold optimization plot from prefix.inputDistribution.txt files."""
    files = sorted(cfg.sample_tss(species, sample).glob("*.inputDistribution.txt"))
    if not files:
        log.info("QC threshold: no *.inputDistribution.txt for %s/%s", species, sample); return

    n = len(files)
    fig, axes = plt.subplots(n, 1, figsize=(8, 5 * n), squeeze=False)

    for ax, f in zip(axes[:, 0], files):
        df = pd.read_csv(f, sep="\t", header=0)
        df.columns = [c.strip() for c in df.columns]

        x_col, tss_col, exon_col, diff_col = (
            "csRNA/input log2 ratio", "TSS CDF", "Exon CDF", "Difference")

        ax.plot(df[x_col], df[tss_col],  color="steelblue",  lw=2, label="TSS CDF")
        ax.plot(df[x_col], df[exon_col], color="darkorange", lw=2, label="Exon CDF")
        ax.plot(df[x_col], df[diff_col], color="gray", lw=1.5, ls="--", label="Difference")

        idx = df[diff_col].idxmax()
        thresh = df.loc[idx, x_col]
        ax.axvline(thresh, color="black", ls=":", lw=1)
        ax.text(thresh + 0.1, 0.05, f"threshold = {thresh:.2f}", fontsize=8)

        s = f.name.replace(".inputDistribution.txt", "")
        ax.set_title(f"Threshold Optimization — {s}")
        ax.set_xlabel("Log2 Ratio of csRNAseq/control")
        ax.set_ylabel("Cumulative Distribution")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(qc_dir / "threshold_optimization.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: threshold_optimization.png")


def qc_tss_nucleotide_freq(cfg, species, sample, qc_dir) -> None:
    """Nucleotide frequency at primary TSS from *.freq.tsv files."""
    files = sorted(cfg.sample_tss(species, sample).glob("*.freq.tsv"))
    if not files:
        log.info("QC TSS nt-freq: no *.freq.tsv for %s/%s", species, sample); return

    n = len(files)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 5), squeeze=False)

    for ax, f in zip(axes[0], files):
        df = pd.read_csv(f, sep="\t", index_col=0)
        nt_cols = {"A frequency": ("A", "steelblue"),
                   "C frequency": ("C", "darkorange"),
                   "G frequency": ("G", "gray"),
                   "T frequency": ("T", "gold")}
        for col, (label, color) in nt_cols.items():
            if col in df.columns:
                ax.plot(df.index, df[col], label=label, color=color, lw=1.5)

        s = f.name.split(".tss.txt")[0]
        ax.set_title(s)
        ax.set_xlabel("Distance from TSS")
        ax.set_ylabel("Nucleotide Frequency")
        ax.set_xlim(-100, 100)
        ax.legend(fontsize=9)

    plt.suptitle("Nucleotide Frequencies at Primary TSS", y=1.0, fontsize=12)
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.savefig(qc_dir / "tss_nucleotide_frequency.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: tss_nucleotide_frequency.png")


def qc_tsr_summary(cfg, species, sample, qc_dir) -> None:
    """Parse findcsRNATSS stats files into a summary table PNG."""
    files = sorted(cfg.sample_tss(species, sample).glob("*.stats.txt"))
    if not files:
        log.info("QC TSR summary: no *.stats.txt for %s/%s", species, sample); return

    has_rna = cfg.combo_tagdir(species, sample, "totalRNA").is_dir()

    rows = []
    for f in files:
        txt = f.read_text()
        s = f.name.replace(".stats.txt", "")

        def _grab(pattern, default="NA"):
            import re
            m = re.search(pattern, txt)
            return m.group(1).strip() if m else default

        row = {
            "Sample":               s,
            "Total csRNA reads":    _grab(r"Total csRNA reads:\s+([\d.]+)"),
            "Total input reads":    _grab(r"Total input reads:\s+([\d.]+)"),
            "Putative TSS":         _grab(r"total putative TSS clusters\s+(\d+)"),
            "Valid TSS":            _grab(r"Valid TSS clusters\s+(\d+)"),
            "% Distal":             _grab(r"Fraction Promoter-Distal.*?:\s+([\d.]+%)"),
            "% Bidirectional":      _grab(r"Fraction of bidirectional.*?:\s+([\d.]+%)"),
            "Log2 vs Input":        _grab(r"log2 fold vs\. input:\s+([\d.\-]+)"),
        }
        if has_rna:
            row.update({
                "Total RNA reads":  _grab(r"Total rna reads:\s+([\d.]+)"),
                "% Stable":         _grab(r"Fraction of stable.*?:\s+([\d.]+%)"),
                "SS":               _grab(r"SS:\s+\d+\s+\(([\d.]+%)"),
                "SU":               _grab(r"SU:\s+\d+\s+\(([\d.]+%)"),
                "S":                _grab(r"\tS:\s+\d+\s+\(([\d.]+%)"),
                "US":               _grab(r"US:\s+\d+\s+\(([\d.]+%)"),
                "UU":               _grab(r"UU:\s+\d+\s+\(([\d.]+%)"),
                "U":                _grab(r"\tU:\s+\d+\s+\(([\d.]+%)"),
                "Log2 vs RNA":      _grab(r"log2 fold vs\. rna:\s+([\d.\-]+)"),
            })
        rows.append(row)

    df = pd.DataFrame(rows).set_index("Sample").T
    df = df[~(df == "NA").all(axis=1)]

    fig, ax = plt.subplots(figsize=(max(6, 3 * len(rows)), len(df) * 0.5 + 1))
    ax.axis("off")
    tbl = ax.table(cellText=df.values, rowLabels=df.index,
                   colLabels=df.columns, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9)
    tbl.auto_set_column_width(col=list(range(len(df.columns) + 1)))
    plt.title("TSR Summary", fontsize=11, pad=10)
    plt.tight_layout()
    plt.savefig(qc_dir / "tsr_summary.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: tsr_summary.png")


def qc_tsr_annotation(cfg, species, sample, qc_dir) -> None:
    """Bar plot of TSR annotation categories from *.tss.txt files."""
    files = sorted(cfg.sample_tss(species, sample).glob("*.tss.txt"))
    if not files:
        log.info("QC TSR annotation: no *.tss.txt for %s/%s", species, sample); return

    cats = ["tss", "firstExon", "singleExon", "tssAntisense",
            "otherExon", "otherExonBidirectional", "other"]
    rows = []
    for f in files:
        df = pd.read_csv(f, sep="\t", low_memory=False)
        df = df[df["chr"].str.startswith("chr", na=False)]
        counts = df["annotation"].value_counts()
        row = {"Sample": f.name.replace(".tss.txt", "")}
        for c in cats:
            row[c] = counts.get(c, 0)
        rows.append(row)

    mdf = pd.DataFrame(rows).set_index("Sample")
    mdf = mdf.loc[:, (mdf > 0).any(axis=0)]
    if mdf.empty:
        log.info("QC TSR annotation: no annotation counts found — skipping"); return
    ax = mdf.plot(kind="bar", stacked=True, figsize=(max(6, 2 * len(rows)), 5),
                  colormap="tab10")
    ax.set_ylabel("Number of TSR clusters")
    ax.set_title("TSR Annotation Categories")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(qc_dir / "tsr_annotation.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: tsr_annotation.png")


def _tagdir_stats_rows(tagdirs_with_labels) -> list[dict]:
    """One row of tagInfo.txt stats per (label, tagdir_path) pair — shared by
    the combo-level table (qc_tagdir_stats) and the per-replicate table
    (qc_tagdir_stats_per_replicate) so both read the exact same fields the
    exact same way."""
    rows = []
    for label, d in tagdirs_with_labels:
        info = d / "tagInfo.txt"
        if not info.exists():
            continue
        txt = info.read_text()

        def _val(key, txt=txt):
            for line in txt.splitlines():
                if line.startswith(key):
                    return line.split("=")[-1].strip()
            return "NA"

        genome_line = next((l for l in txt.splitlines() if l.startswith("genome=")), "")
        parts = genome_line.split("\t")
        rows.append({
            "TagDir":                    label,
            "Total Tags":                parts[2].strip() if len(parts) > 2 else "NA",
            "Unique Positions":          parts[1].strip() if len(parts) > 1 else "NA",
            "Tags per BP":               _val("tagsPerBP"),
            "Avg Tags/Position":         _val("averageTagsPerPosition"),
            "Median Tags/Position":      _val("medianTagsPerPosition"),
            "Avg Read Length":           _val("averageTagLength"),
            "Avg Fragment GC":           _val("averageFragmentGCcontent"),
        })
    return rows


def _save_stats_table(rows: list[dict], title: str, out_path, transpose: bool = True) -> bool:
    """Render tagdir-stats rows as a table PNG. Returns False (and writes
    nothing) if there are no rows, so callers can log why and skip.

    transpose=True (metrics as rows, one column per TagDir) reads well for a
    handful of TagDirs (the combo case: at most csRNA/sRNA/totalRNA). With
    many TagDirs (the per-replicate case, potentially 50+) that layout grows
    the image absurdly WIDE instead of tall, so transpose=False keeps
    TagDirs as rows instead — the image grows vertically, which scrolls
    normally in an HTML report instead of requiring horizontal scrolling.
    """
    if not rows:
        return False
    df = pd.DataFrame(rows).set_index("TagDir")
    n = len(df)
    if transpose:
        df = df.T
        figsize = (max(6, 3 * n), len(df) * 0.5 + 1)
    else:
        figsize = (max(8, 1.3 * len(df.columns)), max(3, 0.35 * n + 1.5))
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")
    tbl = ax.table(cellText=df.values, rowLabels=df.index,
                   colLabels=df.columns, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9)
    tbl.auto_set_column_width(col=list(range(len(df.columns) + 1)))
    plt.title(title, fontsize=11, pad=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight"); plt.close()
    return True


def qc_tagdir_stats(cfg, species, sample, qc_dir) -> None:
    """Table of key stats from tagInfo.txt for each combo tag directory.
    Sample-level / combined-only — unchanged behavior, now built on the
    shared row/table helpers used by the new per-replicate version too."""
    combos = _all_combos(cfg, species, sample)
    rows = _tagdir_stats_rows([(_label(d), d) for d in combos])
    if not rows:
        log.info("QC tagdir stats: no tagInfo.txt found for %s/%s", species, sample); return
    _save_stats_table(rows, "Tag Directory Stats", qc_dir / "tagdir_stats.png")
    log.info("QC: tagdir_stats.png")


def qc_tagdir_stats_per_replicate(cfg, species, sample, qc_dir) -> None:
    """Same table as qc_tagdir_stats, but one column per INDIVIDUAL replicate
    TagDir instead of per combo — so a replicate with an unusually low tag
    count, GC content, or median-tags-per-position is visible even though
    it'll be averaged away once merged into the combo TagDir."""
    leaves = _leaf_tagdirs_for_sample(cfg, species, sample)
    rows = _tagdir_stats_rows(leaves)
    if not rows:
        log.info("QC tagdir stats (per-replicate): no tagInfo.txt found for %s/%s",
                 species, sample)
        return
    _save_stats_table(rows, "Tag Directory Stats — per replicate",
                      qc_dir / "tagdir_stats_per_replicate.png", transpose=False)
    log.info("QC: tagdir_stats_per_replicate.png (%d replicate(s))", len(rows))


def _tags_vs_frac_combined(cfg, species, sample, qc_dir) -> None:
    """Log-log scatter of tags-per-position, csRNA + sRNA overlaid on same axes."""
    frames = []
    for assay in ("csRNA", "sRNA"):
        d = _combo(cfg, species, sample, assay)
        if not d or not (d / "tagCountDistribution.txt").exists():
            continue
        df = pd.read_csv(d / "tagCountDistribution.txt", sep="\t")
        df.columns = ["Tags per tag position", "Fraction of Positions"]
        df = df.replace(0, np.nan).dropna()
        df["Tags per tag position"] = np.log(df["Tags per tag position"])
        df["Fraction of Positions"] = np.log(df["Fraction of Positions"])
        df["Library"] = _label(d)
        frames.append(df)

    if not frames:
        log.info("QC tags-vs-frac combined: no data found."); return

    combined = pd.concat(frames, ignore_index=True)
    plt.figure(figsize=(10, 6))
    ax = sns.scatterplot(data=combined, x="Tags per tag position",
                         y="Fraction of Positions", hue="Library", alpha=0.5)
    ax.set_title("csRNA + sRNA — tags vs fraction of positions")
    sns.move_legend(ax, "upper right")
    plt.tight_layout()
    plt.savefig(qc_dir / "tagsPer_Vs_FracofPos.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: tagsPer_Vs_FracofPos.png")


def _a_plot_combined(cfg, species, sample, qc_dir) -> None:
    """A-frequency plot overlaying csRNA and sRNA on same axes."""
    frames = []
    for assay in ("csRNA", "sRNA"):
        d = _combo(cfg, species, sample, assay)
        if not d or not (d / "tagFreqUniq.txt").exists():
            continue
        df = pd.read_csv(d / "tagFreqUniq.txt", sep="\t")
        df = df.iloc[:, :5].set_index("Offset")
        stacked = df.stack().reset_index()
        stacked.columns = ["Distance to TSS", "nt", "pct"]
        stacked["Library"] = _label(d)
        frames.append(stacked[stacked["nt"] == "A"])

    if not frames:
        log.info("QC A-plot combined: no data found."); return

    combined = pd.concat(frames, ignore_index=True)
    plt.figure(figsize=(12, 5))
    ax = sns.lineplot(data=combined, x="Distance to TSS", y="pct",
                      hue="Library", linewidth=2, alpha=0.7)
    ax.set_ylabel("A [%]"); ax.set_title("csRNA + sRNA — A-plot")
    sns.move_legend(ax, "upper right"); plt.tight_layout()
    plt.savefig(qc_dir / "Aplot.png", dpi=150, bbox_inches="tight"); plt.close()
    log.info("QC: Aplot.png")


# ── Per-replicate plots (always heatmap-based, any replicate count) ──────────
# These read straight from each individual (non-combo) replicate's own TagDir,
# so a problem specific to one replicate — a bad library, an alignment issue,
# an unusual length/nucleotide profile — is visible even after replicates get
# merged into the combo TagDir the sample-level plots above are built from.
# Always rendered the same way regardless of replicate count (one row per
# replicate in a heatmap, assay-grouped/sorted) rather than switching
# presentation style at different scales, so a report looks the same whether
# a sample has 2 replicates or 200.

def qc_read_length_per_replicate(cfg, species, sample, qc_dir) -> None:
    leaves = [(n, d) for n, d in _leaf_tagdirs_for_sample(cfg, species, sample)
              if (d / "tagLengthDistribution.txt").exists()]
    if not leaves:
        log.info("QC read-length (per-replicate): no tagLengthDistribution.txt yet.")
        return
    _read_length_grid(leaves, qc_dir)


def _read_length_grid(leaves, qc_dir) -> None:
    """One read-length-distribution line plot per replicate, laid out in a
    fixed grid (wraps at 4 columns) instead of a single heatmap — each
    replicate keeps its own fully visible curve, same as the sample-level
    read_length_distribution.png but split one-subplot-per-replicate."""
    n = len(leaves)
    fig, axes = _replicate_grid(n)
    for ax, (leaf_name, d) in zip(axes.flat, leaves):
        df = pd.read_csv(d / "tagLengthDistribution.txt", sep="\t", index_col=0)
        ax.plot(df.index, df.iloc[:, 0], color="steelblue", linewidth=1.2)
        ax.set_title(leaf_name, fontsize=9)
        ax.tick_params(labelsize=7)

    fig.supxlabel("Read length (nt)")
    fig.supylabel("Fraction of reads")
    fig.suptitle(f"Read length distribution — {n} replicate(s)")
    plt.tight_layout()
    plt.savefig(qc_dir / "read_length_distribution_per_replicate.png",
               dpi=150, bbox_inches="tight")
    plt.close()
    log.info("QC: read_length_distribution_per_replicate.png (%d replicate(s), grid)", n)


def qc_nucleotide_freq_per_replicate(cfg, species, sample, qc_dir) -> None:
    leaves = [(n, d) for n, d in _leaf_tagdirs_for_sample(cfg, species, sample)
              if (d / "tagFreqUniq.txt").exists()]
    if not leaves:
        log.info("QC nt-freq (per-replicate): no tagFreqUniq.txt yet.")
        return
    _nucleotide_freq_grid(leaves, qc_dir)


def _nucleotide_freq_grid(leaves, qc_dir) -> None:
    """One nucleotide-frequency line plot (A/C/G/T) per replicate, laid out
    in a fixed grid instead of a single A-content-only heatmap — same
    all-4-nucleotide lines as the sample-level nucleotide_frequency.png,
    split one-subplot-per-replicate."""
    n = len(leaves)
    fig, axes = _replicate_grid(n)
    for i, (ax, (leaf_name, d)) in enumerate(zip(axes.flat, leaves)):
        df = pd.read_csv(d / "tagFreqUniq.txt", index_col=0, sep="\t")
        cols = [c for c in ("A", "C", "G", "T") if c in df.columns]
        for col in cols:
            ax.plot(df.index, df[col], label=col, linewidth=1)
        ax.set_title(leaf_name, fontsize=9)
        ax.tick_params(labelsize=7)
        if i == 0 and cols:
            ax.legend(fontsize=6, ncol=len(cols), loc="upper right")

    fig.supxlabel("Distance from 5' end of read")
    fig.supylabel("Nucleotide frequency")
    fig.suptitle(f"Nucleotide frequency — {n} replicate(s)")
    plt.tight_layout()
    plt.savefig(qc_dir / "nucleotide_frequency_per_replicate.png",
               dpi=150, bbox_inches="tight")
    plt.close()
    log.info("QC: nucleotide_frequency_per_replicate.png (%d replicate(s), grid)", n)


def qc_autocorrelation_per_replicate(cfg, species, sample, qc_dir) -> None:
    leaves = [(n, d) for n, d in _leaf_tagdirs_for_sample(cfg, species, sample)
              if (d / "tagAutocorrelation.txt").exists()]
    if not leaves:
        log.info("QC autocorrelation (per-replicate): no tagAutocorrelation.txt yet.")
        return
    _autocorrelation_grid(leaves, qc_dir)


def _autocorrelation_grid(leaves, qc_dir) -> None:
    """One same-strand/opposite-strand autocorrelation line plot per
    replicate, laid out in a fixed grid instead of a single heatmap — same
    windowing as the sample-level autocorrelation.png, split
    one-subplot-per-replicate."""
    n = len(leaves)
    fig, axes = _replicate_grid(n)
    for i, (ax, (leaf_name, d)) in enumerate(zip(axes.flat, leaves)):
        df = pd.read_csv(d / "tagAutocorrelation.txt", sep="\t", index_col=0, header=0,
                         names=["Same Strand", "Opposite Strand"]).iloc[1400:2601]
        ax.plot(df.index, df["Same Strand"], label="Same Strand", linewidth=1)
        ax.plot(df.index, df["Opposite Strand"], label="Opposite Strand", linewidth=1)
        ax.set_title(leaf_name, fontsize=9)
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=6, loc="upper right")

    fig.supxlabel("Distance from 5' end of read")
    fig.supylabel("Read counts rel. to 5' end")
    fig.suptitle(f"Autocorrelation — {n} replicate(s)")
    plt.tight_layout()
    plt.savefig(qc_dir / "autocorrelation_per_replicate.png", dpi=150, bbox_inches="tight")
    plt.close()
    log.info("QC: autocorrelation_per_replicate.png (%d replicate(s), grid)", n)


# ── Trim / alignment tool logs ────────────────────────────────────────────────
# These files already exist on disk as a byproduct of trim.py/mapping.py —
# homerTools' .lengths, skewer's -trimmed.log, STAR's Log.final.out, and
# hisat2's _mappingstats.txt — this only reads and summarizes them; it
# doesn't change how trimming/alignment run. Only the compact per-replicate
# summary tables (trim_stats_summary.png, alignment_stats_summary.png) are
# written into QC/ — the raw logs themselves are read in place and left
# where trim.py/mapping.py wrote them, not copied.

def _parse_homer_lengths(path):
    """homerTools trim's <r1>.lengths: a plain 3-col TSV (Length, # reads,
    Fraction). High-confidence format. Returns (total_reads, dimer_pct,
    avg_length_of_nonzero_reads), or (None, None, None) if unreadable/empty.
    dimer_pct is computed from the read counts directly (length==0 row),
    not trusted from HOMER's own '%'-string Fraction column."""
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        log.warning("QC logs: could not read %s: %s", path, exc)
        return None, None, None
    if df.shape[1] < 2 or df.empty:
        return None, None, None
    length_col, reads_col = df.columns[0], df.columns[1]
    total = df[reads_col].sum()
    if total <= 0:
        return None, None, None
    dimer = df.loc[df[length_col] == 0, reads_col].sum()
    dimer_pct = 100.0 * dimer / total
    nonzero = df[df[length_col] != 0]
    avg_len = (float(np.average(nonzero[length_col], weights=nonzero[reads_col]))
              if nonzero[reads_col].sum() > 0 else float("nan"))
    return int(total), dimer_pct, avg_len


def _parse_skewer_log(path):
    """skewer's <prefix>-trimmed.log. Tries several label patterns seen
    across skewer versions and returns None for whatever doesn't match
    rather than risk silently mis-parsing.

    Returns (total_reads, avail_pct) where avail_pct is the % of reads
    retained/available after trimming."""
    import re
    try:
        txt = path.read_text(errors="replace")
    except Exception as exc:
        log.warning("QC logs: could not read %s: %s", path, exc)
        return None, None

    total = None
    for pattern in (
        r"([\d,]+)\s+reads?\s*(?:pairs?)?\s+processed",
        r"reads?\s+processed:\s+([\d,]+)",
    ):
        m = re.search(pattern, txt, re.I)
        if m:
            total = int(m.group(1).replace(",", ""))
            break

    avail_pct = None
    for pattern in (
        r"([\d.]+)\s*%\)?\s+(?:reads?|pairs?)\s+available",
        r"available;\s*suitable for further analysis:\s*[\d,]+\s*\(([\d.]+)%\)",
        r"\(([\d.]+)%\)\s+trimmed reads? available",
    ):
        m = re.search(pattern, txt, re.I)
        if m:
            avail_pct = float(m.group(1))
            break

    return total, avail_pct


def _parse_star_log(path) -> dict:
    """STAR's <prefix>.Log.final.out — a stable, well-documented 'label |
    value' format. High confidence."""
    import re
    try:
        txt = path.read_text(errors="replace")
    except Exception as exc:
        log.warning("QC logs: could not read %s: %s", path, exc)
        return {}
    patterns = {
        "Input Reads":            (r"Number of input reads \|\s+(\d+)", int),
        "Uniquely Mapped %":      (r"Uniquely mapped reads % \|\s+([\d.]+)%", float),
        "Multi-Mapped %":         (r"% of reads mapped to multiple loci \|\s+([\d.]+)%", float),
        "Too-Many-Loci %":        (r"% of reads mapped to too many loci \|\s+([\d.]+)%", float),
        "Unmapped (mismatch) %":  (r"% of reads unmapped: too many mismatches \|\s+([\d.]+)%", float),
        "Unmapped (short) %":     (r"% of reads unmapped: too short \|\s+([\d.]+)%", float),
        "Unmapped (other) %":     (r"% of reads unmapped: other \|\s+([\d.]+)%", float),
    }
    out = {}
    for key, (pattern, cast) in patterns.items():
        m = re.search(pattern, txt)
        if m:
            out[key] = cast(m.group(1))
    return out


def _parse_hisat2_stats(path) -> dict:
    """hisat2's <prefix>_mappingstats.txt (its own stderr, already captured
    by mapping.py's `-S ... 2> {stats}` redirect) — the standard,
    well-documented hisat2 summary format. High confidence."""
    import re
    try:
        txt = path.read_text(errors="replace")
    except Exception as exc:
        log.warning("QC logs: could not read %s: %s", path, exc)
        return {}
    out = {}
    m = re.search(r"(\d+)\s+reads?;\s+of these", txt)
    if m:
        out["Input Reads"] = int(m.group(1))
    m = re.search(r"([\d.]+)\s*%\)\s+aligned 0 times", txt)
    if m:
        out["Unmapped %"] = float(m.group(1))
    m = re.search(r"([\d.]+)\s*%\)\s+aligned exactly 1 time", txt)
    if m:
        out["Uniquely Mapped %"] = float(m.group(1))
    m = re.search(r"([\d.]+)\s*%\)\s+aligned >1 times", txt)
    if m:
        out["Multi-Mapped %"] = float(m.group(1))
    m = re.search(r"([\d.]+)\s*%\s+overall alignment rate", txt)
    if m:
        out["Overall Aligned %"] = float(m.group(1))
    return out


def _render_log_table(df: pd.DataFrame, title: str, out_path) -> None:
    """Row-per-replicate table (same vertical-growth orientation as the
    per-replicate tagdir-stats table) so it stays readable/scrollable even
    with hundreds of replicates, instead of becoming absurdly wide."""
    fig, ax = plt.subplots(figsize=(max(8, 1.3 * len(df.columns)),
                                    max(3, 0.35 * len(df) + 1.5)))
    ax.axis("off")
    tbl = ax.table(cellText=df.values, rowLabels=df.index,
                   colLabels=df.columns, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9)
    tbl.auto_set_column_width(col=list(range(len(df.columns) + 1)))
    plt.title(title, fontsize=11, pad=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight"); plt.close()


def qc_trim_align_summary(cfg, species, sample, qc_dir) -> None:
    """Trim/alignment tool logs, one row per replicate, covering whichever
    tools actually ran: homerTools trim (csRNA/sRNA) or skewer (totalRNA)
    for trimming; STAR or hisat2 for alignment. Writes two summary tables
    only (trim_stats_summary.png, alignment_stats_summary.png) — the
    underlying log files are read in place and are not copied into QC/.

    Both tools' rows share an "Overall Aligned %" column so alignment
    quality is directly comparable regardless of which aligner ran: hisat2
    reports it natively, STAR's is derived as Uniquely + Multi-Mapped %.
    """
    leaves = _leaf_tagdirs_for_sample(cfg, species, sample)
    if not leaves:
        log.info("QC trim/align logs: no replicates for %s/%s", species, sample)
        return

    r1_by_leaf = {ln: r1 for sp, sa, ln, r1 in iter_leaf_dirs(cfg)
                  if sp == species and sa == sample}

    trim_rows, align_rows = [], []
    for leaf_name, _tagdir in leaves:
        r1 = r1_by_leaf.get(leaf_name)
        assay = assay_of_leaf(leaf_name)
        if r1 is None or assay is None:
            continue
        trimmed_dir = cfg.trimmed_dir(species, sample)
        aligned_dir = cfg.aligned_dir(species, sample)
        prefix = r1.name.split("_R1")[0]

        # ── trimming ──
        if assay == "totalRNA":
            skewer_log = trimmed_dir / f"{prefix}-trimmed.log"
            if skewer_log.exists():
                total, avail_pct = _parse_skewer_log(skewer_log)
                trim_rows.append({
                    "Replicate": leaf_name, "Tool": "skewer",
                    "Input Reads": total if total is not None else "NA",
                    "% Retained": avail_pct if avail_pct is not None else "NA",
                    "% Removed": (round(100 - avail_pct, 2)
                                 if avail_pct is not None else "NA"),
                })
            else:
                log.info("QC trim: no skewer log found for %s (%s)", leaf_name, skewer_log)
        else:
            lengths_file = trimmed_dir / f"{r1.name}.lengths"
            if lengths_file.exists():
                total, dimer_pct, _avg_len = _parse_homer_lengths(lengths_file)
                if total is not None:
                    trim_rows.append({
                        "Replicate": leaf_name, "Tool": "homerTools",
                        "Input Reads": total,
                        "% Retained": round(100 - dimer_pct, 2),
                        "% Removed": round(dimer_pct, 2),
                    })

        # ── alignment ──
        star_log = aligned_dir / f"{prefix}.Log.final.out"
        hisat_log = aligned_dir / f"{prefix}_mappingstats.txt"
        if star_log.exists():
            f = _parse_star_log(star_log)
            if f:
                unmapped = [f[k] for k in f if k.startswith("Unmapped") and f[k] is not None]
                uniq = f.get("Uniquely Mapped %")
                multi = f.get("Multi-Mapped %")
                overall = (round(uniq + multi, 2)
                          if uniq is not None and multi is not None else "NA")
                align_rows.append({
                    "Replicate": leaf_name, "Tool": "STAR",
                    "Input Reads": f.get("Input Reads", "NA"),
                    "Uniquely Mapped %": f.get("Uniquely Mapped %", "NA"),
                    "Multi-Mapped %": f.get("Multi-Mapped %", "NA"),
                    "Unmapped %": round(sum(unmapped), 2) if unmapped else "NA",
                    "Overall Aligned %": overall,
                })
        elif hisat_log.exists():
            f = _parse_hisat2_stats(hisat_log)
            if f:
                align_rows.append({
                    "Replicate": leaf_name, "Tool": "hisat2",
                    "Input Reads": f.get("Input Reads", "NA"),
                    "Uniquely Mapped %": f.get("Uniquely Mapped %", "NA"),
                    "Multi-Mapped %": f.get("Multi-Mapped %", "NA"),
                    "Unmapped %": f.get("Unmapped %", "NA"),
                    "Overall Aligned %": f.get("Overall Aligned %", "NA"),
                })

    if trim_rows:
        _render_log_table(pd.DataFrame(trim_rows).set_index("Replicate"),
                          "Trim Summary", qc_dir / "trim_stats_summary.png")
        log.info("QC: trim_stats_summary.png (%d replicate(s))", len(trim_rows))
    else:
        log.info("QC trim summary: no trim logs found for %s/%s", species, sample)

    if align_rows:
        _render_log_table(pd.DataFrame(align_rows).set_index("Replicate"),
                          "Alignment Summary", qc_dir / "alignment_stats_summary.png")
        log.info("QC: alignment_stats_summary.png (%d replicate(s))", len(align_rows))
    else:
        log.info("QC alignment summary: no alignment logs found for %s/%s", species, sample)



# ── Distal vs. proximal TSS pie chart (combined/sample-level) ────────────────

def qc_distal_proximal_pie(cfg, species, sample, qc_dir) -> None:
    """Promoter-proximal vs. distal TSS proportions, straight from the
    combined *.tss.txt's location column (same detection stability.py's
    _location() uses). Always generated as its own standing QC check,
    independent of whether stability/total-RNA info is available — unlike
    stability.py's combined stable/unstable+location pie, which only ever
    falls back to showing location when there's no total RNA, so a csRNA/
    sRNA-only sample (no total RNA) still gets this chart."""
    files = sorted(cfg.sample_tss(species, sample).glob("*.tss.txt"))
    if not files:
        log.info("QC distal/proximal: no *.tss.txt for %s/%s", species, sample)
        return

    locs = []
    for f in files:
        df = _read_homer_tss(f)
        loc, _how = _location(cfg, df)
        if loc is not None:
            locs.append(loc)
    if not locs:
        log.info("QC distal/proximal: no location column found for %s/%s — skipping.",
                 species, sample)
        return

    all_loc = pd.concat(locs, ignore_index=True).dropna()
    counts = all_loc.value_counts().reindex(["proximal", "distal"]).fillna(0)
    total = counts.sum()
    if total <= 0:
        log.info("QC distal/proximal: zero valid TSSs for %s/%s — skipping.", species, sample)
        return

    plt.figure(figsize=(5, 5))
    plt.pie(counts.values,
            labels=[f"{i.capitalize()} ({int(v)}, {v / total * 100:.1f}%)"
                   for i, v in counts.items()],
            colors=[PROX_C, DISTAL_C], startangle=90)
    plt.title(f"{species}/{sample} — Promoter-proximal vs. distal TSSs")
    plt.tight_layout()
    plt.savefig(qc_dir / "distal_proximal_pie.png", dpi=150, bbox_inches="tight")
    plt.close()
    log.info("QC: distal_proximal_pie.png (%s/%s)", species, sample)


# ── Intermediate-file cleanup (NEW) ───────────────────────────────────────────
# Trimming/alignment intermediates (Trimmed/, Aligned/) are large and, by the
# time QC runs, already fully consumed: tag directories were built from them
# earlier in the pipeline, and the numbers useful for QC purposes have
# already been extracted into the summary tables by qc_trim_align_summary
# (run above, before this ever runs) — note the underlying raw logs
# themselves are NOT copied anywhere, so deleting these directories also
# deletes those logs; only the summary tables survive. Opt-in via
# cfg.cleanup_intermediates so existing pipelines keep their current
# behavior unless a caller asks for it.

def qc_cleanup_intermediates(cfg, species, sample, qc_dir) -> None:
    """Remove this sample's intermediate Trimmed/ and Aligned/ directories
    after QC generation, keeping only the final QC tables/reports. No-op
    unless cfg.cleanup_intermediates is truthy. Missing directories,
    missing Config accessors, or filesystem errors are logged and skipped
    rather than raised, so a cleanup failure never breaks QC generation."""
    if not getattr(cfg, "cleanup_intermediates", False):
        return

    import shutil as _shutil

    targets = (
        ("trimmed_dir", "Trimmed"),
        ("aligned_dir", "Aligned"),
    )
    for accessor_name, label in targets:
        accessor = getattr(cfg, accessor_name, None)
        if accessor is None:
            log.warning("QC cleanup: Config has no %s(); skipping %s cleanup for %s/%s",
                       accessor_name, label, species, sample)
            continue
        try:
            d = accessor(species, sample)
        except Exception as exc:
            log.warning("QC cleanup: could not resolve %s dir for %s/%s: %s",
                       label, species, sample, exc)
            continue
        if d is None or not d.is_dir():
            continue
        try:
            _shutil.rmtree(d)
            log.info("QC cleanup: removed intermediate %s dir %s", label, d)
        except Exception as exc:
            log.warning("QC cleanup: failed to remove %s (%s): %s", d, label, exc)


# ── Main entry point ──────────────────────────────────────────────────────────

def _run_qc_one(cfg, species, sample) -> None:
    qc_dir = cfg.sample_qc(species, sample)
    qc_dir.mkdir(parents=True, exist_ok=True)

    # ── Pipeline logs (trim/align tool output) — summary tables only,
    #    rendered first ───────────────────────────────────────────────────
    qc_trim_align_summary(cfg, species, sample, qc_dir)

    # ── Sample-level QC (from combined/-combo TagDirs) — unchanged ──────────
    qc_read_length(cfg, species, sample, qc_dir)
    qc_nucleotide_freq(cfg, species, sample, qc_dir)
    qc_autocorrelation(cfg, species, sample, qc_dir)

    _median_tags_bar_merged(cfg, species, sample, qc_dir)

    _tags_vs_frac_combined(cfg, species, sample, qc_dir)

    _a_plot_combined(cfg, species, sample, qc_dir)

    qc_threshold_optimization(cfg, species, sample, qc_dir)
    qc_tss_nucleotide_freq(cfg, species, sample, qc_dir)

    qc_tsr_summary(cfg, species, sample, qc_dir)
    qc_tsr_annotation(cfg, species, sample, qc_dir)
    qc_tagdir_stats(cfg, species, sample, qc_dir)
    qc_distal_proximal_pie(cfg, species, sample, qc_dir)

    # ── Per-replicate QC (from each individual leaf TagDir) ─────────────────
    qc_read_length_per_replicate(cfg, species, sample, qc_dir)
    qc_nucleotide_freq_per_replicate(cfg, species, sample, qc_dir)
    qc_autocorrelation_per_replicate(cfg, species, sample, qc_dir)
    qc_tagdir_stats_per_replicate(cfg, species, sample, qc_dir)

    # ── Cleanup — last, so every plot/table above has already read
    #    whatever it needed from Trimmed/ and Aligned/ ───────────────────────
    qc_cleanup_intermediates(cfg, species, sample, qc_dir)


def run_qc(cfg) -> None:
    samples = list(iter_samples(cfg))
    if not samples:
        log.info("QC: no Species/Sample dirs found under %s", cfg.project)
        return
    for species, sample in samples:
        log.info("QC: %s/%s", species, sample)
        _run_qc_one(cfg, species, sample)