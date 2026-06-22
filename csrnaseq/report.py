"""Step — HTML report generation.

Writes self-contained HTML files into <project>/Reports/:
  • <sample>.html          — one per TSS/*.tss.txt
  • alltss.html            — TSS/alltss.txt as a table
  • tss_stats.html         — all TSS/*.stats.txt
  • tagdir_<name>.html     — one per TagDir (tagCountDistribution, tagAutocorrelation, tagLengthDistribution)
  • qc_report.html         — QC/*.png (base64-embedded) + data files

Run standalone:
    csrnaseq --steps report
"""
from __future__ import annotations

import base64
import datetime
from pathlib import Path

from .utils import log

# ── Shared CSS ────────────────────────────────────────────────────────────────

_COMMON_CSS = """
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500&family=JetBrains+Mono:wght@400&display=swap');
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: 'Inter', sans-serif; font-size: 13px; color: #111;
         background: #fff; padding: 60px 48px; max-width: 1400px; margin: 0 auto; }
  header { border-bottom: 1px solid #111; padding-bottom: 24px; margin-bottom: 48px;
           display: flex; justify-content: space-between; align-items: center; }
  header h1 { font-size: 22px; font-weight: 500; letter-spacing: -0.02em; }
  .meta { font-size: 11px; color: #888; }
  footer { border-top: 1px solid #e5e5e5; padding-top: 24px; margin-top: 48px;
           font-size: 11px; color: #aaa; }
  .section-label { font-size: 10px; font-weight: 500; text-transform: uppercase;
                   letter-spacing: 0.1em; color: #aaa; margin-bottom: 24px; margin-top: 48px;
                   border-bottom: 1px solid #e5e5e5; padding-bottom: 8px; }
  .table-wrap { overflow-x: auto; }
  table { width: 100%; border-collapse: collapse; }
  th { font-size: 10px; font-weight: 500; text-transform: uppercase; letter-spacing: 0.08em;
       color: #888; border-bottom: 1px solid #e5e5e5; padding: 8px 12px; text-align: left;
       white-space: nowrap; background: #fff; position: sticky; top: 0; }
  td { padding: 7px 12px; border-bottom: 1px solid #f0f0f0; white-space: nowrap; color: #222;
       font-family: 'JetBrains Mono', monospace; font-size: 11px; }
  tr:last-child td { border-bottom: none; }
  tr:hover td { background: #f9f9f9; }
  pre { font-family: 'JetBrains Mono', monospace; font-size: 11px; line-height: 1.7; color: #333;
        background: #fafafa; border: 1px solid #e5e5e5; border-radius: 4px; padding: 16px;
        overflow-x: auto; white-space: pre-wrap; margin-bottom: 32px; }
"""

# ── TSS CSS ───────────────────────────────────────────────────────────────────

_TSS_CSS = _COMMON_CSS + """
  .tss-section { margin-bottom: 56px; }
  .tss-header { display: flex; align-items: baseline; gap: 14px;
                border-bottom: 1px solid #e5e5e5; padding-bottom: 10px; margin-bottom: 0; }
  .tss-header .fname { font-family: 'JetBrains Mono', monospace; font-size: 13px; font-weight: 500; }
  .tss-header .row-count { font-size: 11px; color: #888; }
  .val-tss   { color: #1a6a1a; font-weight: 500; }
  .val-other { color: #888; }
"""

_KEY_COLS = [
    "tssClusterID", "chr", "start", "end", "strand", "score", "annotation",
    "gene associated with annotation", "Closest TSS Symbol",
    "TSS status", "Promoter Proximal/Distal", "Stable/Unstable", "Bidirectional",
]

# ── TagDir CSS ────────────────────────────────────────────────────────────────

_TAGDIR_CSS = _COMMON_CSS + """
  .tagdir-section { margin-bottom: 56px; }
  .tagdir-section h2 { font-size: 14px; font-weight: 500; margin-bottom: 16px;
                       padding-bottom: 10px; border-bottom: 1px solid #e5e5e5; }
"""

# ── QC CSS ────────────────────────────────────────────────────────────────────

_QC_CSS = _COMMON_CSS + """
  .img-grid { display: flex; flex-direction: column; gap: 56px; margin-top: 8px; }
  .img-item { page-break-inside: avoid; }
  .img-item .img-name { font-family: 'JetBrains Mono', monospace; font-size: 11px;
                color: #888; margin-bottom: 10px; }
  .img-item img { max-width: 100%; width: auto; height: auto;
                  border: 1px solid #e5e5e5; border-radius: 2px;
                  display: block; padding: 8px; background: #fff; }
  .img-caption { margin-top: 16px; border-left: 3px solid #e5e5e5;
                 padding-left: 16px; display: flex; flex-direction: column; gap: 8px; }
  .img-caption .cap-row { display: flex; flex-direction: column; gap: 2px; }
  .img-caption .cap-label { font-size: 10px; font-weight: 600; text-transform: uppercase;
                             letter-spacing: 0.07em; color: #aaa; }
  .img-caption .cap-text  { font-size: 12px; color: #444; line-height: 1.55; }
  .cap-desc { font-size: 12px; color: #333; line-height: 1.6; margin-bottom: 4px; }
  .data-section h2 { font-family: 'JetBrains Mono', monospace; font-size: 11px; font-weight: 500;
                     color: #888; margin-bottom: 8px; text-transform: uppercase;
                     letter-spacing: 0.08em; }
"""

_IMG_ORDER = [
    "tagdir_stats", "read_length_distribution", "nucleotide_frequency",
    "median_tags_per_position",
    "csRNA_Aplot", "sRNA_Aplot", "totalRNA_Aplot", "combined_Aplot",
    "tss_nucleotide_frequency",
    "csRNA_tagsPer_Vs_FracofPos", "sRNA_tagsPer_Vs_FracofPos",
    "totalRNA_tagsPer_Vs_FracofPos", "combined_tagsPer_Vs_FracofPos",
    "autocorrelation", "threshold_optimization",
    "tsr_summary", "tsr_annotation",
    "stability_by_location_stacked_bar", "location_stacked_bar", "tsr_pie",
]

# Descriptions keyed by filename stem prefix (matched with startswith)
_IMG_DESCRIPTIONS: dict[str, tuple[str, str, str, str]] = {
    # key: (what it shows, x-axis, y-axis, what to look for)
    "tagdir_stats": (
        "A summary table of basic sequencing statistics for each sample.",
        "Sample",
        "Statistic",
        "Total tag counts should match your expected sequencing depth. "
        "Median tags per position should be close to 1.",
    ),
    "read_length_distribution": (
        "Shows how long the sequenced reads are across all samples.",
        "Read length (nucleotides)",
        "Fraction of reads",
        "Most reads should fall between 20–55 nt with a peak around 30 nt. "
        "Unexpected spikes may indicate contamination from abundant RNA species.",
    ),
    "nucleotide_frequency": (
        "Shows the frequency of each DNA base (A, C, G, T) along the read, "
        "relative to where sequencing started.",
        "Distance from the start of the read (nt)",
        "Nucleotide frequency",
        "csRNA libraries should show an A enrichment at position +1 and "
        "a slight C just before it. A flat or random pattern suggests the "
        "library may not be well-enriched for transcription start sites.",
    ),
    "median_tags_per_position": (
        "Shows how many reads pile up at any single position in the genome, "
        "on average, for each sample.",
        "Median reads per genomic position",
        "Sample",
        "Values should be close to 1, meaning most positions are covered by "
        "a single read. Higher values indicate PCR duplication.",
    ),
    "csRNA_Aplot": (
        "Shows how often the base Adenine (A) appears near transcription "
        "start sites in the csRNA library.",
        "Distance from the transcription start site (nt)",
        "A frequency (%)",
        "There should be a clear A enrichment right at position +1. "
        "A flat line suggests the library is not well-enriched for "
        "true transcription start sites.",
    ),
    "sRNA_Aplot": (
        "Shows Adenine frequency near transcription start sites "
        "in the input control (sRNA) library.",
        "Distance from the transcription start site (nt)",
        "A frequency (%)",
        "The input control should show a weaker or flatter signal than "
        "the csRNA library at position +1.",
    ),
    "totalRNA_Aplot": (
        "Shows Adenine frequency near transcription start sites "
        "in the total RNA library.",
        "Distance from the transcription start site (nt)",
        "A frequency (%)",
        "Total RNA is not enriched for initiating transcripts so "
        "the profile should be relatively flat near position +1.",
    ),
    "combined_Aplot": (
        "csRNA and sRNA A-frequency profiles shown together "
        "for direct comparison.",
        "Distance from the transcription start site (nt)",
        "A frequency (%)",
        "The csRNA library should show a clear +1 A peak while the "
        "sRNA control should be flatter, confirming the csRNA library "
        "is specifically capturing transcription start sites.",
    ),
    "tss_nucleotide_frequency": (
        "Shows the frequency of all four DNA bases around the "
        "primary transcription start site for each sample.",
        "Distance from the primary TSS (nt)",
        "Nucleotide frequency",
        "A clean pattern with A enrichment at +1 and C just upstream "
        "confirms the pipeline accurately identified true start sites.",
    ),
    "csRNA_tagsPer_Vs_FracofPos": (
        "Shows the relationship between how many reads land at a position "
        "and what fraction of positions have that many reads — csRNA libraries.",
        "ln(Reads per position)",
        "ln(Fraction of positions)",
        "Most positions should have very few reads. A steep downward slope "
        "is healthy. A flattened curve at high read counts suggests duplication.",
    ),
    "sRNA_tagsPer_Vs_FracofPos": (
        "Same read depth distribution as above but for the sRNA input "
        "control libraries.",
        "ln(Reads per position)",
        "ln(Fraction of positions)",
        "Input controls often show slightly more high-depth positions than "
        "csRNA due to abundant small RNAs. Overall the shape should still "
        "slope steeply downward.",
    ),
    "totalRNA_tagsPer_Vs_FracofPos": (
        "Same read depth distribution as above but for total RNA libraries.",
        "ln(Reads per position)",
        "ln(Fraction of positions)",
        "Total RNA tends to have a shallower slope than csRNA because its "
        "reads are spread more broadly across the genome rather than "
        "concentrated at start sites.",
    ),
    "combined_tagsPer_Vs_FracofPos": (
        "csRNA and sRNA read depth distributions overlaid for comparison.",
        "ln(Reads per position)",
        "ln(Fraction of positions)",
        "Both library types should follow a similar steep downward slope. "
        "Any sample that looks very different from the others may have "
        "quality or duplication issues.",
    ),
    "autocorrelation": (
        "Shows how reads relate to each other across the genome — "
        "whether reads on the same or opposite strand cluster nearby.",
        "Distance from the start of each read (nt)",
        "Relative read count",
        "csRNA libraries should show reads clustering together on the same "
        "strand and some signal on the opposite strand nearby, reflecting "
        "the bidirectional nature of transcription. Input controls should "
        "show much weaker clustering.",
    ),
    "threshold_optimization": (
        "Shows how the pipeline chose the cutoff score to distinguish real "
        "transcription start sites from background signal in the input control.",
        "Log2 ratio of csRNA signal vs. input control",
        "Cumulative fraction of sites",
        "The two curves should separate clearly, with true TSS sites rising "
        "at higher ratios than exon background. The dotted line marks the "
        "threshold the pipeline selected. Curves that overlap closely "
        "suggest the sample may have high background contamination.",
    ),
    "tsr_summary": (
        "A table summarizing the key numbers from TSS calling for each sample: "
        "read counts, how many start sites were found, and their breakdown "
        "by stability and location.",
        "Sample",
        "Metric",
        "Valid TSS counts should be in the tens of thousands for a typical "
        "mammalian sample. Enrichment over input should be clearly positive.",
    ),
    "tsr_annotation": (
        "Shows what types of genomic regions the identified transcription "
        "start sites fall in, broken down by sample.",
        "Sample",
        "Number of TSR clusters",
        "Most start sites should fall near annotated gene promoters. "
        "A large fraction landing on gene bodies may indicate background "
        "contamination from highly expressed genes.",
    ),
    "stability_by_location_stacked_bar": (
        "Breaks down transcription start sites by whether they produce "
        "stable or unstable transcripts, and whether they are near or "
        "far from annotated gene promoters. Requires total RNA data.",
        "Stability class",
        "Number of TSR clusters",
        "Stable start sites reflect active gene promoters. Unstable ones "
        "are more transient and often correspond to enhancers or divergent "
        "transcription. The ratio will vary depending on the biology "
        "of your samples.",
    ),
    "location_stacked_bar": (
        "Shows how many transcription start sites fall near annotated "
        "promoters versus further away, per sample. Shown when total "
        "RNA is not available.",
        "Sample",
        "Number of TSR clusters",
        "Start sites near promoters are likely from active genes. Those "
        "further away may represent enhancers or unannotated promoters.",
    ),
    "tsr_pie": (
        "A pooled summary across all samples showing the overall "
        "proportion of stable vs. unstable start sites, or proximal "
        "vs. distal if total RNA is not available.",
        "N/A",
        "N/A",
        "Gives a quick at-a-glance view of the dominant class of "
        "transcription start sites across the experiment.",
    ),
}


def _get_description(stem: str) -> tuple[str, str, str, str] | None:
    """Return (description, x, y, look_for) for a plot stem, or None."""
    for prefix, desc in _IMG_DESCRIPTIONS.items():
        if stem == prefix or stem.startswith(prefix):
            return desc
    return None

# ── Helpers ───────────────────────────────────────────────────────────────────

def _html_page(title: str, project_name: str, now: str, body: str, css: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8">
<title>{title} — {project_name}</title>
<style>{css}</style>
</head>
<body>
<header>
  <div>
    <h1>{title}</h1>
    <div class="meta">Homerun Pipeline &mdash; {project_name}</div>
  </div>
  <div class="meta">Generated {now}</div>
</header>
{body}
<footer>Homerun Pipeline &mdash; {project_name} &mdash; exported {now}</footer>
</body></html>"""


def _tsv_to_html_table(path: Path) -> str:
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        return "<p style='color:#888'>Empty file.</p>"
    header = lines[0].lstrip("#").split("\t")
    head = "".join(f"<th>{c}</th>" for c in header)
    body = ""
    for line in lines[1:]:
        cells = line.split("\t")
        body += "<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>"
    return (f'<div class="table-wrap"><table>'
            f"<thead><tr>{head}</tr></thead>"
            f"<tbody>{body}</tbody></table></div>")


def _cell_class(col: str, val: str) -> str:
    if col == "annotation":
        if val == "tss":   return " class=\"val-tss\""
        if val == "other": return " class=\"val-other\""
    return ""


def _read_tss(path: Path) -> tuple[list[str], list[dict]]:
    with open(path) as fh:
        lines = fh.readlines()
    if not lines:
        return [], []
    header = lines[0].lstrip("#").strip().split("\t")
    rows = [dict(zip(header, line.strip().split("\t")))
            for line in lines[1:] if line.strip()]
    return header, rows


def _tss_table(header: list[str], rows: list[dict]) -> str:
    cols = [c for c in _KEY_COLS if c in header] or header
    head = "".join(f"<th>{c}</th>" for c in cols)
    body = ""
    for row in rows:
        cells = "".join(
            f"<td{_cell_class(c, row.get(c,''))}>{row.get(c,'')}</td>"
            for c in cols
        )
        body += f"<tr>{cells}</tr>"
    return (f'<div class="table-wrap"><table>'
            f"<thead><tr>{head}</tr></thead>"
            f"<tbody>{body}</tbody></table></div>")


def _img_sort_key(name: str) -> tuple:
    stem = Path(name).stem
    for i, prefix in enumerate(_IMG_ORDER):
        if stem == prefix or stem.startswith(prefix):
            return (i, stem)
    return (len(_IMG_ORDER), stem)

# ── TSS builders ─────────────────────────────────────────────────────────────

def _build_sample_tss_html(cfg, f: Path, now: str) -> str:
    header, rows = _read_tss(f)
    if not rows:
        body = '<p style="color:#888">No rows found.</p>'
    else:
        body = f"""
  <section class="tss-section">
    <div class="tss-header">
      <span class="fname">{f.name}</span>
      <span class="row-count">{len(rows):,} clusters</span>
    </div>
    {_tss_table(header, rows)}
  </section>"""
    return _html_page(f.stem, cfg.project.name, now, body, _TSS_CSS)


def _build_alltss_html(cfg, alltss: Path, now: str) -> str:
    body = f"""
  <section class="tss-section">
    <div class="tss-header">
      <span class="fname">alltss.txt</span>
      <span class="row-count">All samples combined</span>
    </div>
    {_tsv_to_html_table(alltss)}
  </section>"""
    return _html_page("All TSS (combined)", cfg.project.name, now, body, _TSS_CSS)


def _build_stats_html(cfg, stats_files: list, now: str) -> str:
    sections = ""
    for f in stats_files:
        content = f.read_text(errors="replace")
        safe = content.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        sections += f"""
  <section class="tss-section">
    <div class="section-label">{f.name}</div>
    <pre>{safe}</pre>
  </section>"""
    body = sections or '<p style="color:#888">No stats files found.</p>'
    return _html_page("TSS Stats", cfg.project.name, now, body, _TSS_CSS)

# ── TagDir builders ───────────────────────────────────────────────────────────

_TAGDIR_FILES = [
    "tagCountDistribution.txt",
    "tagAutocorrelation.txt",
    "tagLengthDistribution.txt",
]


def _build_tagdir_html(cfg, tagdir: Path, now: str) -> str:
    sections = ""
    for fname in _TAGDIR_FILES:
        fpath = tagdir / fname
        if not fpath.exists():
            sections += f"""
  <div class="tagdir-section">
    <h2>{fname}</h2>
    <p style="color:#888">File not found.</p>
  </div>"""
            continue
        # try rendering as table first; fall back to preformatted text
        try:
            table_html = _tsv_to_html_table(fpath)
            sections += f"""
  <div class="tagdir-section">
    <h2>{fname}</h2>
    {table_html}
  </div>"""
        except Exception as exc:
            log.warning("report: could not parse %s: %s", fpath, exc)
            safe = fpath.read_text(errors="replace").replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")
            sections += f"""
  <div class="tagdir-section">
    <h2>{fname}</h2>
    <pre>{safe}</pre>
  </div>"""

    body = sections or '<p style="color:#888">No tagdir files found.</p>'
    return _html_page(tagdir.name, cfg.project.name, now, body, _TAGDIR_CSS)

# ── QC builder ────────────────────────────────────────────────────────────────

def _build_qc_html(cfg, now: str) -> str:
    imgs = sorted(
        [f for f in cfg.qc.glob("*.png") if f.is_file()],
        key=lambda f: _img_sort_key(f.name),
    )
    txts = sorted(
        f for f in cfg.qc.iterdir()
        if f.is_file() and f.suffix in (".txt", ".tsv", ".csv")
    )

    img_section = ""
    if imgs:
        items = ""
        for f in imgs:
            try:
                b64 = base64.b64encode(f.read_bytes()).decode()
                desc = _get_description(f.stem)
                caption_html = ""
                if desc:
                    what, xax, yax, look = desc
                    caption_html = f"""
                <div class="img-caption">
                  <div class="cap-desc">{what}</div>
                  <div class="cap-row"><span class="cap-label">X-axis</span><span class="cap-text">{xax}</span></div>
                  <div class="cap-row"><span class="cap-label">Y-axis</span><span class="cap-text">{yax}</span></div>
                  <div class="cap-row"><span class="cap-label">Look for</span><span class="cap-text">{look}</span></div>
                </div>"""
                items += (f'<div class="img-item">'
                          f'<div class="img-name">{f.name}</div>'
                          f'<img src="data:image/png;base64,{b64}" alt="{f.name}">'
                          f'{caption_html}'
                          f"</div>")
            except Exception as exc:
                log.warning("report: could not embed %s: %s", f.name, exc)
        if items:
            img_section = (f'<div class="section-label">Images</div>'
                           f'<div class="img-grid">{items}</div>')

    txt_section = ""
    if txts:
        parts = ""
        for f in txts:
            try:
                content = f.read_text(errors="replace")
                safe = content.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")
                parts += f'<div class="data-section"><h2>{f.name}</h2><pre>{safe}</pre></div>'
            except Exception as exc:
                log.warning("report: could not read %s: %s", f.name, exc)
        if parts:
            txt_section = f'<div class="section-label">Data Files</div>{parts}'

    body = (img_section + txt_section) or '<p style="color:#888">No QC files found.</p>'
    return _html_page("QC Report", cfg.project.name, now, body, _QC_CSS)

# ── Entry point ───────────────────────────────────────────────────────────────

def run_report(cfg) -> None:
    reports_dir = cfg.project / "Reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    # ── per-sample TSS reports ────────────────────────────────────────────────
    tss_files = sorted(cfg.tss.glob("*.tss.txt"))
    for f in tss_files:
        out = reports_dir / f"{f.stem}.html"
        out.write_text(_build_sample_tss_html(cfg, f, now), encoding="utf-8")
        log.info("report: %s", out.name)

    # ── alltss  (*.alltss.txt) ────────────────────────────────────────────────
    alltss_files = sorted(cfg.tss.glob("*.alltss.txt"))
    if alltss_files:
        sections = ""
        for f in alltss_files:
            sections += f"""
  <section class="tss-section">
    <div class="tss-header">
      <span class="fname">{f.name}</span>
    </div>
    {_tsv_to_html_table(f)}
  </section>"""
        body = sections
        html = _html_page("All TSS Candidates", cfg.project.name, now, body, _TSS_CSS)
        out = reports_dir / "alltss.html"
        out.write_text(html, encoding="utf-8")
        log.info("report: alltss.html (%d file(s))", len(alltss_files))

    # ── stats ─────────────────────────────────────────────────────────────────
    stats_files = sorted(cfg.tss.glob("*.stats.txt"))
    if stats_files:
        out = reports_dir / "tss_stats.html"
        out.write_text(_build_stats_html(cfg, stats_files, now), encoding="utf-8")
        log.info("report: tss_stats.html (%d file(s))", len(stats_files))

    if not tss_files and not alltss_files and not stats_files:
        log.info("report: no TSS files found in %s — skipping TSS reports", cfg.tss)

    # ── per-tagdir reports ────────────────────────────────────────────────────
    tagdirs = sorted(d for d in cfg.tagdirs.iterdir() if d.is_dir())
    for tagdir in tagdirs:
        has_any = any((tagdir / f).exists() for f in _TAGDIR_FILES)
        if not has_any:
            continue
        out = reports_dir / f"tagdir_{tagdir.name}.html"
        out.write_text(_build_tagdir_html(cfg, tagdir, now), encoding="utf-8")
        log.info("report: tagdir_%s.html", tagdir.name)

    # ── QC report ─────────────────────────────────────────────────────────────
    qc_imgs = list(cfg.qc.glob("*.png"))
    qc_txts = [f for f in cfg.qc.iterdir()
               if f.is_file() and f.suffix in (".txt", ".tsv", ".csv")]
    if qc_imgs or qc_txts:
        out = reports_dir / "qc_report.html"
        out.write_text(_build_qc_html(cfg, now), encoding="utf-8")
        log.info("report: qc_report.html (%d image(s), %d data file(s))",
                 len(qc_imgs), len(qc_txts))
    else:
        log.info("report: no QC files found in %s — skipping QC report", cfg.qc)
