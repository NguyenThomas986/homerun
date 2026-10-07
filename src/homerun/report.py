"""Step — HTML + PDF report generation.

Writes one self-contained qc_report.html AND one qc_report.pdf per sample,
both saved under Species/QC/<sample>/ — same content, section order, and image
captions in both. The HTML embeds PNGs as base64; the PDF embeds them as
native reportlab Image flowables built straight from the same PNG files
and the same _IMG_ORDER*/_IMG_DESCRIPTIONS tables, so the two never drift
apart from each other.

PDF generation needs the `reportlab` package (pure-Python, no system
libraries like Cairo/Pango required — deliberately chosen over an
HTML-to-PDF renderer for that reason, since this pipeline runs on HPC
compute nodes where installing system packages usually isn't an option).
If reportlab isn't installed, the PDF is skipped with a warning and the
HTML report is still written — a missing PDF dependency never blocks the
report step.

Run standalone:
    csrnaseq --steps report
"""
from __future__ import annotations

import base64
import datetime
from pathlib import Path

from .utils import log, iter_samples

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
  pre { font-family: 'JetBrains Mono', monospace; font-size: 11px; line-height: 1.7; color: #333;
        background: #fafafa; border: 1px solid #e5e5e5; border-radius: 4px; padding: 16px;
        overflow-x: auto; white-space: pre-wrap; margin-bottom: 32px; }
"""

_QC_CSS = _COMMON_CSS + """
  .img-grid { display: flex; flex-direction: column; gap: 56px; margin-top: 8px; }
  .img-item { page-break-inside: avoid; }
  .img-item .img-name { font-family: 'JetBrains Mono', monospace; font-size: 11px;
                color: #888; margin-bottom: 10px; }
  .img-item img { max-width: 100%; width: auto; height: auto;
                  border: 1px solid #e5e5e5; border-radius: 2px;
                  display: block; padding: 8px; background: #fff; }
  .img-caption { margin-top: 12px; border-left: 3px solid #e5e5e5; padding-left: 16px; }
  .img-caption .cap-source { font-size: 10px; font-weight: 600; text-transform: uppercase;
                             letter-spacing: 0.07em; color: #aaa; margin-bottom: 4px; }
  .img-caption .cap-desc { font-size: 12px; color: #333; line-height: 1.55; }
  .data-section h2 { font-family: 'JetBrains Mono', monospace; font-size: 11px; font-weight: 500;
                     color: #888; margin-bottom: 8px; text-transform: uppercase;
                     letter-spacing: 0.08em; }
  .section-note { font-size: 12px; color: #888; margin-top: -16px; margin-bottom: 24px; }
"""

_IMG_ORDER = [
    "tagdir_stats", "read_length_distribution", "nucleotide_frequency",
    "median_tags_per_position",
    "Aplot",
    "tss_nucleotide_frequency",
    "tagsPer_Vs_FracofPos",
    "autocorrelation", "threshold_optimization",
    "tsr_summary", "tsr_annotation", "ritrie",
    "stability_by_location_stacked_bar", "location_stacked_bar", "tsr_pie",
    "distal_proximal_pie",
]

# Per-replicate images (from individual leaf TagDirs, not the combo) get
# their own report section — see _PER_REPLICATE_STEMS below. Ordered
# separately from _IMG_ORDER so the two sections don't need to interleave.
_IMG_ORDER_PER_REPLICATE = [
    "read_length_distribution_per_replicate",
    "nucleotide_frequency_per_replicate",
    "autocorrelation_per_replicate",
    "tagdir_stats_per_replicate",
]

_PER_REPLICATE_STEMS = frozenset(_IMG_ORDER_PER_REPLICATE)

# Trim/alignment tool log summaries get their own section, rendered FIRST
# (above Sample-Level QC) — this is the "what actually happened during
# trim/align" answer, so it's the first thing in the report rather than
# buried after everything else.
_IMG_ORDER_LOGS = ["trim_stats_summary", "alignment_stats_summary"]

_LOG_STEMS = frozenset(_IMG_ORDER_LOGS)

# key: filename stem prefix (matched with startswith) -> (source, short description)
_IMG_DESCRIPTIONS: dict[str, tuple[str, str]] = {
    "tagdir_stats": (
        "tagInfo.txt (combo TagDirs)",
        "Summary of basic sequencing stats for each combo library.",
    ),
    "read_length_distribution": (
        "tagLengthDistribution.txt (combo TagDirs)",
        "Read length distribution for each combo library.",
    ),
    "nucleotide_frequency": (
        "tagFreqUniq.txt (combo TagDirs)",
        "Base composition near the start of each read, per library — a "
        "pre-TSS-calling check on library quality.",
    ),
    "median_tags_per_position": (
        "tagCountDistribution.txt (combo TagDirs)",
        "Median reads stacked at the same genomic position, per library. "
        "Should be close to 1 — higher values mean many reads are just "
        "PCR copies of the same original molecule rather than independent "
        "biological signal.",
    ),
    "Aplot": (
        "tagFreqUniq.txt (csRNA + sRNA combo TagDirs)",
        "A-frequency near the TSS, csRNA vs. sRNA overlaid.",
    ),
    "tss_nucleotide_frequency": (
        "*.freq.tsv (TSS/)",
        "Base composition around the called primary TSS — checks whether "
        "the called TSS positions look biologically real.",
    ),
    "tagsPer_Vs_FracofPos": (
        "tagCountDistribution.txt (csRNA + sRNA combo TagDirs)",
        "Read-depth distribution across genomic positions, csRNA vs. sRNA.",
    ),
    "autocorrelation": (
        "tagAutocorrelation.txt (combo TagDirs)",
        "How reads cluster relative to each other and to strand.",
    ),
    "threshold_optimization": (
        "*.inputDistribution.txt (TSS/)",
        "How the csRNA-vs-input cutoff score was chosen.",
    ),
    "tsr_summary": (
        "*.stats.txt (TSS/)",
        "Summary table of key TSS-calling numbers per sample.",
    ),
    "tsr_annotation": (
        "*.tss.txt (TSS/)",
        "Genomic annotation breakdown of the called TSR clusters.",
    ),
    "stability_by_location_stacked_bar": (
        "*.tss.txt (TSS/)",
        "Stable vs. unstable TSRs, split by genomic location. Requires total RNA.",
    ),
    "location_stacked_bar": (
        "*.tss.txt (TSS/)",
        "TSRs by genomic location, shown when total RNA isn't available.",
    ),
    "tsr_pie": (
        "*.tss.txt (TSS/)",
        "Pooled stable/unstable split across the sample.",
    ),
    "distal_proximal_pie": (
        "*.tss.txt (TSS/)",
        "Promoter-proximal vs. distal split across every called TSR, "
        "generated regardless of whether total RNA/stability is available.",
    ),
    "read_length_distribution_per_replicate": (
        "tagLengthDistribution.txt (each individual leaf TagDir)",
        "Read length distribution for each individual replicate — an "
        "unusual profile here can be masked once replicates are merged "
        "into the combo library above.",
    ),
    "nucleotide_frequency_per_replicate": (
        "tagFreqUniq.txt (each individual leaf TagDir)",
        "Base composition near the start of each read, per individual "
        "replicate rather than the merged combo library.",
    ),
    "autocorrelation_per_replicate": (
        "tagAutocorrelation.txt (each individual leaf TagDir)",
        "How reads cluster relative to each other and to strand, per "
        "individual replicate rather than the merged combo library.",
    ),
    "tagdir_stats_per_replicate": (
        "tagInfo.txt (each individual leaf TagDir)",
        "Same summary stats as the combo table above, but one row per "
        "individual replicate — a low tag count, unusual GC content, or "
        "poor median-tags-per-position on one replicate is visible here "
        "even though it would be averaged away once merged into the combo.",
    ),
    "trim_stats_summary": (
        "homerTools .lengths / skewer -trimmed.log (each replicate's Trimmed/)",
        "Input reads and % removed per replicate — adapter-dimer rate for "
        "homerTools (csRNA/sRNA), best-effort parse for skewer (totalRNA).",
    ),
    "alignment_stats_summary": (
        "STAR Log.final.out / hisat2 _mappingstats.txt (each replicate's Aligned/)",
        "Input reads and the uniquely-mapped / multi-mapped / unmapped "
        "breakdown per replicate, whichever aligner actually ran.",
    ),
    "ritrie": (
        "RITRIE/ritrie.tsv (per csRNA replicate)",
        "Higher RIT/RIE values generally indicate stronger enrichment of "
        "transcription initiation signal relative to exonic RNA, while lower "
        "values may indicate increased RNA degradation. Thresholds have not "
        "yet been formally established.",
),
}


def _get_description(stem: str) -> tuple[str, str] | None:
    """Return (source, description) for a plot stem, or None. Exact match is
    tried first — otherwise e.g. 'read_length_distribution_per_replicate'
    would match the shorter 'read_length_distribution' prefix before ever
    reaching its own, more specific entry."""
    if stem in _IMG_DESCRIPTIONS:
        return _IMG_DESCRIPTIONS[stem]
    for prefix, desc in _IMG_DESCRIPTIONS.items():
        if stem.startswith(prefix):
            return desc
    return None


def _img_sort_key(name: str, order: list[str]) -> tuple:
    stem = Path(name).stem
    if stem in order:
        return (order.index(stem), stem)
    for i, prefix in enumerate(order):
        if stem.startswith(prefix):
            return (i, stem)
    return (len(order), stem)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _html_page(title: str, subtitle: str, now: str, body: str, css: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8">
<title>{title}</title>
<style>{css}</style>
</head>
<body>
<header>
  <div>
    <h1>QC Report</h1>
    <div class="meta">Homerun Pipeline &mdash; {subtitle}</div>
  </div>
  <div class="meta">Generated {now}</div>
</header>
{body}
<footer>Homerun Pipeline &mdash; {subtitle} &mdash; exported {now}</footer>
</body></html>"""


# ── QC builder ────────────────────────────────────────────────────────────────

def _render_img_grid(imgs: list[Path], heading: str, note: str = "") -> str:
    """Render one section-labeled grid of images, or '' if imgs is empty."""
    if not imgs:
        return ""
    items = ""
    for f in imgs:
        try:
            b64 = base64.b64encode(f.read_bytes()).decode()
            desc = _get_description(f.stem)
            caption_html = ""
            if desc:
                source, text = desc
                caption_html = (f'<div class="img-caption">'
                                f'<div class="cap-source">Source: {source}</div>'
                                f'<div class="cap-desc">{text}</div>'
                                f'</div>')
            items += (f'<div class="img-item">'
                      f'<div class="img-name">{f.name}</div>'
                      f'<img src="data:image/png;base64,{b64}" alt="{f.name}">'
                      f'{caption_html}'
                      f"</div>")
        except Exception as exc:
            log.warning("report: could not embed %s: %s", f.name, exc)
    if not items:
        return ""
    note_html = f'<div class="section-note">{note}</div>' if note else ""
    return (f'<div class="section-label">{heading}</div>'
           f'{note_html}'
           f'<div class="img-grid">{items}</div>')


def _build_qc_html(species: str, sample: str, qc_dir: Path, now: str) -> str:
    all_imgs = [f for f in qc_dir.glob("*.png") if f.is_file()]

    # Three sections, sorted independently against their own order list so
    # they never interleave: pipeline logs (first), sample-level (from the
    # combined/-combo TagDirs), per-replicate (from each individual leaf
    # TagDir).
    log_imgs = sorted(
        (f for f in all_imgs if f.stem in _LOG_STEMS),
        key=lambda f: _img_sort_key(f.name, _IMG_ORDER_LOGS),
    )
    per_rep_imgs = sorted(
        (f for f in all_imgs if f.stem in _PER_REPLICATE_STEMS),
        key=lambda f: _img_sort_key(f.name, _IMG_ORDER_PER_REPLICATE),
    )
    sample_imgs = sorted(
        (f for f in all_imgs
        if f.stem not in _PER_REPLICATE_STEMS and f.stem not in _LOG_STEMS),
        key=lambda f: _img_sort_key(f.name, _IMG_ORDER),
    )

    txts = sorted(
        f for f in qc_dir.iterdir()
        if f.is_file() and f.suffix in (".txt", ".tsv", ".csv")
    )

    img_section = _render_img_grid(
        log_imgs, "Pipeline Logs",
        note="Trim and alignment tool summaries, one row per replicate.",
    )
    img_section += _render_img_grid(sample_imgs, "Sample-Level QC")
    img_section += _render_img_grid(
        per_rep_imgs, "Per-Replicate QC",
        note=("Generated from each individual replicate's own TagDir, not the "
              "merged combo library above — useful for spotting a problem "
              "replicate before it gets averaged away."),
    )

    txt_section = ""
    if txts:
        parts = ""
        for f in txts:
            try:
                content = f.read_text(errors="replace")
                safe = content.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
                parts += f'<div class="data-section"><h2>{f.name}</h2><pre>{safe}</pre></div>'
            except Exception as exc:
                log.warning("report: could not read %s: %s", f.name, exc)
        if parts:
            txt_section = f'<div class="section-label">Data Files</div>{parts}'

    body = (img_section + txt_section) or '<p style="color:#888">No QC files found.</p>'
    return _html_page(f"QC Report — {species}/{sample}", f"{species}/{sample}", now, body, _QC_CSS)


# ── PDF builder ───────────────────────────────────────────────────────────────
# Mirrors _build_qc_html's structure section-for-section (Pipeline Logs,
# Sample-Level QC, Per-Replicate QC, Data Files) using the same
# _IMG_ORDER*/_IMG_DESCRIPTIONS tables, so the PDF and HTML reports never
# show different content or a different order — only the rendering differs.

def _pdf_import_error_hint() -> str:
    return ("PDF report skipped — the 'reportlab' package isn't installed. "
           "Install it with `pip install reportlab` to also get a "
           "qc_report.pdf alongside qc_report.html.")


def _pdf_styles():
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="RptTitle", fontName="Helvetica-Bold", fontSize=18, leading=22,
        spaceAfter=2, alignment=TA_LEFT))
    styles.add(ParagraphStyle(
        name="RptMeta", fontName="Helvetica", fontSize=9,
        textColor=colors.HexColor("#888888"), spaceAfter=4, alignment=TA_LEFT))
    styles.add(ParagraphStyle(
        name="SectionLabel", fontName="Helvetica-Bold", fontSize=10,
        textColor=colors.HexColor("#999999"), spaceBefore=20, spaceAfter=8,
        alignment=TA_LEFT))
    styles.add(ParagraphStyle(
        name="SectionNote", fontName="Helvetica-Oblique", fontSize=9,
        textColor=colors.HexColor("#888888"), spaceAfter=10, alignment=TA_LEFT))
    styles.add(ParagraphStyle(
        name="ImgName", fontName="Courier", fontSize=8,
        textColor=colors.HexColor("#888888"), spaceAfter=4, alignment=TA_LEFT))
    styles.add(ParagraphStyle(
        name="CapSource", fontName="Helvetica-Bold", fontSize=8,
        textColor=colors.HexColor("#aaaaaa"), spaceBefore=6, spaceAfter=2,
        alignment=TA_LEFT))
    styles.add(ParagraphStyle(
        name="CapDesc", fontName="Helvetica", fontSize=9, leading=13,
        textColor=colors.HexColor("#333333"), spaceAfter=16, alignment=TA_LEFT))
    styles.add(ParagraphStyle(
        name="DataHeading", fontName="Courier-Bold", fontSize=9,
        textColor=colors.HexColor("#888888"), spaceBefore=16, spaceAfter=4,
        alignment=TA_LEFT))
    styles.add(ParagraphStyle(
        name="DataBody", fontName="Courier", fontSize=6.5, leading=8.5,
        textColor=colors.HexColor("#333333"), spaceAfter=14, alignment=TA_LEFT))
    return styles


def _pdf_image_flowable(path: Path, avail_width: float, avail_height: float):
    """A reportlab Image scaled to fit within (avail_width, avail_height) —
    constrained by whichever dimension is tighter, so a single image can
    never be too large for one page (which would otherwise raise a
    LayoutError and abort the whole PDF). Very tall per-replicate grid PNGs
    (many replicates -> many rows, see qc.py's _replicate_grid) get scaled
    down a lot as a result — legible detail trades off against guaranteed
    single-page placement; this mirrors the same "many replicates" scaling
    tradeoff qc.py already makes for the PNGs themselves."""
    from reportlab.lib.utils import ImageReader
    from reportlab.platypus import Image as RLImage

    img_w, img_h = ImageReader(str(path)).getSize()
    scale = min(avail_width / img_w, avail_height / img_h, 1.0)
    return RLImage(str(path), width=img_w * scale, height=img_h * scale)


def _pdf_img_block(f: Path, styles, avail_width: float, avail_height: float):
    """Flowables for one image: filename label, the image itself, and its
    source/description caption if one exists — kept together so they never
    split across a page break (safe because the image is already scaled to
    fit a single page by _pdf_image_flowable)."""
    from reportlab.platypus import Paragraph, Spacer, KeepTogether

    block = [
        Paragraph(f.name, styles["ImgName"]),
        _pdf_image_flowable(f, avail_width, avail_height * 0.85),
    ]
    desc = _get_description(f.stem)
    if desc:
        source, text = desc
        block.append(Paragraph(f"Source: {source}", styles["CapSource"]))
        block.append(Paragraph(text, styles["CapDesc"]))
    else:
        block.append(Spacer(1, 16))
    return KeepTogether(block)


def _pdf_img_section(story: list, imgs: list[Path], heading: str, styles,
                     avail_width: float, avail_height: float, note: str = "") -> None:
    if not imgs:
        return
    from reportlab.platypus import Paragraph

    story.append(Paragraph(heading, styles["SectionLabel"]))
    if note:
        story.append(Paragraph(note, styles["SectionNote"]))
    for f in imgs:
        try:
            story.append(_pdf_img_block(f, styles, avail_width, avail_height))
        except Exception as exc:
            log.warning("report (pdf): could not embed %s: %s", f.name, exc)


def _build_qc_pdf(species: str, sample: str, qc_dir: Path, now: str, out_path: Path) -> None:
    from reportlab.lib.pagesizes import letter
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Preformatted, PageBreak

    all_imgs = [f for f in qc_dir.glob("*.png") if f.is_file()]
    log_imgs = sorted(
        (f for f in all_imgs if f.stem in _LOG_STEMS),
        key=lambda f: _img_sort_key(f.name, _IMG_ORDER_LOGS),
    )
    per_rep_imgs = sorted(
        (f for f in all_imgs if f.stem in _PER_REPLICATE_STEMS),
        key=lambda f: _img_sort_key(f.name, _IMG_ORDER_PER_REPLICATE),
    )
    sample_imgs = sorted(
        (f for f in all_imgs
        if f.stem not in _PER_REPLICATE_STEMS and f.stem not in _LOG_STEMS),
        key=lambda f: _img_sort_key(f.name, _IMG_ORDER),
    )
    txts = sorted(
        f for f in qc_dir.iterdir()
        if f.is_file() and f.suffix in (".txt", ".tsv", ".csv")
    )

    margin = 48
    doc = SimpleDocTemplate(
        str(out_path), pagesize=letter,
        leftMargin=margin, rightMargin=margin, topMargin=margin, bottomMargin=margin,
    )
    avail_width = letter[0] - 2 * margin
    avail_height = letter[1] - 2 * margin

    styles = _pdf_styles()
    story = [
        Paragraph("QC Report", styles["RptTitle"]),
        Paragraph(f"Homerun Pipeline &mdash; {species}/{sample}", styles["RptMeta"]),
        Paragraph(f"Generated {now}", styles["RptMeta"]),
        Spacer(1, 12),
    ]

    _pdf_img_section(story, log_imgs, "Pipeline Logs", styles, avail_width, avail_height,
                     note="Trim and alignment tool summaries, one row per replicate.")
    _pdf_img_section(story, sample_imgs, "Sample-Level QC", styles, avail_width, avail_height)
    _pdf_img_section(
        story, per_rep_imgs, "Per-Replicate QC", styles, avail_width, avail_height,
        note=("Generated from each individual replicate's own TagDir, not the "
              "merged combo library above — useful for spotting a problem "
              "replicate before it gets averaged away."),
    )

    if txts:
        story.append(Paragraph("Data Files", styles["SectionLabel"]))
        for f in txts:
            try:
                content = f.read_text(errors="replace")
            except Exception as exc:
                log.warning("report (pdf): could not read %s: %s", f.name, exc)
                continue
            # Preformatted's built-in fonts render a literal tab character as
            # a missing-glyph black box (same failure mode as unsupported
            # unicode chars) — expand to spaces first so TSV/CSV content is
            # actually legible instead of a row of black squares.
            content = content.expandtabs(4)
            story.append(Paragraph(f.name, styles["DataHeading"]))
            story.append(Preformatted(content, styles["DataBody"]))

    if len(story) <= 4:
        story.append(Paragraph("No QC files found.", styles["RptMeta"]))

    doc.build(story)


# ── Entry point ───────────────────────────────────────────────────────────────

def run_report(cfg) -> None:
    samples = list(iter_samples(cfg))
    if not samples:
        log.info("report: no Species/Sample dirs found under %s", cfg.project)
        return

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    warned_missing_reportlab = False

    for species, sample in samples:
        qc_dir = cfg.sample_qc(species, sample)
        qc_imgs = list(qc_dir.glob("*.png")) if qc_dir.is_dir() else []
        qc_txts = ([f for f in qc_dir.iterdir()
                    if f.is_file() and f.suffix in (".txt", ".tsv", ".csv")]
                   if qc_dir.is_dir() else [])
        if not qc_imgs and not qc_txts:
            log.info("report: no QC files for %s/%s — skipping", species, sample)
            continue

        html_out = qc_dir / "qc_report.html"
        html_out.write_text(_build_qc_html(species, sample, qc_dir, now), encoding="utf-8")
        log.info("report: %s/%s/QC/qc_report.html (%d image(s), %d data file(s))",
                 species, sample, len(qc_imgs), len(qc_txts))

        pdf_out = qc_dir / "qc_report.pdf"
        try:
            _build_qc_pdf(species, sample, qc_dir, now, pdf_out)
            log.info("report: %s/%s/QC/qc_report.pdf", species, sample)
        except ImportError:
            if not warned_missing_reportlab:
                log.warning(_pdf_import_error_hint())
                warned_missing_reportlab = True
        except Exception as exc:
            log.warning("report: could not write qc_report.pdf for %s/%s: %s",
                       species, sample, exc)
