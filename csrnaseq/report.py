"""Step — HTML report generation.

Writes one self-contained qc_report.html per Species/Sample, saved directly
into that sample's own QC/ folder (Species/Sample/QC/qc_report.html) —
embedding that sample's PNGs (base64) and any QC data files.

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
]

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
    "ritrie": (
        "RITRIE/ritrie.tsv (per csRNA replicate)",
        "Higher RIT/RIE values generally indicate stronger enrichment of transcription initiation signal relative to exonic RNA, while lower values may indicate increased RNA degradation. Thresholds have not yet  
        been formally established.",
    ),
}


def _get_description(stem: str) -> tuple[str, str] | None:
    """Return (source, description) for a plot stem, or None."""
    for prefix, desc in _IMG_DESCRIPTIONS.items():
        if stem == prefix or stem.startswith(prefix):
            return desc
    return None


def _img_sort_key(name: str) -> tuple:
    stem = Path(name).stem
    for i, prefix in enumerate(_IMG_ORDER):
        if stem == prefix or stem.startswith(prefix):
            return (i, stem)
    return (len(_IMG_ORDER), stem)


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

def _build_qc_html(species: str, sample: str, qc_dir: Path, now: str) -> str:
    imgs = sorted(
        [f for f in qc_dir.glob("*.png") if f.is_file()],
        key=lambda f: _img_sort_key(f.name),
    )
    txts = sorted(
        f for f in qc_dir.iterdir()
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
        if items:
            img_section = (f'<div class="section-label">Images</div>'
                           f'<div class="img-grid">{items}</div>')

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


# ── Entry point ───────────────────────────────────────────────────────────────

def run_report(cfg) -> None:
    samples = list(iter_samples(cfg))
    if not samples:
        log.info("report: no Species/Sample dirs found under %s", cfg.project)
        return

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    for species, sample in samples:
        qc_dir = cfg.sample_qc(species, sample)
        qc_imgs = list(qc_dir.glob("*.png")) if qc_dir.is_dir() else []
        qc_txts = ([f for f in qc_dir.iterdir()
                    if f.is_file() and f.suffix in (".txt", ".tsv", ".csv")]
                   if qc_dir.is_dir() else [])
        if not qc_imgs and not qc_txts:
            log.info("report: no QC files for %s/%s — skipping", species, sample)
            continue

        out = qc_dir / "qc_report.html"
        out.write_text(_build_qc_html(species, sample, qc_dir, now), encoding="utf-8")
        log.info("report: %s/%s/QC/qc_report.html (%d image(s), %d data file(s))",
                 species, sample, len(qc_imgs), len(qc_txts))