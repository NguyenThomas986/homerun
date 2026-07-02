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