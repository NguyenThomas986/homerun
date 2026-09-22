"""Create a project-level Excel workbook containing machine-readable QC.

The workbook complements the HTML/PDF report. It keeps one row per replicate
so lane-merged tag directories and replicate-specific peak calls stay visible.
"""
from __future__ import annotations

from collections import Counter
import csv
import re

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill

from .utils import assay_of_leaf, iter_leaf_dirs, log, replicate_of_leaf


HEADER_FILL = PatternFill("solid", fgColor="1F4E78")
HEADER_FONT = Font(color="FFFFFF", bold=True)
TITLE_FILL = PatternFill("solid", fgColor="D9EAF7")


def _taginfo_value(lines, key):
    for line in lines:
        if line.startswith(key):
            return line.split("=")[-1].strip()
    return "NA"


def _tagdir_rows(cfg):
    source_counts = Counter()
    identities = set()

    for species, sample, leaf_name, _r1 in iter_leaf_dirs(cfg):
        key = (species, sample, leaf_name)
        source_counts[key] += 1
        identities.add(key)

    rows = []

    for species, sample, leaf_name in sorted(identities):
        tagdir = cfg.leaf_tagdir(species, sample, leaf_name)

        if not tagdir.is_dir():
            continue

        row = {
            "Species": species,
            "Sample": sample,
            "Assay": assay_of_leaf(leaf_name) or "NA",
            "Replicate": replicate_of_leaf(leaf_name) or "NA",
            "Library": leaf_name,
            "Source files": source_counts[(species, sample, leaf_name)],
            "TagDir": str(tagdir),
            "Total Tags": "NA",
            "Unique Positions": "NA",
            "Tags per BP": "NA",
            "Average Tags/Position": "NA",
            "Median Tags/Position": "NA",
            "Average Read Length": "NA",
            "Average Fragment GC": "NA",
        }

        info = tagdir / "tagInfo.txt"

        if info.exists():
            lines = info.read_text(errors="replace").splitlines()

            genome_line = next(
                (line for line in lines if line.startswith("genome=")),
                "",
            )

            parts = genome_line.split("\t")

            row.update({
                "Total Tags": parts[2].strip() if len(parts) > 2 else "NA",
                "Unique Positions": parts[1].strip() if len(parts) > 1 else "NA",
                "Tags per BP": _taginfo_value(lines, "tagsPerBP"),
                "Average Tags/Position": _taginfo_value(
                    lines,
                    "averageTagsPerPosition",
                ),
                "Median Tags/Position": _taginfo_value(
                    lines,
                    "medianTagsPerPosition",
                ),
                "Average Read Length": _taginfo_value(
                    lines,
                    "averageTagLength",
                ),
                "Average Fragment GC": _taginfo_value(
                    lines,
                    "averageFragmentGCcontent",
                ),
            })

        rows.append(row)

    return rows


def _grab(text, pattern, cast=str, default="NA"):
    match = re.search(pattern, text)

    if not match:
        return default

    value = match.group(1).strip().rstrip("%")

    try:
        return cast(value)
    except (TypeError, ValueError):
        return default


def _peak_rows(cfg):
    """Collect TSR / findcsRNATSS statistics from Species/TSS/*.stats.txt."""
    rows = []

    for species_dir in sorted(
        path for path in cfg.project.iterdir()
        if path.is_dir()
    ):
        tss_dir = species_dir / "TSS"

        if not tss_dir.is_dir():
            continue

        for stats in sorted(tss_dir.glob("*.stats.txt")):
            text = stats.read_text(errors="replace")

            library = stats.name.removesuffix(".stats.txt")
            sample = library.split("_", 1)[0]

            rows.append({
                "Species": species_dir.name,
                "Sample": sample,
                "Peak call": library,
                "Total csRNA reads": _grab(
                    text,
                    r"Total csRNA reads:\s+([\d.]+)",
                    float,
                ),
                "Total input reads": _grab(
                    text,
                    r"Total input reads:\s+([\d.]+)",
                    float,
                ),
                "Putative TSS": _grab(
                    text,
                    r"total putative TSS clusters\s+(\d+)",
                    int,
                ),
                "Valid TSS": _grab(
                    text,
                    r"Valid TSS clusters\s+(\d+)",
                    int,
                ),
                "Distal %": _grab(
                    text,
                    r"Fraction Promoter-Distal.*?:\s+([\d.]+%)",
                    float,
                ),
                "Bidirectional %": _grab(
                    text,
                    r"Fraction of bidirectional.*?:\s+([\d.]+%)",
                    float,
                ),
                "Stable %": _grab(
                    text,
                    r"Fraction of stable.*?:\s+([\d.]+%)",
                    float,
                ),
                "Log2 vs Input": _grab(
                    text,
                    r"log2 fold vs\. input:\s+([\d.\-]+)",
                    float,
                ),
                "Log2 vs RNA": _grab(
                    text,
                    r"log2 fold vs\. rna:\s+([\d.\-]+)",
                    float,
                ),
                "Stats file": str(stats),
            })

    return rows


def _parse_lengths(path):
    """Return total, retained, adapter/dimer and weighted retained length."""
    try:
        with path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            next(reader, None)

            values = []

            for row in reader:
                if len(row) < 2:
                    continue

                try:
                    values.append(
                        (float(row[0]), float(row[1]))
                    )
                except ValueError:
                    continue

    except OSError:
        return None

    if not values:
        return None

    total = sum(
        count
        for _length, count in values
    )

    if total <= 0:
        return None

    adapters = sum(
        count
        for length, count in values
        if length == 0
    )

    retained = total - adapters

    average = (
        sum(
            length * count
            for length, count in values
            if length != 0
        ) / retained
        if retained > 0
        else "NA"
    )

    return {
        "Input Reads": int(total),
        "Retained Reads": int(retained),
        "Retained %": round(
            100 * retained / total,
            2,
        ),
        "Adapter/Dimer Reads": int(adapters),
        "Adapter/Dimer %": round(
            100 * adapters / total,
            2,
        ),
        "Average Retained Length": (
            round(average, 2)
            if average != "NA"
            else "NA"
        ),
    }


def _parse_star(path):
    text = path.read_text(errors="replace")

    patterns = {
        "Input Reads": (
            r"Number of input reads \|\s+(\d+)",
            int,
        ),
        "Uniquely Mapped %": (
            r"Uniquely mapped reads % \|\s+([\d.]+)%",
            float,
        ),
        "Multi-Mapped %": (
            r"% of reads mapped to multiple loci \|\s+([\d.]+)%",
            float,
        ),
        "Too-Many-Loci %": (
            r"% of reads mapped to too many loci \|\s+([\d.]+)%",
            float,
        ),
        "Unmapped (mismatch) %": (
            r"% of reads unmapped: too many mismatches \|\s+([\d.]+)%",
            float,
        ),
        "Unmapped (short) %": (
            r"% of reads unmapped: too short \|\s+([\d.]+)%",
            float,
        ),
        "Unmapped (other) %": (
            r"% of reads unmapped: other \|\s+([\d.]+)%",
            float,
        ),
    }

    result = {}

    for name, (pattern, cast) in patterns.items():
        match = re.search(pattern, text)

        if match:
            result[name] = cast(match.group(1))

    unmapped = [
        value
        for key, value in result.items()
        if key.startswith("Unmapped (")
    ]

    if unmapped:
        result["Unmapped %"] = round(
            sum(unmapped),
            2,
        )

    uniquely = result.get("Uniquely Mapped %")
    multi = result.get("Multi-Mapped %")
    too_many = result.get("Too-Many-Loci %")

    if all(
        value is not None
        for value in (uniquely, multi, too_many)
    ):
        result["Overall Aligned %"] = round(
            uniquely + multi + too_many,
            2,
        )

    return result


def _parse_hisat2(path):
    text = path.read_text(errors="replace")

    result = {}

    patterns = {
        "Input Reads": (
            r"(\d+)\s+reads?;\s+of these",
            int,
        ),
        "Unmapped %": (
            r"([\d.]+)\s*%\)\s+aligned 0 times",
            float,
        ),
        "Uniquely Mapped %": (
            r"([\d.]+)\s*%\)\s+aligned exactly 1 time",
            float,
        ),
        "Multi-Mapped %": (
            r"([\d.]+)\s*%\)\s+aligned >1 times",
            float,
        ),
        "Overall Aligned %": (
            r"([\d.]+)\s*%\s+overall alignment rate",
            float,
        ),
    }

    for name, (pattern, cast) in patterns.items():
        match = re.search(pattern, text)

        if match:
            result[name] = cast(match.group(1))

    return result


def _trim_alignment_rows(cfg):
    trim_rows = []
    alignment_rows = []

    for species, sample, leaf_name, r1 in iter_leaf_dirs(cfg):
        assay = assay_of_leaf(leaf_name) or "NA"
        replicate = replicate_of_leaf(leaf_name) or "NA"

        common = {
            "Species": species,
            "Sample": sample,
            "Assay": assay,
            "Replicate": replicate,
            "Library": leaf_name,
            "FASTQ": r1.name,
        }

        # Trimming statistics
        lengths = (
            cfg.trimmed_dir(species, sample)
            / f"{r1.name}.lengths"
        )

        if lengths.exists():
            parsed = _parse_lengths(lengths)

            if parsed:
                trim_rows.append({
                    **common,
                    "Tool": "homerTools",
                    **parsed,
                    "Source log": str(lengths),
                })

        # Alignment statistics
        prefix = r1.name.split("_R1")[0]

        aligned_dir = cfg.aligned_dir(
            species,
            sample,
        )

        star_log = (
            aligned_dir
            / f"{prefix}.Log.final.out"
        )

        hisat2_log = (
            aligned_dir
            / f"{prefix}_mappingstats.txt"
        )

        if star_log.exists():
            parsed = _parse_star(star_log)
            tool = "STAR"
            source = star_log

        elif hisat2_log.exists():
            parsed = _parse_hisat2(hisat2_log)
            tool = "HISAT2"
            source = hisat2_log

        else:
            parsed = None

        if parsed:
            alignment_rows.append({
                **common,
                "Aligner": tool,
                "Input Reads": parsed.get(
                    "Input Reads",
                    "NA",
                ),
                "Uniquely Mapped %": parsed.get(
                    "Uniquely Mapped %",
                    "NA",
                ),
                "Multi-Mapped %": parsed.get(
                    "Multi-Mapped %",
                    "NA",
                ),
                "Unmapped %": parsed.get(
                    "Unmapped %",
                    "NA",
                ),
                "Overall Aligned %": parsed.get(
                    "Overall Aligned %",
                    "NA",
                ),
                "Source log": str(source),
            })

    return trim_rows, alignment_rows


def _write_table_sheet(workbook, name, rows):
    """Write a normal worksheet with headers and AutoFilter.

    This intentionally does not create an openpyxl Table object. Excel was
    repairing/removing those table definitions on open, so regular worksheet
    filters are used instead.
    """
    ws = workbook.create_sheet(name)

    if not rows:
        ws.append([
            f"No {name.lower()} data were found."
        ])

        ws["A1"].fill = TITLE_FILL
        ws["A1"].font = Font(bold=True)

        return ws

    headers = list(rows[0])
    ws.append(headers)

    for row in rows:
        ws.append([
            row.get(header, "NA")
            for header in headers
        ])

    # Header formatting
    for cell in ws[1]:
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(
            horizontal="center",
            vertical="center",
            wrap_text=True,
        )

    # Freeze the header row
    ws.freeze_panes = "A2"

    # Native worksheet filter without an Excel Table object
    ws.auto_filter.ref = ws.dimensions

    # Make columns readable without letting paths create huge widths
    for column_cells in ws.columns:
        values = [
            str(cell.value or "")
            for cell in column_cells
        ]

        width = min(
            max(
                max(map(len, values)) + 2,
                11,
            ),
            45,
        )

        ws.column_dimensions[
            column_cells[0].column_letter
        ].width = width

    # Numeric percentage columns are already stored as values like 97.2,
    # so display them as ordinary numbers rather than Excel's 0.972 format.
    for cell in ws[1]:
        if str(cell.value).endswith("%"):
            for data_cell in ws.iter_cols(
                min_col=cell.column,
                max_col=cell.column,
                min_row=2,
            ):
                data_cell[0].number_format = "0.00"

    return ws


def _write_overview(
    workbook,
    cfg,
    tag_rows,
    peak_rows,
    trim_rows,
    alignment_rows,
):
    ws = workbook.create_sheet("Overview")

    ws.append([
        "HOMERun QC Summary"
    ])

    ws["A1"].fill = HEADER_FILL
    ws["A1"].font = Font(
        color="FFFFFF",
        bold=True,
        size=14,
    )

    ws.merge_cells("A1:B1")

    ws.append([
        "Project",
        str(cfg.project),
    ])

    ws.append([
        "Replicate TagDirs",
        len(tag_rows),
    ])

    ws.append([
        "TSR / Peak calls",
        len(peak_rows),
    ])

    ws.append([
        "Trimming records",
        len(trim_rows),
    ])

    ws.append([
        "Alignment records",
        len(alignment_rows),
    ])

    ws.append([
        "Species",
        len({
            row["Species"]
            for row in tag_rows
        }),
    ])

    ws.append([
        "Samples",
        len({
            (
                row["Species"],
                row["Sample"],
            )
            for row in tag_rows
        }),
    ])

    ws.column_dimensions["A"].width = 24
    ws.column_dimensions["B"].width = 70

    for row in range(
        2,
        ws.max_row + 1,
    ):
        ws.cell(
            row,
            1,
        ).font = Font(
            bold=True
        )


def write_qc_workbook(cfg):
    """Write ``QC_summary.xlsx`` at the project root and return its path."""
    tag_rows = _tagdir_rows(cfg)
    peak_rows = _peak_rows(cfg)

    trim_rows, alignment_rows = (
        _trim_alignment_rows(cfg)
    )

    log.info(
        "QC Excel: %d trimming, %d alignment, %d TagDir, %d TSR stat row(s)",
        len(trim_rows),
        len(alignment_rows),
        len(tag_rows),
        len(peak_rows),
    )

    workbook = Workbook()

    # Remove the blank default sheet.
    workbook.remove(
        workbook.active
    )

    _write_overview(
        workbook,
        cfg,
        tag_rows,
        peak_rows,
        trim_rows,
        alignment_rows,
    )

    _write_table_sheet(
        workbook,
        "Trimming",
        trim_rows,
    )

    _write_table_sheet(
        workbook,
        "Alignment",
        alignment_rows,
    )

    _write_table_sheet(
        workbook,
        "Tag Directories",
        tag_rows,
    )

    _write_table_sheet(
        workbook,
        "TSR Stats",
        peak_rows,
    )

    output = (
        cfg.project
        / "QC_summary.xlsx"
    )

    workbook.save(output)

    log.info(
        "QC: wrote Excel workbook %s",
        output,
    )

    return output
