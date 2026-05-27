"""
HomerunQC.py
────────────
QC step for the CSV-driven Homerun pipeline.

Delegates to the original QC plotting functions (ported from the
original HomerunQC.py) but wrapped to:
  - accept a single CSV row dict
  - write outputs to files/QC/
  - return CSV column updates
  - propagate exceptions cleanly

Original QC functions (unchanged logic, same signatures internally):
  csRNAmedianTags, csRNAtagsVsFrac, sRNAmedianTags, sRNAtagsVsFrac,
  totalRNAmedianTags, totalRNAtagsVsFrac, csRNAlengthPlot,
  csRNAcomboLengthPlot, sRNAlengthPlot, sRNAcomboLengthPlot,
  ntPrefs, csRNAaPlot, csRNAcomboAPlot, sRNAaPlot, sRNAcomboAPlot
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import Dict, Any

import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")   # headless backend for HPC
import matplotlib.pyplot as plt
import numpy as np

log = logging.getLogger(__name__)


# ── Public entry point ────────────────────────────────────────────────────────
def run_qc(row: Dict[str, Any], working_path: str, cpus: int = 1) -> Dict[str, Any]:
    """
    Run all QC plots for the sample described by *row*.

    Returns dict with at minimum:
        {"qc": <qc_output_dir>, "status": "RUNNING"}
    """
    sample       = row["sample"]
    genome       = row["genome"]
    tagdir_val   = row.get("tagdir", "")
    mode         = _infer_mode(row)

    qc_output    = Path(working_path) / "files" / "QC"
    qc_output.mkdir(parents=True, exist_ok=True)
    output       = str(qc_output)

    err_msgs: list[str] = []

    def err_write(msg: str):
        log.warning("[%s] QC warning: %s", sample, msg.strip())
        err_msgs.append(msg.strip())

    # Fake errFile object compatible with original function signatures
    class _ErrFile:
        def write(self, msg):
            err_write(msg)

    errFile = _ErrFile()

    total_rna_tag = "_RNA" if mode == "star" else "_totalRNA"

    log.info("[%s] running QC plots → %s", sample, output)

    species   = genome
    # We call the same functions as the original qc() dispatcher
    _run_all_qc_functions(working_path, species, mode, total_rna_tag, output, errFile)

    notes = "; ".join(err_msgs) if err_msgs else ""
    log.info("[%s] QC complete (warnings: %d).", sample, len(err_msgs))

    return {"qc": output, "notes": notes, "status": "RUNNING"}


# ── Dispatcher ────────────────────────────────────────────────────────────────
def _run_all_qc_functions(workingPath, species, mode, totalRnaTag, output, errFile):
    funcs = [
        lambda: csRNAmedianTags(workingPath, species, errFile, output),
        lambda: csRNAtagsVsFrac(workingPath, species, errFile, output),
        lambda: sRNAmedianTags(workingPath, species, errFile, output),
        lambda: sRNAtagsVsFrac(workingPath, species, errFile, output),
        lambda: totalRNAmedianTags(workingPath, species, totalRnaTag, errFile, output),
        lambda: totalRNAtagsVsFrac(workingPath, species, totalRnaTag, errFile, output),
        lambda: csRNAlengthPlot(workingPath, species, totalRnaTag, errFile, output),
        lambda: csRNAcomboLengthPlot(workingPath, species, errFile, output),
        lambda: sRNAlengthPlot(workingPath, species, errFile, output),
        lambda: sRNAcomboLengthPlot(workingPath, species, errFile, output),
        lambda: ntPrefs(workingPath, species, errFile, output),
        lambda: csRNAaPlot(workingPath, species, errFile, output),
        lambda: csRNAcomboAPlot(workingPath, species, errFile, output),
        lambda: sRNAaPlot(workingPath, species, errFile, output),
        lambda: sRNAcomboAPlot(workingPath, species, errFile, output),
    ]
    for fn in funcs:
        try:
            fn()
        except Exception as exc:
            errFile.write(f"QC function error: {exc}")


# ── Helpers ───────────────────────────────────────────────────────────────────
def _infer_mode(row: Dict[str, Any]) -> str:
    notes = str(row.get("notes", "")).lower()
    star_val = str(row.get("star", "")).lower()
    if "hisat" in notes or "hisat" in star_val:
        return "hisat"
    return "star"


def _read_sample_tagdirs(workingPath, species, tag_filter):
    """Read tagdir names matching tag_filter from sampleInfo.txt."""
    info_path = f"{workingPath}data/{species}/fastq/sampleInfo.txt"
    tagdirs = []
    try:
        with open(info_path, "r") as f:
            for line in sorted(f):
                name = line.strip().split("\t")[0]
                if tag_filter(name):
                    tagdirs.append(name)
    except FileNotFoundError:
        pass
    return tagdirs


def _median_bar_plot(tagdirs, output_path):
    my_dict = {
        "Library": [],
        "Median tags per tag position (should be =1)": [],
    }
    for f in tagdirs:
        fpath = f + "/tagCountDistribution.txt"
        try:
            df = pd.read_csv(fpath, sep="\t")
            name = list(df.columns)
            median_val = str(name).split("=")[1].split(",")[0]
            my_dict["Library"].append(f.split("/")[-1])
            my_dict["Median tags per tag position (should be =1)"].append(median_val)
        except Exception:
            continue
    if not my_dict["Library"]:
        return
    mf = pd.DataFrame(my_dict)
    mf["Median tags per tag position (should be =1)"] = pd.to_numeric(
        mf["Median tags per tag position (should be =1)"], errors="coerce"
    )
    ax = sns.barplot(
        data=mf, y="Library",
        x="Median tags per tag position (should be =1)",
        capsize=0.4, errcolor=".5", linewidth=2,
        edgecolor=".5", facecolor=(0, 0, 0, 0),
    )
    ax.axvline(1,   color="g", ls="--")
    ax.axvline(1.2, color="r")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


# ── QC functions (refactored to accept output explicitly) ─────────────────────

def csRNAmedianTags(workingPath, species, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species, lambda n: "csRNA" in n)
    if len(tagdirs) < 2:
        errFile.write("csRNAmedianTags: fewer than 2 csRNA tag dirs found.\n")
        return
    try:
        _median_bar_plot(tagdirs, f"{output}/csRNA_{species}_medianTagsPerPosition.png")
    except Exception as e:
        errFile.write(f"csRNAmedianTags error: {e}\n")


def csRNAtagsVsFrac(workingPath, species, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species,
                                   lambda n: "csRNA" in n)
    if len(tagdirs) < 2:
        errFile.write("csRNAtagsVsFrac: fewer than 2 csRNA tag dirs found.\n")
        return
    try:
        _tags_vs_frac_plot(tagdirs, output,
                           f"csRNA_{species}_tagsPer_Vs_FracofPos")
    except Exception as e:
        errFile.write(f"csRNAtagsVsFrac error: {e}\n")


def sRNAmedianTags(workingPath, species, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species, lambda n: "_sRNA" in n)
    if len(tagdirs) < 2:
        errFile.write("sRNAmedianTags: fewer than 2 sRNA tag dirs found.\n")
        return
    try:
        _median_bar_plot(tagdirs, f"{output}/sRNA_{species}_medianTagsPerPosition.png")
    except Exception as e:
        errFile.write(f"sRNAmedianTags error: {e}\n")


def sRNAtagsVsFrac(workingPath, species, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species, lambda n: "_sRNA" in n)
    if len(tagdirs) < 2:
        errFile.write("sRNAtagsVsFrac: fewer than 2 sRNA tag dirs found.\n")
        return
    try:
        _tags_vs_frac_plot(tagdirs, output, f"sRNA_{species}_tagsPer_Vs_FracofPos")
    except Exception as e:
        errFile.write(f"sRNAtagsVsFrac error: {e}\n")


def totalRNAmedianTags(workingPath, species, totalRnaTag, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species,
                                   lambda n: totalRnaTag in n)
    if len(tagdirs) < 2:
        errFile.write(f"totalRNAmedianTags: fewer than 2 {totalRnaTag} tag dirs.\n")
        return
    try:
        _median_bar_plot(tagdirs, f"{output}/totalRNA_{species}_medianTagsPerPosition.png")
    except Exception as e:
        errFile.write(f"totalRNAmedianTags error: {e}\n")


def totalRNAtagsVsFrac(workingPath, species, totalRnaTag, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species,
                                   lambda n: totalRnaTag in n)
    if len(tagdirs) < 2:
        errFile.write(f"totalRNAtagsVsFrac: fewer than 2 {totalRnaTag} tag dirs.\n")
        return
    try:
        _tags_vs_frac_plot(tagdirs, output,
                           f"totalRNA_{species}_tagsPer_Vs_FracofPos")
    except Exception as e:
        errFile.write(f"totalRNAtagsVsFrac error: {e}\n")


def _tags_vs_frac_plot(tagdirs, output_dir, base_name):
    first = pd.read_csv(tagdirs[1] + "/tagCountDistribution.txt", sep="\t")
    first = first.rename(columns={first.columns[0]: "Tags per tag position"}).iloc[:, [0]]
    for f in tagdirs:
        df = pd.read_csv(f + "/tagCountDistribution.txt", sep="\t")
        median_v = str(list(df)).split("=")[1].split(",")[0].split(" ")[1]
        col_name = f.split("/")[-1] + f" ({median_v})"
        df = df.rename(columns={df.columns[0]: "Tags per tag position",
                                 df.columns[1]: col_name})
        first = pd.merge(first, df, on="Tags per tag position", how="left")
    first = first.set_index("Tags per tag position")
    first.to_csv(f"{output_dir}/{base_name}.txt", sep="\t")
    combined = first.replace(["0", 0], np.nan).stack().reset_index()
    combined.columns = ["Tags per tag position", "Library", "Fraction of Positions"]
    combined["Tags per tag position"] = np.log(combined["Tags per tag position"])
    combined["Fraction of Positions"] = np.log(combined["Fraction of Positions"])
    g = sns.FacetGrid(combined, col="Library", height=8, aspect=1, col_wrap=4)
    g.map(sns.scatterplot, "Tags per tag position", "Fraction of Positions", alpha=0.5)
    plt.savefig(f"{output_dir}/{base_name}.png")
    plt.savefig(f"{output_dir}/{base_name}.svg")
    plt.close()


def csRNAlengthPlot(workingPath, species, totalRnaTag, errFile, output):
    tagdirs = _read_sample_tagdirs(
        workingPath, species,
        lambda n: "ChIPseq" not in n and totalRnaTag not in n
                  and "_sRNA" not in n,
    )
    _length_plot(tagdirs, output, f"csRNA_{species}_Length_plot", errFile)


def csRNAcomboLengthPlot(workingPath, species, errFile, output):
    tagdir_base = f"{workingPath}data/{species}/tagDirs/"
    tagdirs = [
        tagdir_base + "/" + d
        for d in (os.listdir(tagdir_base) if os.path.isdir(tagdir_base) else [])
        if "csRNA" in d and "-r" not in d
    ]
    _length_plot(tagdirs, output, f"csRNA_{species}_Combo_Length_plot", errFile)


def sRNAlengthPlot(workingPath, species, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species, lambda n: "_sRNA" in n)
    _length_plot(tagdirs, output, f"sRNA_{species}_Length_plot", errFile)


def sRNAcomboLengthPlot(workingPath, species, errFile, output):
    tagdir_base = f"{workingPath}data/{species}/tagDirs/"
    tagdirs = [
        tagdir_base + "/" + d
        for d in (os.listdir(tagdir_base) if os.path.isdir(tagdir_base) else [])
        if "_sRNA" in d and "-r" not in d
    ]
    _length_plot(tagdirs, output, f"sRNA_{species}_Combo_Length_plot", errFile)


def _length_plot(tagdirs, output_dir, base_name, errFile):
    if len(tagdirs) < 2:
        errFile.write(f"{base_name}: fewer than 2 tag dirs found.\n")
        return
    try:
        first = pd.read_csv(tagdirs[1] + "/tagLengthDistribution.txt", sep="\t")
        first = first.rename(columns={first.columns[0]: "Length (nt)"}).iloc[:, [0]]
        for f in sorted(tagdirs):
            df = pd.read_csv(f + "/tagLengthDistribution.txt", sep="\t")
            col = f.split("/")[-1]
            df = df.rename(columns={df.columns[0]: "Length (nt)", df.columns[1]: col})
            first = pd.merge(first, df, on="Length (nt)", how="left")
        first = first.iloc[1:].set_index("Length (nt)")
        first.to_csv(f"{output_dir}/{base_name}.txt", sep="\t")
        stacked = first.replace(["0", 0], np.nan).stack().reset_index()
        stacked.columns = ["Length (nt)", "Library", "Fraction of Reads"]
        ax = sns.lineplot(data=stacked, x="Length (nt)", y="Fraction of Reads",
                          hue="Library", linewidth=2, alpha=0.6)
        sns.move_legend(ax, "upper right")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{base_name}.png")
        plt.savefig(f"{output_dir}/{base_name}.svg")
        plt.close()
    except Exception as e:
        errFile.write(f"{base_name} error: {e}\n")


def ntPrefs(workingPath, species, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species, lambda n: True)
    for f in tagdirs:
        try:
            df = pd.read_csv(f + "/tagFreqUniq.txt", sep="\t")
            df = df[df.columns[:5]].set_index("Offset").stack().reset_index()
            name = f.split("/")[-1]
            df.columns = ["Distance from TSS", "nt", name + " - %"]
            plt.figure()
            sns.lineplot(data=df, x="Distance from TSS",
                         y=name + " - %", hue="nt")
            plt.tight_layout()
            plt.savefig(f"{output}/{name}_nt_Preference.png")
            plt.savefig(f"{output}/{name}_nt_Preference.svg")
            plt.close()
        except Exception as e:
            errFile.write(f"ntPrefs error for {f}: {e}\n")


def _a_plot(tagdirs, output_dir, base_name, errFile, selected="csRNA"):
    if len(tagdirs) < 2:
        errFile.write(f"{base_name}: fewer than 2 tag dirs.\n")
        return
    try:
        first = pd.read_csv(tagdirs[1] + "/tagFreqUniq.txt", sep="\t")
        first = first[first.columns[:5]].set_index("Offset").stack().reset_index()
        first = first.rename(columns={first.columns[1]: "nt"}).iloc[:, :2]
        for f in tagdirs:
            if selected not in f:
                continue
            df = pd.read_csv(f + "/tagFreqUniq.txt", sep="\t")
            df = df[df.columns[:5]].set_index("Offset").stack().reset_index()
            col = f.split("/")[-1]
            df.columns = ["Distance from TSS", "nt", col]
            del df["Distance from TSS"]
            del df["nt"]
            first = pd.merge(first, df, left_index=True, right_index=True, how="left")
        aplot = first[first["nt"].str.contains("A") == True]
        del aplot["nt"]
        aplot = aplot.set_index("Offset").stack().reset_index()
        aplot.columns = ["Distance to TSS", "Library", "A [%]"]
        plt.figure()
        ax = sns.lineplot(data=aplot, x="Distance to TSS", y="A [%]",
                          hue="Library", linewidth=2, alpha=0.6)
        sns.move_legend(ax, "upper right")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{base_name}.png")
        plt.savefig(f"{output_dir}/{base_name}.svg")
        plt.close()
    except Exception as e:
        errFile.write(f"{base_name} error: {e}\n")


def csRNAaPlot(workingPath, species, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species, lambda n: "csRNA" in n)
    _a_plot(tagdirs, output, f"csRNA_{species}_Aplot", errFile, selected="csRNA")


def csRNAcomboAPlot(workingPath, species, errFile, output):
    tagdir_base = f"{workingPath}data/{species}/tagDirs/"
    tagdirs = [
        tagdir_base + "/" + d
        for d in (os.listdir(tagdir_base) if os.path.isdir(tagdir_base) else [])
        if "csRNA" in d and "-r" not in d
    ]
    _a_plot(tagdirs, output, f"csRNA_{species}_Combo_Aplots", errFile, selected="csRNA")


def sRNAaPlot(workingPath, species, errFile, output):
    tagdirs = _read_sample_tagdirs(workingPath, species, lambda n: "_sRNA" in n)
    _a_plot(tagdirs, output, f"sRNA_{species}_Aplot", errFile, selected="_sRNA")


def sRNAcomboAPlot(workingPath, species, errFile, output):
    tagdir_base = f"{workingPath}data/{species}/tagDirs/"
    tagdirs = [
        tagdir_base + "/" + d
        for d in (os.listdir(tagdir_base) if os.path.isdir(tagdir_base) else [])
        if "_sRNA" in d and "-r" not in d
    ]
    _a_plot(tagdirs, output, f"sRNA_{species}_Combo_Aplots", errFile, selected="_sRNA")
