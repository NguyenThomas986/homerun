"""Step — TSR characterization: stability (needs total RNA) + genomic location.

Total RNA is OPTIONAL. Without it the stable/unstable split is skipped automatically;
the distal/proximal breakdown still runs. Outputs PNGs into QC/.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt   # noqa: E402
import pandas as pd                     # noqa: E402

from .utils import log                  # noqa: E402

STABLE_C, UNSTABLE_C = "#2c7fb8", "#de2d26"
DISTAL_C, PROX_C     = "#756bb1", "#31a354"


def _read_homer_tss(path):
    df = pd.read_csv(path, sep="\t")
    return df.rename(columns={df.columns[0]: df.columns[0].lstrip("#").strip()})


def _find_col(df, override, *needles, exclude=()):
    if override and override in df.columns:
        return override
    for c in df.columns:
        lc = c.lower()
        if any(n in lc for n in needles) and not any(x in lc for x in exclude):
            return c
    return None


def _stability(cfg, df):
    col = _find_col(df, cfg.stability_col, "stabl")
    if col:
        s = df[col].astype(str).str.lower()
        out = pd.Series(pd.NA, index=df.index, dtype="object")
        out[s.str.contains("unstab")] = "unstable"
        out[s.str.contains("stab") & ~s.str.contains("unstab")] = "stable"
        if out.notna().any():
            return out, f"column '{col}'"
    for c in df.columns:                                       # detect by values
        v = df[c].astype(str).str.lower()
        if v.isin(["stable", "unstable"]).mean() > 0.5:
            return v.where(v.isin(["stable", "unstable"])), f"column '{c}' (values)"
    rcol = _find_col(df, cfg.rna_col, "rna", exclude=("csrna",))   # numeric fallback
    if rcol:
        rna = pd.to_numeric(df[rcol], errors="coerce")
        out = pd.Series("unstable", index=df.index)
        out[rna >= cfg.rna_stable_threshold] = "stable"
        return out, f"RNA '{rcol}' >= {cfg.rna_stable_threshold}"
    return None, "none (no total RNA?)"


def _location(cfg, df):
    col = _find_col(df, cfg.distal_col, "distal", "proximal", "promoter")
    if not col:
        for c in df.columns:
            v = df[c].astype(str).str.lower()
            if v.str.contains("distal|proximal|promoter").mean() > 0.3:
                col = c; break
    if not col:
        return None, "none"
    v = df[col].astype(str).str.lower()
    out = pd.Series(pd.NA, index=df.index, dtype="object")
    out[v.str.contains("distal")] = "distal"
    out[v.str.contains("proximal")] = "proximal"
    if out.isna().all() and v.isin(["true", "1", "yes"]).any():   # boolean "is distal"
        out = pd.Series("proximal", index=df.index)
        out[v.isin(["true", "1", "yes"])] = "distal"
    return out, f"column '{col}'"


def run_stability(cfg) -> None:
    tss_files = sorted(cfg.tss.glob("*.tss.txt"))
    if not tss_files:
        log.info("stability: no *.tss.txt in %s — run findcsRNATSS first.", cfg.tss)
        return

    recs = []
    for f in tss_files:
        df = _read_homer_tss(f)
        sample = f.name.replace(".tss.txt", "")
        stab, how_s = _stability(cfg, df)
        loc,  how_l = _location(cfg, df)
        r = pd.DataFrame(index=df.index)
        r["sample"]    = sample
        r["stability"] = stab if stab is not None else pd.NA
        r["location"]  = loc  if loc  is not None else pd.NA
        recs.append(r)
        log.info("%s: stability=%s; location=%s", sample, how_s, how_l)

    ALL = pd.concat(recs, ignore_index=True)
    has_stability = ALL["stability"].notna().any()
    has_location  = ALL["location"].notna().any()
    log.info("%d TSRs | stability=%s | location=%s",
             len(ALL), has_stability, has_location)

    _plot_stacked_bar(cfg, ALL, has_stability, has_location)
    _plot_pie(cfg, ALL, has_stability, has_location)


def _plot_stacked_bar(cfg, ALL, has_stability, has_location):
    if has_stability and has_location:
        ct = (ALL.dropna(subset=["stability", "location"])
                 .groupby(["stability", "location"]).size().unstack(fill_value=0)
                 .reindex(index=["stable", "unstable"]).fillna(0))
        ct = ct.reindex(columns=[c for c in ["distal", "proximal"] if c in ct.columns])
        ax = ct.plot(kind="bar", stacked=True, figsize=(6, 5),
                     color={"distal": DISTAL_C, "proximal": PROX_C})
        ax.set_ylabel("Number of TSRs"); ax.set_xlabel(""); plt.xticks(rotation=0)
        ax.set_title("Stable vs unstable TSRs, split by genomic location")
        plt.tight_layout()
        plt.savefig(cfg.qc / "stability_by_location_stacked_bar.png", dpi=150); plt.close()
        log.info("plot: stability_by_location_stacked_bar.png")
    elif has_location:
        ct = (ALL.dropna(subset=["location"]).groupby(["sample", "location"]).size()
                 .unstack(fill_value=0))
        ax = ct.plot(kind="bar", stacked=True, figsize=(max(6, 1.6 * len(ct)), 5),
                     color={"distal": DISTAL_C, "proximal": PROX_C})
        ax.set_ylabel("Number of TSRs"); ax.set_xlabel("")
        ax.set_title("TSRs by genomic location  (no total RNA → stability N/A)")
        plt.xticks(rotation=30, ha="right"); plt.tight_layout()
        plt.savefig(cfg.qc / "location_stacked_bar.png", dpi=150); plt.close()
        log.info("plot: location_stacked_bar.png")
    else:
        log.info("stacked bar: no location annotation found.")


def _plot_pie(cfg, ALL, has_stability, has_location):
    if has_stability:
        tot = ALL["stability"].value_counts().reindex(["stable", "unstable"]).fillna(0)
        cols, title = [STABLE_C, UNSTABLE_C], "Stable vs unstable (pooled)"
    elif has_location:
        tot = ALL["location"].value_counts().reindex(["distal", "proximal"]).fillna(0)
        cols, title = [DISTAL_C, PROX_C], "Distal vs proximal (pooled; no total RNA)"
    else:
        return
    if tot.sum() <= 0:
        return
    plt.figure(figsize=(5, 5))
    plt.pie(tot.values, labels=[f"{i} ({int(v)})" for i, v in tot.items()],
            autopct="%1.1f%%", startangle=90, colors=cols)
    plt.title(title); plt.tight_layout()
    plt.savefig(cfg.qc / "tsr_pie.png", dpi=150); plt.close()
    log.info("plot: tsr_pie.png")
