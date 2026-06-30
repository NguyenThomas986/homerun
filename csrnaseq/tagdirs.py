"""Step 3 — HOMER tag directories.

Builds two kinds of tag directories per sample:
  • <sample>-combo     — all replicates merged together (existing behavior)
  • <sample>_r<N>       — one tag dir per individual replicate SAM

Both are built from the same aligned SAM files; the combo dir is unaffected
by the addition of per-replicate dirs.
"""
from __future__ import annotations
import re
from .utils import run, log, seq_type, done


def _make_tagdir(cmd_input_sams: str, tagdir, base: str, label: str, cfg) -> None:
    if done(tagdir):
        log.info("  skip (done): %s", tagdir.name)
        return
    st = seq_type(base)
    if st in ("csRNA", "sRNA"):
        cmd = (f"makeTagDirectory {tagdir} {cmd_input_sams} "
               f"-genome {cfg.genome} -checkGC -fragLength 150 -omitSN")
    elif st == "totalRNA":
        cmd = (f"makeTagDirectory {tagdir} {cmd_input_sams} "
               f"-genome {cfg.genome} -checkGC -fragLength 150 -read2")
    else:
        log.warning("tagdir: skipping untyped %s", base)
        return
    run(cmd, label=label)


def run_tagdirs(cfg) -> None:
    rep_sams = sorted(cfg.aligned.glob("*[_-]r1*.Aligned.out.sam"))
    if not rep_sams:
        log.info("tagdir: no *[_-]r1*.Aligned.out.sam in %s", cfg.aligned)
        return

    for sam in rep_sams:
        base = re.split(r"[_-]r1", sam.name)[0]

        # ── Combo: all replicates merged (existing behavior) ──────────────────
        combo_dir  = cfg.tagdirs / f"{base}-combo"
        combo_sams = f"{cfg.aligned}/{base}*.sam"
        _make_tagdir(combo_sams, combo_dir, base, f"tagdir {base}-combo", cfg)

        # ── Per-replicate: one tag dir per individual SAM ──────────────────────
        rep_sams_for_base = sorted(cfg.aligned.glob(f"{base}[_-]r*.Aligned.out.sam"))
        for rep_sam in rep_sams_for_base:
            m = re.search(r"[_-](r\d+)", rep_sam.name)
            if not m:
                log.warning("tagdir: could not parse replicate label from %s", rep_sam.name)
                continue
            rep_label = m.group(1)
            rep_dir   = cfg.tagdirs / f"{base}_{rep_label}"
            _make_tagdir(str(rep_sam), rep_dir, base,
                         f"tagdir {base}_{rep_label}", cfg)
