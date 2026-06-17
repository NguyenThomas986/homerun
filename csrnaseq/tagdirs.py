"""Step 3 — HOMER tag directories. Replicates merged into <sample>-combo dirs."""
from __future__ import annotations
import re
from .utils import run, log, seq_type, done

def run_tagdirs(cfg) -> None:
    r1_sams = sorted(cfg.aligned.glob("*[_-]r1*.Aligned.out.sam"))
    if not r1_sams:
        log.info("tagdir: no *[_-]r1*.Aligned.out.sam in %s", cfg.aligned)
        return
    for sam in r1_sams:
        base       = re.split(r"[_-]r1", sam.name)[0]
        tagdir     = cfg.tagdirs / f"{base}-combo"
        input_sams = f"{cfg.aligned}/{base}*.sam"
        if done(tagdir):
            log.info("  skip (done): %s", tagdir.name); continue
        st = seq_type(base)
        if st in ("csRNA", "sRNA"):
            cmd = (f"makeTagDirectory {tagdir} {input_sams} "
                   f"-genome {cfg.genome} -checkGC -fragLength 150 -omitSN")
        elif st == "totalRNA":
            cmd = (f"makeTagDirectory {tagdir} {input_sams} "
                   f"-genome {cfg.genome} -checkGC -fragLength 150 -read2")
        else:
            log.warning("tagdir: skipping untyped %s", base)
            continue
        run(cmd, label=f"tagdir {base}-combo")