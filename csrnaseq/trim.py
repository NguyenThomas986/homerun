"""Step 1 — Trim. homerTools for everything now: csRNA/sRNA single-end via
the plain trim command, totalRNA paired-end via homerTools trim's own -pe
mate-pairing. Previously totalRNA went through skewer instead; homerTools
-pe replaces that so trimming is done by one tool everywhere, and so QC
(qc.py) can read every replicate's stats from the same .lengths format
rather than a separate skewer log parser.

R1's R2 mate is located via utils.find_r2_for_r1() — matched by PARSED
identity (species/sample/leaf_name), not by swapping "_R1" for "_R2" in the
filename. That substitution assumes both mates share the same basename
apart from the R1/R2 tag, which breaks for ENCODE-style downloads where
each mate has its own independent accession (e.g. an R1 file's accession
number is completely unrelated to its R2 mate's) — swapping _R1->_R2 on
such a name silently builds a path to a file that was never going to
exist, and homerTools -pe then fails to open it.

NOTE: -pe's mate-synchronized read removal (so R1/R2 stay paired 1:1 after
-min drops short reads) is only confirmed from HOMER's own usage text, not
independently verified here against a real paired-end run — after your
first real totalRNA run, it's worth a quick sanity check that the R1/R2
'.trimmed' outputs have the same read count (e.g. `zcat X_R1... | wc -l`
vs the R2 equivalent, both divisible by 4) before trusting it in production.

Operates per R1 file so it can run inside a SLURM array (one task per sample).
homerTools writes outputs next to the input (i.e. into the sample's own
nested RawData/), so each call moves ONLY its own sample's *.trimmed/*.lengths
into that same sample's Trimmed/ — safe under concurrency.
"""
from __future__ import annotations

import shutil

from .utils import run, log, seq_type, done, list_r1, leaf_dir, find_r2_for_r1


def trim_one(cfg, r1) -> None:
    st = seq_type(r1.name)
    trimmed_dir = leaf_dir(r1) / "Trimmed"
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    if st in ("csRNA", "sRNA"):                              # single-end
        out = trimmed_dir / f"{r1.name}.trimmed"
        if done(out):
            log.info("  skip (done): %s", r1.name); return
        run(f"homerTools trim -3 {cfg.trim_adapter} -mis {cfg.trim_mis} "
            f"-minMatchLength {cfg.trim_minmatch} -min {cfg.trim_min} "
            f"-max {cfg.trim_max} {r1}", label=f"trim SE {r1.name}")
        for suffix in (".trimmed", ".lengths"):              # move only THIS sample's outputs
            src = r1.parent / f"{r1.name}{suffix}"           # homerTools writes next to r1 (RawData/)
            if src.exists():
                shutil.move(str(src), str(trimmed_dir / src.name))
    elif st == "totalRNA":                                   # paired-end, via homerTools -pe
        r2 = find_r2_for_r1(r1)
        if r2 is None:
            log.warning("trim: could not uniquely identify R2 mate for %s in %s "
                       "(same-basename substitution doesn't hold for e.g. ENCODE "
                       "downloads, where mates have different accessions) — skipping.",
                       r1.name, r1.parent)
            return
        out1 = trimmed_dir / f"{r1.name}.trimmed"
        if done(out1):
            log.info("  skip (done): %s", r1.name); return
        run(f"homerTools trim -3 {cfg.trim_adapter} -mis {cfg.trim_mis} "
            f"-minMatchLength {cfg.trim_minmatch} -min {cfg.totalrna_trim_min} "
            f"-max {cfg.totalrna_trim_max} -pe {r1}", label=f"trim PE {r1.name}")
        for src_name in (r1.name, r2.name):                   # move BOTH mates' outputs
            for suffix in (".trimmed", ".lengths"):
                src = r1.parent / f"{src_name}{suffix}"       # homerTools writes next to r1/r2 (RawData/)
                if src.exists():
                    shutil.move(str(src), str(trimmed_dir / src.name))
    else:
        log.warning("trim: skipping untyped file %s", r1.name)


def run_trim(cfg, sample_index=None) -> None:
    r1s = list_r1(cfg)
    if not r1s:
        log.info("trim: no *_R1*.fastq[.gz] under nested RawData/ dirs in %s", cfg.project); return
    if sample_index is not None:
        if not (0 <= sample_index < len(r1s)):
            raise IndexError(f"sample_index {sample_index} out of range (0-{len(r1s)-1})")
        r1s = [r1s[sample_index]]
    for r1 in r1s:
        trim_one(cfg, r1)