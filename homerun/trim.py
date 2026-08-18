"""Step 1 — Trim. homerTools trim for csRNA/sRNA (single-end). totalRNA
(paired-end) goes through skewer instead — homerTools trim's -pe mode needs
Config attributes (totalrna_trim_min/max) that don't exist here, so totalRNA
stays on skewer as before rather than switching tools.

R1's R2 mate is located via utils.find_r2_for_r1() — matched by PARSED
identity (species/sample/leaf_name), not by swapping "_R1" for "_R2" in the
filename. That substitution assumes both mates share the same basename
apart from the R1/R2 tag, which breaks for ENCODE-style downloads where
each mate has its own independent accession (e.g. an R1 file's accession
number is completely unrelated to its R2 mate's) — swapping _R1->_R2 on
such a name silently builds a path to a file that was never going to
exist, and skewer then fails to open it.

Operates per R1 file so it can run inside a SLURM array (one task per sample).
homerTools writes outputs next to the input (i.e. into the sample's own
nested RawData/), so each call moves ONLY its own sample's *.trimmed/*.lengths
into that same sample's Trimmed/ — safe under concurrency. skewer is pointed
at Trimmed/ directly via -o, and its paired outputs are renamed to the same
"<matefile>.trimmed" convention homerTools uses, so mapping.py and qc.py
don't need to know which tool produced a given sample's trimmed reads.
"""
from __future__ import annotations

import shutil

from .utils import run, log, seq_type, done, list_r1, leaf_dir, find_r2_for_r1


def _skewer_cmd(cfg, r1, r2, out_prefix: str) -> str:
    """skewer paired-end trim command for totalRNA. Adapter is shared with
    the homerTools SE trim (cfg.trim_adapter); quality/min-length are
    skewer-specific and read from cfg with safe defaults so a Config that
    doesn't define them still runs, rather than crashing the way referencing
    the homerTools-only cfg.totalrna_trim_min/max did."""
    quality = getattr(cfg, "skewer_quality", 20)
    min_len = getattr(cfg, "skewer_min_length", 18)
    return (f"skewer -m pe -x {cfg.trim_adapter} "
            f"-q {quality} -l {min_len} "
            f"-t {cfg.threads} "
            f"-o {out_prefix} {r1} {r2}")


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
    elif st == "totalRNA":                                   # paired-end, via skewer
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
        prefix = r1.name.split("_R1")[0]
        out_prefix = trimmed_dir / prefix
        run(_skewer_cmd(cfg, r1, r2, str(out_prefix)), label=f"trim PE {r1.name}")

        # skewer (-o out_prefix) writes "<prefix>-trimmed-pair1.fastq[.gz]",
        # "<prefix>-trimmed-pair2.fastq[.gz]", and "<prefix>-trimmed.log"
        # directly into trimmed_dir already — nothing to move dirs for, only
        # the pair outputs need renaming so mapping.py/qc.py can find them
        # under the same "<matefile>.trimmed" name homerTools trim uses.
        pair_map = (
            (r1, "pair1", out1),
            (r2, "pair2", trimmed_dir / f"{r2.name}.trimmed"),
        )
        for mate, tag, dest in pair_map:
            for ext in (".fastq", ".fastq.gz"):
                src = trimmed_dir / f"{prefix}-trimmed-{tag}{ext}"
                if src.exists():
                    shutil.move(str(src), str(dest))
                    break
            else:
                log.warning("trim: expected skewer output for %s (%s-trimmed-%s.fastq[.gz]) "
                           "not found in %s", mate.name, prefix, tag, trimmed_dir)
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