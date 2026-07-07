"""Step 5 — Find TSRs with findcsRNATSS.pl, per Species/Sample.

For each sample's csRNA-combo TagDir, uses the matching sRNA-combo as input
control and, if present, the totalRNA-combo as the -rna reference. The -rna
reference is optional: without total RNA the TSRs are still called, just
without stable/unstable assignment. Output lands in Species/Sample/TSS/,
one TSS per sample (not per assay).
"""
from __future__ import annotations

from .utils import run, log, done, iter_samples


def run_tss(cfg) -> None:
    found_any = False

    for species, sample in iter_samples(cfg):
        cs_dir = cfg.combo_tagdir(species, sample, "csRNA")
        if not cs_dir.is_dir():
            continue
        found_any = True

        sRNA_dir = cfg.combo_tagdir(species, sample, "sRNA")
        if not sRNA_dir.is_dir():
            log.warning("TSS: %s/%s has csRNA but no sRNA-combo — skipping (need an input control).",
                        species, sample)
            continue

        tss_dir = cfg.sample_tss(species, sample)
        tss_dir.mkdir(parents=True, exist_ok=True)
        out = tss_dir / sample
        if done(f"{out}.tss.txt"):
            log.info("  skip (done): %s/%s/TSS/%s.tss.txt", species, sample, sample)
            continue

        rna_dir = cfg.combo_tagdir(species, sample, "totalRNA")
        if not rna_dir.is_dir():
            rna_dir = None

        cmd = (f"findcsRNATSS.pl {cs_dir} -o {out} -genome {cfg.genome} "
               f"-ntagThreshold {cfg.ntag_threshold} -i {sRNA_dir}")
        if rna_dir:
            cmd += f" -rna {rna_dir}"
            log.info("  using total-RNA reference: %s/%s/totalRNA-combo", species, sample)
        else:
            log.info("  no total-RNA reference for %s/%s (stability will be unavailable)",
                     species, sample)
        run(cmd, label=f"findcsRNATSS {species}/{sample}", cwd=tss_dir)

    if not found_any:
        log.info("TSS: no csRNA-combo TagDirs under %s", cfg.project)