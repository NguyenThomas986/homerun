"""Step 5 — call TSRs for matched individual and combined TagDirs.

Each csRNA replicate uses the sRNA replicate with the same replicate marker
and condition context. A matching totalRNA replicate is included when it
exists. Outputs are written under Species/TSS/ with replicate-specific names,
for example ``K562_csRNA_r1.tss.txt``. Each condition-preserving csRNA combo
is also called against the corresponding sRNA/totalRNA combos, producing an
output such as ``Testicle_D2_Post_BDAD_NO_RA_csRNA-combo.tss.txt``.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import re

from .utils import (
    assay_of_leaf,
    done,
    iter_leaf_dirs,
    iter_samples,
    log,
    replicate_of_leaf,
    run,
)


def _condition_key(leaf_name: str) -> tuple[str, ...]:
    """Return non-assay, non-replicate tokens used to match controls."""
    tokens = leaf_name.split("_")
    kept = []

    for token in tokens:
        low = token.lower()

        if re.fullmatch(r"r(ep)?\d+", token):
            continue

        if (
            low.startswith("csrna")
            or low.startswith("srna")
            or low.startswith("totalrna")
            or low == "rna"
        ):
            continue

        kept.append(token)

    return tuple(kept)


def _replicate_tagdirs(cfg, species, sample):
    """Return one existing TagDir per distinct leaf identity."""
    found = {}

    for sp, sa, leaf_name, _r1 in iter_leaf_dirs(cfg):
        if sp != species or sa != sample:
            continue

        tagdir = cfg.leaf_tagdir(species, sample, leaf_name)

        if tagdir.is_dir():
            found[leaf_name] = tagdir

    return found


def _combo_tagdirs(cfg, species, sample):
    """Return ``combo_leaf -> path`` for this sample's combo TagDirs."""
    found = {}
    prefix = f"{sample}_"

    for tagdir in sorted(
        (cfg.species_dir(species) / "TagDirs").glob("*-combo")
    ):
        name = tagdir.name

        if not tagdir.is_dir() or not name.startswith(prefix):
            continue

        combo_leaf = name[len(prefix):-len("-combo")]
        found[combo_leaf] = tagdir

    return found


def _matching_leaf(leaves, source_leaf, assay):
    rep = replicate_of_leaf(source_leaf)
    context = _condition_key(source_leaf)

    matches = [
        (leaf, path)
        for leaf, path in leaves.items()
        if assay_of_leaf(leaf) == assay
        and replicate_of_leaf(leaf) == rep
        and _condition_key(leaf) == context
    ]

    if len(matches) == 1:
        return matches[0]

    return None


def _matching_combo(combos, source_leaf, assay):
    context = _condition_key(source_leaf)

    matches = [
        (leaf, path)
        for leaf, path in combos.items()
        if assay_of_leaf(leaf) == assay
        and _condition_key(leaf) == context
    ]

    if len(matches) == 1:
        return matches[0]

    return None


def run_tss(cfg, group=None) -> None:
    """Call TSRs for each csRNA replicate and condition-preserving combo."""
    found_csrna = False
    jobs = []

    for species, sample in iter_samples(cfg):
        if group is not None and (species, sample) != group:
            continue

        leaves = _replicate_tagdirs(cfg, species, sample)

        cs_leaves = sorted(
            (leaf, path)
            for leaf, path in leaves.items()
            if assay_of_leaf(leaf) == "csRNA"
        )

        if not cs_leaves:
            continue

        found_csrna = True

        tss_dir = cfg.sample_tss(species, sample)
        tss_dir.mkdir(parents=True, exist_ok=True)

        # Individual replicate TSS calls
        for cs_leaf, cs_dir in cs_leaves:
            srna_match = _matching_leaf(
                leaves,
                cs_leaf,
                "sRNA",
            )

            if srna_match is None:
                log.warning(
                    "TSS: %s/%s/%s has no unique matching sRNA replicate — skipping",
                    species,
                    sample,
                    cs_leaf,
                )
                continue

            srna_leaf, srna_dir = srna_match

            out = tss_dir / f"{sample}_{cs_leaf}"

            if done(f"{out}.tss.txt"):
                log.info(
                    "  skip (done): %s.tss.txt",
                    out,
                )
                continue

            cmd = (
                f"findcsRNATSS.pl {cs_dir} "
                f"-o {out} "
                f"-genome {cfg.genome} "
                f"-ntagThreshold {cfg.ntag_threshold} "
                f"-i {srna_dir}"
            )

            rna_match = _matching_leaf(
                leaves,
                cs_leaf,
                "totalRNA",
            )

            if rna_match is not None:
                rna_leaf, rna_dir = rna_match
                cmd += f" -rna {rna_dir}"

                log.info(
                    "  %s uses controls %s and %s",
                    cs_leaf,
                    srna_leaf,
                    rna_leaf,
                )
            else:
                log.info(
                    "  %s uses control %s; no matching totalRNA replicate",
                    cs_leaf,
                    srna_leaf,
                )

            jobs.append(
                (
                    cmd,
                    f"findcsRNATSS {species}/{sample}/{cs_leaf}",
                    tss_dir,
                )
            )

        # Combined TagDir TSS calls
        combos = _combo_tagdirs(
            cfg,
            species,
            sample,
        )

        cs_combos = sorted(
            (leaf, path)
            for leaf, path in combos.items()
            if assay_of_leaf(leaf) == "csRNA"
        )

        for cs_leaf, cs_dir in cs_combos:
            srna_match = _matching_combo(
                combos,
                cs_leaf,
                "sRNA",
            )

            if srna_match is None:
                log.warning(
                    "TSS combo: %s/%s/%s has no unique matching sRNA combo — skipping",
                    species,
                    sample,
                    cs_leaf,
                )
                continue

            srna_leaf, srna_dir = srna_match

            out = tss_dir / f"{sample}_{cs_leaf}-combo"

            if done(f"{out}.tss.txt"):
                log.info(
                    "  skip (done): %s.tss.txt",
                    out,
                )
                continue

            cmd = (
                f"findcsRNATSS.pl {cs_dir} "
                f"-o {out} "
                f"-genome {cfg.genome} "
                f"-ntagThreshold {cfg.ntag_threshold} "
                f"-i {srna_dir}"
            )

            rna_match = _matching_combo(
                combos,
                cs_leaf,
                "totalRNA",
            )

            if rna_match is not None:
                rna_leaf, rna_dir = rna_match
                cmd += f" -rna {rna_dir}"

                log.info(
                    "  combo %s uses controls %s and %s",
                    cs_leaf,
                    srna_leaf,
                    rna_leaf,
                )
            else:
                log.info(
                    "  combo %s uses control %s; no matching totalRNA combo",
                    cs_leaf,
                    srna_leaf,
                )

            jobs.append(
                (
                    cmd,
                    f"findcsRNATSS {species}/{sample}/{cs_leaf}-combo",
                    tss_dir,
                )
            )

    if jobs:
        workers = max(1, cfg.threads)

        log.info(
            "TSS: running %d job(s) with %d worker(s)",
            len(jobs),
            workers,
        )

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(
                    run,
                    cmd,
                    label=label,
                    cwd=cwd,
                )
                for cmd, label, cwd in jobs
            ]

            for future in as_completed(futures):
                future.result()

    if not found_csrna:
        log.info(
            "TSS: no csRNA TagDirs under %s",
            cfg.project,
        )