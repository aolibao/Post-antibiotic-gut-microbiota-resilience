#!/usr/bin/env python3
"""Subject-aware longitudinal validation for Fig. 3 HDR/LDR dynamics.

The inference unit is the participant.  Sample rows are never treated as
independent replicates.  For each taxon and participant, this module estimates
four prespecified summaries on log10(CPM + 1) relative abundance:

``cip_vs_pre``
    Mean during the fixed ciprofloxacin interval (0.5 < day <= 5.5) minus the
    pre-exposure mean (day <= 0.5).
``late_recovery_vs_pre``
    Mean at C4-C6 minus the pre-exposure mean.
``late_recovery_vs_cip``
    Mean at C4-C6 minus the fixed-interval mean.
``recovery_slope_per_10_days``
    Within-participant OLS slope over all day > 5.5 observations, scaled to ten
    days.

HDR-minus-LDR contrasts are accompanied by stratified participant bootstrap
intervals, participant-label permutation p-values, Hedges' g, and BH FDR.  The
late-recovery grouping was itself defined from C4-C6 diversity, so tests using
C4-C6 are explicitly labelled descriptive/selection-coupled rather than
independent confirmation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

try:  # script execution and package-style imports are both supported
    from .summarize_posteriors import (
        CANONICAL_STAGES,
        DEFAULT_ANALYSIS_ROOT,
        bh_fdr,
        load_abundance_inputs,
    )
except ImportError:  # pragma: no cover - exercised by CLI invocation
    from summarize_posteriors import (
        CANONICAL_STAGES,
        DEFAULT_ANALYSIS_ROOT,
        bh_fdr,
        load_abundance_inputs,
    )


DEFAULT_TABLE_DIR = DEFAULT_ANALYSIS_ROOT / "tables"
DEFAULT_SEED = 20_261_009
CIP_START = 0.5
CIP_END = 5.5
LATE_RECOVERY_STAGES = ("C4", "C5", "C6")


def hedges_g(hdr: np.ndarray, ldr: np.ndarray) -> float:
    """Bias-corrected standardized HDR-minus-LDR mean difference."""

    x = np.asarray(hdr, dtype=float)
    y = np.asarray(ldr, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return np.nan
    degrees = len(x) + len(y) - 2
    pooled_variance = (
        (len(x) - 1) * np.var(x, ddof=1) + (len(y) - 1) * np.var(y, ddof=1)
    ) / degrees
    if pooled_variance <= 0:
        difference = float(np.mean(x) - np.mean(y))
        return 0.0 if difference == 0 else float(np.sign(difference) * np.inf)
    correction = 1.0 - 3.0 / (4.0 * degrees - 1.0)
    return float(correction * (np.mean(x) - np.mean(y)) / np.sqrt(pooled_variance))


def _stratified_bootstrap_difference(
    hdr: np.ndarray,
    ldr: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    x = np.asarray(hdr, dtype=float)
    y = np.asarray(ldr, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2 or n_bootstrap <= 0:
        return np.nan, np.nan
    x_index = rng.integers(0, len(x), size=(n_bootstrap, len(x)))
    y_index = rng.integers(0, len(y), size=(n_bootstrap, len(y)))
    differences = x[x_index].mean(axis=1) - y[y_index].mean(axis=1)
    low, high = np.quantile(differences, [0.025, 0.975])
    return float(low), float(high)


def _label_permutation_pvalue(
    hdr: np.ndarray,
    ldr: np.ndarray,
    n_permutations: int,
    rng: np.random.Generator,
) -> float:
    x = np.asarray(hdr, dtype=float)
    y = np.asarray(ldr, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2 or n_permutations <= 0:
        return np.nan
    pooled = np.concatenate([x, y])
    observed = abs(float(np.mean(x) - np.mean(y)))
    exceedances = 0
    n_hdr = len(x)
    for _ in range(n_permutations):
        permuted = rng.permutation(len(pooled))
        difference = abs(
            float(
                pooled[permuted[:n_hdr]].mean()
                - pooled[permuted[n_hdr:]].mean()
            )
        )
        exceedances += difference >= observed - 1e-15
    return float((exceedances + 1) / (n_permutations + 1))


def two_group_subject_inference(
    subject_contrasts: pd.DataFrame,
    *,
    n_bootstrap: int = 2_000,
    n_permutations: int = 9_999,
    seed: int = DEFAULT_SEED,
) -> pd.DataFrame:
    """Compare HDR/LDR taxon contrasts using participants as independent units."""

    required = {"subject", "group", "taxon", "contrast", "value"}
    missing = sorted(required - set(subject_contrasts.columns))
    if missing:
        raise ValueError(f"subject_contrasts missing columns: {missing}")
    duplicate = subject_contrasts.duplicated(["subject", "taxon", "contrast"])
    if duplicate.any():
        raise ValueError("Expected one value per subject, taxon and contrast")
    rng = np.random.default_rng(seed)
    records: list[dict[str, Any]] = []
    # Resampling indices are shared across taxa within a contrast.  This both
    # preserves the cross-taxon dependence structure and avoids millions of
    # Python-level permutation loops.
    for contrast, block in subject_contrasts.groupby(
        "contrast", sort=False, observed=True
    ):
        taxon_order = block["taxon"].drop_duplicates().tolist()
        matrix = block.pivot(
            index=["group", "subject"], columns="taxon", values="value"
        ).reindex(columns=taxon_order)
        hdr_matrix = matrix.loc["HDR"].to_numpy(dtype=float)
        ldr_matrix = matrix.loc["LDR"].to_numpy(dtype=float)
        n_hdr, n_taxa = hdr_matrix.shape
        n_ldr = ldr_matrix.shape[0]
        hdr_n = np.isfinite(hdr_matrix).sum(axis=0)
        ldr_n = np.isfinite(ldr_matrix).sum(axis=0)
        with np.errstate(invalid="ignore"):
            hdr_mean = np.nanmean(hdr_matrix, axis=0)
            ldr_mean = np.nanmean(ldr_matrix, axis=0)
        observed = hdr_mean - ldr_mean

        if n_bootstrap > 0 and n_hdr >= 2 and n_ldr >= 2:
            hdr_indices = rng.integers(
                0, n_hdr, size=(n_bootstrap, n_hdr)
            )
            ldr_indices = rng.integers(
                0, n_ldr, size=(n_bootstrap, n_ldr)
            )
            with np.errstate(invalid="ignore"):
                bootstrap = (
                    np.nanmean(hdr_matrix[hdr_indices], axis=1)
                    - np.nanmean(ldr_matrix[ldr_indices], axis=1)
                )
            ci_low, ci_high = np.nanquantile(
                bootstrap, [0.025, 0.975], axis=0
            )
        else:
            ci_low = np.full(n_taxa, np.nan)
            ci_high = np.full(n_taxa, np.nan)

        if n_permutations > 0 and n_hdr >= 2 and n_ldr >= 2:
            pooled = np.vstack([hdr_matrix, ldr_matrix])
            n_total = len(pooled)
            random_keys = rng.random((n_permutations, n_total))
            chosen = np.argpartition(
                random_keys, kth=n_hdr - 1, axis=1
            )[:, :n_hdr]
            hdr_selector = np.zeros(
                (n_permutations, n_total), dtype=float
            )
            np.put_along_axis(hdr_selector, chosen, 1.0, axis=1)
            ldr_selector = 1.0 - hdr_selector
            finite = np.isfinite(pooled).astype(float)
            filled = np.nan_to_num(pooled, nan=0.0)
            hdr_denominator = hdr_selector @ finite
            ldr_denominator = ldr_selector @ finite
            with np.errstate(divide="ignore", invalid="ignore"):
                permuted = (
                    (hdr_selector @ filled) / hdr_denominator
                    - (ldr_selector @ filled) / ldr_denominator
                )
            exceedances = np.sum(
                np.abs(permuted) >= np.abs(observed)[None, :] - 1e-15,
                axis=0,
            )
            pvalues = (exceedances + 1.0) / (n_permutations + 1.0)
        else:
            pvalues = np.full(n_taxa, np.nan)

        for index, taxon in enumerate(taxon_order):
            hdr = hdr_matrix[:, index]
            ldr = ldr_matrix[:, index]
            records.append(
                {
                    "taxon": taxon,
                    "contrast": contrast,
                    "estimate_HDR_minus_LDR": observed[index],
                    "ci025_subject_bootstrap": ci_low[index],
                    "ci975_subject_bootstrap": ci_high[index],
                    "p_subject_label_permutation": pvalues[index],
                    "hedges_g_HDR_minus_LDR": hedges_g(hdr, ldr),
                    "HDR_mean": hdr_mean[index],
                    "LDR_mean": ldr_mean[index],
                    "HDR_n_subjects": int(hdr_n[index]),
                    "LDR_n_subjects": int(ldr_n[index]),
                    "bootstrap_replicates": n_bootstrap,
                    "permutation_replicates": n_permutations,
                    "independent_unit": "subject",
                    "group_is_outcome_defined": True,
                    "selection_coupled_to_group_definition": True,
                    "contrast_directly_overlaps_C4_C6_grouping_window": str(
                        contrast
                    ).startswith("late_recovery"),
                }
            )
    results = pd.DataFrame(records)
    if not results.empty:
        results["q_BH_within_contrast"] = np.nan
        for contrast, indices in results.groupby("contrast", sort=False).groups.items():
            results.loc[indices, "q_BH_within_contrast"] = bh_fdr(
                results.loc[indices, "p_subject_label_permutation"].to_numpy()
            )
        results["significant_q05"] = results["q_BH_within_contrast"] < 0.05
    return results


def _phase_means(
    log_abundance: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Return one mean per participant/taxon for each prespecified phase."""

    sample_by_taxon = log_abundance.T
    joined = sample_by_taxon.join(
        metadata[["subject", "time", "canon", "recovery_group"]], how="inner"
    )
    taxon_columns = list(log_abundance.index)
    masks = {
        "pre": joined["time"] <= CIP_START,
        "cip": (joined["time"] > CIP_START) & (joined["time"] <= CIP_END),
        "recovery": joined["time"] > CIP_END,
        "late_recovery": joined["canon"].isin(LATE_RECOVERY_STAGES),
    }
    pieces: list[pd.DataFrame] = []
    for phase, mask in masks.items():
        selected = joined.loc[mask]
        if selected.empty:
            raise ValueError(f"No samples available for phase {phase}")
        summary = selected.groupby(
            ["recovery_group", "subject"], observed=True
        )[taxon_columns].mean()
        summary["phase"] = phase
        pieces.append(summary.reset_index())
    return pd.concat(pieces, ignore_index=True)


def calculate_subject_contrasts(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate prespecified within-participant changes for every taxon."""

    # CPM is scale-independent after within-sample normalization.  A +1 CPM
    # offset is documented and fixed before any group comparison.
    log_abundance = np.log10(abundance * 1_000_000.0 + 1.0)
    phase_long = _phase_means(log_abundance, metadata)
    taxa = list(log_abundance.index)
    phase_wide = phase_long.set_index(
        ["recovery_group", "subject", "phase"]
    )[taxa].unstack("phase")
    required_phases = {"pre", "cip", "late_recovery"}
    present = set(phase_wide.columns.get_level_values("phase"))
    if not required_phases.issubset(present):
        raise ValueError(f"Missing phase means: {sorted(required_phases - present)}")

    contrast_arrays = {
        "cip_vs_pre": phase_wide.xs("cip", level="phase", axis=1)
        - phase_wide.xs("pre", level="phase", axis=1),
        "late_recovery_vs_pre": phase_wide.xs(
            "late_recovery", level="phase", axis=1
        ) - phase_wide.xs("pre", level="phase", axis=1),
        "late_recovery_vs_cip": phase_wide.xs(
            "late_recovery", level="phase", axis=1
        ) - phase_wide.xs("cip", level="phase", axis=1),
    }
    records: list[pd.DataFrame] = []
    for name, matrix in contrast_arrays.items():
        frame = matrix.reset_index().melt(
            id_vars=["recovery_group", "subject"],
            var_name="taxon",
            value_name="value",
        )
        frame = frame.rename(columns={"recovery_group": "group"})
        frame["contrast"] = name
        records.append(frame)

    # Within-subject recovery slope; no participant contributes more than one
    # value to the subsequent group comparison.
    sample_by_taxon = log_abundance.T.join(
        metadata[["subject", "time", "canon", "recovery_group"]], how="inner"
    )
    recovery = sample_by_taxon.loc[sample_by_taxon["time"] > CIP_END]
    slope_rows: list[dict[str, Any]] = []
    for (group, subject), block in recovery.groupby(
        ["recovery_group", "subject"], observed=True
    ):
        x = block["time"].to_numpy(dtype=float)
        x = x - np.mean(x)
        denominator = float(np.sum(x**2))
        for taxon in taxa:
            y = block[taxon].to_numpy(dtype=float)
            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() < 3 or denominator <= 0:
                slope = np.nan
            else:
                centred_x = x[valid]
                centred_y = y[valid] - np.mean(y[valid])
                slope = float(np.sum(centred_x * centred_y) / np.sum(centred_x**2) * 10.0)
            slope_rows.append(
                {
                    "group": group,
                    "subject": subject,
                    "taxon": taxon,
                    "value": slope,
                    "contrast": "recovery_slope_per_10_days",
                }
            )
    records.append(pd.DataFrame(slope_rows))
    result = pd.concat(records, ignore_index=True)
    return result[["subject", "group", "taxon", "contrast", "value"]]


def canonical_stage_summary(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    n_bootstrap: int = 2_000,
    seed: int = DEFAULT_SEED,
) -> pd.DataFrame:
    """Subject-level stage means with stratified participant bootstrap CIs."""

    log_abundance = np.log10(abundance * 1_000_000.0 + 1.0)
    joined = log_abundance.T.join(
        metadata[["subject", "canon", "recovery_group"]], how="inner"
    )
    taxa = list(log_abundance.index)
    subject_stage = joined.groupby(
        ["recovery_group", "subject", "canon"], observed=True
    )[taxa].mean()
    rng = np.random.default_rng(seed)
    records: list[dict[str, Any]] = []
    for (group, canon), block in subject_stage.groupby(
        ["recovery_group", "canon"], observed=True
    ):
        values = block.to_numpy(dtype=float)
        if n_bootstrap > 0 and len(values) >= 2:
            indices = rng.integers(0, len(values), size=(n_bootstrap, len(values)))
            boot = values[indices].mean(axis=1)
            low, high = np.quantile(boot, [0.025, 0.975], axis=0)
        else:
            low = np.full(len(taxa), np.nan)
            high = np.full(len(taxa), np.nan)
        means = np.nanmean(values, axis=0)
        medians = np.nanmedian(values, axis=0)
        for index, taxon in enumerate(taxa):
            records.append(
                {
                    "group": group,
                    "canon": canon,
                    "stage_index": CANONICAL_STAGES.index(canon)
                    if canon in CANONICAL_STAGES
                    else np.nan,
                    "taxon": taxon,
                    "mean_log10_CPM_plus1": means[index],
                    "median_log10_CPM_plus1": medians[index],
                    "ci025_subject_bootstrap": low[index],
                    "ci975_subject_bootstrap": high[index],
                    "n_subjects": len(values),
                    "independent_unit": "subject",
                }
            )
    return pd.DataFrame(records).sort_values(
        ["group", "stage_index", "taxon"], kind="mergesort"
    )


def run_longitudinal_analysis(
    ldr_reads: Path,
    ldr_metadata: Path,
    hdr_reads: Path,
    hdr_metadata: Path,
    output_dir: Path = DEFAULT_TABLE_DIR,
    *,
    ldr_taxonomy: Path | None = None,
    hdr_taxonomy: Path | None = None,
    n_bootstrap: int = 2_000,
    n_permutations: int = 9_999,
    seed: int = DEFAULT_SEED,
) -> dict[str, Any]:
    """Run all subject-aware longitudinal analyses and write audit tables."""

    ldr_abundance, ldr_meta = load_abundance_inputs(ldr_reads, ldr_metadata, "LDR")
    hdr_abundance, hdr_meta = load_abundance_inputs(hdr_reads, hdr_metadata, "HDR")
    if list(ldr_abundance.index) != list(hdr_abundance.index):
        raise ValueError("HDR and LDR reads do not have identical taxon order")
    abundance = pd.concat([ldr_abundance, hdr_abundance], axis=1)
    metadata = pd.concat([ldr_meta, hdr_meta], axis=0)
    if metadata.index.duplicated().any():
        raise ValueError("Sample identifiers overlap between HDR and LDR inputs")

    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    subject_contrasts = calculate_subject_contrasts(abundance, metadata)
    group_contrasts = two_group_subject_inference(
        subject_contrasts,
        n_bootstrap=n_bootstrap,
        n_permutations=n_permutations,
        seed=seed,
    )
    stage_summary = canonical_stage_summary(
        abundance, metadata, n_bootstrap=n_bootstrap, seed=seed + 1
    )
    subject_contrasts.to_csv(
        output_dir / "longitudinal_subject_contrasts.csv", index=False
    )
    group_contrasts.to_csv(
        output_dir / "longitudinal_group_contrasts.csv", index=False
    )
    stage_summary.to_csv(
        output_dir / "longitudinal_stage_summary.csv", index=False
    )

    genus_manifest: dict[str, Any] = {"status": "not_requested"}
    if ldr_taxonomy is not None and hdr_taxonomy is not None:
        ldr_taxonomy_frame = pd.read_csv(ldr_taxonomy, sep="\t", index_col=0)
        hdr_taxonomy_frame = pd.read_csv(hdr_taxonomy, sep="\t", index_col=0)
        if "Genus" not in ldr_taxonomy_frame or "Genus" not in hdr_taxonomy_frame:
            raise ValueError("Taxonomy inputs must contain a Genus column")
        ldr_mapping = ldr_taxonomy_frame["Genus"].astype(str)
        hdr_mapping = hdr_taxonomy_frame["Genus"].astype(str)
        ldr_mapping = ldr_mapping.reindex(ldr_abundance.index)
        hdr_mapping = hdr_mapping.reindex(hdr_abundance.index)
        if ldr_mapping.isna().any() or hdr_mapping.isna().any():
            raise ValueError("Taxonomy inputs do not cover every fixed species")
        if not ldr_mapping.equals(hdr_mapping):
            raise ValueError("HDR and LDR genus mappings differ")
        genus_abundance = abundance.groupby(ldr_mapping, sort=True).sum()
        genus_subject = calculate_subject_contrasts(genus_abundance, metadata)
        genus_group = two_group_subject_inference(
            genus_subject,
            n_bootstrap=n_bootstrap,
            n_permutations=n_permutations,
            seed=seed + 2,
        )
        genus_stage = canonical_stage_summary(
            genus_abundance, metadata, n_bootstrap=n_bootstrap, seed=seed + 3
        )
        genus_subject.to_csv(
            output_dir / "longitudinal_genus_subject_contrasts.csv", index=False
        )
        genus_group.to_csv(
            output_dir / "longitudinal_genus_group_contrasts.csv", index=False
        )
        genus_stage.to_csv(
            output_dir / "longitudinal_genus_stage_summary.csv", index=False
        )
        focus = ["Faecalibacterium", "Bacteroides", "Alistipes"]
        pd.concat(
            [
                genus_group.loc[genus_group["taxon"].isin(focus)].assign(
                    table="group_contrast"
                ),
                genus_stage.loc[genus_stage["taxon"].isin(focus)].assign(
                    table="stage_summary"
                ),
            ],
            ignore_index=True,
            sort=False,
        ).to_csv(output_dir / "longitudinal_focus_genera.csv", index=False)
        genus_manifest = {
            "status": "written",
            "n_genera": int(genus_abundance.shape[0]),
            "focus_genera": focus,
            "taxonomy_LDR": str(Path(ldr_taxonomy).resolve()),
            "taxonomy_HDR": str(Path(hdr_taxonomy).resolve()),
        }

    manifest = {
        "schema_version": 1,
        "estimand": "HDR minus LDR mean of within-subject log10(CPM+1) changes",
        "independent_unit": "subject",
        "n_subjects": metadata.groupby("recovery_group")["subject"].nunique().to_dict(),
        "taxa": int(abundance.shape[0]),
        "samples": int(abundance.shape[1]),
        "ciprofloxacin_interval_days": [CIP_START, CIP_END],
        "pre_definition": f"actual day <= {CIP_START}",
        "cip_definition": f"{CIP_START} < actual day <= {CIP_END}",
        "recovery_definition": f"actual day > {CIP_END}",
        "late_recovery_definition": list(LATE_RECOVERY_STAGES),
        "transform": "log10(relative abundance * 1e6 + 1)",
        "bootstrap": {
            "type": "stratified nonparametric subject bootstrap",
            "replicates": n_bootstrap,
            "interval": "percentile 95%",
        },
        "permutation": {
            "type": "two-sided subject-label permutation with fixed group sizes",
            "replicates": n_permutations,
            "plus_one_correction": True,
        },
        "multiplicity": "Benjamini-Hochberg separately within each prespecified contrast",
        "seed": seed,
        "genus_analysis": genus_manifest,
        "selection_warning": (
            "HDR/LDR were selected from C4-C6 diversity, so every between-group result "
            "is conditional on an outcome-derived grouping and is descriptive rather "
            "than independent validation. Late-recovery contrasts additionally overlap "
            "the C4-C6 grouping window directly."
        ),
        "qPCR_warning": (
            "No measured qPCR is available; abundance statistics and MDSINE2 parameters "
            "must not be interpreted as absolute biomass."
        ),
    }
    (output_dir / "longitudinal_methods.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ldr-reads", type=Path,
        default=DEFAULT_ANALYSIS_ROOT / "inputs" / "LDR" / "reads.tsv",
    )
    parser.add_argument(
        "--ldr-metadata", type=Path,
        default=DEFAULT_ANALYSIS_ROOT / "inputs" / "LDR" / "metadata.tsv",
    )
    parser.add_argument(
        "--hdr-reads", type=Path,
        default=DEFAULT_ANALYSIS_ROOT / "inputs" / "HDR" / "reads.tsv",
    )
    parser.add_argument(
        "--hdr-metadata", type=Path,
        default=DEFAULT_ANALYSIS_ROOT / "inputs" / "HDR" / "metadata.tsv",
    )
    parser.add_argument(
        "--ldr-taxonomy", type=Path,
        default=DEFAULT_ANALYSIS_ROOT / "inputs" / "LDR" / "taxonomy.tsv",
    )
    parser.add_argument(
        "--hdr-taxonomy", type=Path,
        default=DEFAULT_ANALYSIS_ROOT / "inputs" / "HDR" / "taxonomy.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_TABLE_DIR)
    parser.add_argument("--bootstrap", type=int, default=2_000)
    parser.add_argument("--permutations", type=int, default=9_999)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    manifest = run_longitudinal_analysis(
        args.ldr_reads,
        args.ldr_metadata,
        args.hdr_reads,
        args.hdr_metadata,
        args.output_dir,
        ldr_taxonomy=args.ldr_taxonomy,
        hdr_taxonomy=args.hdr_taxonomy,
        n_bootstrap=args.bootstrap,
        n_permutations=args.permutations,
        seed=args.seed,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
