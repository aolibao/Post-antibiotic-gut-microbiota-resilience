#!/usr/bin/env python3
"""Reconstruct manuscript Figure 2F with the installed LEfSe implementation.

Raster reverse-calibration establishes that the panel has 45 (not 30) species,
that all 45 x 16 cells retain an LDA value, and that colour sign follows the
raw LDR-minus-HDR mean rather than LEfSe's occasionally unstable class label.
LEfSe therefore runs without inferential filtering on the complete candidate
universe.  Exact Mann-Whitney P values are BH-adjusted over the 45 displayed
species at each time and form an independent star layer.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import run_all as audit


PREVALENCE_MIN = 0.20
MEAN_ABUNDANCE_MIN_PERCENT = 0.01
LDA_THRESHOLD = 2.0
N_BOOTSTRAPS = 30
# LEfSe applies the Kruskal-Wallis gate with a strict ``p < alpha`` test.
# Using alpha=1.0 therefore still drops mathematically valid p=1 features,
# despite this reconstruction's intended no-inferential-filter LDA run.  A
# value infinitesimally above one disables that gate without changing any
# observed P value or the separate BH-controlled annotation layer.
KW_ALPHA_NO_FILTER = 1.0000001
SOURCE_SPECIES_ORDER = [
    "Flavonifractor plautii",
    "Enterocloster bolteae",
    "Blautia wexlerae",
    "Bifidobacterium longum",
    "Agathobacter rectalis",
    "Bifidobacterium pseudocatenulatum",
    "Mediterraneibacter gnavus",
    "Roseburia hominis",
    "Lawsonibacter asaccharolyticus",
    "Anaerostipes hadrus",
    "Lachnospira eligens",
    "Ruminococcus bicirculans (ex Wegman et al. 2014)",
    "Simiaoa sunii",
    "Akkermansia muciniphila",
    "Faecalibacterium prausnitzii",
    "Blautia obeum",
    "Faecalibacterium sp. I4-1-79",
    "Segatella hominis",
    "Collinsella aerofaciens",
    "Blautia sp. KLE 1732 HM 1032",
    "Dorea longicatena",
    "Faecalibacterium sp. I3-3-33",
    "Faecalibacterium duncaniae",
    "Faecalibacterium sp. IP-1-18",
    "Vescimonas coprocola",
    "Vescimonas fastidiosa",
    "Faecalibacterium sp. I3-3-89",
    "uncultured Alistipes sp.",
    "Anaerobutyricum hallii",
    "Alistipes shahii",
    "Coprococcus sp. ART55/1",
    "Wujia chipingensis",
    "Coprococcus eutactus",
    "Alistipes communis",
    "Eshraghiella crossota",
    "Alistipes finegoldii",
    "Akkermansia massiliensis",
    "Faecalibacterium sp. I2-3-92",
    "Faecalibacterium sp. IP-3-29",
    "Ruminococcus sp. FMBCY1",
    "Dorea formicigenerans",
    "Alistipes onderdonkii",
    "Schaalia turicensis",
    "Segatella copri",
    "[Ruminococcus] torques",
]
N_DISPLAY = len(SOURCE_SPECIES_ORDER)


def _tool_paths() -> tuple[Path, Path, dict[str, str]]:
    formatter = shutil.which("lefse_format_input.py")
    runner = shutil.which("lefse_run.py")
    if formatter is None or runner is None:
        raise RuntimeError("LEfSe executables were not found on PATH")
    formatter_path = Path(formatter).resolve()
    runner_path = Path(runner).resolve()
    prefix = runner_path.parent.parent
    r_home = prefix / "lib" / "R"
    if not r_home.exists():
        raise RuntimeError(f"LEfSe R runtime was not found at {r_home}")
    env = os.environ.copy()
    env["R_HOME"] = str(r_home)
    r_lib = str(r_home / "lib")
    prefix_lib = str(prefix / "lib")
    prior = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = ":".join(x for x in [r_lib, prefix_lib, prior] if x)
    return formatter_path, runner_path, env


def _run_lda(
    canon: str,
    values: np.ndarray,
    safe_ids: list[str],
    formatter: Path,
    runner: Path,
    env: dict[str, str],
) -> pd.DataFrame:
    if values.shape[0] == 0:
        return pd.DataFrame(columns=["safe_id", "enriched_class", "lda_abs"])
    with tempfile.TemporaryDirectory(prefix=f"lefse_{canon}_") as temporary:
        root = Path(temporary)
        source = root / "input.tsv"
        formatted = root / "input.in"
        result = root / "result.tsv"
        sample_classes = ["LDR"] * 15 + ["HDR"] * 15
        with source.open("w", encoding="utf-8") as handle:
            handle.write("class\t" + "\t".join(sample_classes) + "\n")
            # LEfSe 1.1.2's formatter expects a subject row in this build.
            handle.write("subject\t" + "\t".join(f"S{i:02d}" for i in range(30)) + "\n")
            for safe_id, row in zip(safe_ids, values, strict=True):
                handle.write(safe_id + "\t" + "\t".join(f"{value:.12g}" for value in row) + "\n")
        subprocess.run(
            [
                str(formatter),
                str(source),
                str(formatted),
                "-c",
                "1",
                "-u",
                "2",
                "-o",
                "1000000",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            [
                str(runner),
                str(formatted),
                str(result),
                "-a",
                str(KW_ALPHA_NO_FILTER),
                "--wilc",
                "0",
                "-l",
                "0.0",
                "-b",
                str(N_BOOTSTRAPS),
            ],
            check=True,
            capture_output=True,
            text=True,
            env=env,
        )
        parsed = pd.read_csv(
            result,
            sep="\t",
            header=None,
            names=["safe_id", "log10_max_mean", "enriched_class", "lda_abs", "wilcoxon"],
        )
    return parsed[["safe_id", "enriched_class", "lda_abs"]]


def eligible_feature_universe(ds: audit.Dataset) -> pd.Index:
    """Return the LEfSe candidate universe for the dataset's current groups."""
    selected = ds.metrics.loc[ds.metrics["recovery_group"].isin(["LDR", "HDR"])].copy()
    selected_ids = selected["ids"].tolist()
    prevalence = (ds.species_named[selected_ids] > 0).mean(axis=1)
    mean_abundance = ds.species_named[selected_ids].mean(axis=1)
    eligible = (prevalence >= PREVALENCE_MIN) & (mean_abundance >= MEAN_ABUNDANCE_MIN_PERCENT)
    return pd.Index(ds.species_named.index[eligible], name="species")


def lefse_all_timepoints(
    ds: audit.Dataset,
    tables: Path,
    *,
    frozen_features: pd.Index | list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    selected = ds.metrics.loc[ds.metrics["recovery_group"].isin(["LDR", "HDR"])].copy()
    if frozen_features is None:
        features = eligible_feature_universe(ds)
        feature_rule = "current-group prevalence/mean-abundance filter"
    else:
        features = pd.Index([str(value) for value in frozen_features], name="species")
        if features.duplicated().any():
            raise ValueError("Frozen LEfSe candidate universe contains duplicate species")
        missing = features.difference(ds.species_named.index)
        if len(missing):
            raise ValueError(
                "Frozen LEfSe candidate universe contains species absent from the dataset: "
                f"{missing[:10].tolist()}"
            )
        feature_rule = "frozen D0 candidate universe"
    pd.DataFrame({"species": features, "feature_rule": feature_rule}).to_csv(
        tables / "Fig2F_LEfSe_candidate_universe.csv", index=False
    )
    formatter, runner, env = _tool_paths()

    statistics: dict[str, pd.DataFrame] = {}
    jobs: dict[object, tuple[str, list[str]]] = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for canon in audit.TIME_ORDER:
            rows = selected.loc[selected["canon"] == canon]
            ldr_ids = rows.loc[rows["recovery_group"] == "LDR", "ids"].tolist()
            hdr_ids = rows.loc[rows["recovery_group"] == "HDR", "ids"].tolist()
            if len(ldr_ids) != 15 or len(hdr_ids) != 15:
                raise ValueError(f"{canon}: expected 15 LDR and 15 HDR samples")
            ldr = ds.species_named.loc[features, ldr_ids].to_numpy(dtype=float)
            hdr = ds.species_named.loc[features, hdr_ids].to_numpy(dtype=float)
            test = mannwhitneyu(ldr, hdr, axis=1, alternative="two-sided", method="exact")
            pvalues = np.asarray(test.pvalue, dtype=float)
            identical = (
                (np.ptp(ldr, axis=1) == 0)
                & (np.ptp(hdr, axis=1) == 0)
                & (ldr[:, 0] == hdr[:, 0])
            )
            pvalues[identical] = 1.0
            qvalues = audit.bh_adjust(pvalues)
            frame = pd.DataFrame(
                {
                    "canon": canon,
                    "species": features,
                    "mean_LDR_percent": ldr.mean(axis=1),
                    "mean_HDR_percent": hdr.mean(axis=1),
                    "p_exact_mannwhitney": pvalues,
                    "q_BH_candidate_universe": qvalues,
                }
            )
            frame["passes_candidate_q"] = frame["q_BH_candidate_universe"] < 0.05
            frame["safe_id"] = [f"F{i:05d}" for i in range(len(frame))]
            statistics[canon] = frame
            # Compute the LDA colour value for every prefiltered feature.  The
            # manuscript matrix is coloured in nearly every cell and uses the
            # q/LDA gates for star annotations; limiting LDA itself to q hits
            # incorrectly creates a sparse white checkerboard.
            values = np.concatenate([ldr, hdr], axis=1)
            safe_ids = frame["safe_id"].tolist()
            future = pool.submit(_run_lda, canon, values, safe_ids, formatter, runner, env)
            jobs[future] = (canon, safe_ids)

        lda_results: dict[str, pd.DataFrame] = {}
        for future in as_completed(jobs):
            canon, _ = jobs[future]
            lda_results[canon] = future.result()

    output = []
    for canon in audit.TIME_ORDER:
        frame = statistics[canon].merge(lda_results[canon], on="safe_id", how="left")
        frame = frame.rename(columns={"enriched_class": "lefse_enriched_class"})
        mean_difference = frame["mean_LDR_percent"] - frame["mean_HDR_percent"]
        frame["mean_direction_class"] = np.where(
            mean_difference > 0,
            "LDR",
            np.where(mean_difference < 0, "HDR", "Tie"),
        )
        # The source raster follows raw mean direction in all 720 cells.  LEfSe
        # supplies only the absolute effect magnitude here.
        frame["lda_score_LDR_over_HDR"] = frame["lda_abs"] * np.sign(mean_difference)
        frame["passes_caption_lda"] = frame["lda_abs"].ge(LDA_THRESHOLD).fillna(False)
        output.append(frame)
    full = pd.concat(output, ignore_index=True)
    full.to_csv(tables / "Fig2F_LEfSe_all_eligible_species.csv", index=False)
    pd.DataFrame(
        {
            "parameter": [
                "prevalence_min",
                "mean_abundance_min_percent",
                "BH_q_threshold_within_time",
                "BH_adjustment_universe_for_stars",
                "LDA_run_threshold",
                "LEfSe_KW_alpha_for_no_filter",
                "LEfSe_KW_gate_note",
                "caption_LDA_abs_threshold_audit_only",
                "LDA_bootstraps",
                "displayed_species",
                "LDA_sign_rule",
                "cell_mask_rule",
                "candidate_universe_rule",
            ],
            "value": [
                PREVALENCE_MIN,
                MEAN_ABUNDANCE_MIN_PERCENT,
                0.05,
                "45 displayed source species per time point",
                0.0,
                KW_ALPHA_NO_FILTER,
                "alpha > 1 because LEfSe uses strict p < alpha; preserves exact-p=1 features",
                LDA_THRESHOLD,
                N_BOOTSTRAPS,
                N_DISPLAY,
                "sign(mean_LDR_percent - mean_HDR_percent)",
                "none; every source species/time cell retains its LDA score",
                feature_rule,
            ],
        }
    ).to_csv(tables / "Fig2F_LEfSe_parameters.csv", index=False)
    return full


def _star(qvalue: float) -> str:
    if not np.isfinite(qvalue) or qvalue >= 0.05:
        return ""
    if qvalue < 0.001:
        return "***"
    if qvalue < 0.01:
        return "**"
    if qvalue < 0.05:
        return "*"
    return ""


def _digitize_source_lda(source_tiff: Path) -> pd.DataFrame:
    """Invert source cell-centre RGB values on the linear -5..5 RdBu scale."""
    from matplotlib import colormaps
    from PIL import Image

    with Image.open(source_tiff) as image:
        rgb_image = np.asarray(image.convert("RGB"))
    if rgb_image.shape[:2] != (2834, 2441):
        raise ValueError(f"Unexpected source Figure 2 raster geometry: {rgb_image.shape[:2]}")
    x_centres = np.linspace(164, 2356, N_DISPLAY).round().astype(int)
    y_centres = np.linspace(1733, 2500, len(audit.TIME_ORDER)).round().astype(int)
    observed_rgb = rgb_image[np.ix_(y_centres, x_centres)]
    coordinates = np.linspace(0.0, 1.0, 10_001)
    lookup = (colormaps["RdBu_r"](coordinates)[:, :3] * 255).round().astype(int)
    values = np.empty((len(y_centres), len(x_centres)), dtype=float)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            distance = ((lookup - observed_rgb[row, column]) ** 2).sum(axis=1)
            values[row, column] = -5.0 + 10.0 * coordinates[int(distance.argmin())]
    return pd.DataFrame(values, index=audit.TIME_ORDER, columns=SOURCE_SPECIES_ORDER)


def _write_source_digitization_audit(
    ds: audit.Dataset,
    reconstructed: pd.DataFrame,
    tables: Path,
) -> None:
    """Export source-vs-reconstruction measurements without calibrating data."""
    source_tiff = ds.root / "抗生素稿件 isme j" / "Figures" / "Figure2.tif"
    if not source_tiff.exists():
        return
    source = _digitize_source_lda(source_tiff)
    source.to_csv(tables / "Fig2F_source_digitized_LDA_matrix.csv")
    comparison = (
        source.rename_axis(index="canon", columns="species")
        .stack()
        .rename("source_digitized_lda")
        .to_frame()
        .join(
            reconstructed.rename_axis(index="canon", columns="species")
            .stack()
            .rename("reconstructed_lda")
        )
        .reset_index()
    )
    comparison["difference_reconstructed_minus_source"] = (
        comparison["reconstructed_lda"] - comparison["source_digitized_lda"]
    )
    comparison["direction_match"] = np.sign(comparison["reconstructed_lda"]) == np.sign(
        comparison["source_digitized_lda"]
    )
    comparison.to_csv(tables / "Fig2F_source_digitized_comparison.csv", index=False)

    x = comparison["reconstructed_lda"].to_numpy(dtype=float)
    y = comparison["source_digitized_lda"].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(x)), x])
    intercept, slope = np.linalg.lstsq(design, y, rcond=None)[0]
    calibrated = intercept + slope * x
    pd.DataFrame(
        {
            "metric": [
                "cells",
                "direction_match_fraction",
                "pearson_r",
                "raw_rmse",
                "source_on_reconstructed_intercept",
                "source_on_reconstructed_slope",
                "affine_residual_rmse",
            ],
            "value": [
                len(x),
                comparison["direction_match"].mean(),
                np.corrcoef(x, y)[0, 1],
                np.sqrt(np.mean((x - y) ** 2)),
                intercept,
                slope,
                np.sqrt(np.mean((calibrated - y) ** 2)),
            ],
        }
    ).to_csv(tables / "Fig2F_source_digitized_audit_summary.csv", index=False)


def figure2f_tables(
    ds: audit.Dataset,
    tables: Path,
    *,
    source_digitization_audit: bool = True,
    frozen_features: pd.Index | list[str] | tuple[str, ...] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop-in replacement for ``render_selected_panels.figure2f_tables``."""
    full = lefse_all_timepoints(ds, tables, frozen_features=frozen_features)
    available = set(full["species"])
    missing = [species for species in SOURCE_SPECIES_ORDER if species not in available]
    if missing:
        raise RuntimeError(f"Source Figure 2F species absent from Bracken candidate universe: {missing}")
    order = pd.Index(SOURCE_SPECIES_ORDER, name="species")
    matrix = (
        full.loc[full["species"].isin(order)]
        .pivot(index="canon", columns="species", values="lda_score_LDR_over_HDR")
        .reindex(index=audit.TIME_ORDER, columns=order)
    )
    if matrix.isna().any().any():
        bad = np.argwhere(matrix.isna().to_numpy())
        raise RuntimeError(f"LEfSe did not return all source Figure 2F cells: {bad[:10].tolist()}")
    annotations = full.loc[
        full["species"].isin(order),
        [
            "canon",
            "species",
            "lda_score_LDR_over_HDR",
            "lda_abs",
            "lefse_enriched_class",
            "mean_direction_class",
            "mean_LDR_percent",
            "mean_HDR_percent",
            "p_exact_mannwhitney",
            "passes_caption_lda",
        ],
    ].copy()
    annotations["q_BH_displayed_45"] = np.nan
    for canon in audit.TIME_ORDER:
        selected = annotations["canon"].eq(canon)
        annotations.loc[selected, "q_BH_displayed_45"] = audit.bh_adjust(
            annotations.loc[selected, "p_exact_mannwhitney"].to_numpy(dtype=float)
        )
    annotations["displayed_significant"] = annotations["q_BH_displayed_45"] < 0.05
    annotations["symbol"] = annotations["q_BH_displayed_45"].map(_star)
    annotations["method"] = (
        "Unfiltered LEfSe 1.1.2/MASS LDA magnitude on the 268-feature candidate universe; "
        "sign from raw mean LDR-HDR difference; stars from exact Mann-Whitney P values "
        "BH-adjusted within the 45 displayed species and graded by q"
    )
    matrix.to_csv(tables / "Fig2F_LEfSe_LDA_matrix.csv")
    annotations.to_csv(tables / "Fig2F_LEfSe_annotations.csv", index=False)
    if source_digitization_audit:
        _write_source_digitization_audit(ds, matrix, tables)
    return matrix, annotations
