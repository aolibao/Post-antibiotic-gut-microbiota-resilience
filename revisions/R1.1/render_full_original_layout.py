#!/usr/bin/env python3
"""Re-run complete manuscript-layout Figures 2 and 4 for R1.1 definitions.

The alternative recovery definitions are fixed by
``R1_1_group_manifest_D0_D6.csv``.  Group labels come from the primary
Methods-faithful rarefied branch, whereas all downstream panels are rebuilt
with the source-faithful Bracken branch used by the manuscript restoration.
This deliberate branch boundary is recorded in every definition manifest.

Figure 2B is the sole visual substitution: it displays the definition-specific
ranking and 15/30/15 assignment.  Every other Figure 2/4 panel is recalculated
after replacing the subject labels; no D0 panel-level result is copied.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Callable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-r11")

import numpy as np
import pandas as pd
from PIL import Image


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
REPRODUCTION = PROJECT_ROOT / "reproduction"
sys.path.insert(0, str(REPRODUCTION))

import lefse_reconstruction as lefse  # noqa: E402
import fastspar_fig4 as fastspar_pipeline  # noqa: E402
import load_raw_bracken as raw_loader  # noqa: E402
import render_fig4_sfig4 as fig4_renderer  # noqa: E402
import render_high_fidelity as hifi  # noqa: E402
import render_selected_panels as selected_panels  # noqa: E402
import run_all as audit  # noqa: E402


GROUP_MANIFEST = HERE / "R1_1_group_manifest_D0_D6.csv"
DEFAULT_OUTPUT = HERE / "full_original_layout"
DEFINITION_SLUGS = {
    "D1": "D1_endpoint_abs_shannon",
    "D2": "D2_endpoint_bray",
    "D3": "D3_endpoint_jaccard",
    "D4": "D4_auc_abs_shannon",
    "D5": "D5_auc_bray",
    "D6": "D6_auc_jaccard",
}
PANEL_B_LABELS = {
    "D1": "Endpoint |ΔShannon|",
    "D2": "Endpoint Bray–Curtis",
    "D3": "Endpoint Jaccard",
    "D4": "AUC |ΔShannon|",
    "D5": "AUC Bray–Curtis",
    "D6": "AUC Jaccard",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_dump(payload: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _mapping_sha256(frame: pd.DataFrame) -> str:
    """Hash a normalized subject/group mapping rather than file formatting."""
    normalized = (
        frame.loc[:, ["subject", "recovery_group"]]
        .assign(
            subject=lambda value: value["subject"].astype(str),
            recovery_group=lambda value: value["recovery_group"].astype(str),
        )
        .sort_values("subject", kind="mergesort")
        .reset_index(drop=True)
    )
    return hashlib.sha256(
        normalized.to_csv(index=False, lineterminator="\n").encode("utf-8")
    ).hexdigest()


def validate_d0_grouping_regression(
    base_dataset: audit.Dataset,
    group_manifest: Path,
    frozen_lefse_features: pd.Index,
    output_path: Path,
) -> dict[str, object]:
    """Fail before rendering if the frozen reference grouping has drifted."""
    checks: list[dict[str, object]] = []

    def record(name: str, passed: bool, details: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "details": details})

    d0 = load_definition_table(group_manifest, "D0")
    manifest_mapping = d0.loc[:, ["subject", "recovery_group"]].copy()
    base_mapping = (
        base_dataset.group_table.reset_index()
        .loc[:, ["subject", "recovery_group"]]
        .copy()
    )
    manifest_mapping["subject"] = manifest_mapping["subject"].astype(str)
    base_mapping["subject"] = base_mapping["subject"].astype(str)
    manifest_mapping["recovery_group"] = manifest_mapping["recovery_group"].astype(str)
    base_mapping["recovery_group"] = base_mapping["recovery_group"].astype(str)

    manifest_subjects = set(manifest_mapping["subject"])
    base_subjects = set(base_mapping["subject"])
    record(
        "D0_subject_set_exact",
        manifest_subjects == base_subjects and len(manifest_subjects) == 60,
        {
            "manifest_n": len(manifest_subjects),
            "base_n": len(base_subjects),
            "missing_from_manifest": sorted(base_subjects - manifest_subjects),
            "extra_in_manifest": sorted(manifest_subjects - base_subjects),
        },
    )
    comparison = manifest_mapping.merge(
        base_mapping,
        on="subject",
        how="outer",
        suffixes=("_manifest", "_base"),
        indicator=True,
    )
    mismatches = comparison.loc[
        comparison["recovery_group_manifest"].ne(comparison["recovery_group_base"])
        | comparison["_merge"].ne("both")
    ]
    record(
        "D0_subject_to_group_exact",
        mismatches.empty,
        {
            "n_mismatches": int(len(mismatches)),
            "mismatches": mismatches.astype(object).where(pd.notna(mismatches), None).to_dict("records"),
        },
    )
    expected_counts = {"HDR": 15, "LDR": 15, "Middle": 30}
    manifest_counts = manifest_mapping["recovery_group"].value_counts().to_dict()
    base_counts = base_mapping["recovery_group"].value_counts().to_dict()
    record(
        "D0_group_counts_15_15_30",
        manifest_counts == expected_counts and base_counts == expected_counts,
        {"expected": expected_counts, "manifest": manifest_counts, "base": base_counts},
    )
    feature_values = pd.Index(frozen_lefse_features).astype(str)
    feature_ok = (
        len(feature_values) == 268
        and feature_values.nunique() == 268
        and not feature_values.isna().any()
    )
    record(
        "D0_frozen_LEfSe_candidate_universe_268",
        feature_ok,
        {
            "n_features": int(len(feature_values)),
            "n_unique": int(feature_values.nunique()),
        },
    )

    passed = all(bool(item["passed"]) for item in checks)
    payload = {
        "status": "PASS" if passed else "FAIL",
        "purpose": (
            "Pre-render regression gate: the raw-Bracken reconstruction must retain "
            "the exact D0 participant grouping encoded in the fixed R1.1 manifest."
        ),
        "group_manifest": str(group_manifest.resolve()),
        "group_manifest_sha256": sha256(group_manifest),
        "manifest_mapping_sha256": _mapping_sha256(manifest_mapping),
        "base_dataset_mapping_sha256": _mapping_sha256(base_mapping),
        "checks": checks,
    }
    json_dump(payload, output_path)
    if not passed:
        failed = [str(item["check"]) for item in checks if not item["passed"]]
        raise AssertionError(f"D0 grouping regression failed: {failed}; see {output_path}")
    return payload


def load_definition_table(path: Path, definition_id: str) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype={"subject": str, "definition_id": str})
    required = {
        "subject",
        "definition_id",
        "definition_short",
        "definition_full",
        "metric_column",
        "high_recovery_direction",
        "metric_value",
        "metric_rank_ascending",
        "recovery_group",
        "dataset_branch",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Grouping manifest lacks columns: {sorted(missing)}")
    selected = frame.loc[frame["definition_id"].eq(definition_id)].copy()
    if len(selected) != 60 or selected["subject"].nunique() != 60:
        raise ValueError(f"{definition_id}: expected 60 unique subjects, found {len(selected)} rows")
    if selected["metric_value"].isna().any():
        raise ValueError(f"{definition_id}: grouping metric contains missing values")
    counts = selected["recovery_group"].value_counts().to_dict()
    if counts != {"Middle": 30, "HDR": 15, "LDR": 15}:
        raise ValueError(f"{definition_id}: invalid group sizes {counts}")
    if selected["high_recovery_direction"].nunique() != 1:
        raise ValueError(f"{definition_id}: recovery direction is not unique")
    direction = str(selected["high_recovery_direction"].iloc[0])
    if direction not in {"larger", "smaller"}:
        raise ValueError(f"{definition_id}: invalid recovery direction {direction!r}")

    ordered = selected.sort_values(["metric_value", "subject"], kind="mergesort")
    values = ordered["metric_value"].to_numpy(dtype=float)
    lower_tie = values[14] == values[15]
    upper_tie = values[44] == values[45]
    if lower_tie or upper_tie:
        raise ValueError(
            f"{definition_id}: exact boundary tie detected "
            f"(lower={lower_tie}, upper={upper_tie}); fixed n=15 groups are not unique"
        )
    expected_hdr = set(
        ordered.iloc[:15]["subject"]
        if direction == "smaller"
        else ordered.iloc[45:]["subject"]
    )
    expected_ldr = set(
        ordered.iloc[45:]["subject"]
        if direction == "smaller"
        else ordered.iloc[:15]["subject"]
    )
    observed_hdr = set(selected.loc[selected["recovery_group"].eq("HDR"), "subject"])
    observed_ldr = set(selected.loc[selected["recovery_group"].eq("LDR"), "subject"])
    if expected_hdr != observed_hdr or expected_ldr != observed_ldr:
        raise ValueError(f"{definition_id}: manifest labels do not match its full-precision ranking")
    return selected


def inject_definition(ds: audit.Dataset, definition: pd.DataFrame) -> audit.Dataset:
    definition_id = str(definition["definition_id"].iloc[0])
    lookup = definition.set_index("subject")
    group_table = ds.group_table.copy()
    group_table.index = group_table.index.astype(str)
    missing = sorted(set(group_table.index).difference(lookup.index))
    extra = sorted(set(lookup.index).difference(group_table.index))
    if missing or extra:
        raise ValueError(
            f"{definition_id}: subject mismatch; missing={missing[:5]}, extra={extra[:5]}"
        )
    aligned = lookup.reindex(group_table.index)
    group_table["recovery_group"] = aligned["recovery_group"].astype(str)
    group_table["alternative_definition_id"] = definition_id
    group_table["alternative_metric_label"] = PANEL_B_LABELS[definition_id]
    group_table["alternative_high_recovery_direction"] = aligned[
        "high_recovery_direction"
    ].astype(str)
    group_table["alternative_metric_value"] = aligned["metric_value"].astype(float)
    group_table["alternative_metric_rank_ascending"] = aligned[
        "metric_rank_ascending"
    ].astype(int)
    group_table["subject_sort_key"] = group_table.index.astype(str)

    metrics = ds.metrics.drop(columns=["recovery_group"], errors="ignore").copy()
    metrics["subject"] = metrics["subject"].astype(str)
    metrics = metrics.merge(
        group_table[["recovery_group"]],
        left_on="subject",
        right_index=True,
        how="left",
        validate="many_to_one",
    )
    metrics["recovery_group"] = metrics["recovery_group"].fillna("Control")
    selected_counts = (
        metrics.loc[metrics["recovery_group"].isin(["HDR", "LDR"])]
        .groupby("recovery_group")["subject"]
        .nunique()
        .to_dict()
    )
    if selected_counts != {"HDR": 15, "LDR": 15}:
        raise ValueError(f"{definition_id}: injected sample table has invalid counts {selected_counts}")
    return replace(ds, metrics=metrics, group_table=group_table)


def make_lefse_without_source_calibration(
    frozen_features: pd.Index,
) -> Callable[[audit.Dataset, Path], tuple[pd.DataFrame, pd.DataFrame]]:
    def run(
        ds: audit.Dataset, tables: Path
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        return lefse.figure2f_tables(
            ds,
            tables,
            source_digitization_audit=False,
            frozen_features=frozen_features,
        )

    return run


def make_fig2_export(definition_id: str) -> Callable[..., dict[str, object]]:
    def export(figure, stem: Path, dpi: int = 600) -> dict[str, object]:
        if stem.name == "Figure2_high_fidelity":
            target = stem.parent / f"R1_1_{definition_id}_Fig2"
        else:
            target = stem
        record = selected_panels.export_pdf_tiff(figure, target, dpi=dpi)
        record["canvas_policy"] = "fixed manuscript page; bbox_inches=tight disabled"
        return record

    return export


def definition_audit(definition: pd.DataFrame) -> dict[str, object]:
    row = definition.iloc[0]
    ordered = definition.sort_values(["metric_value", "subject"], kind="mergesort")
    values = ordered["metric_value"].to_numpy(dtype=float)
    return {
        "definition_id": str(row["definition_id"]),
        "definition_short": str(row["definition_short"]),
        "definition_full": str(row["definition_full"]),
        "metric_column": str(row["metric_column"]),
        "panel_B_axis_label": PANEL_B_LABELS[str(row["definition_id"])],
        "high_recovery_direction": str(row["high_recovery_direction"]),
        "group_rule": "all 60 ranked at full precision; best 15 HDR, worst 15 LDR, middle 30 not analysed",
        "group_counts": definition["recovery_group"].value_counts().to_dict(),
        "boundary_values_ascending": {
            "rank15": float(values[14]),
            "rank16": float(values[15]),
            "rank45": float(values[44]),
            "rank46": float(values[45]),
        },
        "exact_boundary_tie": bool(values[14] == values[15] or values[44] == values[45]),
        "grouping_branch": str(row["dataset_branch"]),
        "downstream_panel_branch": "source-faithful unrarefied bacterial Bracken leaves",
        "interpretation_boundary": (
            "Figure 2B is definitional. Robustness is judged from the remaining "
            "cross-domain panels and recurrence across independently defined groups."
        ),
    }


def build_fastspar_argv(
    data_root: Path,
    group_manifest: Path,
    definition_id: str,
    output_dir: Path,
    panel_jobs: int,
    jobs: int,
) -> list[str]:
    """Build a parser-testable FastSpar command without shell interpolation."""
    return [
        "--data-root",
        str(data_root),
        "--group-manifest",
        str(group_manifest),
        "--definition-id",
        definition_id,
        "--output-dir",
        str(output_dir),
        "--features",
        "275",
        "--bootstraps",
        "500",
        "--iterations",
        "50",
        "--exclusion-iterations",
        "10",
        "--panel-jobs",
        str(panel_jobs),
        "--jobs",
        str(jobs),
        "--threads-per-fit",
        "1",
        "--bootstrap-threads",
        "8",
        "--pvalue-threads",
        "8",
        "--fastspar-bin-dir",
        "/home/qqt/.conda/envs/fastspar_env/bin",
    ]


def _inspect_figure_export(
    definition_id: str,
    figure_name: str,
    figures_dir: Path,
    qa_dir: Path,
    checks: list[dict[str, object]],
) -> dict[str, object]:
    """Validate a PDF/TIFF pair and create a Python/PIL QA thumbnail."""

    def record(name: str, passed: bool, details: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "details": details})

    stem = figures_dir / f"R1_1_{definition_id}_{figure_name}"
    pdf_path = stem.with_suffix(".pdf")
    tif_path = stem.with_suffix(".tif")
    pdf_ok = pdf_path.is_file() and pdf_path.stat().st_size > 0
    tif_ok = tif_path.is_file() and tif_path.stat().st_size > 0
    record(
        f"{figure_name}_PDF_exists_nonempty",
        pdf_ok,
        {"path": str(pdf_path.resolve()), "bytes": pdf_path.stat().st_size if pdf_ok else 0},
    )
    record(
        f"{figure_name}_TIFF_exists_nonempty",
        tif_ok,
        {"path": str(tif_path.resolve()), "bytes": tif_path.stat().st_size if tif_ok else 0},
    )
    if not pdf_ok or not tif_ok:
        return {
            "figure": figure_name,
            "pdf_path": str(pdf_path.resolve()),
            "tif_path": str(tif_path.resolve()),
            "status": "FAIL",
        }

    with pdf_path.open("rb") as handle:
        pdf_signature = handle.read(5)
    record(
        f"{figure_name}_PDF_signature",
        pdf_signature == b"%PDF-",
        {"signature": pdf_signature.decode("latin-1", errors="replace")},
    )

    expected_pixels = {"Fig2": (4882, 5668), "Fig4": (4722, 6795)}
    with Image.open(tif_path) as image:
        image.load()
        width, height = image.size
        mode = image.mode
        compression = str(image.info.get("compression", ""))
        dpi_value = image.info.get("dpi", (float("nan"), float("nan")))
        if isinstance(dpi_value, (int, float)):
            dpi_x = dpi_y = float(dpi_value)
        else:
            dpi_x = float(dpi_value[0])
            dpi_y = float(dpi_value[1])
        preview = image.convert("RGB")

    record(f"{figure_name}_TIFF_RGB", mode == "RGB", {"mode": mode})
    record(
        f"{figure_name}_TIFF_600dpi",
        abs(dpi_x - 600.0) <= 1.0 and abs(dpi_y - 600.0) <= 1.0,
        {"dpi_x": dpi_x, "dpi_y": dpi_y},
    )
    record(
        f"{figure_name}_TIFF_LZW",
        compression.lower() in {"tiff_lzw", "lzw"},
        {"compression": compression},
    )
    expected_width, expected_height = expected_pixels[figure_name]
    dimensions_ok = (
        abs(width - expected_width) / expected_width <= 0.02
        and abs(height - expected_height) / expected_height <= 0.02
        and 3_000 <= min(width, height)
        and max(width, height) <= 10_000
    )
    record(
        f"{figure_name}_TIFF_reasonable_fixed_canvas",
        dimensions_ok,
        {
            "width_px": int(width),
            "height_px": int(height),
            "nominal_width_px": expected_width,
            "nominal_height_px": expected_height,
            "tolerance": "2% per dimension",
        },
    )

    preview.thumbnail((1_400, 1_600), Image.Resampling.LANCZOS)
    thumbnail = qa_dir / f"R1_1_{definition_id}_{figure_name}_QA_thumbnail.png"
    preview.save(thumbnail, format="PNG", optimize=True)
    record(
        f"{figure_name}_PIL_thumbnail_written",
        thumbnail.is_file() and thumbnail.stat().st_size > 0,
        {
            "path": str(thumbnail.resolve()),
            "width_px": int(preview.width),
            "height_px": int(preview.height),
        },
    )
    return {
        "figure": figure_name,
        "status": "PASS",
        "pdf_path": str(pdf_path.resolve()),
        "pdf_bytes": int(pdf_path.stat().st_size),
        "pdf_sha256": sha256(pdf_path),
        "tif_path": str(tif_path.resolve()),
        "tif_bytes": int(tif_path.stat().st_size),
        "tif_sha256": sha256(tif_path),
        "mode": mode,
        "compression": compression,
        "dpi_x": dpi_x,
        "dpi_y": dpi_y,
        "width_px": int(width),
        "height_px": int(height),
        "thumbnail_path": str(thumbnail.resolve()),
        "thumbnail_bytes": int(thumbnail.stat().st_size),
    }


def validate_definition_outputs(
    definition_id: str,
    definition: pd.DataFrame,
    dataset: audit.Dataset,
    definition_dir: Path,
    fig2f_mode: str,
    fig4c_mode: str,
    completed_fastspar: Any | None,
) -> dict[str, object]:
    """Run machine-readable export, matrix, grouping, and provenance gates."""
    figures = definition_dir / "figures"
    tables = definition_dir / "tables"
    qa_dir = definition_dir / "qa"
    qa_dir.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, object]] = []

    def record(name: str, passed: bool, details: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "details": details})

    figure_exports = [
        _inspect_figure_export(definition_id, figure_name, figures, qa_dir, checks)
        for figure_name in ("Fig2", "Fig4")
    ]

    matrix_path = (
        tables / "Fig2F_LEfSe_LDA_matrix.csv"
        if fig2f_mode == "lefse"
        else tables / "figure2_effect_heatmap_high_fidelity.csv"
    )
    matrix_exists = matrix_path.is_file() and matrix_path.stat().st_size > 0
    record(
        "Fig2F_matrix_exists_nonempty",
        matrix_exists,
        {"path": str(matrix_path.resolve()), "bytes": matrix_path.stat().st_size if matrix_exists else 0},
    )
    if matrix_exists:
        matrix_frame = pd.read_csv(matrix_path)
        matrix_values = matrix_frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        finite = bool(np.isfinite(matrix_values.to_numpy(dtype=float)).all())
        matrix_ok = matrix_values.shape == (16, 45) and finite
        record(
            "Fig2F_matrix_16_by_45_no_missing_or_nonfinite",
            matrix_ok,
            {
                "matrix_rows": int(matrix_values.shape[0]),
                "matrix_columns": int(matrix_values.shape[1]),
                "missing_cells": int(matrix_values.isna().sum().sum()),
                "all_finite": finite,
            },
        )
    else:
        record(
            "Fig2F_matrix_16_by_45_no_missing_or_nonfinite",
            False,
            "matrix unavailable",
        )

    if fig2f_mode == "lefse":
        candidate_path = tables / "Fig2F_LEfSe_candidate_universe.csv"
        candidate_exists = candidate_path.is_file() and candidate_path.stat().st_size > 0
        if candidate_exists:
            candidate_frame = pd.read_csv(candidate_path)
            candidate_columns_ok = {"species", "feature_rule"}.issubset(
                candidate_frame.columns
            )
            if candidate_columns_ok:
                species = candidate_frame["species"].astype("string")
                feature_rules = candidate_frame["feature_rule"].astype("string")
                candidate_ok = (
                    len(candidate_frame) == 268
                    and species.notna().all()
                    and species.nunique(dropna=True) == 268
                    and feature_rules.notna().all()
                    and feature_rules.eq("frozen D0 candidate universe").all()
                )
                candidate_details: object = {
                    "path": str(candidate_path.resolve()),
                    "n_rows": int(len(candidate_frame)),
                    "n_unique_species": int(species.nunique(dropna=True)),
                    "n_missing_species": int(species.isna().sum()),
                    "feature_rule_values": sorted(
                        feature_rules.dropna().astype(str).unique().tolist()
                    ),
                }
            else:
                candidate_ok = False
                candidate_details = {
                    "path": str(candidate_path.resolve()),
                    "columns": candidate_frame.columns.astype(str).tolist(),
                    "missing_columns": sorted(
                        {"species", "feature_rule"}.difference(candidate_frame.columns)
                    ),
                }
        else:
            candidate_ok = False
            candidate_details = {
                "path": str(candidate_path.resolve()),
                "reason": "missing or empty",
            }
        record(
            "Fig2F_LEfSe_frozen_D0_candidate_universe_exact_268",
            candidate_ok,
            candidate_details,
        )

        parameter_path = tables / "Fig2F_LEfSe_parameters.csv"
        parameter_exists = parameter_path.is_file() and parameter_path.stat().st_size > 0
        if parameter_exists:
            parameter_frame = pd.read_csv(parameter_path, dtype=str)
            parameter_columns_ok = {"parameter", "value"}.issubset(
                parameter_frame.columns
            )
            required_parameters = {
                "LEfSe_KW_alpha_for_no_filter",
                "LDA_run_threshold",
                "cell_mask_rule",
            }
            if parameter_columns_ok:
                parameter_counts = parameter_frame["parameter"].value_counts()
                parameters_unique = all(
                    int(parameter_counts.get(name, 0)) == 1
                    for name in required_parameters
                )
                values = (
                    parameter_frame.set_index("parameter")["value"]
                    if parameters_unique
                    else pd.Series(dtype=str)
                )
                try:
                    kw_alpha = float(values["LEfSe_KW_alpha_for_no_filter"])
                    lda_threshold = float(values["LDA_run_threshold"])
                    cell_mask_raw = str(values["cell_mask_rule"]).strip()
                    # The exported value may retain an explanatory clause after a
                    # semicolon; the machine-readable rule token itself must be none.
                    cell_mask_token = cell_mask_raw.split(";", maxsplit=1)[0].strip().casefold()
                    parameter_ok = (
                        np.isfinite(kw_alpha)
                        and kw_alpha > 1.0
                        and np.isfinite(lda_threshold)
                        and lda_threshold == 0.0
                        and cell_mask_token == "none"
                    )
                    parameter_details: object = {
                        "path": str(parameter_path.resolve()),
                        "LEfSe_KW_alpha_for_no_filter": kw_alpha,
                        "LDA_run_threshold": lda_threshold,
                        "cell_mask_rule_raw": cell_mask_raw,
                        "cell_mask_rule_token": cell_mask_token,
                        "required_parameter_rows_unique": parameters_unique,
                    }
                except (KeyError, TypeError, ValueError):
                    parameter_ok = False
                    parameter_details = {
                        "path": str(parameter_path.resolve()),
                        "required_parameter_rows_unique": parameters_unique,
                        "available_parameters": sorted(
                            parameter_frame["parameter"].dropna().astype(str).tolist()
                        ),
                    }
            else:
                parameter_ok = False
                parameter_details = {
                    "path": str(parameter_path.resolve()),
                    "columns": parameter_frame.columns.astype(str).tolist(),
                    "missing_columns": sorted(
                        {"parameter", "value"}.difference(parameter_frame.columns)
                    ),
                }
        else:
            parameter_ok = False
            parameter_details = {
                "path": str(parameter_path.resolve()),
                "reason": "missing or empty",
            }
        record(
            "Fig2F_LEfSe_no_filter_and_no_cell_mask_parameters",
            parameter_ok,
            parameter_details,
        )

        eligible_path = tables / "Fig2F_LEfSe_all_eligible_species.csv"
        eligible_exists = eligible_path.is_file() and eligible_path.stat().st_size > 0
        if eligible_exists:
            eligible_frame = pd.read_csv(eligible_path)
            eligible_required = {
                "canon",
                "species",
                "lda_abs",
                "lda_score_LDR_over_HDR",
            }
            eligible_columns_ok = eligible_required.issubset(eligible_frame.columns)
            if eligible_columns_ok:
                lda_values = eligible_frame.loc[
                    :, ["lda_abs", "lda_score_LDR_over_HDR"]
                ].apply(pd.to_numeric, errors="coerce")
                pair_duplicates = int(
                    eligible_frame.duplicated(subset=["canon", "species"]).sum()
                )
                rows_per_canon = eligible_frame.groupby("canon", dropna=False).size()
                rows_per_species = eligible_frame.groupby("species", dropna=False).size()
                all_finite = bool(
                    np.isfinite(lda_values.to_numpy(dtype=float)).all()
                )
                eligible_ok = (
                    len(eligible_frame) == 16 * 268
                    and eligible_frame["canon"].nunique(dropna=True) == 16
                    and eligible_frame["species"].nunique(dropna=True) == 268
                    and pair_duplicates == 0
                    and (rows_per_canon == 268).all()
                    and (rows_per_species == 16).all()
                    and all_finite
                )
                eligible_details: object = {
                    "path": str(eligible_path.resolve()),
                    "n_rows": int(len(eligible_frame)),
                    "expected_rows": int(16 * 268),
                    "n_unique_canons": int(
                        eligible_frame["canon"].nunique(dropna=True)
                    ),
                    "n_unique_species": int(
                        eligible_frame["species"].nunique(dropna=True)
                    ),
                    "duplicate_canon_species_pairs": pair_duplicates,
                    "rows_per_canon_values": sorted(
                        {int(value) for value in rows_per_canon.tolist()}
                    ),
                    "rows_per_species_values": sorted(
                        {int(value) for value in rows_per_species.tolist()}
                    ),
                    "missing_or_nonfinite_lda_cells": int(
                        (~np.isfinite(lda_values.to_numpy(dtype=float))).sum()
                    ),
                    "all_lda_values_finite": all_finite,
                }
            else:
                eligible_ok = False
                eligible_details = {
                    "path": str(eligible_path.resolve()),
                    "columns": eligible_frame.columns.astype(str).tolist(),
                    "missing_columns": sorted(
                        eligible_required.difference(eligible_frame.columns)
                    ),
                }
        else:
            eligible_ok = False
            eligible_details = {
                "path": str(eligible_path.resolve()),
                "reason": "missing or empty",
            }
        record(
            "Fig2F_LEfSe_all_eligible_complete_16_by_268_and_finite",
            eligible_ok,
            eligible_details,
        )

    expected_counts = {"HDR": 15, "LDR": 15, "Middle": 30}
    manifest_counts = {
        str(key): int(value)
        for key, value in definition["recovery_group"].value_counts().to_dict().items()
    }
    dataset_mapping = (
        dataset.group_table.reset_index()
        .loc[:, ["subject", "recovery_group"]]
        .assign(subject=lambda value: value["subject"].astype(str))
    )
    dataset_counts = {
        str(key): int(value)
        for key, value in dataset_mapping["recovery_group"].value_counts().to_dict().items()
    }
    record(
        "grouping_60_total_15_HDR_15_LDR_30_Middle",
        len(definition) == 60
        and definition["subject"].nunique() == 60
        and manifest_counts == expected_counts
        and dataset_counts == expected_counts,
        {
            "n_manifest_rows": int(len(definition)),
            "n_manifest_subjects": int(definition["subject"].nunique()),
            "manifest_counts": manifest_counts,
            "injected_dataset_counts": dataset_counts,
            "expected_counts": expected_counts,
        },
    )
    manifest_mapping = definition.loc[:, ["subject", "recovery_group"]].copy()
    record(
        "injected_subject_to_group_exact",
        _mapping_sha256(manifest_mapping) == _mapping_sha256(dataset_mapping),
        {
            "manifest_mapping_sha256": _mapping_sha256(manifest_mapping),
            "injected_mapping_sha256": _mapping_sha256(dataset_mapping),
        },
    )

    if fig4c_mode == "fastspar":
        if completed_fastspar is None:
            record("FastSpar_complete_12_panels", False, "completed result object unavailable")
            record("FastSpar_group_assignment_exact", False, "completed result object unavailable")
        else:
            summary = completed_fastspar.summary.copy()
            observed_panels = {
                f"{str(row.group)}__{str(row.stage)}"
                for row in summary.itertuples(index=False)
            }
            expected_panels = {
                f"{group}__{stage}"
                for group in fastspar_pipeline.GROUP_ORDER
                for stage in fastspar_pipeline.NETWORK_STAGES
            }
            bootstrap_values = sorted(
                {int(value) for value in summary["bootstraps"].tolist()}
            )
            complete = (
                len(summary) == 12
                and observed_panels == expected_panels
                and len(completed_fastspar.panel_summaries) == 12
                and bootstrap_values == [500]
            )
            record(
                "FastSpar_complete_12_panels",
                complete,
                {
                    "n_summary_rows": int(len(summary)),
                    "n_panel_summaries": int(len(completed_fastspar.panel_summaries)),
                    "missing_panels": sorted(expected_panels - observed_panels),
                    "extra_panels": sorted(observed_panels - expected_panels),
                    "bootstrap_values": bootstrap_values,
                },
            )
            assignment_path = definition_dir / "fastspar" / "group_assignments.tsv"
            if assignment_path.is_file():
                fastspar_mapping = pd.read_csv(
                    assignment_path, sep="\t", dtype={"subject": str}
                ).loc[:, ["subject", "recovery_group"]]
                assignment_exact = _mapping_sha256(fastspar_mapping) == _mapping_sha256(
                    manifest_mapping
                )
                recorded_hash = str(
                    completed_fastspar.manifest.get("group_assignment_table_sha256", "")
                )
                file_hash = sha256(assignment_path)
                manifest_identity_ok = (
                    str(completed_fastspar.manifest.get("grouping_definition_id"))
                    == definition_id
                )
                hash_ok = recorded_hash == file_hash
            else:
                assignment_exact = False
                manifest_identity_ok = False
                hash_ok = False
                recorded_hash = ""
                file_hash = ""
            record(
                "FastSpar_group_assignment_exact",
                assignment_exact and manifest_identity_ok and hash_ok,
                {
                    "path": str(assignment_path.resolve()),
                    "mapping_matches_definition": assignment_exact,
                    "manifest_definition_matches": manifest_identity_ok,
                    "recorded_file_sha256": recorded_hash,
                    "observed_file_sha256": file_hash,
                    "file_hash_matches_manifest": hash_ok,
                },
            )
    else:
        record(
            "FastSpar_complete_12_panels",
            True,
            "SKIPPED: --fig4c-mode proxy is diagnostic-only; finalized runs require fastspar.",
        )
        record(
            "FastSpar_group_assignment_exact",
            True,
            "SKIPPED: --fig4c-mode proxy is diagnostic-only; injected grouping was checked above.",
        )

    passed = all(bool(item["passed"]) for item in checks)
    payload = {
        "definition_id": definition_id,
        "status": "PASS" if passed else "FAIL",
        "fig2f_mode": fig2f_mode,
        "fig4c_mode": fig4c_mode,
        "figure_exports": figure_exports,
        "checks": checks,
    }
    qa_path = qa_dir / "figure_export_QA.json"
    json_dump(payload, qa_path)
    if not passed:
        failed = [str(item["check"]) for item in checks if not item["passed"]]
        raise AssertionError(f"{definition_id}: output QA failed: {failed}; see {qa_path}")
    return payload


def render_definition(
    base_dataset: audit.Dataset,
    definition: pd.DataFrame,
    project_root: Path,
    output_root: Path,
    fig2f_mode: str,
    frozen_lefse_features: pd.Index,
    group_manifest: Path,
    fig4c_mode: str,
    fastspar_panel_jobs: int,
    fastspar_jobs: int,
) -> dict[str, object]:
    definition_id = str(definition["definition_id"].iloc[0])
    definition_dir = output_root / DEFINITION_SLUGS[definition_id]
    figures = definition_dir / "figures"
    tables = definition_dir / "tables"
    qa = definition_dir / "qa"
    for path in (figures, tables, qa):
        path.mkdir(parents=True, exist_ok=True)

    ds = inject_definition(base_dataset, definition)
    ds.group_table.to_csv(tables / "subject_group_assignments_and_panelB_values.csv")
    ds.metrics.sort_values(["subject", "canons"]).to_csv(
        tables / "sample_metrics_with_alternative_groups.csv", index=False
    )

    hifi.set_hifi_style()
    hifi.export = make_fig2_export(definition_id)
    selected_panels.figure2f_tables = (
        make_lefse_without_source_calibration(frozen_lefse_features)
        if fig2f_mode == "lefse"
        else selected_panels.__dict__["_r11_proxy_figure2f_tables"]
    )
    fig2_record = hifi.figure2_hifi(ds, figures, tables)

    completed_fastspar = None
    if fig4c_mode == "fastspar":
        fastspar_output = definition_dir / "fastspar"
        returncode = fastspar_pipeline.main(
            build_fastspar_argv(
                project_root,
                group_manifest,
                definition_id,
                fastspar_output,
                fastspar_panel_jobs,
                fastspar_jobs,
            )
        )
        if returncode != 0:
            raise RuntimeError(f"{definition_id}: FastSpar returned {returncode}")
        completed_fastspar = fastspar_pipeline.load_completed_results(
            fastspar_output / "results"
        )

    fig4_renderer.set_style()
    definition_full = str(definition["definition_full"].iloc[0])
    fig4_record = fig4_renderer.draw_figure4(
        ds,
        figures,
        tables,
        fastspar_results=completed_fastspar,
        group_definition=(
            f"{definition_id}: {definition_full}; full-precision ranking, "
            "15 HDR and 15 LDR"
        ),
        summary_uses_full_eligible_network=True,
    )
    fig4_source = figures / "Fig4_restored"
    fig4_targets = {}
    for suffix in (".pdf", ".tif"):
        source = fig4_source.with_suffix(suffix)
        target = figures / f"R1_1_{definition_id}_Fig4{suffix}"
        source.replace(target)
        fig4_targets[suffix] = target

    qa_payload = validate_definition_outputs(
        definition_id,
        definition,
        ds,
        definition_dir,
        fig2f_mode,
        fig4c_mode,
        completed_fastspar,
    )

    audit_payload = definition_audit(definition)
    audit_payload.update(
        {
            "figure2f_mode": fig2f_mode,
            "figure2f_method": (
                "LEfSe 1.1.2/MASS LDA reconstruction on the fixed source species panel"
                if fig2f_mode == "lefse"
                else "rank-biserial/Mann-Whitney fallback"
            ),
            "figure4c_method": (
                "FastSpar with 500 bootstraps on the frozen global top-275 genus universe; "
                "recomputed for this definition, not copied from D0"
                if fig4c_mode == "fastspar"
                else "original renderer's fast subject-mean Spearman network fallback; "
                "recomputed for this definition, not copied from D0"
            ),
            "figure4c_mode": fig4c_mode,
            "quality_assurance": qa_payload,
            "figure_records": [fig2_record, fig4_record],
            "outputs": [
                {
                    "path": str(path.resolve()),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in sorted(figures.glob("R1_1_*_Fig[24].*"))
                if path.suffix in {".pdf", ".tif"}
            ],
        }
    )
    json_dump(audit_payload, definition_dir / "definition.json")
    return audit_payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--group-manifest", type=Path, default=GROUP_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--definitions",
        nargs="+",
        choices=tuple(DEFINITION_SLUGS),
        default=list(DEFINITION_SLUGS),
    )
    parser.add_argument(
        "--fig2f-mode",
        choices=("lefse", "proxy"),
        default="lefse",
        help="LEfSe matches the finalized source restoration; proxy is a fast diagnostic fallback.",
    )
    parser.add_argument(
        "--fig4c-mode",
        choices=("fastspar", "proxy"),
        default="fastspar",
    )
    parser.add_argument("--fastspar-panel-jobs", type=int, default=4)
    parser.add_argument("--fastspar-jobs", type=int, default=8)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    project_root = args.project_root.resolve()
    group_manifest = args.group_manifest.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Preserve the genuine proxy function before replacing the public hook.
    if "_r11_proxy_figure2f_tables" not in selected_panels.__dict__:
        selected_panels.__dict__["_r11_proxy_figure2f_tables"] = (
            selected_panels.figure2f_tables
        )
    dataset = raw_loader.load_raw_bracken_dataset(project_root)
    frozen_lefse_features = lefse.eligible_feature_universe(dataset)
    d0_regression = validate_d0_grouping_regression(
        dataset,
        group_manifest,
        frozen_lefse_features,
        output / "D0_grouping_regression.json",
    )
    pd.DataFrame({"species": frozen_lefse_features}).to_csv(
        output / "D0_frozen_Fig2F_LEfSe_candidate_universe.csv", index=False
    )
    records = []
    for definition_id in args.definitions:
        definition = load_definition_table(group_manifest, definition_id)
        records.append(
            render_definition(
                dataset,
                definition,
                project_root,
                output,
                args.fig2f_mode,
                frozen_lefse_features,
                group_manifest,
                args.fig4c_mode,
                args.fastspar_panel_jobs,
                args.fastspar_jobs,
            )
        )

    index_rows = []
    for record in records:
        definition_id = str(record["definition_id"])
        definition_dir = output / DEFINITION_SLUGS[definition_id]
        fig3_pdf = definition_dir / "figures" / "Fig3_restored.pdf"
        fig3_tif = definition_dir / "figures" / "Fig3_restored.tif"
        fig3_marker = definition_dir / "mdsine2" / "finalized.json"
        if fig3_marker.is_file() and fig3_pdf.is_file() and fig3_tif.is_file():
            fig3_status = "validated developmental"
        elif fig3_marker.is_file():
            fig3_status = "finalized marker present but export missing"
        else:
            fig3_status = "pending rapid MDSINE2"
        index_rows.append(
            {
                "definition_id": definition_id,
                "definition_short": record["definition_short"],
                "directory": str(definition_dir.resolve()),
                "fig2_pdf": str((definition_dir / "figures" / f"R1_1_{definition_id}_Fig2.pdf").resolve()),
                "fig2_tif": str((definition_dir / "figures" / f"R1_1_{definition_id}_Fig2.tif").resolve()),
                "fig3_pdf": str(fig3_pdf.resolve()) if fig3_pdf.is_file() else "",
                "fig3_tif": str(fig3_tif.resolve()) if fig3_tif.is_file() else "",
                "fig3_status": fig3_status,
                "fig4_pdf": str((definition_dir / "figures" / f"R1_1_{definition_id}_Fig4.pdf").resolve()),
                "fig4_tif": str((definition_dir / "figures" / f"R1_1_{definition_id}_Fig4.tif").resolve()),
            }
        )
    pd.DataFrame(index_rows).to_csv(output / "reviewer_ready_index.csv", index=False)

    visual_qa_rows = []
    for record in records:
        qa_payload = record["quality_assurance"]
        for figure_export in qa_payload["figure_exports"]:
            visual_qa_rows.append(
                {
                    "definition_id": record["definition_id"],
                    "definition_short": record["definition_short"],
                    "definition_QA_status": qa_payload["status"],
                    "figure": figure_export["figure"],
                    "figure_export_status": figure_export["status"],
                    "pdf_path": figure_export.get("pdf_path", ""),
                    "pdf_bytes": figure_export.get("pdf_bytes", 0),
                    "pdf_sha256": figure_export.get("pdf_sha256", ""),
                    "tif_path": figure_export.get("tif_path", ""),
                    "tif_bytes": figure_export.get("tif_bytes", 0),
                    "tif_sha256": figure_export.get("tif_sha256", ""),
                    "width_px": figure_export.get("width_px", ""),
                    "height_px": figure_export.get("height_px", ""),
                    "mode": figure_export.get("mode", ""),
                    "dpi_x": figure_export.get("dpi_x", ""),
                    "dpi_y": figure_export.get("dpi_y", ""),
                    "compression": figure_export.get("compression", ""),
                    "thumbnail_path": figure_export.get("thumbnail_path", ""),
                }
            )
    pd.DataFrame(visual_qa_rows).to_csv(output / "visual_QA.csv", index=False)
    json_dump(
        {
            "renderer": str(Path(__file__).resolve()),
            "backend": "Python/matplotlib",
            "definitions": args.definitions,
            "fig2f_mode": args.fig2f_mode,
            "fig4c_mode": args.fig4c_mode,
            "D0_grouping_regression": d0_regression,
            "figure_contract": {
                "core_conclusion": (
                    "The definition-specific Fig. 2 and Fig. 4 reruns identify which "
                    "HDR-LDR patterns recur under independent baseline-referenced "
                    "recovery definitions; the integrated Fig. 2-4 conclusion awaits "
                    "the separately supervised rapid Fig. 3 sensitivity runs."
                ),
                "archetype": "definition-specific complete quantitative multi-panel figures",
                "panel_policy": (
                    "Fig. 2B is definition-specific; all remaining panels preserve the "
                    "manuscript layout and are recalculated after group replacement."
                ),
                "reviewer_risk": (
                    "Fig. 2B is definitional and must not be treated as independent "
                    "support; MDSINE2 rapid single-chain results require explicit caveats."
                ),
            },
            "records": records,
        },
        output / "run_manifest_fig2_fig4.json",
    )
    print(f"Rendered {len(records)} complete Fig. 2/Fig. 4 definition sets -> {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
