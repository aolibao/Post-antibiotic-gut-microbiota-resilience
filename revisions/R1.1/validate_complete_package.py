#!/usr/bin/env python3
"""Read-only completion gate for the Reviewer 1.1 response package.

The gate deliberately writes no report or cache files.  It exits with status 0
only when all six alternative definitions have complete, validated Fig. 2--4
exports and the bilingual response files pass their structural checks.

``--allow-draft`` relaxes only the response-text readiness check (INTERNAL
markers and unfinished-work wording).  It never relaxes a figure, summary
CSV, QA, provenance, index, or archive-integrity check.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import sys
import zipfile
from xml.etree import ElementTree

from PIL import Image


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT_ROOT = HERE / "full_original_layout"
DEFINITION_SLUGS = {
    "D1": "D1_endpoint_abs_shannon",
    "D2": "D2_endpoint_bray",
    "D3": "D3_endpoint_jaccard",
    "D4": "D4_auc_abs_shannon",
    "D5": "D5_auc_bray",
    "D6": "D6_auc_jaccard",
}
FIG3_STATUS = "validated developmental"
MIN_TIFF_DPI = 300.0
MAX_TIFF_DPI = 2400.0
MIN_TIFF_DIMENSION = 1000
EXPECTED_DEFINITION_IDS = set(DEFINITION_SLUGS)
EXPECTED_EVIDENCE_TIER = "pilot_single_or_short_chain"
EXPECTED_CHAINS_PER_GROUP = 1
EXPECTED_RETAINED_DRAWS = 100
EXPECTED_RHAT_STATUS = "unavailable_single_chain"
FIG3_SUMMARY_FILES = {
    "audit": "Fig3_across_definitions_audit.csv",
    "direction_consistency": "Fig3_across_definitions_direction_consistency.csv",
    "focus_genera": "Fig3_focus_genera_HDR_minus_LDR.csv",
}
FIG24_RECURRENCE_FILE = "Fig24_D0_vs_D1_D6_recurrence.csv"
FIG2F_AGREEMENT_FILE = "Fig2F_D0_direction_agreement.csv"
EXPECTED_D0_D6_IDS = {"D0", *EXPECTED_DEFINITION_IDS}
EXPECTED_FIG24_ITEMS = 101
EXPECTED_FIG24_ROWS = EXPECTED_FIG24_ITEMS * len(EXPECTED_D0_D6_IDS)


@dataclass(frozen=True)
class Check:
    section: str
    name: str
    passed: bool
    detail: str


class Gate:
    """Collect checks without mutating the package under inspection."""

    def __init__(self) -> None:
        self.checks: list[Check] = []
        self._sha256_cache: dict[Path, str] = {}

    def add(self, section: str, name: str, passed: bool, detail: object) -> bool:
        self.checks.append(Check(section, name, bool(passed), str(detail)))
        return bool(passed)

    def sha256(self, path: Path) -> str:
        canonical = path.resolve()
        if canonical not in self._sha256_cache:
            digest = hashlib.sha256()
            with canonical.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            self._sha256_cache[canonical] = digest.hexdigest()
        return self._sha256_cache[canonical]


def expected_paths(output_root: Path, definition_id: str) -> dict[str, Path]:
    root = output_root / DEFINITION_SLUGS[definition_id]
    figure_dir = root / "figures"
    return {
        "directory": root,
        "fig2_pdf": figure_dir / f"R1_1_{definition_id}_Fig2.pdf",
        "fig2_tif": figure_dir / f"R1_1_{definition_id}_Fig2.tif",
        "fig3_pdf": figure_dir / "Fig3_restored.pdf",
        "fig3_tif": figure_dir / "Fig3_restored.tif",
        "fig4_pdf": figure_dir / f"R1_1_{definition_id}_Fig4.pdf",
        "fig4_tif": figure_dir / f"R1_1_{definition_id}_Fig4.tif",
    }


def same_path(recorded: object, expected: Path) -> bool:
    text = str(recorded or "").strip()
    if not text:
        return False
    return Path(text).expanduser().resolve(strict=False) == expected.resolve(strict=False)


def read_json(gate: Gate, path: Path, section: str, name: str) -> object | None:
    exists = path.is_file() and path.stat().st_size > 0
    if not gate.add(section, f"{name} exists and is non-empty", exists, path):
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # a malformed audit file is a gate failure
        gate.add(section, f"{name} is valid JSON", False, f"{type(exc).__name__}: {exc}")
        return None
    gate.add(section, f"{name} is valid JSON", True, path)
    return payload


def check_pdf(gate: Gate, path: Path, section: str, label: str) -> None:
    nonempty = path.is_file() and path.stat().st_size > 0
    if not gate.add(
        section,
        f"{label} PDF exists and is non-empty",
        nonempty,
        f"{path} ({path.stat().st_size if path.is_file() else 0} bytes)",
    ):
        return
    try:
        with path.open("rb") as handle:
            signature = handle.read(5)
    except OSError as exc:
        gate.add(section, f"{label} PDF is readable", False, exc)
        return
    gate.add(section, f"{label} PDF signature", signature == b"%PDF-", repr(signature))


def check_tiff(gate: Gate, path: Path, section: str, label: str) -> None:
    nonempty = path.is_file() and path.stat().st_size > 0
    if not gate.add(
        section,
        f"{label} TIFF exists and is non-empty",
        nonempty,
        f"{path} ({path.stat().st_size if path.is_file() else 0} bytes)",
    ):
        return
    try:
        with Image.open(path) as image:
            image_format = str(image.format or "")
            mode = str(image.mode)
            width, height = (int(image.width), int(image.height))
            raw_dpi = image.info.get("dpi", (0.0, 0.0))
            try:
                dpi = tuple(float(value) for value in raw_dpi[:2])
            except (TypeError, ValueError, IndexError):
                dpi = (0.0, 0.0)
            compression = str(image.info.get("compression", ""))
            compression_tag = image.tag_v2.get(259)
    except Exception as exc:
        gate.add(section, f"{label} TIFF is readable", False, f"{type(exc).__name__}: {exc}")
        return

    gate.add(section, f"{label} TIFF format", image_format.upper() == "TIFF", image_format)
    gate.add(section, f"{label} TIFF color mode", mode == "RGB", mode)
    is_lzw = "lzw" in compression.lower() or compression_tag == 5
    gate.add(
        section,
        f"{label} TIFF LZW compression",
        is_lzw,
        f"info={compression!r}; tag259={compression_tag!r}",
    )
    dpi_ok = (
        len(dpi) == 2
        and all(MIN_TIFF_DPI <= value <= MAX_TIFF_DPI for value in dpi)
    )
    gate.add(
        section,
        f"{label} TIFF resolution is reasonable",
        dpi_ok,
        f"dpi={dpi}; accepted={MIN_TIFF_DPI:g}--{MAX_TIFF_DPI:g}",
    )
    dimensions_ok = width >= MIN_TIFF_DIMENSION and height >= MIN_TIFF_DIMENSION
    gate.add(
        section,
        f"{label} TIFF dimensions are reasonable",
        dimensions_ok,
        f"{width}x{height}; each dimension >= {MIN_TIFF_DIMENSION}",
    )


def check_figure_files(gate: Gate, output_root: Path, definition_id: str) -> None:
    section = f"{definition_id} figure exports"
    paths = expected_paths(output_root, definition_id)
    for figure in ("fig2", "fig3", "fig4"):
        label = figure.replace("fig", "Fig. ")
        check_pdf(gate, paths[f"{figure}_pdf"], section, label)
        check_tiff(gate, paths[f"{figure}_tif"], section, label)


def check_figure24_qa(gate: Gate, output_root: Path, definition_id: str) -> None:
    section = f"{definition_id} Fig. 2/4 QA"
    root = output_root / DEFINITION_SLUGS[definition_id]
    qa_path = root / "qa" / "figure_export_QA.json"
    payload = read_json(gate, qa_path, section, "figure_export_QA.json")
    if not isinstance(payload, dict):
        if payload is not None:
            gate.add(section, "QA payload is an object", False, type(payload).__name__)
        return
    gate.add(section, "QA overall status is PASS", payload.get("status") == "PASS", payload.get("status"))

    checks = payload.get("checks")
    check_rows_ok = (
        isinstance(checks, list)
        and bool(checks)
        and all(isinstance(row, dict) and row.get("passed") is True for row in checks)
    )
    failed_names = []
    if isinstance(checks, list):
        failed_names = [
            str(row.get("check", "unnamed"))
            for row in checks
            if not isinstance(row, dict) or row.get("passed") is not True
        ]
    gate.add(
        section,
        "all existing definition-level QA checks PASS",
        check_rows_ok,
        f"n={len(checks) if isinstance(checks, list) else 0}; failed={failed_names}",
    )

    exports = payload.get("figure_exports")
    export_map = {
        str(row.get("figure")): row
        for row in exports
        if isinstance(row, dict)
    } if isinstance(exports, list) else {}
    for figure in ("Fig2", "Fig4"):
        row = export_map.get(figure, {})
        expected = expected_paths(output_root, definition_id)
        passed = (
            row.get("status") == "PASS"
            and same_path(row.get("pdf_path"), expected[f"{figure.lower()}_pdf"])
            and same_path(row.get("tif_path"), expected[f"{figure.lower()}_tif"])
        )
        gate.add(
            section,
            f"{figure} existing export QA is PASS with consistent paths",
            passed,
            f"status={row.get('status')}; pdf={row.get('pdf_path')}; tif={row.get('tif_path')}",
        )


def check_visual_qa_index(gate: Gate, output_root: Path) -> None:
    """Verify the consolidated, pre-existing Fig. 2/4 QA table."""

    section = "consolidated Fig. 2/4 QA"
    path = output_root / "visual_QA.csv"
    if not gate.add(section, "visual_QA.csv exists and is non-empty", path.is_file() and path.stat().st_size > 0, path):
        return
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except Exception as exc:
        gate.add(section, "visual_QA.csv is readable", False, f"{type(exc).__name__}: {exc}")
        return
    pairs = {(str(row.get("definition_id")), str(row.get("figure"))) for row in rows}
    expected_pairs = {(definition_id, figure) for definition_id in DEFINITION_SLUGS for figure in ("Fig2", "Fig4")}
    statuses_ok = all(
        row.get("definition_QA_status") == "PASS" and row.get("figure_export_status") == "PASS"
        for row in rows
    )
    gate.add(
        section,
        "visual_QA.csv has exactly the D1--D6 Fig2/Fig4 PASS grid",
        len(rows) == 12 and pairs == expected_pairs and statuses_ok,
        f"rows={len(rows)}; pairs_match={pairs == expected_pairs}; statuses_PASS={statuses_ok}",
    )


def _read_nonempty_csv(
    gate: Gate,
    path: Path,
    section: str,
    label: str,
) -> tuple[list[dict[str, str]], set[str]] | None:
    nonempty = path.is_file() and path.stat().st_size > 0
    if not gate.add(
        section,
        f"{label} exists and is non-empty",
        nonempty,
        f"{path} ({path.stat().st_size if path.is_file() else 0} bytes)",
    ):
        return None
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
            fields = set(reader.fieldnames or [])
    except Exception as exc:
        gate.add(section, f"{label} is readable CSV", False, f"{type(exc).__name__}: {exc}")
        return None
    if not gate.add(
        section,
        f"{label} contains data rows",
        bool(rows),
        f"rows={len(rows)}; columns={sorted(fields)}",
    ):
        return None
    return rows, fields


def _all_numeric_equal(rows: list[dict[str, str]], column: str, expected: float) -> bool:
    try:
        return bool(rows) and all(float(row.get(column, "")) == expected for row in rows)
    except (TypeError, ValueError):
        return False


def _normalized_boundary(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip()).casefold()


def check_fig24_d0_recurrence_summaries(gate: Gate, output_root: Path) -> None:
    """Validate the R-derived D0 versus D1--D6 Fig. 2/4 audit tables."""

    section = "Fig. 2/4 D0 recurrence R summaries"
    recurrence_result = _read_nonempty_csv(
        gate,
        output_root / FIG24_RECURRENCE_FILE,
        section,
        FIG24_RECURRENCE_FILE,
    )
    agreement_result = _read_nonempty_csv(
        gate,
        output_root / FIG2F_AGREEMENT_FILE,
        section,
        FIG2F_AGREEMENT_FILE,
    )

    if recurrence_result is not None:
        rows, fields = recurrence_result
        required = {"item_id", "definition_id"}
        missing = sorted(required - fields)
        gate.add(
            section,
            "recurrence table has item and definition keys",
            not missing,
            f"missing={missing}",
        )
        if not missing:
            keys = [
                (str(row.get("item_id", "")).strip(), str(row.get("definition_id", "")).strip())
                for row in rows
            ]
            unique_keys = set(keys)
            item_ids = {item_id for item_id, _ in keys if item_id}
            definition_ids = {definition_id for _, definition_id in keys if definition_id}
            gate.add(
                section,
                "recurrence table has exactly 707 keyed rows",
                len(rows) == EXPECTED_FIG24_ROWS
                and len(keys) == len(unique_keys)
                and all(item_id and definition_id for item_id, definition_id in keys),
                (
                    f"rows={len(rows)}; expected={EXPECTED_FIG24_ROWS}; "
                    f"unique_keys={len(unique_keys)}; blank_keys="
                    f"{sum(not item_id or not definition_id for item_id, definition_id in keys)}"
                ),
            )
            expected_grid = {
                (item_id, definition_id)
                for item_id in item_ids
                for definition_id in EXPECTED_D0_D6_IDS
            }
            missing_pairs = expected_grid - unique_keys
            extra_pairs = unique_keys - expected_grid
            per_item_complete = all(
                {definition_id for observed_item, definition_id in unique_keys if observed_item == item_id}
                == EXPECTED_D0_D6_IDS
                for item_id in item_ids
            )
            gate.add(
                section,
                "each of 101 predefined items has one complete D0--D6 grid",
                len(item_ids) == EXPECTED_FIG24_ITEMS
                and definition_ids == EXPECTED_D0_D6_IDS
                and per_item_complete
                and not missing_pairs
                and not extra_pairs,
                (
                    f"items={len(item_ids)}; definitions={sorted(definition_ids)}; "
                    f"missing_pairs={len(missing_pairs)}; extra_pairs={len(extra_pairs)}"
                ),
            )

    if agreement_result is not None:
        rows, fields = agreement_result
        required = {"definition_id", "n_cells_significant_in_all_D1_D6"}
        missing = sorted(required - fields)
        gate.add(
            section,
            "Fig2F agreement table has definition and universal-significance fields",
            not missing,
            f"missing={missing}",
        )
        if not missing:
            definition_ids = [str(row.get("definition_id", "")).strip() for row in rows]
            gate.add(
                section,
                "Fig2F agreement table has exactly one D0--D6 row",
                len(rows) == len(EXPECTED_D0_D6_IDS)
                and len(definition_ids) == len(set(definition_ids))
                and set(definition_ids) == EXPECTED_D0_D6_IDS,
                f"rows={len(rows)}; definitions={sorted(definition_ids)}",
            )
            gate.add(
                section,
                "Fig2F universal significant-cell count is zero",
                _all_numeric_equal(rows, "n_cells_significant_in_all_D1_D6", 0.0),
                "expected n_cells_significant_in_all_D1_D6=0 in all D0--D6 rows",
            )


def check_fig3_cross_definition_summaries(gate: Gate, output_root: Path) -> None:
    """Validate the three audit-only R summaries without recomputing results."""

    section = "Fig. 3 cross-definition R summaries"
    loaded: dict[str, tuple[list[dict[str, str]], set[str]]] = {}
    for key, filename in FIG3_SUMMARY_FILES.items():
        result = _read_nonempty_csv(gate, output_root / filename, section, filename)
        if result is not None:
            loaded[key] = result
    if set(loaded) != set(FIG3_SUMMARY_FILES):
        # Missing files have already produced specific failures above.  Do not
        # create misleading schema/boundary cascades from partial inputs.
        return

    audit, audit_fields = loaded["audit"]
    focus, focus_fields = loaded["focus_genera"]
    consistency, consistency_fields = loaded["direction_consistency"]

    common_fields = {
        "definition_id",
        "evidence_tier",
        "chains_per_group",
        "retained_draws_per_chain",
        "rhat_status",
        "interpretation_boundary",
    }
    for label, rows, fields in (
        ("audit", audit, audit_fields),
        ("focus-genera", focus, focus_fields),
    ):
        missing = sorted(common_fields - fields)
        gate.add(
            section,
            f"{label} summary has evidence-boundary schema",
            not missing,
            f"missing={missing}",
        )
        if missing:
            continue
        observed_ids = {str(row.get("definition_id", "")).strip() for row in rows}
        gate.add(
            section,
            f"{label} summary contains exactly D1--D6",
            observed_ids == EXPECTED_DEFINITION_IDS,
            f"ids={sorted(observed_ids)}",
        )
        tiers = {str(row.get("evidence_tier", "")).strip() for row in rows}
        gate.add(
            section,
            f"{label} evidence tier is uniformly developmental",
            tiers == {EXPECTED_EVIDENCE_TIER},
            f"values={sorted(tiers)}",
        )
        gate.add(
            section,
            f"{label} is uniformly one chain per group",
            _all_numeric_equal(rows, "chains_per_group", EXPECTED_CHAINS_PER_GROUP),
            f"expected={EXPECTED_CHAINS_PER_GROUP}",
        )
        gate.add(
            section,
            f"{label} retained-draw boundary is uniform",
            _all_numeric_equal(rows, "retained_draws_per_chain", EXPECTED_RETAINED_DRAWS),
            f"expected={EXPECTED_RETAINED_DRAWS}",
        )
        rhat_values = {str(row.get("rhat_status", "")).strip() for row in rows}
        gate.add(
            section,
            f"{label} uniformly records unavailable single-chain R-hat",
            rhat_values == {EXPECTED_RHAT_STATUS},
            f"values={sorted(rhat_values)}",
        )

    direction_columns = {f"D{index}_direction" for index in range(1, 7)}
    direction_required = direction_columns | {
        "n_definitions",
        "evidence_boundary",
    }
    direction_missing = sorted(direction_required - consistency_fields)
    gate.add(
        section,
        "direction-consistency summary has D1--D6 and boundary schema",
        not direction_missing,
        f"missing={direction_missing}",
    )
    if not direction_missing:
        all_directions_present = all(
            all(str(row.get(column, "")).strip() for column in direction_columns)
            for row in consistency
        )
        gate.add(
            section,
            "direction-consistency rows contain all D1--D6 directions",
            all_directions_present
            and _all_numeric_equal(consistency, "n_definitions", len(EXPECTED_DEFINITION_IDS)),
            f"rows={len(consistency)}; n_definitions_expected={len(EXPECTED_DEFINITION_IDS)}",
        )

    boundaries: dict[str, set[str]] = {}
    if "interpretation_boundary" in audit_fields:
        boundaries["audit"] = {_normalized_boundary(row.get("interpretation_boundary")) for row in audit}
    if "interpretation_boundary" in focus_fields:
        boundaries["focus_genera"] = {_normalized_boundary(row.get("interpretation_boundary")) for row in focus}
    if "evidence_boundary" in consistency_fields:
        boundaries["direction_consistency"] = {
            _normalized_boundary(row.get("evidence_boundary")) for row in consistency
        }
    nonempty_singletons = (
        len(boundaries) == 3
        and all(len(values) == 1 and "" not in values for values in boundaries.values())
    )
    combined = set().union(*boundaries.values()) if boundaries else set()
    gate.add(
        section,
        "all three summaries use one identical evidence boundary",
        nonempty_singletons and len(combined) == 1,
        f"unique_counts={{{', '.join(f'{key}: {len(value)}' for key, value in boundaries.items())}}}; cross_file_unique={len(combined)}",
    )
    boundary = next(iter(combined), "") if len(combined) == 1 else ""
    required_phrases = (
        "pilot/developmental",
        "one-chain-per-group",
        "r-hat is unavailable",
        "no formal multichain convergence",
    )
    missing_phrases = [phrase for phrase in required_phrases if phrase not in boundary]
    gate.add(
        section,
        "shared boundary states developmental tier, one chain, and unavailable R-hat",
        bool(boundary) and not missing_phrases,
        f"missing_phrases={missing_phrases}",
    )


def check_d0_gate(gate: Gate, output_root: Path) -> None:
    section = "D0 grouping regression gate"
    path = output_root / "D0_grouping_regression.json"
    payload = read_json(gate, path, section, "D0_grouping_regression.json")
    if not isinstance(payload, dict):
        if payload is not None:
            gate.add(section, "D0 gate payload is an object", False, type(payload).__name__)
        return
    checks = payload.get("checks")
    all_pass = (
        payload.get("status") == "PASS"
        and isinstance(checks, list)
        and bool(checks)
        and all(isinstance(row, dict) and row.get("passed") is True for row in checks)
    )
    failed_names = []
    if isinstance(checks, list):
        failed_names = [
            str(row.get("check", "unnamed"))
            for row in checks
            if not isinstance(row, dict) or row.get("passed") is not True
        ]
    gate.add(
        section,
        "D0 grouping regression status and all checks PASS",
        all_pass,
        f"status={payload.get('status')}; n={len(checks) if isinstance(checks, list) else 0}; failed={failed_names}",
    )


def check_fig3_finalization(gate: Gate, output_root: Path, definition_id: str) -> None:
    section = f"{definition_id} Fig. 3 finalization"
    root = output_root / DEFINITION_SLUGS[definition_id]
    expected = expected_paths(output_root, definition_id)
    marker_path = root / "mdsine2" / "finalized.json"
    payload = read_json(gate, marker_path, section, "finalized.json")
    if isinstance(payload, dict):
        gate.add(
            section,
            "finalized marker definition matches",
            payload.get("definition_id") == definition_id,
            payload.get("definition_id"),
        )
        outputs = payload.get("outputs")
        output_rows = outputs if isinstance(outputs, list) else []
        for key in ("fig3_pdf", "fig3_tif"):
            artifact = expected[key]
            matches = [
                row for row in output_rows
                if isinstance(row, dict) and same_path(row.get("path"), artifact)
            ]
            if len(matches) != 1:
                gate.add(
                    section,
                    f"finalized marker records {key}",
                    False,
                    f"matching records={len(matches)}; expected={artifact}",
                )
                continue
            row = matches[0]
            file_ready = artifact.is_file() and artifact.stat().st_size > 0
            try:
                recorded_bytes = int(row.get("bytes", -1))
            except (TypeError, ValueError):
                recorded_bytes = -1
            bytes_match = file_ready and recorded_bytes == artifact.stat().st_size
            try:
                hash_match = file_ready and str(row.get("sha256", "")) == gate.sha256(artifact)
            except OSError:
                hash_match = False
            gate.add(
                section,
                f"finalized marker {key} size and SHA-256 match",
                bytes_match and hash_match,
                f"bytes_match={bytes_match}; sha256_match={hash_match}",
            )
    elif payload is not None:
        gate.add(section, "finalized marker payload is an object", False, type(payload).__name__)

    validation_path = root / "qa" / "validation_checks.tsv"
    exists = validation_path.is_file() and validation_path.stat().st_size > 0
    if not gate.add(section, "validation_checks.tsv exists and is non-empty", exists, validation_path):
        return
    try:
        with validation_path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except Exception as exc:
        gate.add(section, "validation_checks.tsv is readable", False, f"{type(exc).__name__}: {exc}")
        return
    required = {"name", "status", "severity"}
    headers = set(rows[0]) if rows else set()
    if not gate.add(
        section,
        "validation_checks.tsv has the required schema",
        bool(rows) and required.issubset(headers),
        f"rows={len(rows)}; columns={sorted(headers)}",
    ):
        return
    error_rows = [row for row in rows if str(row.get("severity", "")).strip().lower() == "error"]
    failed_errors = [str(row.get("name", "unnamed")) for row in error_rows if row.get("status") != "PASS"]
    gate.add(
        section,
        "all error-level Fig. 3 validation checks PASS",
        bool(error_rows) and not failed_errors,
        f"error_checks={len(error_rows)}; failed={failed_errors}",
    )


def check_reviewer_index(gate: Gate, output_root: Path) -> None:
    section = "reviewer-ready index"
    path = output_root / "reviewer_ready_index.csv"
    if not gate.add(section, "reviewer_ready_index.csv exists and is non-empty", path.is_file() and path.stat().st_size > 0, path):
        return
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except Exception as exc:
        gate.add(section, "reviewer_ready_index.csv is readable", False, f"{type(exc).__name__}: {exc}")
        return
    by_definition: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_definition.setdefault(str(row.get("definition_id", "")), []).append(row)
    exact_ids = len(rows) == 6 and set(by_definition) == set(DEFINITION_SLUGS) and all(
        len(value) == 1 for value in by_definition.values()
    )
    gate.add(
        section,
        "index has exactly one D1--D6 row",
        exact_ids,
        f"rows={len(rows)}; ids={sorted(by_definition)}",
    )
    for definition_id in DEFINITION_SLUGS:
        expected = expected_paths(output_root, definition_id)
        candidates = by_definition.get(definition_id, [])
        if len(candidates) != 1:
            gate.add(section, f"{definition_id} index row is unique", False, f"rows={len(candidates)}")
            continue
        row = candidates[0]
        gate.add(
            section,
            f"{definition_id} Fig. 3 status is {FIG3_STATUS!r}",
            row.get("fig3_status") == FIG3_STATUS,
            row.get("fig3_status", ""),
        )
        path_fields = ("directory", "fig2_pdf", "fig2_tif", "fig3_pdf", "fig3_tif", "fig4_pdf", "fig4_tif")
        mismatches = [
            field for field in path_fields
            if not same_path(row.get(field), expected[field])
        ]
        gate.add(
            section,
            f"{definition_id} indexed paths match deterministic outputs",
            not mismatches,
            f"mismatches={mismatches}",
        )


def docx_text_and_integrity(gate: Gate, path: Path, section: str) -> str | None:
    nonempty = path.is_file() and path.stat().st_size > 0
    if not gate.add(
        section,
        "bilingual DOCX exists and is non-empty",
        nonempty,
        f"{path} ({path.stat().st_size if path.is_file() else 0} bytes)",
    ):
        return None
    try:
        with zipfile.ZipFile(path, "r") as archive:
            names = set(archive.namelist())
            required = {"[Content_Types].xml", "word/document.xml"}
            missing = sorted(required - names)
            bad_member = archive.testzip()
            document_xml = archive.read("word/document.xml") if not missing else b""
    except Exception as exc:
        gate.add(section, "DOCX ZIP archive is readable", False, f"{type(exc).__name__}: {exc}")
        return None
    archive_ok = not missing and bad_member is None
    gate.add(
        section,
        "DOCX ZIP archive is complete",
        archive_ok,
        f"missing_required={missing}; first_bad_CRC_member={bad_member}",
    )
    if not archive_ok:
        return None
    try:
        root = ElementTree.fromstring(document_xml)
        return " ".join(piece for piece in root.itertext() if piece)
    except ElementTree.ParseError as exc:
        gate.add(section, "DOCX document.xml is valid XML", False, exc)
        return None


def response_readiness_markers(text: str) -> list[str]:
    marker_patterns = {
        "INTERNAL": r"\binternal\b",
        "draft": r"\bdraft\b",
        "NOT READY FOR SUBMISSION": r"not\s+ready\s+for\s+submission",
        "草稿": r"草稿",
        "暂不可提交": r"暂不可提交",
        "in progress": r"\bin[\s-]+progress\b",
        "pending": r"\bpending\b",
        "awaiting": r"\bawaiting\b",
        "being generated/run/finalized": (
            r"\b(?:is|are|remains?|remain)\s+(?:still\s+)?being\s+"
            r"(?:generated|run|processed|completed|finali[sz]ed)\b"
        ),
        "not yet complete/available": (
            r"\bnot\s+yet\s+(?:complete|completed|available|finali[sz]ed)\b"
        ),
        "to be inserted/completed": (
            r"\bto\s+be\s+(?:inserted|completed|finali[sz]ed|provided|added|updated)\b"
        ),
        "will be provided/completed": (
            r"\bwill\s+be\s+(?:provided|supplied|completed|finali[sz]ed|generated|added|updated|treated|interpreted)\b"
        ),
        "provisional until completion": r"\b(?:remain\s+)?provisional\b[^.;\n]{0,100}\buntil\b",
        "正在运行": r"正在运行|仍在运行|尚在运行",
        "正在生成": r"正在生成|正使用[^。；;\n]{0,40}(?:生成|重跑|运行)",
        "待完成/尚未完成": r"待完成|尚未完成|尚未全部完成|尚待完成",
        "待补充/待插入": r"待[^。；;\n]{0,30}(?:补充|插入|更新|提供)",
        "完成前/暂定": r"完成前|保持暂定|结论暂定",
    }
    return [label for label, pattern in marker_patterns.items() if re.search(pattern, text, flags=re.IGNORECASE)]


def check_response_files(gate: Gate, response_root: Path, allow_draft: bool) -> None:
    section = "bilingual response"
    md_path = response_root / "R1.1_response_bilingual.md"
    docx_path = response_root / "R1.1_response_bilingual.docx"
    md_nonempty = md_path.is_file() and md_path.stat().st_size > 0
    gate.add(
        section,
        "bilingual Markdown exists and is non-empty",
        md_nonempty,
        f"{md_path} ({md_path.stat().st_size if md_path.is_file() else 0} bytes)",
    )
    md_text: str | None = None
    if md_nonempty:
        try:
            md_text = md_path.read_text(encoding="utf-8")
        except Exception as exc:
            gate.add(section, "bilingual Markdown is UTF-8 readable", False, f"{type(exc).__name__}: {exc}")
    docx_text = docx_text_and_integrity(gate, docx_path, section)

    embedded_media: list[str] = []
    if docx_path.is_file():
        try:
            with zipfile.ZipFile(docx_path, "r") as archive:
                embedded_media = sorted(
                    name for name in archive.namelist()
                    if name.startswith("word/media/") and not name.endswith("/")
                )
        except Exception:
            embedded_media = []
    gate.add(
        section,
        "final DOCX embeds exactly 18 response figures",
        len(embedded_media) == 18,
        f"embedded_media={len(embedded_media)}",
    )
    non_png = [name for name in embedded_media if Path(name).suffix.lower() != ".png"]
    gate.add(
        section,
        "all embedded response figures are Word-ready PNG files",
        bool(embedded_media) and not non_png,
        f"non_png={non_png}",
    )
    expected_figure_labels = [
        f"R1.1-D{definition}-Fig. {figure}"
        for definition in range(1, 7)
        for figure in (2, 3, 4)
    ]
    missing_figure_labels = [
        label for label in expected_figure_labels
        if label not in (docx_text or "")
    ]
    gate.add(
        section,
        "final DOCX contains all 18 response-figure labels",
        not missing_figure_labels,
        f"missing={missing_figure_labels}",
    )

    if allow_draft:
        gate.add(
            section,
            "response text is submission-ready",
            True,
            "SKIPPED only because --allow-draft was supplied; figures and summary CSVs remain strict",
        )
        return
    for label, text in (("Markdown", md_text), ("DOCX", docx_text)):
        markers = response_readiness_markers(text or "") if text is not None else ["unreadable response text"]
        gate.add(
            section,
            f"{label} has no draft, unfinished-work, or revision-location placeholder wording",
            not markers,
            f"markers={markers}",
        )


def evaluate(output_root: Path, response_root: Path, allow_draft: bool) -> Gate:
    gate = Gate()
    check_d0_gate(gate, output_root)
    check_visual_qa_index(gate, output_root)
    check_fig24_d0_recurrence_summaries(gate, output_root)
    for definition_id in DEFINITION_SLUGS:
        check_figure_files(gate, output_root, definition_id)
        check_figure24_qa(gate, output_root, definition_id)
        check_fig3_finalization(gate, output_root, definition_id)
    check_fig3_cross_definition_summaries(gate, output_root)
    check_reviewer_index(gate, output_root)
    check_response_files(gate, response_root, allow_draft)
    return gate


def print_report(gate: Gate, verbose: bool) -> None:
    sections: list[str] = []
    for check in gate.checks:
        if check.section not in sections:
            sections.append(check.section)
    for section in sections:
        checks = [check for check in gate.checks if check.section == section]
        failures = [check for check in checks if not check.passed]
        status = "PASS" if not failures else "FAIL"
        print(f"[{status}] {section}: {len(checks) - len(failures)}/{len(checks)} checks passed")
        shown = checks if verbose else failures
        for check in shown:
            flag = "PASS" if check.passed else "FAIL"
            print(f"  [{flag}] {check.name}: {check.detail}")

    failures = [check for check in gate.checks if not check.passed]
    overall = "PASS" if not failures else "FAIL"
    print()
    print(
        f"R1.1 COMPLETE PACKAGE GATE: {overall} "
        f"({len(gate.checks) - len(failures)}/{len(gate.checks)} checks passed; "
        f"{len(failures)} failed)"
    )
    if failures and not verbose:
        print("Only failures are expanded above; use --verbose to print every check.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Directory containing D1--D6 full_original_layout outputs.",
    )
    parser.add_argument(
        "--response-root",
        type=Path,
        default=HERE,
        help="Directory containing R1.1_response_bilingual.md/.docx.",
    )
    parser.add_argument(
        "--allow-draft",
        action="store_true",
        help=(
            "Skip only the response text-readiness check (INTERNAL markers, "
            "unfinished wording, and location placeholders). Figure, summary "
            "CSV, and QA requirements remain strict."
        ),
    )
    parser.add_argument("--verbose", action="store_true", help="Print passing checks as well as failures.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    gate = evaluate(
        output_root=args.output_root.expanduser().resolve(),
        response_root=args.response_root.expanduser().resolve(),
        allow_draft=bool(args.allow_draft),
    )
    print_report(gate, verbose=bool(args.verbose))
    return 0 if all(check.passed for check in gate.checks) else 1


if __name__ == "__main__":
    sys.exit(main())
