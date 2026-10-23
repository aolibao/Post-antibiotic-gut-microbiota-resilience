#!/usr/bin/env python3
"""Validated entry point for the manuscript's official no-qPCR MDSINE2 runner.

This file deliberately does not reimplement the model.  It verifies the two
no-qPCR biomass constants in the pinned manuscript runner, audits the prepared
``Study``, and then executes that runner unchanged with ``runpy``.
"""

from __future__ import annotations

import argparse
import ast
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import sys
import traceback
from typing import Any, Iterable

import numpy as np


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _numeric_ast_value(node: ast.AST) -> float | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if (
        isinstance(node, ast.UnaryOp)
        and isinstance(node.op, (ast.UAdd, ast.USub))
        and isinstance(node.operand, ast.Constant)
        and isinstance(node.operand.value, (int, float))
    ):
        sign = -1.0 if isinstance(node.op, ast.USub) else 1.0
        return sign * float(node.operand.value)
    return None


def inspect_runner_contract(runner: Path) -> dict[str, Any]:
    """Extract the hard-coded no-qPCR biomass assumptions from source AST."""
    source = runner.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(runner))
    means: list[float] = []
    cvs: list[float] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func_name = None
        if isinstance(node.func, ast.Name):
            func_name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            func_name = node.func.attr
        if func_name == "initialize_graph_no_qpcr":
            for keyword in node.keywords:
                if keyword.arg == "biomass_mean":
                    value = _numeric_ast_value(keyword.value)
                    if value is not None:
                        means.append(value)
        if func_name == "FilteringWithoutQPCR":
            for keyword in node.keywords:
                if keyword.arg == "biomass_cv":
                    value = _numeric_ast_value(keyword.value)
                    if value is not None:
                        cvs.append(value)
    return {
        "path": str(runner.resolve()),
        "sha256": sha256_file(runner),
        "biomass_mean_literals": sorted(set(means)),
        "biomass_cv_literals": sorted(set(cvs)),
    }


def _flag_value(arguments: list[str], flag: str) -> str:
    try:
        index = arguments.index(flag)
    except ValueError as exc:
        raise ValueError(f"Official runner argument {flag!r} is required") from exc
    if index + 1 >= len(arguments):
        raise ValueError(f"Official runner argument {flag!r} has no value")
    return arguments[index + 1]


def validate_runner_arguments(arguments: list[str]) -> Path:
    if arguments and arguments[0] == "--":
        arguments.pop(0)
    if "--resume" in arguments or "-r" in arguments:
        raise ValueError(
            "--resume is intentionally forbidden: the manuscript no-qPCR runner "
            "parses it but does not implement trace continuation"
        )
    required = (
        "--input",
        "--negbin",
        "--seed",
        "--burnin",
        "--n-samples",
        "--checkpoint",
        "--basepath",
        "--rename-study",
        "--interaction-ind-prior",
        "--perturbation-ind-prior",
    )
    for flag in required:
        if flag not in arguments:
            raise ValueError(f"Official runner argument {flag!r} is required")
    return Path(_flag_value(arguments, "--input")).resolve()


def _same_times(left: Iterable[float], right: Iterable[float]) -> bool:
    return set(map(float, left)) == set(map(float, right))


def audit_study(
    input_path: Path,
    expected_subjects: int,
    expected_taxa: int,
    perturbation_start: float,
    perturbation_end: float,
) -> dict[str, Any]:
    """Load and structurally audit a trusted, locally generated Study pickle."""
    import mdsine2 as md2

    study = md2.Study.load(str(input_path))
    subjects = list(study)
    taxa_names = [str(taxon.name) for taxon in study.taxa]
    if len(subjects) != expected_subjects:
        raise ValueError(
            f"Expected {expected_subjects} subjects, found {len(subjects)} in {input_path}"
        )
    if len(taxa_names) != expected_taxa:
        raise ValueError(
            f"Expected {expected_taxa} taxa, found {len(taxa_names)} in {input_path}"
        )
    if len(set(taxa_names)) != len(taxa_names):
        raise ValueError("Taxon names are not unique")

    subject_names: list[str] = []
    sample_counts: list[int] = []
    depths: list[int] = []
    all_times: list[float] = []
    for subject in subjects:
        name = str(subject.name)
        subject_names.append(name)
        read_times = sorted(map(float, subject.reads.keys()))
        if not read_times:
            raise ValueError(f"Subject {name!r} has no read observations")
        if len(read_times) != len(set(read_times)):
            raise ValueError(f"Subject {name!r} has duplicate observation times")
        if not _same_times(subject.qpcr.keys(), read_times):
            raise ValueError(
                f"Subject {name!r} structural qPCR keys do not match read times"
            )
        sample_counts.append(len(read_times))
        all_times.extend(read_times)
        for timepoint in read_times:
            values = np.asarray(subject.reads[timepoint])
            if values.shape != (expected_taxa,):
                raise ValueError(
                    f"Subject {name!r}, time {timepoint:g}: expected read vector "
                    f"length {expected_taxa}, found shape {values.shape}"
                )
            if not np.issubdtype(values.dtype, np.integer):
                if not np.all(np.isfinite(values)) or not np.all(values == np.floor(values)):
                    raise ValueError(
                        f"Subject {name!r}, time {timepoint:g}: reads are not integer counts"
                    )
            if np.any(values < 0):
                raise ValueError(
                    f"Subject {name!r}, time {timepoint:g}: reads contain negative counts"
                )
            depth = int(values.sum())
            if depth <= 0:
                raise ValueError(
                    f"Subject {name!r}, time {timepoint:g}: total read count is not positive"
                )
            depths.append(depth)
            qpcr = subject.qpcr[timepoint]
            raw_qpcr = np.asarray(qpcr._raw_data, dtype=float)
            if raw_qpcr.size == 0 or not np.all(np.isfinite(raw_qpcr)) or np.any(raw_qpcr <= 0):
                raise ValueError(
                    f"Subject {name!r}, time {timepoint:g}: invalid structural qPCR placeholder"
                )

    if len(set(subject_names)) != len(subject_names):
        raise ValueError("Subject names are not unique")
    perturbations = list(study.perturbations or [])
    if len(perturbations) != 1:
        raise ValueError(
            f"Expected one ciprofloxacin perturbation, found {len(perturbations)}"
        )
    perturbation = perturbations[0]
    missing_starts = sorted(set(subject_names) - set(perturbation.starts))
    missing_ends = sorted(set(subject_names) - set(perturbation.ends))
    if missing_starts or missing_ends:
        raise ValueError(
            "Perturbation does not cover every subject: "
            f"missing starts={missing_starts}, missing ends={missing_ends}"
        )
    starts = [float(perturbation.starts[name]) for name in subject_names]
    ends = [float(perturbation.ends[name]) for name in subject_names]
    if not all(math.isclose(x, perturbation_start, abs_tol=1e-9) for x in starts):
        raise ValueError(f"Perturbation starts are not uniformly {perturbation_start:g}")
    if not all(math.isclose(x, perturbation_end, abs_tol=1e-9) for x in ends):
        raise ValueError(f"Perturbation ends are not uniformly {perturbation_end:g}")

    return {
        "study_name_before_rename": str(study.name),
        "n_subjects": len(subjects),
        "n_taxa": len(taxa_names),
        "n_samples": len(depths),
        "samples_per_subject_min": min(sample_counts),
        "samples_per_subject_max": max(sample_counts),
        "read_depth_min": min(depths),
        "read_depth_max": max(depths),
        "time_min": min(all_times),
        "time_max": max(all_times),
        "perturbation_name": str(perturbation.name),
        "perturbation_start": perturbation_start,
        "perturbation_end": perturbation_end,
        "qpcr_role": (
            "structural placeholder only; the official no-qPCR runner replaces all values "
            "with its fixed biomass prior before model initialization"
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--official-runner", type=Path, required=True)
    parser.add_argument("--audit-json", type=Path, required=True)
    parser.add_argument("--expected-group", choices=("LDR", "HDR"), required=True)
    parser.add_argument("--expected-subjects", type=int, default=15)
    parser.add_argument("--expected-taxa", type=int, default=259)
    parser.add_argument("--perturbation-start", type=float, default=0.5)
    parser.add_argument("--perturbation-end", type=float, default=5.5)
    parser.add_argument("--expected-biomass-mean", type=float, default=1e11)
    parser.add_argument("--expected-biomass-cv", type=float, default=0.10)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Stop after the runner-contract and Study preflight audit.",
    )
    parser.add_argument("runner_args", nargs=argparse.REMAINDER)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    runner = args.official_runner.resolve()
    audit: dict[str, Any] = {
        "schema_version": 1,
        "audited_at_utc": utc_now(),
        "group": args.expected_group,
        "status": "preflight_failed",
    }
    try:
        if not runner.is_file():
            raise FileNotFoundError(f"Official no-qPCR runner not found: {runner}")
        contract = inspect_runner_contract(runner)
        if contract["biomass_mean_literals"] != [float(args.expected_biomass_mean)]:
            raise RuntimeError(
                "Official runner biomass-mean contract changed: "
                f"{contract['biomass_mean_literals']}"
            )
        if contract["biomass_cv_literals"] != [float(args.expected_biomass_cv)]:
            raise RuntimeError(
                "Official runner biomass-CV contract changed: "
                f"{contract['biomass_cv_literals']}"
            )
        runner_args = list(args.runner_args)
        input_path = validate_runner_arguments(runner_args)
        if not input_path.is_file():
            raise FileNotFoundError(f"MDSINE2 Study input not found: {input_path}")
        study_audit = audit_study(
            input_path=input_path,
            expected_subjects=args.expected_subjects,
            expected_taxa=args.expected_taxa,
            perturbation_start=args.perturbation_start,
            perturbation_end=args.perturbation_end,
        )
        audit.update(
            {
                "status": "preflight_passed",
                "input": str(input_path),
                "input_sha256": sha256_file(input_path),
                "runner_contract": contract,
                "study": study_audit,
            }
        )
        atomic_json(args.audit_json.resolve(), audit)
        if args.validate_only:
            return
        sys.argv = [str(runner), *runner_args]
        runpy.run_path(str(runner), run_name="__main__")
    except BaseException as exc:
        audit.update(
            {
                "status": "preflight_or_run_failed",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(),
            }
        )
        atomic_json(args.audit_json.resolve(), audit)
        raise


if __name__ == "__main__":
    main()
