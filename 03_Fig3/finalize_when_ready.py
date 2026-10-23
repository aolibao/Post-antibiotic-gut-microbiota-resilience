#!/usr/bin/env python3
"""Wait for six formal Fig. 3 chains, then run the formal post-processing once.

This is deliberately a foreground watcher.  It does not launch MDSINE2, does
not repair inference state, and does not daemonize itself.  A caller may place
it in tmux/nohup after reviewing the command documented in README.md.

The post-processing gate is intentionally redundant: the inference supervisor
must report ``succeeded`` and both ``status.tsv`` and each selected immutable
attempt's ``status.json`` must independently prove success for the fixed
2 groups x 3 seeds contract.  Model directories are read from those status
records rather than guessed from an attempt number.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from dataclasses import asdict, dataclass
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time
import traceback
from types import SimpleNamespace
from typing import Any, Callable, Sequence


SCRIPT_PATH = Path(__file__).resolve()
MODULE_DIR = SCRIPT_PATH.parent
PROJECT_ROOT = MODULE_DIR.parents[1]
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "reproduction" / "outputs" / "fig3_mdsine2"
DEFAULT_INFERENCE_ROOT = DEFAULT_OUTPUT_ROOT / "inference"
MDSINE_PYTHON = Path("/home/qqt/.conda/envs/mdsine2/bin/python")

GROUPS = ("LDR", "HDR")
FORMAL_SEEDS = (20261009, 20261010, 20261011)
REQUIRED_MODEL_FILES = ("mcmc.pkl", "metadata.txt", "params.pkl", "traces.hdf5")
FORBIDDEN_FORMAL_FLAGS = ("--development-fallback", "--allow-development")
STATE_SCHEMA_VERSION = 2


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    """Write a JSON document atomically in the destination directory."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _integer(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _zero(value: Any) -> bool:
    return _integer(value) == 0


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


@dataclass(frozen=True)
class FormalChain:
    group: str
    seed: int
    attempt: int
    model_dir: Path
    status_path: Path
    mcmc_path: Path
    semantic_audit: dict[str, Any] | None = None

    def to_json(self) -> dict[str, Any]:
        payload = asdict(self)
        return {key: str(value) if isinstance(value, Path) else value for key, value in payload.items()}


@dataclass(frozen=True)
class GateResult:
    ready: bool
    supervisor_state: str
    reasons: tuple[str, ...]
    chains: tuple[FormalChain, ...]
    checked_at_utc: str

    def to_json(self) -> dict[str, Any]:
        return {
            "ready": self.ready,
            "supervisor_state": self.supervisor_state,
            "reasons": list(self.reasons),
            "chains": [chain.to_json() for chain in self.chains],
            "checked_at_utc": self.checked_at_utc,
            "required_groups": list(GROUPS),
            "required_seeds": list(FORMAL_SEEDS),
        }


@dataclass(frozen=True)
class CommandPlan:
    name: str
    argv: tuple[str, ...]
    cwd: Path
    log_path: Path

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "argv": list(self.argv),
            "shell_escaped_display_only": shlex.join(self.argv),
            "cwd": str(self.cwd),
            "log_path": str(self.log_path),
            "state": "planned",
        }


def _read_status_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "profile",
            "group",
            "seed",
            "attempt",
            "state",
            "returncode",
            "child_returncode",
            "model_dir",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"status.tsv is missing required columns: {', '.join(sorted(missing))}"
            )
        return [dict(row) for row in reader]


def _formal_contract() -> dict[str, Any]:
    return {
        "groups": list(GROUPS),
        "seeds": list(FORMAL_SEEDS),
        "chains_per_group": len(FORMAL_SEEDS),
        "bayes_factor_threshold": 10,
        "evidence_tier": "formal_multichain",
        "development_fallback_forbidden": True,
        "allow_development_forbidden": True,
        "artifact_symlink_escape_forbidden": True,
        "live_semantic_output_reaudit": True,
        "strict_qa_revalidated_on_idempotent_exit": True,
    }


def _live_semantic_audit(chain: FormalChain) -> dict[str, Any]:
    """Re-run the canonical no-qPCR output audit against current bytes.

    ``validate_model_output`` only needs ``model_dir``, ``seed`` and the
    sampling profile from its spec.  A small adapter lets the finalizer audit
    a dynamically discovered model directory without guessing its basename.
    """

    import run_inference

    spec = SimpleNamespace(
        model_dir=chain.model_dir,
        seed=chain.seed,
        profile=run_inference.PROFILES["formal"],
    )
    audit = run_inference.validate_model_output(spec)
    if not isinstance(audit, dict):
        raise TypeError("validate_model_output did not return a JSON object")
    return audit


SemanticAuditor = Callable[[FormalChain], dict[str, Any]]


def evaluate_formal_gate(
    inference_root: Path,
    *,
    semantic_auditor: SemanticAuditor | None = None,
) -> GateResult:
    """Cross-check the supervisor, status table, and six attempt status files."""

    inference_root = Path(inference_root).resolve()
    checked_at = utc_now()
    reasons: list[str] = []
    chains: list[FormalChain] = []

    supervisor_path = inference_root / "supervisor_status.json"
    supervisor_state = "missing"
    if not supervisor_path.is_file():
        reasons.append(f"missing supervisor status: {supervisor_path}")
    else:
        try:
            supervisor = load_json(supervisor_path)
            supervisor_state = str(supervisor.get("state", "unknown"))
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            supervisor_state = "unreadable"
            reasons.append(f"cannot read supervisor status: {exc}")
        if supervisor_state != "succeeded":
            reasons.append(
                f"supervisor state is {supervisor_state!r}, not 'succeeded'"
            )

    table_path = inference_root / "status.tsv"
    if not table_path.is_file():
        reasons.append(f"missing status table: {table_path}")
        return GateResult(False, supervisor_state, tuple(reasons), (), checked_at)
    try:
        rows = _read_status_rows(table_path)
    except (OSError, ValueError, csv.Error) as exc:
        reasons.append(f"cannot read status table: {exc}")
        return GateResult(False, supervisor_state, tuple(reasons), (), checked_at)

    expected = {(group, seed) for group in GROUPS for seed in FORMAL_SEEDS}
    formal_rows = [row for row in rows if row.get("profile") == "formal"]
    successful_rows: dict[tuple[str, int], list[dict[str, str]]] = {}
    for row in formal_rows:
        if row.get("state") != "succeeded":
            continue
        seed = _integer(row.get("seed"))
        if seed is None:
            reasons.append(f"formal succeeded row has invalid seed: {row.get('seed')!r}")
            continue
        key = (str(row.get("group", "")), seed)
        successful_rows.setdefault(key, []).append(row)

    unexpected = sorted(set(successful_rows) - expected)
    if unexpected:
        reasons.append(f"unexpected formal succeeded group/seed rows: {unexpected}")

    for group in GROUPS:
        for seed in FORMAL_SEEDS:
            key = (group, seed)
            matches = successful_rows.get(key, [])
            if len(matches) != 1:
                observed = [
                    row.get("state", "")
                    for row in formal_rows
                    if row.get("group") == group and _integer(row.get("seed")) == seed
                ]
                reasons.append(
                    f"expected exactly one formal succeeded row for {group}/{seed}; "
                    f"found {len(matches)} (all observed states: {observed or ['missing']})"
                )
                continue
            row = matches[0]
            attempt = _integer(row.get("attempt"))
            if attempt is None or attempt < 1:
                reasons.append(f"invalid attempt for {group}/{seed}: {row.get('attempt')!r}")
                continue
            if not _zero(row.get("returncode")) or not _zero(row.get("child_returncode")):
                reasons.append(
                    f"status.tsv lacks zero return codes for {group}/{seed}/attempt_{attempt:03d}"
                )
                continue

            attempt_root = (
                inference_root
                / "formal"
                / group
                / f"seed_{seed}"
                / f"attempt_{attempt:03d}"
            ).resolve()
            if not _is_within(attempt_root, inference_root):
                reasons.append(
                    f"attempt root escapes inference root for {group}/{seed}: {attempt_root}"
                )
                continue
            status_path = attempt_root / "status.json"
            if not status_path.is_file():
                reasons.append(f"missing attempt status for {group}/{seed}: {status_path}")
                continue
            try:
                resolved_status = status_path.resolve(strict=True)
            except OSError as exc:
                reasons.append(f"cannot resolve attempt status {status_path}: {exc}")
                continue
            if not _is_within(resolved_status, attempt_root):
                reasons.append(
                    f"attempt status symlink escapes its immutable attempt root: {status_path}"
                )
                continue
            try:
                status = load_json(resolved_status)
            except (OSError, ValueError, json.JSONDecodeError) as exc:
                reasons.append(f"cannot read attempt status {status_path}: {exc}")
                continue

            expected_fields = {
                "profile": "formal",
                "group": group,
                "seed": seed,
                "attempt": attempt,
                "state": "succeeded",
            }
            mismatch = {
                field: (status.get(field), wanted)
                for field, wanted in expected_fields.items()
                if status.get(field) != wanted
            }
            if mismatch:
                reasons.append(f"attempt status field mismatch at {status_path}: {mismatch}")
                continue
            if not _zero(status.get("returncode")) or not _zero(
                status.get("child_returncode")
            ):
                reasons.append(f"attempt status lacks zero return codes: {status_path}")
                continue
            recorded_audit = status.get("output_audit")
            if not isinstance(recorded_audit, dict) or recorded_audit.get("complete") is not True:
                reasons.append(f"attempt output audit is not complete: {status_path}")
                continue

            row_model = str(row.get("model_dir", "")).strip()
            status_model = str(status.get("model_dir", "")).strip()
            if not row_model or not status_model:
                reasons.append(f"missing model_dir for {group}/{seed}/attempt_{attempt:03d}")
                continue
            model_dir = Path(status_model).expanduser().resolve()
            if Path(row_model).expanduser().resolve() != model_dir:
                reasons.append(
                    f"status.tsv/status.json model_dir mismatch for {group}/{seed}"
                )
                continue
            if not _is_within(model_dir, attempt_root):
                reasons.append(
                    f"model_dir escapes its immutable attempt root for {group}/{seed}: {model_dir}"
                )
                continue

            missing_model_files: list[str] = []
            escaped_model_files: list[str] = []
            resolved_model_files: dict[str, Path] = {}
            for name in REQUIRED_MODEL_FILES:
                artifact = model_dir / name
                try:
                    resolved_artifact = artifact.resolve(strict=True)
                except OSError:
                    missing_model_files.append(name)
                    continue
                if not artifact.is_file() or artifact.stat().st_size <= 0:
                    missing_model_files.append(name)
                    continue
                if not _is_within(resolved_artifact, model_dir) or not _is_within(
                    resolved_artifact, attempt_root
                ):
                    escaped_model_files.append(name)
                    continue
                resolved_model_files[name] = resolved_artifact
            if missing_model_files:
                reasons.append(
                    f"missing/empty formal artifacts for {group}/{seed}: {missing_model_files}"
                )
                continue
            if escaped_model_files:
                reasons.append(
                    f"formal artifact symlink escapes model/attempt root for "
                    f"{group}/{seed}: {escaped_model_files}"
                )
                continue
            chains.append(
                FormalChain(
                    group=group,
                    seed=seed,
                    attempt=attempt,
                    model_dir=model_dir,
                    status_path=resolved_status,
                    mcmc_path=resolved_model_files["mcmc.pkl"],
                )
            )

    chain_keys = {(chain.group, chain.seed) for chain in chains}
    if chain_keys != expected:
        missing = sorted(expected - chain_keys)
        if missing:
            reasons.append(f"formal chain gate remains incomplete: {missing}")

    # The persisted audit proves what was present when the runner completed;
    # this second audit proves that the bytes currently selected for formal
    # post-processing still satisfy the fixed 2000/500/20 sampling contract.
    # Delay the relatively expensive probes until every structural condition
    # and the supervisor terminal state agree, so polling partial runs is cheap.
    if not reasons and len(chains) == len(expected):
        auditor = semantic_auditor or _live_semantic_audit
        audited_chains: list[FormalChain] = []
        for chain in chains:
            try:
                audit = auditor(chain)
            except Exception as exc:
                reasons.append(
                    f"live semantic output audit failed to execute for "
                    f"{chain.group}/{chain.seed}: {type(exc).__name__}: {exc}"
                )
                audited_chains.append(chain)
                continue
            audited_chains.append(
                FormalChain(
                    group=chain.group,
                    seed=chain.seed,
                    attempt=chain.attempt,
                    model_dir=chain.model_dir,
                    status_path=chain.status_path,
                    mcmc_path=chain.mcmc_path,
                    semantic_audit=audit,
                )
            )
            if audit.get("complete") is not True:
                reasons.append(
                    f"live semantic output audit is incomplete for "
                    f"{chain.group}/{chain.seed}: "
                    f"{audit.get('semantic_errors') or audit.get('missing_or_empty') or 'unknown error'}"
                )
        chains = audited_chains

    chains.sort(key=lambda item: (GROUPS.index(item.group), FORMAL_SEEDS.index(item.seed)))
    return GateResult(
        ready=not reasons and len(chains) == len(expected),
        supervisor_state=supervisor_state,
        reasons=tuple(reasons),
        chains=tuple(chains),
        checked_at_utc=checked_at,
    )


def build_formal_commands(
    gate: GateResult,
    output_root: Path,
    log_dir: Path,
    *,
    python: Path = MDSINE_PYTHON,
    project_root: Path = PROJECT_ROOT,
) -> tuple[CommandPlan, ...]:
    """Construct the only permitted formal post-processing sequence."""

    if not gate.ready:
        raise ValueError("Refusing to construct formal commands before the gate is ready")
    expected = {(group, seed) for group in GROUPS for seed in FORMAL_SEEDS}
    observed = {(chain.group, chain.seed) for chain in gate.chains}
    if observed != expected or len(gate.chains) != len(expected):
        raise ValueError("Formal gate does not contain exactly the fixed six-chain contract")

    output_root = Path(output_root).resolve()
    log_dir = Path(log_dir).resolve()
    project_root = Path(project_root).resolve()
    python_text = str(Path(python))
    tables = output_root / "tables"
    inputs = output_root / "inputs"

    ordered = {
        group: sorted(
            (chain for chain in gate.chains if chain.group == group),
            key=lambda item: FORMAL_SEEDS.index(item.seed),
        )
        for group in GROUPS
    }
    summarize: list[str] = [python_text, str(MODULE_DIR / "summarize_posteriors.py")]
    for chain in ordered["LDR"]:
        summarize.extend(("--ldr-mcmc", str(chain.mcmc_path)))
    for chain in ordered["HDR"]:
        summarize.extend(("--hdr-mcmc", str(chain.mcmc_path)))
    summarize.extend(
        (
            "--output-dir",
            str(tables),
            "--ldr-reads",
            str(inputs / "LDR" / "reads.tsv"),
            "--ldr-metadata",
            str(inputs / "LDR" / "metadata.tsv"),
            "--hdr-reads",
            str(inputs / "HDR" / "reads.tsv"),
            "--hdr-metadata",
            str(inputs / "HDR" / "metadata.tsv"),
            "--bf-threshold",
            "10",
            "--evidence-tier",
            "formal_multichain",
        )
    )

    longitudinal = (
        python_text,
        str(MODULE_DIR / "longitudinal_stats.py"),
        "--ldr-reads",
        str(inputs / "LDR" / "reads.tsv"),
        "--ldr-metadata",
        str(inputs / "LDR" / "metadata.tsv"),
        "--hdr-reads",
        str(inputs / "HDR" / "reads.tsv"),
        "--hdr-metadata",
        str(inputs / "HDR" / "metadata.tsv"),
        "--ldr-taxonomy",
        str(inputs / "LDR" / "taxonomy.tsv"),
        "--hdr-taxonomy",
        str(inputs / "HDR" / "taxonomy.tsv"),
        "--output-dir",
        str(tables),
        "--bootstrap",
        "2000",
        "--permutations",
        "9999",
        "--seed",
        "20261009",
    )
    render = (
        python_text,
        str(MODULE_DIR / "render_fig3.py"),
        "--project-root",
        str(project_root),
        "--output",
        str(output_root),
    )
    validate = (
        python_text,
        str(MODULE_DIR / "validate_fig3.py"),
        "--project-root",
        str(project_root),
        "--output",
        str(output_root),
        "--strict",
    )
    plans = (
        CommandPlan("summarize_posteriors", tuple(summarize), project_root, log_dir / "01_summarize_posteriors.log"),
        CommandPlan("longitudinal_stats", longitudinal, project_root, log_dir / "02_longitudinal_stats.log"),
        CommandPlan("render_fig3", render, project_root, log_dir / "03_render_fig3.log"),
        CommandPlan("validate_fig3_strict", validate, project_root, log_dir / "04_validate_fig3_strict.log"),
    )
    for plan in plans:
        forbidden = set(plan.argv) & set(FORBIDDEN_FORMAL_FLAGS)
        if forbidden:
            raise AssertionError(f"Forbidden development flags in formal plan: {forbidden}")
    return plans


def _validate_dependencies(plans: Sequence[CommandPlan]) -> None:
    python = Path(plans[0].argv[0])
    if not python.is_file() or not os.access(python, os.X_OK):
        raise FileNotFoundError(f"MDSINE2 Python is unavailable or not executable: {python}")
    for plan in plans:
        script = Path(plan.argv[1])
        if not script.is_file():
            raise FileNotFoundError(f"Post-processing script is missing: {script}")


def _archive_previous_state(state_path: Path, payload: dict[str, Any]) -> Path:
    history = state_path.parent / "finalization_history"
    history.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    destination = history / f"finalize_{stamp}.json"
    atomic_json(destination, payload)
    return destination


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_strict_qa(output_root: Path) -> tuple[bool, list[str]]:
    """Verify the current formal manifests, strict report, and rendered bytes."""

    output_root = output_root.resolve()
    reasons: list[str] = []
    report_path = output_root / "qa" / "validation_report.json"
    try:
        report = load_json(report_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return False, [f"cannot read strict QA report {report_path}: {exc}"]
    if report.get("strict") is not True:
        reasons.append("QA report was not produced in strict mode")
    if report.get("allow_development") is not False:
        reasons.append("QA report allowed development evidence")
    counts = report.get("status_counts")
    if not isinstance(counts, dict) or _integer(counts.get("FAIL", 0)) != 0:
        reasons.append(f"strict QA report contains failures: {counts!r}")
    checks = report.get("checks")
    if not isinstance(checks, list) or not checks:
        reasons.append("strict QA report has no recorded checks")
    else:
        failed_names: list[str] = []
        for item in checks:
            if not isinstance(item, dict):
                failed_names.append("malformed check")
            elif item.get("status") == "FAIL":
                failed_names.append(str(item.get("name", "unnamed")))
        if failed_names:
            reasons.append(f"strict QA contains failed checks: {failed_names}")
        passed_names = {
            str(item.get("name"))
            for item in checks
            if isinstance(item, dict) and item.get("status") == "PASS"
        }
        required_passes = {
            "LDR posterior is formal multichain",
            "HDR posterior is formal multichain",
            "render uses formal posterior",
        }
        missing_passes = sorted(required_passes - passed_names)
        if missing_passes:
            reasons.append(
                f"strict QA lacks required formal PASS checks: {missing_passes}"
            )

    posterior_path = output_root / "tables" / "posterior_summary_manifest.json"
    try:
        posterior = load_json(posterior_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        reasons.append(f"cannot read posterior manifest {posterior_path}: {exc}")
    else:
        manifest_groups = posterior.get("groups")
        if not isinstance(manifest_groups, dict):
            reasons.append("posterior manifest groups field is malformed")
            manifest_groups = {}
        for group in GROUPS:
            info = manifest_groups.get(group, {})
            if not isinstance(info, dict):
                reasons.append(f"posterior manifest entry is malformed for {group}")
                info = {}
            draws = info.get("posterior_draws_by_chain", [])
            if info.get("evidence_tier") != "formal_multichain":
                reasons.append(
                    f"posterior manifest is not formal_multichain for {group}"
                )
            if len(info.get("chains", [])) != len(FORMAL_SEEDS):
                reasons.append(f"posterior manifest does not contain three {group} chains")
            try:
                enough_draws = len(draws) == len(FORMAL_SEEDS) and all(
                    int(value) >= 1500 for value in draws
                )
            except (TypeError, ValueError):
                enough_draws = False
            if not enough_draws:
                reasons.append(
                    f"posterior manifest lacks >=1500 retained draws in all {group} chains"
                )

    render_path = output_root / "render_manifest.json"
    try:
        render = load_json(render_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        reasons.append(f"cannot read render manifest {render_path}: {exc}")
    else:
        if render.get("development_only") is not False:
            reasons.append("render manifest is development-only")
        expected_outputs = {
            (output_root / "figures" / f"Fig3_restored.{suffix}").resolve()
            for suffix in ("svg", "pdf", "png", "tif")
        }
        observed_outputs: set[Path] = set()
        records = render.get("outputs")
        if not isinstance(records, list):
            reasons.append("render manifest has no output records")
            records = []
        for record in records:
            if not isinstance(record, dict) or not str(record.get("path", "")).strip():
                reasons.append(f"invalid render output record: {record!r}")
                continue
            path = Path(str(record["path"])).expanduser()
            try:
                resolved = path.resolve(strict=True)
            except OSError as exc:
                reasons.append(f"rendered output is missing/unresolvable {path}: {exc}")
                continue
            observed_outputs.add(resolved)
            if not _is_within(resolved, output_root / "figures"):
                reasons.append(f"rendered output escapes figure root: {resolved}")
                continue
            size = resolved.stat().st_size
            if _integer(record.get("bytes")) != size:
                reasons.append(f"rendered output size changed since QA: {resolved}")
            expected_hash = str(record.get("sha256", ""))
            if not expected_hash or _sha256(resolved) != expected_hash:
                reasons.append(f"rendered output hash changed since QA: {resolved}")
        if observed_outputs != expected_outputs:
            reasons.append(
                "render manifest output set differs from the fixed SVG/PDF/PNG/TIFF contract"
            )
    return not reasons, reasons


def _revalidate_succeeded_state(
    previous: dict[str, Any],
    output_root: Path,
    inference_root: Path,
    *,
    semantic_auditor: SemanticAuditor | None = None,
) -> dict[str, Any]:
    """Fail closed before treating a prior ``succeeded`` record as idempotent."""

    reasons: list[str] = []
    output_root = output_root.resolve()
    inference_root = inference_root.resolve()
    if _integer(previous.get("schema_version")) != STATE_SCHEMA_VERSION:
        reasons.append(
            f"state schema is {previous.get('schema_version')!r}, expected {STATE_SCHEMA_VERSION}"
        )
    for key, expected in (
        ("output_root", output_root),
        ("inference_root", inference_root),
    ):
        raw = str(previous.get(key, "")).strip()
        if not raw or Path(raw).expanduser().resolve() != expected:
            reasons.append(f"stored {key} does not match this invocation")
    if previous.get("formal_contract") != _formal_contract():
        reasons.append("stored formal contract differs from the current fixed contract")

    commands = previous.get("commands")
    expected_commands = [
        "summarize_posteriors",
        "longitudinal_stats",
        "render_fig3",
        "validate_fig3_strict",
    ]
    if (
        not isinstance(commands, list)
        or len(commands) != len(expected_commands)
        or not all(isinstance(item, dict) for item in commands)
        or [item.get("name") for item in commands] != expected_commands
    ):
        reasons.append("stored command sequence is not the fixed formal sequence")
    elif any(item.get("state") != "succeeded" or not _zero(item.get("returncode")) for item in commands):
        reasons.append("not every stored formal command succeeded with return code zero")
    else:
        argv_lists = [item.get("argv") for item in commands]
        if not all(
            isinstance(argv, list) and all(isinstance(arg, str) for arg in argv)
            for argv in argv_lists
        ):
            reasons.append("stored command argv records are malformed")
        else:
            final_argv = argv_lists[-1]
            all_argv = {arg for argv in argv_lists for arg in argv}
            invalid = "--strict" not in final_argv or bool(
                set(FORBIDDEN_FORMAL_FLAGS) & all_argv
            )
            if invalid:
                reasons.append("stored command sequence did not enforce strict formal-only QA")

    gate = evaluate_formal_gate(
        inference_root,
        semantic_auditor=semantic_auditor,
    )
    if not gate.ready:
        reasons.extend(f"formal gate revalidation: {reason}" for reason in gate.reasons)
    qa_ok, qa_reasons = _verify_strict_qa(output_root)
    if not qa_ok:
        reasons.extend(f"strict QA revalidation: {reason}" for reason in qa_reasons)
    return {
        "valid": not reasons,
        "checked_at_utc": utc_now(),
        "reasons": reasons,
        "gate": gate.to_json(),
        "strict_qa_valid": qa_ok,
    }


def _initial_payload(
    state_path: Path,
    inference_root: Path,
    output_root: Path,
    poll_seconds: float,
    previous_state: dict[str, Any] | None,
) -> dict[str, Any]:
    return {
        "schema_version": STATE_SCHEMA_VERSION,
        "state": "waiting",
        "run_id": dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        + f"-pid{os.getpid()}",
        "pid": os.getpid(),
        "started_at_utc": utc_now(),
        "updated_at_utc": utc_now(),
        "state_path": str(state_path),
        "inference_root": str(inference_root),
        "output_root": str(output_root),
        "poll_seconds": poll_seconds,
        "formal_contract": _formal_contract(),
        "previous_state_archived": previous_state,
        "commands": [],
    }


def _run_commands(
    plans: Sequence[CommandPlan],
    payload: dict[str, Any],
    state_path: Path,
) -> None:
    _validate_dependencies(plans)
    payload["commands"] = [plan.to_json() for plan in plans]
    payload["state"] = "running"
    payload["postprocessing_started_at_utc"] = utc_now()
    payload["updated_at_utc"] = utc_now()
    atomic_json(state_path, payload)

    environment = os.environ.copy()
    environment.update(
        {
            "MPLCONFIGDIR": "/tmp/mdsine2_mpl",
            "PYTHONDONTWRITEBYTECODE": "1",
        }
    )
    for index, plan in enumerate(plans):
        command_state = payload["commands"][index]
        command_state["state"] = "running"
        command_state["started_at_utc"] = utc_now()
        payload["updated_at_utc"] = utc_now()
        atomic_json(state_path, payload)

        plan.log_path.parent.mkdir(parents=True, exist_ok=True)
        with plan.log_path.open("w", encoding="utf-8") as log_handle:
            log_handle.write(f"started_at_utc={command_state['started_at_utc']}\n")
            log_handle.write(f"cwd={plan.cwd}\n")
            log_handle.write(f"argv={shlex.join(plan.argv)}\n\n")
            log_handle.flush()
            result = subprocess.run(
                list(plan.argv),
                cwd=str(plan.cwd),
                env=environment,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        command_state["returncode"] = result.returncode
        command_state["finished_at_utc"] = utc_now()
        command_state["state"] = "succeeded" if result.returncode == 0 else "failed"
        payload["updated_at_utc"] = utc_now()
        atomic_json(state_path, payload)
        if result.returncode != 0:
            raise RuntimeError(
                f"{plan.name} failed with return code {result.returncode}; see {plan.log_path}"
            )


def run_finalizer(args: argparse.Namespace) -> int:
    output_root = args.output.resolve()
    inference_root = (args.inference_root or (output_root / "inference")).resolve()
    # State and lock locations are intentionally not configurable.  The lock
    # is scoped to the output root so two callers cannot select different
    # inference roots yet concurrently mutate the same tables and figures.
    state_path = (output_root / "formal_finalization_status.json").resolve()
    lock_path = (output_root / ".formal_finalization.lock").resolve()
    if args.poll_seconds <= 0:
        raise ValueError("--poll-seconds must be positive")
    if args.timeout_seconds is not None and args.timeout_seconds <= 0:
        raise ValueError("--timeout-seconds must be positive when supplied")

    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock_handle:
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                f"Another Fig. 3 finalizer already holds the lock: {lock_path}"
            ) from exc

        previous: dict[str, Any] | None = None
        archived: str | None = None
        if state_path.is_file():
            previous = load_json(state_path)
            previous_state = previous.get("state")
            if previous_state == "succeeded":
                revalidation = _revalidate_succeeded_state(
                    previous,
                    output_root,
                    inference_root,
                )
                if revalidation["valid"]:
                    previous["last_idempotent_revalidation"] = revalidation
                    previous["updated_at_utc"] = utc_now()
                    atomic_json(state_path, previous)
                    print(
                        "Formal Fig. 3 post-processing already succeeded and was "
                        f"revalidated: {state_path}"
                    )
                    return 0
                previous.update(
                    {
                        "state": "failed",
                        "failure_phase": "succeeded_state_revalidation",
                        "succeeded_state_revalidation": revalidation,
                        "error_type": "SucceededStateInvalidated",
                        "error": "; ".join(revalidation["reasons"]),
                        "finished_at_utc": utc_now(),
                        "updated_at_utc": utc_now(),
                    }
                )
                atomic_json(state_path, previous)
                previous_state = "failed"
            if previous_state in {"running", "failed"} and not args.retry_failed:
                raise RuntimeError(
                    f"Previous finalization state is {previous_state!r}; inspect {state_path} "
                    "and use --retry-failed only after resolving the recorded failure"
                )
            if args.retry_failed:
                archived = str(_archive_previous_state(state_path, previous))

        payload = _initial_payload(
            state_path,
            inference_root,
            output_root,
            args.poll_seconds,
            {"path": archived, "state": previous.get("state")} if archived and previous else None,
        )
        atomic_json(state_path, payload)
        monotonic_start = time.monotonic()
        try:
            while True:
                gate = evaluate_formal_gate(inference_root)
                payload["gate"] = gate.to_json()
                payload["updated_at_utc"] = utc_now()
                atomic_json(state_path, payload)
                if gate.ready:
                    break
                if args.once:
                    print(json.dumps(gate.to_json(), indent=2, ensure_ascii=False))
                    return 4
                if (
                    args.timeout_seconds is not None
                    and time.monotonic() - monotonic_start >= args.timeout_seconds
                ):
                    raise TimeoutError(
                        f"Formal inference gate did not become ready within {args.timeout_seconds}s"
                    )
                time.sleep(args.poll_seconds)

            log_dir = output_root / "finalize_logs" / payload["run_id"]
            plans = build_formal_commands(
                gate,
                output_root,
                log_dir,
                python=args.python,
                project_root=args.project_root,
            )
            _run_commands(plans, payload, state_path)
            payload["state"] = "succeeded"
            payload["finished_at_utc"] = utc_now()
            payload["updated_at_utc"] = payload["finished_at_utc"]
            atomic_json(state_path, payload)
            print(f"Formal Fig. 3 post-processing succeeded: {state_path}")
            return 0
        except BaseException as exc:
            payload.update(
                {
                    "state": "failed",
                    "finished_at_utc": utc_now(),
                    "updated_at_utc": utc_now(),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "traceback": traceback.format_exc(),
                }
            )
            atomic_json(state_path, payload)
            raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--inference-root",
        type=Path,
        default=None,
        help="Defaults to <output>/inference.",
    )
    parser.add_argument("--python", type=Path, default=MDSINE_PYTHON)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--poll-seconds", type=float, default=60.0)
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=None,
        help="Optional finite wait; the default waits indefinitely.",
    )
    parser.add_argument(
        "--once",
        action="store_true",
        help="Evaluate once and exit 4 without running commands when not ready.",
    )
    parser.add_argument(
        "--retry-failed",
        action="store_true",
        help="Archive a prior failed/stale-running record and explicitly start a new attempt.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    return run_finalizer(build_parser().parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
