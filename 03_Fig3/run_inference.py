#!/usr/bin/env python3
"""Plan, validate, launch, supervise, and audit Fig. 3 MDSINE2 chains.

The launcher uses the manuscript repository's no-qPCR implementation verbatim.
It never passes ``--resume``: a failed chain is restarted in a new immutable
``attempt_NNN`` directory so a partial trace cannot be mistaken for a resume.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from dataclasses import asdict, dataclass
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import pickle
import platform
import re
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from typing import Any, Iterable, Sequence


SCRIPT_PATH = Path(__file__).resolve()
MODULE_DIR = SCRIPT_PATH.parent
PROJECT_ROOT = MODULE_DIR.parents[1]
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "reproduction" / "outputs" / "fig3_mdsine2"
DEFAULT_INPUT_ROOT = DEFAULT_OUTPUT_ROOT / "inputs"
DEFAULT_INFERENCE_ROOT = DEFAULT_OUTPUT_ROOT / "inference"

MDSINE_PYTHON = Path("/home/qqt/.conda/envs/mdsine2/bin/python")
OFFICIAL_RUNNER = Path(
    "/home/qqt/MDSINE2_Server_Package/MDSINE2_Paper/scripts/"
    "semisynthetic2/inference/mdsine2/helpers/mdsine2_ra.py"
)
PAPER_REPOSITORY = Path("/home/qqt/MDSINE2_Server_Package/MDSINE2_Paper")
MDSINE_REPOSITORY = Path("/home/qqt/MDSINE2_Server_Package/MDSINE2")
WRAPPER = MODULE_DIR / "mdsine2_noqpcr_wrapper.py"
NO_QPCR_LOADER = MODULE_DIR / "mdsine2_noqpcr_loader.py"

GROUPS = ("LDR", "HDR")
FORMAL_SEEDS = (20261009, 20261010, 20261011)
SIDECARS = (
    "metadata.tsv",
    "reads.tsv",
    "taxonomy.tsv",
    "perturbations.tsv",
    "qpcr_structural_placeholder.tsv",
)
REQUIRED_MODEL_FILES = (
    "mcmc.pkl",
    "metadata.txt",
    "params.pkl",
    "traces.hdf5",
)
# ``initialize_graph_no_qpcr`` in the pinned manuscript runner does not call
# the standard MDSINE2 helper that serializes ``subjset.pkl``.  Requiring that
# file therefore turns a successful no-qPCR run into a false failure.
INTENTIONALLY_ABSENT_NO_QPCR_FILES = ("subjset.pkl",)
REQUIRED_TRACE_DATASETS = (
    "Growth parameter",
    "Interactions object",
    "ciprofloxacin",
)
THREADS_PER_GROUP = 16
# The host has two 96-core NUMA nodes.  Keeping each group on 16 physical cores
# from a different node makes the allocations disjoint and avoids shared-core
# contention between the two simultaneous fits.
CPU_AFFINITY = {
    "LDR": "0-15",
    "HDR": "96-111",
}
# Three formal chains run concurrently inside each 16-CPU group allocation.
# Giving every process the full group mask together with OMP_PROC_BIND=TRUE
# pins every Python leader to the first CPU in that mask (all LDR leaders to
# CPU 0 and all HDR leaders to CPU 96).  The largely serial clustering update
# then time-slices on one CPU.  Partitioning the group mask by seed preserves
# the user-requested 16 CPUs per group while preventing leader contention and
# BLAS/OpenMP oversubscription between chains.
FORMAL_CHAIN_CPU_AFFINITY = {
    "LDR": {
        20261009: "0-4",
        20261010: "5-9",
        20261011: "10-15",
    },
    "HDR": {
        20261009: "96-100",
        20261010: "101-105",
        20261011: "106-111",
    },
}
NO_QPCR_MULTIPROCESSING_NOTE = (
    "The official no-qPCR runner is kept at --multiprocessing 0 because its "
    "FilteringWithoutQPCR.update() unconditionally indexes self.pool[ridx]; "
    "with --multiprocessing 1 self.pool is PersistentPool, which has no "
    "__getitem__. Thread-capable NumPy/BLAS/Numba components receive the "
    "per-chain thread allowance inside the fixed 16-CPU group allocation."
)


@dataclass(frozen=True)
class Profile:
    name: str
    n_samples: int
    burnin: int
    checkpoint: int
    seeds: tuple[int, ...]


PROFILES = {
    "pilot": Profile("pilot", n_samples=140, burnin=40, checkpoint=20, seeds=(20261009,)),
    "formal": Profile(
        "formal",
        n_samples=2000,
        burnin=500,
        checkpoint=20,
        seeds=FORMAL_SEEDS,
    ),
}


@dataclass(frozen=True)
class ChainSpec:
    group: str
    profile: Profile
    seed: int
    input_root: Path
    inference_root: Path
    attempt: int

    @property
    def input_dir(self) -> Path:
        return self.input_root / self.group

    @property
    def input_path(self) -> Path:
        return self.input_dir / f"{self.group}_study.pkl"

    @property
    def logical_root(self) -> Path:
        return self.inference_root / self.profile.name / self.group / f"seed_{self.seed}"

    @property
    def attempt_root(self) -> Path:
        return self.logical_root / f"attempt_{self.attempt:03d}"

    @property
    def study_name(self) -> str:
        return f"fig3-{self.group.lower()}-{self.profile.name}-seed{self.seed}"

    @property
    def model_dir(self) -> Path:
        return self.attempt_root / self.study_name

    @property
    def status_path(self) -> Path:
        return self.attempt_root / "status.json"

    @property
    def audit_path(self) -> Path:
        return self.attempt_root / "input_and_runner_audit.json"

    @property
    def log_path(self) -> Path:
        return self.attempt_root / "run.log"

    @property
    def command_path(self) -> Path:
        return self.attempt_root / "command.json"

    @property
    def session_name(self) -> str:
        prefix = "p" if self.profile.name == "pilot" else "f"
        return f"fig3m2-{prefix}-{self.group.lower()}-{self.seed}-a{self.attempt:03d}"

    @property
    def cpu_affinity(self) -> str:
        if self.profile.name == "formal":
            return FORMAL_CHAIN_CPU_AFFINITY[self.group][self.seed]
        return CPU_AFFINITY[self.group]

    @property
    def threads_per_chain(self) -> int:
        return len(expand_cpu_list(self.cpu_affinity))


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(value, encoding="utf-8")
    os.replace(temporary, path)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def subprocess_capture(command: Sequence[str], cwd: Path | None = None) -> tuple[int, str]:
    result = subprocess.run(
        list(command),
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return result.returncode, result.stdout.strip()


def git_provenance(repository: Path) -> dict[str, Any]:
    if not repository.is_dir():
        return {"path": str(repository), "available": False}
    code, commit = subprocess_capture(("git", "-C", str(repository), "rev-parse", "HEAD"))
    if code != 0:
        return {
            "path": str(repository.resolve()),
            "available": False,
            "error": commit,
        }
    _, status = subprocess_capture(
        ("git", "-C", str(repository), "status", "--porcelain=v1", "--untracked-files=all")
    )
    diff = subprocess.run(
        ("git", "-C", str(repository), "diff", "--binary", "HEAD"),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    ).stdout
    return {
        "path": str(repository.resolve()),
        "available": True,
        "commit": commit,
        "dirty": bool(status),
        "status_porcelain": status.splitlines(),
        "tracked_diff_sha256": sha256_bytes(diff),
        "tracked_diff_bytes": len(diff),
    }


def environment_provenance() -> dict[str, Any]:
    probe = (
        "import importlib.metadata as m,json,platform,sys;"
        "d=sorted((x.metadata.get('Name') or '',x.version) for x in m.distributions());"
        "print(json.dumps({'python':sys.version,'executable':sys.executable,"
        "'platform':platform.platform(),'packages':d},sort_keys=True))"
    )
    code, output = subprocess_capture((str(MDSINE_PYTHON), "-c", probe))
    if code == 0:
        try:
            return json.loads(output.splitlines()[-1])
        except (IndexError, json.JSONDecodeError):
            pass
    return {"error": output, "executable": str(MDSINE_PYTHON)}


def expand_cpu_list(value: str) -> set[int]:
    """Expand Linux taskset list syntax (for example ``0-3,8``)."""
    cpus: set[int] = set()
    for part in value.split(","):
        part = part.strip()
        if not part:
            raise ValueError(f"Empty CPU-list component in {value!r}")
        if "-" in part:
            start_text, end_text = part.split("-", 1)
            start, end = int(start_text), int(end_text)
            if start < 0 or end < start:
                raise ValueError(f"Invalid CPU range {part!r}")
            cpus.update(range(start, end + 1))
        else:
            cpu = int(part)
            if cpu < 0:
                raise ValueError(f"Invalid CPU {cpu}")
            cpus.add(cpu)
    return cpus


def validate_cpu_contract() -> None:
    allocations = {group: expand_cpu_list(value) for group, value in CPU_AFFINITY.items()}
    for group, cpus in allocations.items():
        if len(cpus) != THREADS_PER_GROUP:
            raise ValueError(
                f"{group} affinity {CPU_AFFINITY[group]!r} contains {len(cpus)} CPUs; "
                f"expected {THREADS_PER_GROUP}"
            )
    overlap = allocations["LDR"] & allocations["HDR"]
    if overlap:
        raise ValueError(f"LDR/HDR CPU allocations overlap: {sorted(overlap)}")
    formal_chain_sets: list[set[int]] = []
    for group, group_cpus in allocations.items():
        chain_sets = {
            seed: expand_cpu_list(value)
            for seed, value in FORMAL_CHAIN_CPU_AFFINITY[group].items()
        }
        formal_chain_sets.extend(chain_sets.values())
        if set(chain_sets) != set(FORMAL_SEEDS):
            raise ValueError(
                f"{group} formal chain affinity seeds do not match {FORMAL_SEEDS}"
            )
        chain_union: set[int] = set()
        for seed, cpus in chain_sets.items():
            if chain_union & cpus:
                raise ValueError(
                    f"{group} formal chain CPU allocation overlaps at seed {seed}"
                )
            chain_union.update(cpus)
        if chain_union != group_cpus:
            raise ValueError(
                f"{group} formal chain allocations do not exactly cover the "
                f"16-CPU group mask"
            )
    if hasattr(os, "sched_getaffinity"):
        allowed = set(os.sched_getaffinity(0))
        requested_union = set().union(*allocations.values())
        # The top-level supervisor sees the host-wide mask.  A chain launcher is
        # itself intentionally taskset-bound and therefore sees exactly one of
        # the two configured group masks.
        running_inside_known_mask = allowed in [*allocations.values(), *formal_chain_sets]
        unavailable = requested_union - allowed
        if unavailable and not running_inside_known_mask:
            raise RuntimeError(
                f"Requested Fig. 3 CPUs are outside this process's allowed mask: "
                f"{sorted(unavailable)}"
            )


def write_provenance(inference_root: Path) -> Path:
    destination = inference_root / "provenance.json"
    payload = {
        "schema_version": 1,
        "captured_at_utc": utc_now(),
        "host": {
            "node": platform.node(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
        },
        "fixed_model_contract": {
            "runner": str(OFFICIAL_RUNNER),
            "no_qpcr": True,
            "biomass_mean": 1e11,
            "biomass_cv": 0.10,
            "negative_binomial_a0": 0.0025,
            "negative_binomial_a1": 0.025,
            "interaction_indicator_prior": "weak-agnostic",
            "perturbation_indicator_prior": "weak-agnostic",
            "resume_supported": False,
            "mdsine2_multiprocessing": 0,
            "mdsine2_multiprocessing_note": NO_QPCR_MULTIPROCESSING_NOTE,
            "threads_per_group": THREADS_PER_GROUP,
            "cpu_affinity": CPU_AFFINITY,
            "formal_chain_cpu_affinity": FORMAL_CHAIN_CPU_AFFINITY,
            "formal_chain_affinity_rationale": (
                "Disjoint per-seed partitions prevent OMP_PROC_BIND from pinning "
                "all concurrently running chain leaders to the same group-first CPU."
            ),
            "cpu_sets_disjoint": not bool(
                expand_cpu_list(CPU_AFFINITY["LDR"])
                & expand_cpu_list(CPU_AFFINITY["HDR"])
            ),
        },
        "files": {
            str(path): sha256_file(path)
            for path in (OFFICIAL_RUNNER, WRAPPER, NO_QPCR_LOADER, SCRIPT_PATH)
            if path.is_file()
        },
        "git": {
            "MDSINE2": git_provenance(MDSINE_REPOSITORY),
            "MDSINE2_Paper": git_provenance(PAPER_REPOSITORY),
        },
        "environment": environment_provenance(),
    }
    atomic_json(destination, payload)
    return destination


def validate_static_dependencies() -> None:
    problems = []
    for label, path in (
        ("lowercase MDSINE2 environment Python", MDSINE_PYTHON),
        ("official no-qPCR runner", OFFICIAL_RUNNER),
        ("validated wrapper", WRAPPER),
    ):
        if not path.is_file():
            problems.append(f"{label} not found: {path}")
    if shutil.which("tmux") is None:
        problems.append("tmux is not available on PATH")
    if shutil.which("taskset") is None:
        problems.append("taskset is not available on PATH")
    if problems:
        raise FileNotFoundError("\n".join(problems))
    validate_cpu_contract()


def validate_input_bundle(group: str, input_root: Path) -> list[Path]:
    input_dir = input_root / group
    required = [input_dir / f"{group}_study.pkl", *[input_dir / x for x in SIDECARS]]
    missing = [path for path in required if not path.is_file()]
    empty = [path for path in required if path.is_file() and path.stat().st_size == 0]
    if missing or empty:
        parts = []
        if missing:
            parts.append("missing=" + ", ".join(map(str, missing)))
        if empty:
            parts.append("empty=" + ", ".join(map(str, empty)))
        raise FileNotFoundError(f"Invalid {group} input bundle: {'; '.join(parts)}")
    return required


def official_runner_arguments(spec: ChainSpec) -> list[str]:
    return [
        "--input",
        str(spec.input_path.resolve()),
        "--negbin",
        "0.0025",
        "0.025",
        "--seed",
        str(spec.seed),
        "--burnin",
        str(spec.profile.burnin),
        "--n-samples",
        str(spec.profile.n_samples),
        "--checkpoint",
        str(spec.profile.checkpoint),
        "--multiprocessing",
        "0",
        "--rename-study",
        spec.study_name,
        "--basepath",
        str(spec.attempt_root.resolve()),
        "--interaction-ind-prior",
        "weak-agnostic",
        "--perturbation-ind-prior",
        "weak-agnostic",
        "--log-every",
        "20",
    ]


def wrapper_command(spec: ChainSpec, validate_only: bool = False) -> list[str]:
    command = [
        str(MDSINE_PYTHON),
        str(WRAPPER),
        "--official-runner",
        str(OFFICIAL_RUNNER),
        "--audit-json",
        str(spec.audit_path.resolve()),
        "--expected-group",
        spec.group,
        "--expected-subjects",
        "15",
        "--expected-taxa",
        "259",
        "--perturbation-start",
        "0.5",
        "--perturbation-end",
        "5.5",
        "--expected-biomass-mean",
        "100000000000",
        "--expected-biomass-cv",
        "0.10",
    ]
    if validate_only:
        command.append("--validate-only")
    command.extend(("--", *official_runner_arguments(spec)))
    return command


def affinity_wrapped_command(spec: ChainSpec, validate_only: bool = False) -> list[str]:
    """Return a command whose entire MDSINE2 process tree is CPU-bound."""
    return [
        "taskset",
        "--cpu-list",
        spec.cpu_affinity,
        *wrapper_command(spec, validate_only=validate_only),
    ]


def execution_environment(spec: ChainSpec) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update(
        {
            "NUMBA_CACHE_DIR": "/tmp/mdsine2_numba",
            "MPLCONFIGDIR": "/tmp/mdsine2_mpl",
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONUNBUFFERED": "1",
            "OMP_NUM_THREADS": str(spec.threads_per_chain),
            "MKL_NUM_THREADS": str(spec.threads_per_chain),
            "OPENBLAS_NUM_THREADS": str(spec.threads_per_chain),
            "NUMEXPR_NUM_THREADS": str(spec.threads_per_chain),
            "NUMBA_NUM_THREADS": str(spec.threads_per_chain),
            "OMP_DYNAMIC": "FALSE",
            "MKL_DYNAMIC": "FALSE",
            "OMP_PROC_BIND": "TRUE",
            "OMP_PLACES": "cores",
            "LOG_FILEPATH": str(spec.log_path.resolve()),
        }
    )
    return environment


def base_status(spec: ChainSpec, state: str) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "updated_at_utc": utc_now(),
        "state": state,
        "group": spec.group,
        "profile": spec.profile.name,
        "seed": spec.seed,
        "attempt": spec.attempt,
        "session": spec.session_name,
        "input": str(spec.input_path.resolve()),
        "attempt_root": str(spec.attempt_root.resolve()),
        "model_dir": str(spec.model_dir.resolve()),
        "log": str(spec.log_path.resolve()),
        "parameters": {
            **asdict(spec.profile),
            "negative_binomial_a0": 0.0025,
            "negative_binomial_a1": 0.025,
            "biomass_mean": 1e11,
            "biomass_cv": 0.10,
            "multiprocessing": 0,
            "multiprocessing_note": NO_QPCR_MULTIPROCESSING_NOTE,
            "threads_per_group": THREADS_PER_GROUP,
            "threads_per_chain": spec.threads_per_chain,
            "cpu_affinity": spec.cpu_affinity,
            "interaction_indicator_prior": "weak-agnostic",
            "perturbation_indicator_prior": "weak-agnostic",
        },
        "environment_overrides": {
            key: value
            for key, value in execution_environment(spec).items()
            if key
            in {
                "NUMBA_CACHE_DIR",
                "MPLCONFIGDIR",
                "PYTHONDONTWRITEBYTECODE",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "NUMBA_NUM_THREADS",
                "OMP_DYNAMIC",
                "MKL_DYNAMIC",
                "OMP_PROC_BIND",
                "OMP_PLACES",
                "LOG_FILEPATH",
            }
        },
        "resume_used": False,
    }


def attempts(logical_root: Path) -> list[int]:
    found: list[int] = []
    if logical_root.is_dir():
        for path in logical_root.glob("attempt_[0-9][0-9][0-9]"):
            try:
                found.append(int(path.name.removeprefix("attempt_")))
            except ValueError:
                continue
    return sorted(found)


def latest_status_for(group: str, profile: Profile, seed: int, input_root: Path, inference_root: Path) -> tuple[ChainSpec | None, dict[str, Any] | None]:
    template = ChainSpec(group, profile, seed, input_root, inference_root, attempt=1)
    numbers = attempts(template.logical_root)
    if not numbers:
        return None, None
    spec = ChainSpec(group, profile, seed, input_root, inference_root, attempt=numbers[-1])
    if not spec.status_path.is_file():
        return spec, None
    return spec, load_json(spec.status_path)


def reserve_attempt(
    group: str,
    profile: Profile,
    seed: int,
    input_root: Path,
    inference_root: Path,
    max_attempts: int,
) -> tuple[ChainSpec, bool]:
    previous_spec, previous_status = latest_status_for(
        group, profile, seed, input_root, inference_root
    )
    if previous_spec is not None:
        state = (previous_status or {}).get("state", "unknown")
        if state == "succeeded":
            return previous_spec, False
        if state in {"planned", "running"}:
            return previous_spec, False
        next_attempt = previous_spec.attempt + 1
    else:
        next_attempt = 1
    if next_attempt > max_attempts:
        raise RuntimeError(
            f"{profile.name}/{group}/seed {seed} exhausted max_attempts={max_attempts}; "
            "partial output was retained and will not be resumed or overwritten"
        )
    spec = ChainSpec(group, profile, seed, input_root, inference_root, next_attempt)
    spec.logical_root.mkdir(parents=True, exist_ok=True)
    try:
        spec.attempt_root.mkdir()
    except FileExistsError as exc:
        raise RuntimeError(f"Attempt reservation race at {spec.attempt_root}") from exc
    status = base_status(spec, "planned")
    status["created_at_utc"] = utc_now()
    atomic_json(spec.status_path, status)
    atomic_json(
        spec.command_path,
        {
            "created_at_utc": utc_now(),
            "argv": affinity_wrapped_command(spec),
            "shell_escaped_display_only": shlex.join(affinity_wrapped_command(spec)),
            "cwd": str(PAPER_REPOSITORY),
            "resume_deliberately_omitted": True,
            "cpu_affinity": spec.cpu_affinity,
            "threads_per_group": THREADS_PER_GROUP,
            "threads_per_chain": spec.threads_per_chain,
        },
    )
    atomic_text(spec.logical_root / "latest_attempt.txt", f"{spec.attempt}\n")
    return spec, True


def _file_edges(path: Path, width: int = 16) -> tuple[bytes, bytes]:
    """Read a bounded prefix/suffix without materialising a multi-GB pickle."""

    with path.open("rb") as handle:
        prefix = handle.read(width)
        handle.seek(max(0, path.stat().st_size - width))
        suffix = handle.read(width)
    return prefix, suffix


def _audit_pickle_envelope(path: Path) -> tuple[dict[str, Any], list[str]]:
    """Check the binary-pickle envelope without loading an untrusted object."""

    details: dict[str, Any] = {"path": str(path.resolve())}
    errors: list[str] = []
    try:
        prefix, suffix = _file_edges(path)
    except OSError as exc:
        return details, [f"{path.name}: cannot read pickle envelope: {exc}"]
    protocol = prefix[1] if len(prefix) >= 2 and prefix[:1] == b"\x80" else None
    details.update(
        {
            "binary_pickle_protocol": protocol,
            "ends_with_stop_opcode": suffix.endswith(b"."),
        }
    )
    if protocol is None or not 2 <= protocol <= pickle.HIGHEST_PROTOCOL:
        errors.append(f"{path.name}: not a supported binary pickle")
    if not suffix.endswith(b"."):
        errors.append(f"{path.name}: pickle has no terminal STOP opcode")
    return details, errors


def _metadata_value(text: str, label: str) -> str | None:
    match = re.search(rf"(?m)^{re.escape(label)}:\s*(.*?)\s*$", text)
    return match.group(1) if match else None


def _audit_metadata(path: Path, spec: ChainSpec) -> tuple[dict[str, Any], list[str]]:
    details: dict[str, Any] = {"path": str(path.resolve())}
    errors: list[str] = []
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        return details, [f"metadata.txt: cannot decode UTF-8 metadata: {exc}"]

    expected_integer_fields = {
        "Random seed": spec.seed,
        "Total number of Gibb steps": spec.profile.n_samples,
        "Number of Gibb steps for burn-in": spec.profile.burnin,
    }
    observed: dict[str, Any] = {}
    for label, expected in expected_integer_fields.items():
        raw = _metadata_value(text, label)
        try:
            value = int(raw) if raw is not None else None
        except ValueError:
            value = None
        observed[label] = value
        if value != expected:
            errors.append(
                f"metadata.txt: {label!r} is {raw!r}; expected {expected}"
            )

    saved_location = _metadata_value(text, "Saved location")
    observed["Saved location"] = saved_location
    if saved_location is None or Path(saved_location).resolve() != spec.model_dir.resolve():
        errors.append(
            "metadata.txt: saved location does not match the immutable model directory"
        )

    inference_time = _metadata_value(text, "Time for inference")
    hours: float | None = None
    if inference_time is not None:
        match = re.fullmatch(
            r"([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?)\s+hours",
            inference_time,
        )
        if match:
            hours = float(match.group(1))
    observed["Time for inference hours"] = hours
    # This line is appended only after md2.run_graph returns successfully.  It
    # distinguishes a completed run from preallocated/partially checkpointed
    # HDF5 and pickle files.
    if hours is None or not math.isfinite(hours) or hours < 0:
        errors.append(
            "metadata.txt: missing or invalid terminal 'Time for inference' marker"
        )
    details["observed"] = observed
    return details, errors


def _audit_params_pickle(path: Path, spec: ChainSpec) -> tuple[dict[str, Any], list[str]]:
    """Load the small trusted config pickle in the pinned MDSINE2 environment."""

    details: dict[str, Any] = {"path": str(path.resolve())}
    probe = r'''
import json
import pickle
import sys

with open(sys.argv[1], "rb") as handle:
    value = pickle.load(handle)

names = (
    "SEED", "N_SAMPLES", "BURNIN", "CHECKPOINT", "MODEL_PATH",
    "OUTPUT_BASEPATH", "NEGBIN_A0", "NEGBIN_A1", "MP_FILTERING",
)
observed = {}
for name in names:
    item = getattr(value, name, None)
    if isinstance(item, (str, int, float, bool)) or item is None:
        observed[name] = item
    else:
        observed[name] = str(item)
print("FIG3_PARAMS_AUDIT=" + json.dumps({
    "class": type(value).__module__ + "." + type(value).__qualname__,
    "observed": observed,
}, sort_keys=True))
'''
    environment = os.environ.copy()
    environment.update(
        {
            "NUMBA_CACHE_DIR": "/tmp/mdsine2_numba",
            "MPLCONFIGDIR": "/tmp/mdsine2_mpl",
            "PYTHONDONTWRITEBYTECODE": "1",
        }
    )
    try:
        result = subprocess.run(
            (str(MDSINE_PYTHON), "-c", probe, str(path.resolve())),
            cwd=str(PAPER_REPOSITORY),
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=120,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return details, [f"params.pkl: semantic probe failed to execute: {exc}"]
    marker = "FIG3_PARAMS_AUDIT="
    payload_line = next(
        (line[len(marker) :] for line in result.stdout.splitlines() if line.startswith(marker)),
        None,
    )
    if result.returncode != 0 or payload_line is None:
        details["probe_returncode"] = result.returncode
        details["probe_output_tail"] = "\n".join(result.stdout.splitlines()[-20:])
        return details, ["params.pkl: MDSINE2 semantic probe did not complete"]
    try:
        payload = json.loads(payload_line)
    except json.JSONDecodeError as exc:  # pragma: no cover - controlled probe output
        return details, [f"params.pkl: semantic probe returned invalid JSON: {exc}"]

    details.update(payload)
    observed = payload["observed"]
    expected = {
        "SEED": spec.seed,
        "N_SAMPLES": spec.profile.n_samples,
        "BURNIN": spec.profile.burnin,
        "CHECKPOINT": spec.profile.checkpoint,
        "MODEL_PATH": str(spec.model_dir.resolve()),
        "OUTPUT_BASEPATH": str(spec.model_dir.resolve()),
        "NEGBIN_A0": 0.0025,
        "NEGBIN_A1": 0.025,
    }
    errors: list[str] = []
    for key, expected_value in expected.items():
        observed_value = observed.get(key)
        if isinstance(expected_value, float) and observed_value is not None:
            try:
                matches = math.isclose(
                    float(observed_value), expected_value, rel_tol=0, abs_tol=1e-15
                )
            except (TypeError, ValueError):
                matches = False
        else:
            matches = observed_value == expected_value
        if not matches:
            errors.append(
                f"params.pkl: {key} is {observed_value!r}; expected {expected_value!r}"
            )
    return details, errors


def _audit_hdf5(path: Path, spec: ChainSpec) -> tuple[dict[str, Any], list[str], list[str]]:
    details: dict[str, Any] = {"path": str(path.resolve())}
    errors: list[str] = []
    warnings: list[str] = []
    try:
        prefix, _ = _file_edges(path, width=8)
    except OSError as exc:
        return details, [f"traces.hdf5: cannot read file signature: {exc}"], warnings
    details["hdf5_signature"] = prefix == b"\x89HDF\r\n\x1a\n"
    if not details["hdf5_signature"]:
        errors.append("traces.hdf5: invalid HDF5 signature")
        return details, errors, warnings

    try:
        import h5py
    except ImportError:  # pragma: no cover - h5py is installed in the workflow env
        details["h5py_semantic_audit"] = "unavailable"
        warnings.append(
            "traces.hdf5: h5py unavailable; only the HDF5 signature was checked"
        )
        return details, errors, warnings

    try:
        with h5py.File(path, "r") as handle:
            attrs = {
                key: int(handle.attrs[key]) if key in handle.attrs else None
                for key in ("n_samples", "burnin", "ckpt")
            }
            datasets: dict[str, list[int]] = {}

            def collect(name: str, obj: Any) -> None:
                if isinstance(obj, h5py.Dataset):
                    datasets[name] = list(obj.shape)

            handle.visititems(collect)
    except (OSError, RuntimeError, ValueError) as exc:
        return details, [f"traces.hdf5: cannot open complete HDF5 structure: {exc}"], warnings

    details.update({"h5py_semantic_audit": "performed", "attrs": attrs, "datasets": datasets})
    expected_attrs = {
        "n_samples": spec.profile.n_samples,
        "burnin": spec.profile.burnin,
        "ckpt": spec.profile.checkpoint,
    }
    for key, expected in expected_attrs.items():
        if attrs.get(key) != expected:
            errors.append(
                f"traces.hdf5: root attr {key!r} is {attrs.get(key)!r}; expected {expected}"
            )
    missing_datasets = sorted(set(REQUIRED_TRACE_DATASETS) - set(datasets))
    if missing_datasets:
        errors.append(
            "traces.hdf5: required posterior datasets absent: "
            + ", ".join(missing_datasets)
        )
    wrong_leading_dimension = {
        name: shape
        for name, shape in datasets.items()
        if shape and shape[0] != spec.profile.n_samples
    }
    if wrong_leading_dimension:
        errors.append(
            "traces.hdf5: trace datasets have an unexpected leading dimension: "
            + json.dumps(wrong_leading_dimension, sort_keys=True)
        )
    if not datasets:
        errors.append("traces.hdf5: contains no trace datasets")
    return details, errors, warnings


def validate_model_output(spec: ChainSpec) -> dict[str, Any]:
    """Audit official no-qPCR artifacts without expecting ``subjset.pkl``.

    File presence alone is insufficient because MDSINE2 preallocates several
    artifacts before sampling completes.  The gate therefore also checks the
    terminal metadata marker, the configured seed/sampling contract, pickle
    envelopes, and HDF5 attributes/dataset dimensions.
    """

    files: dict[str, dict[str, Any]] = {}
    missing: list[str] = []
    for filename in REQUIRED_MODEL_FILES:
        path = spec.model_dir / filename
        if path.is_file() and path.stat().st_size > 0:
            files[filename] = {"bytes": path.stat().st_size, "path": str(path.resolve())}
        else:
            missing.append(filename)

    checks: dict[str, Any] = {}
    errors: list[str] = []
    warnings: list[str] = []
    if not missing:
        for filename in ("mcmc.pkl", "params.pkl"):
            check, found_errors = _audit_pickle_envelope(spec.model_dir / filename)
            checks[f"{filename}_envelope"] = check
            errors.extend(found_errors)
        params_check, params_errors = _audit_params_pickle(
            spec.model_dir / "params.pkl", spec
        )
        checks["params"] = params_check
        errors.extend(params_errors)
        metadata_check, metadata_errors = _audit_metadata(
            spec.model_dir / "metadata.txt", spec
        )
        checks["metadata"] = metadata_check
        errors.extend(metadata_errors)
        hdf5_check, hdf5_errors, hdf5_warnings = _audit_hdf5(
            spec.model_dir / "traces.hdf5", spec
        )
        checks["traces_hdf5"] = hdf5_check
        errors.extend(hdf5_errors)
        warnings.extend(hdf5_warnings)

    return {
        "schema_version": 2,
        "complete": not missing and not errors,
        "required_files": list(REQUIRED_MODEL_FILES),
        "intentionally_not_required": list(INTENTIONALLY_ABSENT_NO_QPCR_FILES),
        "files": files,
        "missing_or_empty": missing,
        "semantic_checks": checks,
        "semantic_errors": errors,
        "warnings": warnings,
    }


def tail_text(path: Path, lines: int = 50) -> str:
    if not path.is_file():
        return ""
    with path.open("rb") as handle:
        try:
            handle.seek(-min(path.stat().st_size, 128 * 1024), os.SEEK_END)
        except OSError:
            handle.seek(0)
        text = handle.read().decode("utf-8", errors="replace")
    return "\n".join(text.splitlines()[-lines:])


def update_status_table(inference_root: Path) -> Path:
    rows: list[dict[str, Any]] = []
    for path in sorted(inference_root.glob("*/*/seed_*/attempt_*/status.json")):
        try:
            status = load_json(path)
        except (OSError, json.JSONDecodeError):
            continue
        rows.append(
            {
                "profile": status.get("profile", ""),
                "group": status.get("group", ""),
                "seed": status.get("seed", ""),
                "attempt": status.get("attempt", ""),
                "state": status.get("state", ""),
                "returncode": status.get("returncode", ""),
                "child_returncode": status.get("child_returncode", ""),
                "started_at_utc": status.get("started_at_utc", ""),
                "finished_at_utc": status.get("finished_at_utc", ""),
                "pid": status.get("pid", ""),
                "session": status.get("session", ""),
                "model_dir": status.get("model_dir", ""),
                "log": status.get("log", ""),
            }
        )
    destination = inference_root / "status.tsv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    columns = (
        "profile",
        "group",
        "seed",
        "attempt",
        "state",
        "returncode",
        "child_returncode",
        "started_at_utc",
        "finished_at_utc",
        "pid",
        "session",
        "model_dir",
        "log",
    )
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, destination)
    return destination


def execute(spec: ChainSpec) -> int:
    validate_static_dependencies()
    validate_input_bundle(spec.group, spec.input_root)
    if not spec.attempt_root.is_dir() or not spec.status_path.is_file():
        raise RuntimeError(
            f"Attempt {spec.attempt} was not reserved at {spec.attempt_root}; use launch"
        )
    status = load_json(spec.status_path)
    if status.get("state") not in {"planned", "running"}:
        raise RuntimeError(
            f"Refusing to execute attempt in state {status.get('state')!r}: {spec.status_path}"
        )
    status.update(
        {
            "state": "running",
            "updated_at_utc": utc_now(),
            "started_at_utc": status.get("started_at_utc", utc_now()),
            "pid": os.getpid(),
            "command": affinity_wrapped_command(spec),
        }
    )
    atomic_json(spec.status_path, status)
    update_status_table(spec.inference_root)
    command = affinity_wrapped_command(spec)
    returncode = 1
    child_returncode: int | None = None
    child_pid = None
    try:
        spec.log_path.parent.mkdir(parents=True, exist_ok=True)
        with spec.log_path.open("ab", buffering=0) as log_handle:
            header = (
                f"\n[{utc_now()}] starting immutable attempt {spec.attempt}\n"
                f"cwd={PAPER_REPOSITORY}\n"
                f"command={shlex.join(command)}\n"
            ).encode("utf-8")
            log_handle.write(header)
            process = subprocess.Popen(
                command,
                cwd=str(PAPER_REPOSITORY),
                env=execution_environment(spec),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
            )
            child_pid = process.pid
            status["child_pid"] = child_pid
            status["updated_at_utc"] = utc_now()
            atomic_json(spec.status_path, status)
            child_returncode = process.wait()
            returncode = child_returncode
        output_audit = validate_model_output(spec)
        state = "succeeded" if returncode == 0 and output_audit["complete"] else "failed"
        if returncode == 0 and not output_audit["complete"]:
            returncode = 2
        status.update(
            {
                "state": state,
                "returncode": returncode,
                "child_returncode": child_returncode,
                "updated_at_utc": utc_now(),
                "finished_at_utc": utc_now(),
                "output_audit": output_audit,
            }
        )
        if state == "failed":
            status["log_tail"] = tail_text(spec.log_path)
    except BaseException as exc:
        status.update(
            {
                "state": "failed",
                "returncode": returncode,
                "child_returncode": child_returncode,
                "child_pid": child_pid,
                "updated_at_utc": utc_now(),
                "finished_at_utc": utc_now(),
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(),
                "log_tail": tail_text(spec.log_path),
            }
        )
        atomic_json(spec.status_path, status)
        update_status_table(spec.inference_root)
        raise
    atomic_json(spec.status_path, status)
    update_status_table(spec.inference_root)
    return returncode


def tmux_session_exists(name: str) -> bool:
    result = subprocess.run(
        ("tmux", "has-session", "-t", name),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def execute_cli_command(spec: ChainSpec) -> list[str]:
    return [
        sys.executable,
        str(SCRIPT_PATH),
        "execute",
        "--profile",
        spec.profile.name,
        "--group",
        spec.group,
        "--seed",
        str(spec.seed),
        "--attempt",
        str(spec.attempt),
        "--input-root",
        str(spec.input_root.resolve()),
        "--inference-root",
        str(spec.inference_root.resolve()),
    ]


def launch_tmux(spec: ChainSpec) -> None:
    if tmux_session_exists(spec.session_name):
        raise RuntimeError(f"tmux session already exists: {spec.session_name}")
    command = (
        "tmux",
        "new-session",
        "-d",
        "-s",
        spec.session_name,
        "taskset",
        "--cpu-list",
        spec.cpu_affinity,
        *execute_cli_command(spec),
    )
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        status = load_json(spec.status_path)
        status.update(
            {
                "state": "failed",
                "updated_at_utc": utc_now(),
                "finished_at_utc": utc_now(),
                "failure_reason": "tmux launch failed",
                "tmux_error": (result.stderr or result.stdout).strip(),
            }
        )
        atomic_json(spec.status_path, status)
        update_status_table(spec.inference_root)
        raise RuntimeError(status["tmux_error"])


def profile_specs(profile: Profile, groups: Iterable[str], input_root: Path, inference_root: Path) -> list[ChainSpec]:
    return [
        ChainSpec(group, profile, seed, input_root, inference_root, attempt=1)
        for group in groups
        for seed in profile.seeds
    ]


def launch_profile(
    profile: Profile,
    groups: Sequence[str],
    input_root: Path,
    inference_root: Path,
    mode: str,
    max_attempts: int,
) -> list[dict[str, Any]]:
    validate_static_dependencies()
    for group in groups:
        validate_input_bundle(group, input_root)
    write_provenance(inference_root)
    outcomes: list[dict[str, Any]] = []
    for template in profile_specs(profile, groups, input_root, inference_root):
        spec, newly_reserved = reserve_attempt(
            template.group,
            profile,
            template.seed,
            input_root,
            inference_root,
            max_attempts=max_attempts,
        )
        prior_state = load_json(spec.status_path).get("state") if spec.status_path.is_file() else None
        if not newly_reserved:
            outcomes.append(
                {
                    "group": spec.group,
                    "seed": spec.seed,
                    "attempt": spec.attempt,
                    "action": "unchanged",
                    "state": prior_state,
                }
            )
            continue
        if mode == "tmux":
            launch_tmux(spec)
            action = "launched"
            returncode = None
        else:
            action = "completed_foreground"
            returncode = execute(spec)
        outcomes.append(
            {
                "group": spec.group,
                "seed": spec.seed,
                "attempt": spec.attempt,
                "action": action,
                "returncode": returncode,
                "session": spec.session_name,
                "status": str(spec.status_path),
            }
        )
    update_status_table(inference_root)
    return outcomes


def logical_states(profile: Profile, groups: Sequence[str], input_root: Path, inference_root: Path) -> list[dict[str, Any]]:
    states = []
    for template in profile_specs(profile, groups, input_root, inference_root):
        spec, status = latest_status_for(
            template.group, profile, template.seed, input_root, inference_root
        )
        if spec is None:
            states.append(
                {
                    "group": template.group,
                    "seed": template.seed,
                    "state": "not_started",
                }
            )
        else:
            states.append(
                {
                    "group": spec.group,
                    "seed": spec.seed,
                    "attempt": spec.attempt,
                    "state": (status or {}).get("state", "unknown"),
                    "status": str(spec.status_path),
                }
            )
    return states


def process_is_alive(pid: Any) -> bool:
    try:
        numeric_pid = int(pid)
    except (TypeError, ValueError):
        return False
    if numeric_pid <= 0:
        return False
    try:
        os.kill(numeric_pid, 0)
    except (ProcessLookupError, ValueError):
        return False
    except PermissionError:
        return True
    return True


def mark_latest_attempts_stopped(
    profile: Profile,
    groups: Sequence[str],
    input_root: Path,
    inference_root: Path,
    reason: str,
) -> list[dict[str, Any]]:
    """Close dead planned/running attempts without deleting or overwriting them."""
    results: list[dict[str, Any]] = []
    for template in profile_specs(profile, groups, input_root, inference_root):
        spec, status = latest_status_for(
            template.group, profile, template.seed, input_root, inference_root
        )
        if spec is None or status is None:
            results.append(
                {"group": template.group, "seed": template.seed, "action": "not_found"}
            )
            continue
        state = status.get("state")
        if state not in {"planned", "running"}:
            results.append(
                {
                    "group": spec.group,
                    "seed": spec.seed,
                    "attempt": spec.attempt,
                    "action": "unchanged",
                    "state": state,
                }
            )
            continue
        live_pids = [
            int(pid)
            for pid in (status.get("pid"), status.get("child_pid"))
            if process_is_alive(pid)
        ]
        if live_pids:
            raise RuntimeError(
                f"Refusing to mark {spec.status_path} stopped; live PIDs remain: {live_pids}"
            )
        status.update(
            {
                "state": "stopped",
                "updated_at_utc": utc_now(),
                "finished_at_utc": utc_now(),
                "stop_reason": reason,
                "partial_outputs_preserved": True,
                "resume_used": False,
            }
        )
        atomic_json(spec.status_path, status)
        results.append(
            {
                "group": spec.group,
                "seed": spec.seed,
                "attempt": spec.attempt,
                "action": "marked_stopped",
                "status": str(spec.status_path),
            }
        )
    update_status_table(inference_root)
    return results


def _successful_child_evidence(status: dict[str, Any]) -> tuple[bool, str]:
    """Return whether a status proves that the original child exited zero.

    New statuses retain ``child_returncode`` explicitly.  Before that field was
    added, this launcher used orchestrator return code 2 *only* when the child
    returned zero but the five-file gate failed.  We accept that legacy
    inference only for the exact historical false-negative signature: the four
    real official artifacts were recorded and ``subjset.pkl`` alone was absent.
    """

    if "child_returncode" in status:
        return status.get("child_returncode") == 0, "explicit child_returncode"

    old_audit = status.get("output_audit")
    legacy_gate_only_failure = (
        status.get("state") == "failed"
        and status.get("returncode") == 2
        and isinstance(old_audit, dict)
        and old_audit.get("complete") is False
        and set(old_audit.get("missing_or_empty", []))
        == set(INTENTIONALLY_ABSENT_NO_QPCR_FILES)
        and set(old_audit.get("files", {})) == set(REQUIRED_MODEL_FILES)
        and not status.get("error")
        and not status.get("error_type")
        and not status.get("failure_reason")
    )
    if legacy_gate_only_failure:
        return True, (
            "legacy launcher returncode=2 with only subjset.pkl absent; under "
            "the historical execute() branch this proves child returncode=0"
        )
    return False, "no unambiguous zero child-returncode evidence"


def reconcile_latest_outputs(
    profile: Profile,
    groups: Sequence[str],
    input_root: Path,
    inference_root: Path,
    apply: bool = False,
) -> list[dict[str, Any]]:
    """Re-audit terminal attempts and safely repair legacy gate-only failures.

    A status changes to ``succeeded`` only if (1) the original child exit zero
    is explicit or follows the exact old ``subjset.pkl`` false-failure branch,
    and (2) the current four-file semantic gate passes.  Running, stopped, or
    ambiguously failed attempts are never promoted.
    """

    outcomes: list[dict[str, Any]] = []
    changed = False
    for template in profile_specs(profile, groups, input_root, inference_root):
        spec, status = latest_status_for(
            template.group, profile, template.seed, input_root, inference_root
        )
        if spec is None or status is None:
            outcomes.append(
                {"group": template.group, "seed": template.seed, "action": "not_found"}
            )
            continue
        if status.get("state") == "succeeded":
            outcomes.append(
                {
                    "group": spec.group,
                    "seed": spec.seed,
                    "attempt": spec.attempt,
                    "action": "unchanged",
                    "state": "succeeded",
                }
            )
            continue
        if status.get("state") != "failed":
            outcomes.append(
                {
                    "group": spec.group,
                    "seed": spec.seed,
                    "attempt": spec.attempt,
                    "action": "ineligible_state",
                    "state": status.get("state", "unknown"),
                }
            )
            continue

        child_ok, child_evidence = _successful_child_evidence(status)
        if not child_ok:
            outcomes.append(
                {
                    "group": spec.group,
                    "seed": spec.seed,
                    "attempt": spec.attempt,
                    "action": "refused",
                    "reason": child_evidence,
                }
            )
            continue
        output_audit = validate_model_output(spec)
        if not output_audit["complete"]:
            outcomes.append(
                {
                    "group": spec.group,
                    "seed": spec.seed,
                    "attempt": spec.attempt,
                    "action": "refused",
                    "reason": "current four-file semantic output gate failed",
                    "child_returncode_evidence": child_evidence,
                    "output_audit": output_audit,
                }
            )
            continue

        outcome = {
            "group": spec.group,
            "seed": spec.seed,
            "attempt": spec.attempt,
            "action": "reconciled" if apply else "would_reconcile",
            "child_returncode_evidence": child_evidence,
            "status": str(spec.status_path),
        }
        if apply:
            prior_state = status.get("state")
            prior_returncode = status.get("returncode")
            reconciliation = {
                "reconciled_at_utc": utc_now(),
                "from_state": prior_state,
                "from_returncode": prior_returncode,
                "prior_status_sha256": sha256_file(spec.status_path),
                "child_returncode_evidence": child_evidence,
                "gate_schema_version": output_audit["schema_version"],
                "reason": (
                    "official initialize_graph_no_qpcr does not emit subjset.pkl; "
                    "the corrected four-file semantic gate passed"
                ),
            }
            status.update(
                {
                    "state": "succeeded",
                    "returncode": 0,
                    "child_returncode": 0,
                    "updated_at_utc": reconciliation["reconciled_at_utc"],
                    "output_audit": output_audit,
                    "status_reconciliation": reconciliation,
                }
            )
            atomic_json(spec.status_path, status)
            changed = True
        outcomes.append(outcome)
    if changed:
        update_status_table(inference_root)
    return outcomes


def wait_for_profile(
    profile: Profile,
    groups: Sequence[str],
    input_root: Path,
    inference_root: Path,
    poll_seconds: int,
) -> tuple[bool, list[dict[str, Any]]]:
    while True:
        states = logical_states(profile, groups, input_root, inference_root)
        if all(item["state"] == "succeeded" for item in states):
            return True, states
        if any(item["state"] in {"failed", "unknown", "not_started"} for item in states):
            # A legacy launcher that was running on entry can atomically write
            # its old gate-only ``failed`` status at any later poll.  Reconcile
            # *after* observing that terminal state, then refresh before
            # deciding whether the profile really failed.
            reconcile_latest_outputs(
                profile, groups, input_root, inference_root, apply=True
            )
            states = logical_states(profile, groups, input_root, inference_root)
            if all(item["state"] == "succeeded" for item in states):
                return True, states
            if any(
                item["state"] in {"failed", "unknown", "not_started"}
                for item in states
            ):
                return False, states
        time.sleep(poll_seconds)


def supervisor_worker(args: argparse.Namespace) -> int:
    inference_root = args.inference_root.resolve()
    supervisor_status = inference_root / "supervisor_status.json"
    inference_root.mkdir(parents=True, exist_ok=True)
    if supervisor_status.is_file():
        history_dir = inference_root / "supervisor_history"
        history_dir.mkdir(parents=True, exist_ok=True)
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        history_path = history_dir / f"supervisor_{stamp}.json"
        history_path.write_bytes(supervisor_status.read_bytes())
    lock_path = inference_root / ".supervisor.lock"
    with lock_path.open("a+") as lock_handle:
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another Fig. 3 MDSINE2 supervisor is active") from exc
        payload = {
            "schema_version": 1,
            "state": "launching_pilots",
            "started_at_utc": utc_now(),
            "pid": os.getpid(),
            "argv": sys.argv,
            "thread_contract": {
                "threads_per_group": THREADS_PER_GROUP,
                "cpu_affinity": CPU_AFFINITY,
                "formal_chain_cpu_affinity": FORMAL_CHAIN_CPU_AFFINITY,
                "disjoint": True,
                "mdsine2_multiprocessing": 0,
                "note": NO_QPCR_MULTIPROCESSING_NOTE,
            },
            "pilot_acceptance": (
                "both no-qPCR child processes exit 0 and each model directory passes the "
                "four-file metadata/params/HDF5 semantic gate: "
                + ", ".join(REQUIRED_MODEL_FILES)
            ),
        }
        atomic_json(supervisor_status, payload)
        try:
            # On supervisor restart, repair only the exact historical
            # subjset.pkl gate false-negative before reserving attempts.  This
            # prevents a completed attempt_002 from becoming attempt_003.
            payload["pilot_output_reconciliation"] = reconcile_latest_outputs(
                PROFILES["pilot"],
                GROUPS,
                args.input_root.resolve(),
                inference_root,
                apply=True,
            )
            payload["pilot_launch"] = launch_profile(
                PROFILES["pilot"],
                GROUPS,
                args.input_root.resolve(),
                inference_root,
                mode="tmux",
                max_attempts=args.max_attempts,
            )
            payload["state"] = "waiting_for_pilots"
            payload["updated_at_utc"] = utc_now()
            atomic_json(supervisor_status, payload)
            pilots_ok, pilot_states = wait_for_profile(
                PROFILES["pilot"],
                GROUPS,
                args.input_root.resolve(),
                inference_root,
                args.poll_seconds,
            )
            payload["pilot_states"] = pilot_states
            if not pilots_ok:
                payload.update(
                    {
                        "state": "failed",
                        "failure_reason": "at least one pilot did not pass output validation",
                        "finished_at_utc": utc_now(),
                    }
                )
                atomic_json(supervisor_status, payload)
                return 2
            payload["state"] = "launching_formal"
            payload["updated_at_utc"] = utc_now()
            atomic_json(supervisor_status, payload)
            payload["formal_launch"] = launch_profile(
                PROFILES["formal"],
                GROUPS,
                args.input_root.resolve(),
                inference_root,
                mode="tmux",
                max_attempts=args.max_attempts,
            )
            payload["state"] = "waiting_for_formal"
            payload["updated_at_utc"] = utc_now()
            atomic_json(supervisor_status, payload)
            formal_ok, formal_states = wait_for_profile(
                PROFILES["formal"],
                GROUPS,
                args.input_root.resolve(),
                inference_root,
                args.poll_seconds,
            )
            payload["formal_states"] = formal_states
            payload["state"] = "succeeded" if formal_ok else "failed"
            if not formal_ok:
                payload["failure_reason"] = "at least one formal chain failed output validation"
            payload["finished_at_utc"] = utc_now()
            atomic_json(supervisor_status, payload)
            return 0 if formal_ok else 3
        except BaseException as exc:
            payload.update(
                {
                    "state": "failed",
                    "finished_at_utc": utc_now(),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "traceback": traceback.format_exc(),
                }
            )
            atomic_json(supervisor_status, payload)
            raise


def detach_supervisor(args: argparse.Namespace) -> dict[str, Any]:
    validate_static_dependencies()
    session = "fig3m2-supervisor"
    if tmux_session_exists(session):
        raise RuntimeError(f"Supervisor tmux session already exists: {session}")
    command = [
        sys.executable,
        str(SCRIPT_PATH),
        "supervisor-worker",
        "--input-root",
        str(args.input_root.resolve()),
        "--inference-root",
        str(args.inference_root.resolve()),
        "--poll-seconds",
        str(args.poll_seconds),
        "--max-attempts",
        str(args.max_attempts),
    ]
    result = subprocess.run(
        ("tmux", "new-session", "-d", "-s", session, *command),
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError((result.stderr or result.stdout).strip())
    return {
        "session": session,
        "state": "launched",
        "supervisor_status": str(args.inference_root.resolve() / "supervisor_status.json"),
        "command": command,
    }


def dry_run_plan(profile: Profile, groups: Sequence[str], input_root: Path, inference_root: Path) -> list[dict[str, Any]]:
    plan = []
    for template in profile_specs(profile, groups, input_root, inference_root):
        spec = ChainSpec(
            template.group,
            profile,
            template.seed,
            input_root,
            inference_root,
            attempt=1,
        )
        plan.append(
            {
                "group": spec.group,
                "profile": profile.name,
                "seed": spec.seed,
                "input": str(spec.input_path.resolve()),
                "attempt_root_template": str(spec.attempt_root.resolve()),
                "model_dir_template": str(spec.model_dir.resolve()),
                "session_template": spec.session_name,
                "command": affinity_wrapped_command(spec),
                "cpu_affinity": spec.cpu_affinity,
                "threads_per_group": THREADS_PER_GROUP,
                "threads_per_chain": spec.threads_per_chain,
                "resume_used": False,
            }
        )
    return plan


def add_paths(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--inference-root", type=Path, default=DEFAULT_INFERENCE_ROOT)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)

    plan_parser = subparsers.add_parser("plan", help="Print deterministic commands without writing")
    plan_parser.add_argument("--profile", choices=tuple(PROFILES), required=True)
    plan_parser.add_argument("--groups", nargs="+", choices=GROUPS, default=list(GROUPS))
    plan_parser.add_argument("--validate-files", action="store_true")
    add_paths(plan_parser)

    validate_parser = subparsers.add_parser(
        "validate-inputs", help="Deep-audit both Study pickles without fitting"
    )
    validate_parser.add_argument("--groups", nargs="+", choices=GROUPS, default=list(GROUPS))
    add_paths(validate_parser)

    launch_parser = subparsers.add_parser("launch", help="Launch one fixed profile")
    launch_parser.add_argument("--profile", choices=tuple(PROFILES), required=True)
    launch_parser.add_argument("--groups", nargs="+", choices=GROUPS, default=list(GROUPS))
    launch_parser.add_argument("--mode", choices=("tmux", "foreground"), default="tmux")
    launch_parser.add_argument("--max-attempts", type=int, default=1)
    add_paths(launch_parser)

    execute_parser = subparsers.add_parser("execute", help=argparse.SUPPRESS)
    execute_parser.add_argument("--profile", choices=tuple(PROFILES), required=True)
    execute_parser.add_argument("--group", choices=GROUPS, required=True)
    execute_parser.add_argument("--seed", type=int, required=True)
    execute_parser.add_argument("--attempt", type=int, required=True)
    add_paths(execute_parser)

    status_parser = subparsers.add_parser("status", help="Refresh and print chain status")
    add_paths(status_parser)

    stopped_parser = subparsers.add_parser(
        "mark-stopped",
        help="Mark dead planned/running attempts stopped while retaining all files",
    )
    stopped_parser.add_argument("--profile", choices=tuple(PROFILES), required=True)
    stopped_parser.add_argument("--groups", nargs="+", choices=GROUPS, default=list(GROUPS))
    stopped_parser.add_argument("--reason", required=True)
    add_paths(stopped_parser)

    reconcile_parser = subparsers.add_parser(
        "reconcile-outputs",
        help=(
            "Re-audit terminal attempts; dry-run unless --apply, and promote only "
            "proven zero-exit legacy gate-only failures"
        ),
    )
    reconcile_parser.add_argument("--profile", choices=tuple(PROFILES), required=True)
    reconcile_parser.add_argument(
        "--groups", nargs="+", choices=GROUPS, default=list(GROUPS)
    )
    reconcile_parser.add_argument("--apply", action="store_true")
    add_paths(reconcile_parser)

    supervise_parser = subparsers.add_parser(
        "supervise", help="Detach pilot-to-formal supervisor in tmux"
    )
    supervise_parser.add_argument("--poll-seconds", type=int, default=60)
    supervise_parser.add_argument("--max-attempts", type=int, default=1)
    add_paths(supervise_parser)

    worker_parser = subparsers.add_parser("supervisor-worker", help=argparse.SUPPRESS)
    worker_parser.add_argument("--poll-seconds", type=int, default=60)
    worker_parser.add_argument("--max-attempts", type=int, default=1)
    add_paths(worker_parser)
    return parser.parse_args()


def validate_profile_seed(profile: Profile, seed: int) -> None:
    if seed not in profile.seeds:
        raise ValueError(
            f"Seed {seed} is not allowed for profile {profile.name}; expected {profile.seeds}"
        )


def deep_validate_inputs(args: argparse.Namespace) -> list[dict[str, Any]]:
    validate_static_dependencies()
    results = []
    validation_root = args.inference_root.resolve() / "input_validation"
    validation_root.mkdir(parents=True, exist_ok=True)
    for group in args.groups:
        validate_input_bundle(group, args.input_root.resolve())
        profile = PROFILES["pilot"]
        spec = ChainSpec(
            group,
            profile,
            profile.seeds[0],
            args.input_root.resolve(),
            validation_root / "command_scaffold",
            attempt=1,
        )
        audit_path = validation_root / f"{group}_input_audit.json"
        command = wrapper_command(spec, validate_only=True)
        audit_index = command.index("--audit-json") + 1
        command[audit_index] = str(audit_path)
        result = subprocess.run(
            command,
            cwd=str(PAPER_REPOSITORY),
            env=execution_environment(spec),
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"Deep input validation failed for {group}:\n{result.stdout}\n{result.stderr}"
            )
        results.append(load_json(audit_path))
    return results


def main() -> int:
    args = parse_args()
    if hasattr(args, "poll_seconds") and args.poll_seconds < 1:
        raise ValueError("--poll-seconds must be >= 1")
    if hasattr(args, "max_attempts") and args.max_attempts < 1:
        raise ValueError("--max-attempts must be >= 1")
    if args.action == "plan":
        if args.validate_files:
            validate_static_dependencies()
            for group in args.groups:
                validate_input_bundle(group, args.input_root.resolve())
        print(
            json.dumps(
                dry_run_plan(
                    PROFILES[args.profile],
                    args.groups,
                    args.input_root.resolve(),
                    args.inference_root.resolve(),
                ),
                indent=2,
                ensure_ascii=False,
            )
        )
        return 0
    if args.action == "validate-inputs":
        print(json.dumps(deep_validate_inputs(args), indent=2, ensure_ascii=False))
        return 0
    if args.action == "launch":
        outcomes = launch_profile(
            PROFILES[args.profile],
            args.groups,
            args.input_root.resolve(),
            args.inference_root.resolve(),
            mode=args.mode,
            max_attempts=args.max_attempts,
        )
        print(json.dumps(outcomes, indent=2, ensure_ascii=False))
        return 0
    if args.action == "execute":
        profile = PROFILES[args.profile]
        validate_profile_seed(profile, args.seed)
        return execute(
            ChainSpec(
                args.group,
                profile,
                args.seed,
                args.input_root.resolve(),
                args.inference_root.resolve(),
                args.attempt,
            )
        )
    if args.action == "status":
        table = update_status_table(args.inference_root.resolve())
        payload = {
            profile: logical_states(
                PROFILES[profile], GROUPS, args.input_root.resolve(), args.inference_root.resolve()
            )
            for profile in PROFILES
        }
        payload["status_table"] = str(table)
        supervisor = args.inference_root.resolve() / "supervisor_status.json"
        if supervisor.is_file():
            payload["supervisor"] = load_json(supervisor)
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 0
    if args.action == "mark-stopped":
        print(
            json.dumps(
                mark_latest_attempts_stopped(
                    PROFILES[args.profile],
                    args.groups,
                    args.input_root.resolve(),
                    args.inference_root.resolve(),
                    args.reason,
                ),
                indent=2,
                ensure_ascii=False,
            )
        )
        return 0
    if args.action == "reconcile-outputs":
        print(
            json.dumps(
                reconcile_latest_outputs(
                    PROFILES[args.profile],
                    args.groups,
                    args.input_root.resolve(),
                    args.inference_root.resolve(),
                    apply=args.apply,
                ),
                indent=2,
                ensure_ascii=False,
            )
        )
        return 0
    if args.action == "supervise":
        print(json.dumps(detach_supervisor(args), indent=2, ensure_ascii=False))
        return 0
    if args.action == "supervisor-worker":
        return supervisor_worker(args)
    raise RuntimeError(f"Unhandled action: {args.action}")


if __name__ == "__main__":
    raise SystemExit(main())
