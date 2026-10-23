#!/usr/bin/env python3
"""Prepare, launch, supervise and finalize R1.1 fast MDSINE2 reruns.

Each recovery definition receives one independent no-qPCR chain per group with
the repository's previously validated pilot settings (140 total iterations,
40 burn-in, 100 retained draws, checkpoint every 20 iterations).  These runs
are sensitivity/developmental evidence, not formal multichain convergence runs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time
import traceback

import pandas as pd


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
FIG3_CODE = PROJECT_ROOT / "reproduction" / "fig3_mdsine2"
sys.path.insert(0, str(FIG3_CODE))

# The pinned MDSINE2 source tree is read-only in this workspace.  Keep all
# runtime caches in writable temporary locations before importing the package.
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/mdsine2_numba")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mdsine2_mpl")
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

import prepare_inputs as prepare  # noqa: E402
import run_inference as runner  # noqa: E402


OUTPUT_ROOT = HERE / "full_original_layout"
GROUP_MANIFEST = HERE / "R1_1_group_manifest_D0_D6.csv"
DEFINITION_SLUGS = {
    "D1": "D1_endpoint_abs_shannon",
    "D2": "D2_endpoint_bray",
    "D3": "D3_endpoint_jaccard",
    "D4": "D4_auc_abs_shannon",
    "D5": "D5_auc_bray",
    "D6": "D6_auc_jaccard",
}
CPU_PAIRS = {
    "D1": {"LDR": "16-31", "HDR": "112-127"},
    "D2": {"LDR": "32-47", "HDR": "128-143"},
    "D3": {"LDR": "48-63", "HDR": "144-159"},
    "D4": {"LDR": "64-79", "HDR": "160-175"},
    "D5": {"LDR": "80-95", "HDR": "176-191"},
    # D6 is queued and reuses the first completed definition's pair.
    "D6": {"LDR": "16-31", "HDR": "112-127"},
}
MDSINE_PYTHON = Path("/home/qqt/.conda/envs/mdsine2/bin/python")
PROFILE = runner.Profile(
    "r11_fast",
    n_samples=140,
    burnin=40,
    checkpoint=20,
    seeds=(20261009,),
)
SEED = PROFILE.seeds[0]


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class RapidSpec(runner.ChainSpec):
    definition_id: str
    requested_cpu_affinity: str

    @property
    def study_name(self) -> str:
        return (
            f"r11-{self.definition_id.lower()}-{self.group.lower()}-"
            f"fast-seed{self.seed}"
        )

    @property
    def session_name(self) -> str:
        return f"r11m2-{self.definition_id.lower()}-{self.group.lower()}-a{self.attempt:03d}"

    @property
    def cpu_affinity(self) -> str:
        return self.requested_cpu_affinity


def definition_root(output_root: Path, definition_id: str) -> Path:
    return output_root / DEFINITION_SLUGS[definition_id]


def assignment_path(output_root: Path, definition_id: str) -> Path:
    return definition_root(output_root, definition_id) / "mdsine2" / "assignments.csv"


def input_root(output_root: Path, definition_id: str) -> Path:
    return definition_root(output_root, definition_id) / "inputs"


def inference_root(output_root: Path, definition_id: str) -> Path:
    return definition_root(output_root, definition_id) / "mdsine2" / "inference"


def make_spec(
    output_root: Path,
    definition_id: str,
    group: str,
    attempt: int = 1,
    cpu_affinity: str | None = None,
) -> RapidSpec:
    return RapidSpec(
        group=group,
        profile=PROFILE,
        seed=SEED,
        input_root=input_root(output_root, definition_id),
        inference_root=inference_root(output_root, definition_id),
        attempt=attempt,
        definition_id=definition_id,
        requested_cpu_affinity=(
            cpu_affinity if cpu_affinity is not None else CPU_PAIRS[definition_id][group]
        ),
    )


def prepare_definition(
    output_root: Path,
    group_manifest: Path,
    definition_id: str,
) -> dict[str, object]:
    manifest = pd.read_csv(
        group_manifest,
        dtype={"subject": str, "definition_id": str},
    )
    selected = manifest.loc[manifest["definition_id"].eq(definition_id)].copy()
    if len(selected) != 60 or selected["subject"].nunique() != 60:
        raise ValueError(f"{definition_id}: expected 60 unique subjects")
    counts = selected["recovery_group"].value_counts().to_dict()
    if counts != {"Middle": 30, "HDR": 15, "LDR": 15}:
        raise ValueError(f"{definition_id}: invalid group counts {counts}")
    assignments = selected[
        ["subject", "metric_rank_ascending", "recovery_group"]
    ].rename(columns={"metric_rank_ascending": "rank"})
    destination = assignment_path(output_root, definition_id)
    destination.parent.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(destination, index=False)

    prepared = prepare.prepare_inputs(
        mpa_path=PROJECT_ROOT / "raw" / "combine_brk.mpa",
        metadata_path=PROJECT_ROOT / "raw" / "meta.csv",
        assignments_path=destination,
        common_table_path=prepare.DEFAULT_COMMON_TABLE,
        taxonomy_source_path=prepare.DEFAULT_TAXONOMY_SOURCE,
        output_dir=input_root(output_root, definition_id),
        depth=prepare.RAREFACTION_DEPTH,
        seed=prepare.RAREFACTION_SEED,
        make_study_pickles=True,
        sample_stable_rarefaction=True,
    )
    return {
        "definition_id": definition_id,
        "assignments": str(destination.resolve()),
        "assignments_sha256": sha256(destination),
        "inputs": str(input_root(output_root, definition_id).resolve()),
        "subjects": prepared["selection"]["observed_subjects_per_group"],
        "sample_stable_rarefaction": True,
    }


def reserve_spec(spec: RapidSpec) -> bool:
    if spec.status_path.is_file():
        state = runner.load_json(spec.status_path).get("state")
        if state in {"planned", "running", "succeeded"}:
            return False
        raise RuntimeError(
            f"{spec.definition_id}/{spec.group}: prior attempt is {state}; "
            "files were preserved and automatic destructive retry is disabled"
        )
    spec.logical_root.mkdir(parents=True, exist_ok=True)
    spec.attempt_root.mkdir()
    status = runner.base_status(spec, "planned")
    status.update(
        {
            "created_at_utc": utc_now(),
            "definition_id": spec.definition_id,
            "evidence_tier": "pilot_single_or_short_chain",
            "interpretation": "fast sensitivity/developmental run; not formal multichain convergence evidence",
        }
    )
    runner.atomic_json(spec.status_path, status)
    runner.atomic_json(
        spec.command_path,
        {
            "created_at_utc": utc_now(),
            "worker_argv": worker_command(spec),
            "official_runner_argv": runner.affinity_wrapped_command(spec),
            "resume_deliberately_omitted": True,
            "definition_id": spec.definition_id,
            "cpu_affinity": spec.cpu_affinity,
            "threads_per_chain": spec.threads_per_chain,
        },
    )
    runner.atomic_text(spec.logical_root / "latest_attempt.txt", f"{spec.attempt}\n")
    runner.update_status_table(spec.inference_root)
    return True


def worker_command(spec: RapidSpec) -> list[str]:
    return [
        str(MDSINE_PYTHON),
        str(Path(__file__).resolve()),
        "worker",
        "--definition",
        spec.definition_id,
        "--group",
        spec.group,
        "--output-root",
        str(spec.input_root.parent.parent.resolve()),
        "--cpu-affinity",
        spec.cpu_affinity,
    ]


def launch_worker(spec: RapidSpec) -> dict[str, object]:
    runner.validate_input_bundle(spec.group, spec.input_root)
    newly_reserved = reserve_spec(spec)
    state = runner.load_json(spec.status_path).get("state")
    if not newly_reserved:
        return {"session": spec.session_name, "state": state, "action": "unchanged"}
    if runner.tmux_session_exists(spec.session_name):
        raise RuntimeError(f"tmux session already exists: {spec.session_name}")
    command = [
        "tmux",
        "new-session",
        "-d",
        "-s",
        spec.session_name,
        "taskset",
        "--cpu-list",
        spec.cpu_affinity,
        *worker_command(spec),
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode:
        status = runner.load_json(spec.status_path)
        status.update(
            {
                "state": "failed",
                "finished_at_utc": utc_now(),
                "failure_reason": "tmux launch failed",
                "tmux_error": (completed.stderr or completed.stdout).strip(),
            }
        )
        runner.atomic_json(spec.status_path, status)
        raise RuntimeError(status["tmux_error"])
    return {"session": spec.session_name, "state": "launched", "action": "launched"}


def execute_worker(spec: RapidSpec) -> int:
    """Execute one immutable chain using the audited runner's gates."""
    runner.validate_input_bundle(spec.group, spec.input_root)
    if not spec.status_path.is_file():
        raise FileNotFoundError(f"Reserved status is missing: {spec.status_path}")
    status = runner.load_json(spec.status_path)
    if status.get("state") not in {"planned", "running"}:
        raise RuntimeError(f"Refusing worker in state {status.get('state')!r}")
    status.update(
        {
            "state": "running",
            "updated_at_utc": utc_now(),
            "started_at_utc": status.get("started_at_utc", utc_now()),
            "pid": os.getpid(),
            "command": runner.affinity_wrapped_command(spec),
        }
    )
    runner.atomic_json(spec.status_path, status)
    runner.update_status_table(spec.inference_root)
    returncode = 1
    child_returncode = None
    child_pid = None
    try:
        spec.log_path.parent.mkdir(parents=True, exist_ok=True)
        with spec.log_path.open("ab", buffering=0) as log_handle:
            log_handle.write(
                (
                    f"\n[{utc_now()}] starting R1.1 fast immutable attempt\n"
                    f"definition={spec.definition_id}\n"
                    f"cwd={runner.PAPER_REPOSITORY}\n"
                    f"command={shlex.join(runner.affinity_wrapped_command(spec))}\n"
                ).encode("utf-8")
            )
            process = subprocess.Popen(
                runner.affinity_wrapped_command(spec),
                cwd=str(runner.PAPER_REPOSITORY),
                env=runner.execution_environment(spec),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
            )
            child_pid = process.pid
            status["child_pid"] = child_pid
            status["updated_at_utc"] = utc_now()
            runner.atomic_json(spec.status_path, status)
            child_returncode = process.wait()
            returncode = child_returncode
        output_audit = runner.validate_model_output(spec)
        state = (
            "succeeded"
            if child_returncode == 0 and output_audit["complete"]
            else "failed"
        )
        if child_returncode == 0 and not output_audit["complete"]:
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
            status["log_tail"] = runner.tail_text(spec.log_path)
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
                "log_tail": runner.tail_text(spec.log_path),
            }
        )
        runner.atomic_json(spec.status_path, status)
        runner.update_status_table(spec.inference_root)
        raise
    runner.atomic_json(spec.status_path, status)
    runner.update_status_table(spec.inference_root)
    return int(returncode)


def chain_state(output_root: Path, definition_id: str, group: str) -> str:
    spec = make_spec(output_root, definition_id, group)
    if not spec.status_path.is_file():
        return "not_started"
    return str(runner.load_json(spec.status_path).get("state", "unknown"))


def definition_state(output_root: Path, definition_id: str) -> dict[str, str]:
    return {
        group: chain_state(output_root, definition_id, group)
        for group in ("LDR", "HDR")
    }


def run_checked(command: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(f"[{utc_now()}] {shlex.join(command)}\n")
        handle.flush()
        completed = subprocess.run(
            command,
            check=False,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if completed.returncode:
        raise RuntimeError(
            f"Post-processing command returned {completed.returncode}: {shlex.join(command)}"
        )


def _definition_short_labels() -> dict[str, str]:
    """Return one stable human-readable label for every alternative definition."""

    labels = {definition_id: DEFINITION_SLUGS[definition_id] for definition_id in DEFINITION_SLUGS}
    if GROUP_MANIFEST.is_file():
        manifest = pd.read_csv(
            GROUP_MANIFEST,
            usecols=["definition_id", "definition_short"],
            dtype=str,
        ).drop_duplicates("definition_id")
        labels.update(
            manifest.set_index("definition_id")["definition_short"].dropna().to_dict()
        )
    return labels


def update_reviewer_index_fig3(
    output_root: Path,
    definition_id: str,
    pdf_path: Path,
    tif_path: Path,
) -> Path:
    """Atomically register one validated developmental Fig. 3 export.

    The Fig. 2/4 renderer and the MDSINE2 supervisor can finish in either
    order.  If the shared index does not yet exist, initialise all six rows
    with their deterministic output paths so a later definition-level update
    never discards another completed Fig. 3.
    """

    missing = [
        str(path)
        for path in (pdf_path, tif_path)
        if not path.is_file() or path.stat().st_size <= 0
    ]
    if missing:
        raise FileNotFoundError(
            f"Cannot register {definition_id} Fig. 3; missing/empty exports: {missing}"
        )
    index_path = output_root / "reviewer_ready_index.csv"
    labels = _definition_short_labels()
    if index_path.is_file() and index_path.stat().st_size:
        index = pd.read_csv(index_path, dtype=str, keep_default_na=False)
    else:
        rows: list[dict[str, str]] = []
        for current_id, slug in DEFINITION_SLUGS.items():
            root = definition_root(output_root, current_id)
            rows.append(
                {
                    "definition_id": current_id,
                    "definition_short": labels[current_id],
                    "directory": str(root.resolve()),
                    "fig2_pdf": str((root / "figures" / f"R1_1_{current_id}_Fig2.pdf").resolve()),
                    "fig2_tif": str((root / "figures" / f"R1_1_{current_id}_Fig2.tif").resolve()),
                    "fig3_pdf": "",
                    "fig3_tif": "",
                    "fig3_status": "pending rapid MDSINE2",
                    "fig4_pdf": str((root / "figures" / f"R1_1_{current_id}_Fig4.pdf").resolve()),
                    "fig4_tif": str((root / "figures" / f"R1_1_{current_id}_Fig4.tif").resolve()),
                }
            )
        index = pd.DataFrame(rows)

    if "definition_id" not in index.columns:
        raise ValueError(f"Reviewer index lacks definition_id: {index_path}")
    matches = index["definition_id"].astype(str).eq(definition_id)
    if int(matches.sum()) != 1:
        raise ValueError(
            f"Reviewer index must contain exactly one {definition_id} row; "
            f"observed {int(matches.sum())}"
        )
    for column in ("fig3_pdf", "fig3_tif", "fig3_status"):
        if column not in index.columns:
            index[column] = ""
    index.loc[matches, "fig3_pdf"] = str(pdf_path.resolve())
    index.loc[matches, "fig3_tif"] = str(tif_path.resolve())
    index.loc[matches, "fig3_status"] = "validated developmental"

    temporary = index_path.with_name(f".{index_path.name}.{os.getpid()}.tmp")
    index_path.parent.mkdir(parents=True, exist_ok=True)
    index.to_csv(temporary, index=False, lineterminator="\n")
    os.replace(temporary, index_path)
    return index_path


def finalize_definition(output_root: Path, definition_id: str) -> dict[str, object]:
    root = definition_root(output_root, definition_id)
    marker = root / "mdsine2" / "finalized.json"
    if marker.is_file():
        payload = json.loads(marker.read_text(encoding="utf-8"))
        update_reviewer_index_fig3(
            output_root,
            definition_id,
            root / "figures" / "Fig3_restored.pdf",
            root / "figures" / "Fig3_restored.tif",
        )
        return payload
    states = definition_state(output_root, definition_id)
    if set(states.values()) != {"succeeded"}:
        raise RuntimeError(f"{definition_id}: cannot finalize chain states {states}")
    specs = {group: make_spec(output_root, definition_id, group) for group in ("LDR", "HDR")}
    mcmc = {group: spec.model_dir / "mcmc.pkl" for group, spec in specs.items()}
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    log_path = root / "mdsine2" / "finalize.log"
    inputs = input_root(output_root, definition_id)

    run_checked(
        [
            str(MDSINE_PYTHON),
            str(FIG3_CODE / "summarize_posteriors.py"),
            "--ldr-mcmc",
            str(mcmc["LDR"]),
            "--hdr-mcmc",
            str(mcmc["HDR"]),
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
            "--evidence-tier",
            "pilot_single_or_short_chain",
        ],
        log_path,
    )
    run_checked(
        [
            str(MDSINE_PYTHON),
            str(FIG3_CODE / "longitudinal_stats.py"),
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
        ],
        log_path,
    )
    run_checked(
        [
            str(MDSINE_PYTHON),
            str(FIG3_CODE / "render_fig3.py"),
            "--project-root",
            str(PROJECT_ROOT),
            "--output",
            str(root),
        ],
        log_path,
    )
    run_checked(
        [
            str(MDSINE_PYTHON),
            str(FIG3_CODE / "validate_fig3.py"),
            "--project-root",
            str(PROJECT_ROOT),
            "--output",
            str(root),
            "--strict",
            "--allow-development",
        ],
        log_path,
    )
    outputs = [root / "figures" / f"Fig3_restored.{suffix}" for suffix in ("pdf", "tif")]
    payload = {
        "definition_id": definition_id,
        "finalized_at_utc": utc_now(),
        "profile": {
            "name": PROFILE.name,
            "chains_per_group": 1,
            "seed": SEED,
            "n_samples": PROFILE.n_samples,
            "burnin": PROFILE.burnin,
            "retained_draws": PROFILE.n_samples - PROFILE.burnin,
            "checkpoint": PROFILE.checkpoint,
            "multiprocessing": 0,
            "evidence_tier": "pilot_single_or_short_chain",
        },
        "interpretation_boundary": (
            "Fast single-chain sensitivity analysis; R-hat is unavailable and the "
            "result is not a substitute for formal multichain convergence."
        ),
        "outputs": [
            {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in outputs
        ],
    }
    runner.atomic_json(marker, payload)
    update_reviewer_index_fig3(
        output_root,
        definition_id,
        outputs[0],
        outputs[1],
    )
    return payload


def write_supervisor_status(output_root: Path, payload: dict[str, object]) -> None:
    path = output_root / "mdsine2_fast_supervisor_status.json"
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def supervise(output_root: Path, poll_seconds: int = 60) -> int:
    active = ["D1", "D2", "D3", "D4", "D5"]
    for definition_id in active:
        for group in ("LDR", "HDR"):
            launch_worker(make_spec(output_root, definition_id, group))
    d6_launched = any(
        chain_state(output_root, "D6", group) != "not_started"
        for group in ("LDR", "HDR")
    )
    finalized: set[str] = {
        definition_id
        for definition_id in DEFINITION_SLUGS
        if (definition_root(output_root, definition_id) / "mdsine2" / "finalized.json").is_file()
    }
    finalization_errors: dict[str, dict[str, str]] = {}
    while True:
        states = {
            definition_id: definition_state(output_root, definition_id)
            for definition_id in DEFINITION_SLUGS
        }
        if not d6_launched:
            free_definition = next(
                (
                    definition_id
                    for definition_id in active
                    if set(states[definition_id].values()) == {"succeeded"}
                ),
                None,
            )
            if free_definition is not None:
                CPU_PAIRS["D6"] = CPU_PAIRS[free_definition].copy()
                for group in ("LDR", "HDR"):
                    launch_worker(make_spec(output_root, "D6", group))
                d6_launched = True
                states["D6"] = definition_state(output_root, "D6")

        for definition_id, group_states in states.items():
            if definition_id not in finalized and set(group_states.values()) == {"succeeded"}:
                try:
                    finalize_definition(output_root, definition_id)
                except Exception as exc:
                    finalization_errors[definition_id] = {
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                        "traceback": traceback.format_exc(),
                    }
                else:
                    finalized.add(definition_id)

        failed = {
            definition_id: group_states
            for definition_id, group_states in states.items()
            if "failed" in group_states.values()
        }
        for definition_id, error in finalization_errors.items():
            failed[definition_id] = {
                **states[definition_id],
                "finalize": "failed",
                "finalization_error": error["error"],
            }
        payload = {
            "updated_at_utc": utc_now(),
            "profile": {
                "n_samples": PROFILE.n_samples,
                "burnin": PROFILE.burnin,
                "checkpoint": PROFILE.checkpoint,
                "chains_per_group": 1,
                "seed": SEED,
            },
            "states": states,
            "finalized": sorted(finalized),
            "failed": failed,
            "finalization_error": finalization_errors,
            "d6_launched": d6_launched,
            "status": "failed" if failed else "running",
        }
        write_supervisor_status(output_root, payload)
        if failed:
            return 1
        if d6_launched and len(finalized) == len(DEFINITION_SLUGS):
            payload["status"] = "complete"
            payload["completed_at_utc"] = utc_now()
            write_supervisor_status(output_root, payload)
            return 0
        time.sleep(poll_seconds)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=("prepare", "launch", "supervise", "worker", "status", "finalize"),
    )
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--group-manifest", type=Path, default=GROUP_MANIFEST)
    parser.add_argument("--definition", choices=tuple(DEFINITION_SLUGS))
    parser.add_argument("--group", choices=("LDR", "HDR"))
    parser.add_argument("--cpu-affinity")
    parser.add_argument("--poll-seconds", type=int, default=60)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    if args.action == "prepare":
        definitions = [args.definition] if args.definition else list(DEFINITION_SLUGS)
        records = [
            prepare_definition(output_root, args.group_manifest.resolve(), definition_id)
            for definition_id in definitions
        ]
        path = output_root / "mdsine2_fast_input_manifest.json"
        path.write_text(json.dumps(records, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps(records, indent=2, ensure_ascii=False))
        return 0
    if args.action == "worker":
        if not args.definition or not args.group:
            raise ValueError("worker requires --definition and --group")
        return execute_worker(
            make_spec(
                output_root,
                args.definition,
                args.group,
                cpu_affinity=args.cpu_affinity,
            )
        )
    if args.action == "supervise":
        return supervise(output_root, args.poll_seconds)
    if args.action == "launch":
        # Input preparation is synchronous so no chain can start from a partial bundle.
        records = [
            prepare_definition(output_root, args.group_manifest.resolve(), definition_id)
            for definition_id in DEFINITION_SLUGS
        ]
        (output_root / "mdsine2_fast_input_manifest.json").write_text(
            json.dumps(records, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        session = "r11m2-fast-supervisor"
        if runner.tmux_session_exists(session):
            print(json.dumps({"session": session, "action": "unchanged", "state": "running"}))
            return 0
        command = [
            "tmux",
            "new-session",
            "-d",
            "-s",
            session,
            str(MDSINE_PYTHON),
            str(Path(__file__).resolve()),
            "supervise",
            "--output-root",
            str(output_root),
            "--poll-seconds",
            str(args.poll_seconds),
        ]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode:
            raise RuntimeError((completed.stderr or completed.stdout).strip())
        print(json.dumps({"session": session, "action": "launched", "definitions": list(DEFINITION_SLUGS)}))
        return 0
    if args.action == "status":
        payload = {
            definition_id: definition_state(output_root, definition_id)
            for definition_id in DEFINITION_SLUGS
        }
        status_path = output_root / "mdsine2_fast_supervisor_status.json"
        if status_path.is_file():
            payload["supervisor"] = json.loads(status_path.read_text(encoding="utf-8"))
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 0
    if args.action == "finalize":
        definitions = [args.definition] if args.definition else list(DEFINITION_SLUGS)
        results = [finalize_definition(output_root, definition_id) for definition_id in definitions]
        print(json.dumps(results, indent=2, ensure_ascii=False))
        return 0
    raise AssertionError(args.action)


if __name__ == "__main__":
    raise SystemExit(main())
