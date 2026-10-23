#!/usr/bin/env python3
"""Prepare auditable, no-qPCR MDSINE2 inputs for restored Figure 3.

The Bracken MPA table is a hierarchical table.  Only rows whose lineage starts
with ``k__Bacteria`` and whose terminal rank is ``s__`` are treated as bacterial
species leaves.  Rarefaction is deliberately applied to that *complete* species
matrix before selecting the fixed 259-species panel used by the source P4 table.

The generated qPCR table is only a structural API placeholder.  The companion
``mdsine2_ra.py`` runner replaces every value with a biomass prior (mean 1e11,
CV 10%); the placeholder must never be interpreted as measured qPCR or absolute
biomass.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MPA = PROJECT_ROOT / "raw" / "combine_brk.mpa"
DEFAULT_METADATA = PROJECT_ROOT / "raw" / "meta.csv"
DEFAULT_ASSIGNMENTS = (
    PROJECT_ROOT
    / "reproduction"
    / "outputs"
    / "integrated_fig124_raw"
    / "tables"
    / "raw_bracken_hdr_ldr_assignments.csv"
)
DEFAULT_COMMON_TABLE = Path(
    "/home/qqt/MDSINE2_Server_Package/AbxData/P4_OTU_table.csv"
)
DEFAULT_TAXONOMY_SOURCE = Path(
    "/home/qqt/MDSINE2_Server_Package/AbxData/P4_taxonomy.csv"
)
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "reproduction" / "outputs" / "fig3_mdsine2" / "inputs"
)

RAREFACTION_DEPTH = 2_693_000
RAREFACTION_SEED = 20_261_009
PERTURBATION_START = 0.5
PERTURBATION_END = 5.5
GROUP_ALIASES = {"LDR": "P1", "HDR": "P4"}
GROUP_ORDER = {"LDR": 0, "HDR": 1}
QPCR_PLACEHOLDER_VALUE = 1.0
SEQUENCE_PLACEHOLDER = "A" * 50


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest without modifying ``path``."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _require_columns(frame: pd.DataFrame, required: Iterable[str], label: str) -> None:
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def _jsonable(value: Any) -> Any:
    """Convert NumPy/Pandas scalar containers into JSON-safe Python objects."""

    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if pd.isna(value):
        return None
    return value


def write_json(payload: Mapping[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def load_common_taxa(common_table_path: Path) -> list[str]:
    """Load the source-calibrated taxon order from P4_OTU_table.csv."""

    table = pd.read_csv(common_table_path, usecols=[0])
    taxa = table.iloc[:, 0].astype(str).tolist()
    if len(taxa) != len(set(taxa)):
        raise ValueError("The fixed common-taxon table contains duplicate row names")
    if not taxa or not all(name.startswith("s__") for name in taxa):
        raise ValueError("The fixed common-taxon table is not a species-level s__ panel")
    return taxa


def load_selected_metadata(
    metadata_path: Path,
    assignments_path: Path,
    expected_subjects_per_group: int = 15,
    expected_samples_per_subject: int = 16,
) -> pd.DataFrame:
    """Join HDR/LDR assignments to sample metadata and impose stable ordering."""

    metadata = pd.read_csv(metadata_path)
    assignments = pd.read_csv(assignments_path)
    _require_columns(
        metadata,
        ["ids", "subject", "day", "canon", "canons"],
        "sample metadata",
    )
    _require_columns(
        assignments,
        ["subject", "rank", "recovery_group"],
        "HDR/LDR assignments",
    )

    if metadata["ids"].astype(str).duplicated().any():
        raise ValueError("Sample IDs are not unique in metadata")
    selected_assignments = assignments.loc[
        assignments["recovery_group"].isin(GROUP_ORDER),
        ["subject", "rank", "recovery_group"],
    ].copy()
    if selected_assignments["subject"].astype(str).duplicated().any():
        raise ValueError("A selected subject has more than one HDR/LDR assignment")

    group_sizes = selected_assignments.groupby("recovery_group")["subject"].nunique()
    for group in GROUP_ORDER:
        observed = int(group_sizes.get(group, 0))
        if observed != expected_subjects_per_group:
            raise ValueError(
                f"Expected {expected_subjects_per_group} {group} subjects, found {observed}"
            )

    selected = metadata.merge(
        selected_assignments,
        on="subject",
        how="inner",
        validate="many_to_one",
    )
    selected["ids"] = selected["ids"].astype(str)
    selected["subject"] = selected["subject"].astype(str)
    selected["day"] = pd.to_numeric(selected["day"], errors="raise")
    selected["rank"] = pd.to_numeric(selected["rank"], errors="raise").astype(int)

    counts = selected.groupby("subject")["ids"].size()
    bad_counts = counts[counts != expected_samples_per_subject]
    if not bad_counts.empty:
        raise ValueError(
            "Selected subjects do not all have the expected number of samples: "
            + bad_counts.to_dict().__repr__()
        )
    if selected.duplicated(["subject", "day"]).any():
        duplicate = selected.loc[
            selected.duplicated(["subject", "day"], keep=False),
            ["subject", "day", "ids"],
        ]
        raise ValueError(
            "MDSINE2 requires a unique time per subject; duplicates found: "
            + duplicate.head().to_dict("records").__repr__()
        )

    selected["_group_order"] = selected["recovery_group"].map(GROUP_ORDER)
    selected = selected.sort_values(
        ["_group_order", "rank", "day", "ids"], kind="stable"
    ).drop(columns="_group_order")
    return selected.reset_index(drop=True)


def read_bacterial_species_leaves(
    mpa_path: Path,
    sample_ids: Iterable[str],
) -> tuple[pd.DataFrame, pd.Series]:
    """Read selected samples and return all bacterial terminal species rows.

    Returns
    -------
    counts
        Integer species-by-sample table indexed by terminal ``s__`` names.
    lineages
        Mapping from terminal species name to full MPA lineage.
    """

    sample_ids = [str(value) for value in sample_ids]
    header = pd.read_csv(mpa_path, sep="\t", nrows=0).columns.tolist()
    if not header:
        raise ValueError("The Bracken MPA table has no columns")
    lineage_column = header[0]
    source_by_id: dict[str, str] = {}
    for column in header[1:]:
        sample_id = (
            column[: -len(".brk.report")]
            if column.endswith(".brk.report")
            else column
        )
        if sample_id in source_by_id:
            raise ValueError(f"Duplicate Bracken sample after suffix removal: {sample_id}")
        source_by_id[sample_id] = column

    missing = sorted(set(sample_ids) - set(source_by_id))
    if missing:
        raise ValueError(f"Selected sample IDs absent from Bracken table: {missing[:10]}")

    source_columns = [source_by_id[sample_id] for sample_id in sample_ids]
    raw = pd.read_csv(
        mpa_path,
        sep="\t",
        usecols=[lineage_column, *source_columns],
    )
    lineage = raw[lineage_column].astype(str)
    terminal = lineage.str.rsplit("|", n=1).str[-1]
    mask = lineage.str.startswith("k__Bacteria") & terminal.str.startswith("s__")
    species = raw.loc[mask, source_columns].copy()
    species.columns = sample_ids
    species.index = terminal.loc[mask].to_numpy()
    species.index.name = "name"

    if species.index.duplicated().any():
        duplicate = species.index[species.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate terminal bacterial species labels: {duplicate[:10]}")
    values = species.to_numpy()
    if not np.issubdtype(values.dtype, np.integer):
        if np.any(~np.isfinite(values)) or np.any(values != np.floor(values)):
            raise ValueError("Bracken species counts must be finite integers")
        species = species.astype(np.int64)
        values = species.to_numpy()
    if np.any(values < 0):
        raise ValueError("Bracken species counts contain negative values")

    lineage_map = pd.Series(
        lineage.loc[mask].to_numpy(), index=species.index, name="source_lineage"
    )
    return species.astype(np.int64, copy=False), lineage_map


def rarefy_counts(
    counts: pd.DataFrame,
    depth: int = RAREFACTION_DEPTH,
    seed: int = RAREFACTION_SEED,
    sample_stable: bool = False,
) -> pd.DataFrame:
    """Rarefy columns without replacement using a stable sequential RNG stream."""

    if depth <= 0:
        raise ValueError("Rarefaction depth must be positive")
    values = counts.to_numpy(dtype=np.int64, copy=False)
    if values.ndim != 2 or values.shape[0] == 0 or values.shape[1] == 0:
        raise ValueError("Rarefaction input must be a non-empty 2D count table")
    if np.any(values < 0):
        raise ValueError("Rarefaction input contains negative values")
    observed_depths = values.sum(axis=0, dtype=np.int64)
    bad = np.flatnonzero(observed_depths < depth)
    if len(bad):
        detail = {
            str(counts.columns[i]): int(observed_depths[i]) for i in bad[:10]
        }
        raise ValueError(f"Samples below rarefaction depth {depth}: {detail}")

    rarefied = np.empty(values.shape, dtype=np.int64)
    for column_index in range(values.shape[1]):
        if sample_stable:
            sample_id = str(counts.columns[column_index])
            digest = hashlib.sha256(f"{seed}\0{sample_id}".encode("utf-8")).digest()
            sample_seed = int.from_bytes(digest[:8], byteorder="little", signed=False)
            rng = np.random.default_rng(sample_seed)
        elif column_index == 0:
            rng = np.random.default_rng(seed)
        rarefied[:, column_index] = rng.multivariate_hypergeometric(
            values[:, column_index], depth
        )
    return pd.DataFrame(rarefied, index=counts.index, columns=counts.columns)


def build_taxonomy_table(
    common_taxa: list[str], taxonomy_source_path: Path
) -> pd.DataFrame:
    """Build the MDSINE2 taxonomy table in the exact fixed-taxon order."""

    source = pd.read_csv(taxonomy_source_path)
    _require_columns(
        source,
        [
            "Kingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
            "Taxa_Clean",
        ],
        "taxonomy source",
    )
    if source["Taxa_Clean"].astype(str).duplicated().any():
        raise ValueError("Taxonomy source has duplicate Taxa_Clean identifiers")
    taxonomy = source.set_index(source["Taxa_Clean"].astype(str)).reindex(common_taxa)
    if taxonomy["Taxa_Clean"].isna().any():
        missing = taxonomy.index[taxonomy["Taxa_Clean"].isna()].tolist()
        raise ValueError(f"Fixed taxa absent from taxonomy source: {missing[:10]}")
    taxonomy = taxonomy[
        ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    ].copy()
    taxonomy = taxonomy.fillna("NA").astype(str)
    taxonomy["sequence"] = SEQUENCE_PLACEHOLDER
    taxonomy.index.name = "OTU_ID"
    return taxonomy


def create_study_pickle(
    group: str,
    group_directory: Path,
    pickle_path: Path,
) -> None:
    """Parse exported tables with MDSINE2 and serialize a validated Study."""

    try:
        import mdsine2 as md2
    except ImportError as exc:  # pragma: no cover - depends on execution environment
        raise RuntimeError(
            "MDSINE2 is required to create study pickles. Run this script with "
            "/home/qqt/.conda/envs/mdsine2/bin/python or pass --skip-study-pickle."
        ) from exc

    study = md2.dataset.parse(
        name=f"Fig3_{group}_{GROUP_ALIASES[group]}_no_qPCR",
        metadata=str(group_directory / "metadata.tsv"),
        taxonomy=str(group_directory / "taxonomy.tsv"),
        reads=str(group_directory / "reads.tsv"),
        qpcr=str(group_directory / "qpcr_structural_placeholder.tsv"),
        perturbations=str(group_directory / "perturbations.tsv"),
        sep="\t",
    )
    if len(study.taxa) != 259:
        raise RuntimeError(f"Parsed Study has {len(study.taxa)} taxa instead of 259")
    if len(study) != 15:
        raise RuntimeError(f"Parsed Study has {len(study)} subjects instead of 15")
    study.save(pickle_path)


def prepare_inputs(
    *,
    mpa_path: Path = DEFAULT_MPA,
    metadata_path: Path = DEFAULT_METADATA,
    assignments_path: Path = DEFAULT_ASSIGNMENTS,
    common_table_path: Path = DEFAULT_COMMON_TABLE,
    taxonomy_source_path: Path = DEFAULT_TAXONOMY_SOURCE,
    output_dir: Path = DEFAULT_OUTPUT,
    depth: int = RAREFACTION_DEPTH,
    seed: int = RAREFACTION_SEED,
    expected_subjects_per_group: int = 15,
    expected_samples_per_subject: int = 16,
    make_study_pickles: bool = True,
    sample_stable_rarefaction: bool = False,
) -> dict[str, Any]:
    """Prepare both group inputs and return the global audit manifest."""

    paths = {
        "bracken_mpa": Path(mpa_path).resolve(),
        "sample_metadata": Path(metadata_path).resolve(),
        "hdr_ldr_assignments": Path(assignments_path).resolve(),
        "fixed_259_table": Path(common_table_path).resolve(),
        "taxonomy_source": Path(taxonomy_source_path).resolve(),
    }
    for label, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    selected = load_selected_metadata(
        paths["sample_metadata"],
        paths["hdr_ldr_assignments"],
        expected_subjects_per_group=expected_subjects_per_group,
        expected_samples_per_subject=expected_samples_per_subject,
    )
    sample_ids = selected["ids"].tolist()
    common_taxa = load_common_taxa(paths["fixed_259_table"])
    full_counts, lineage_map = read_bacterial_species_leaves(
        paths["bracken_mpa"], sample_ids
    )
    missing_taxa = sorted(set(common_taxa) - set(full_counts.index))
    if missing_taxa:
        raise ValueError(f"Fixed taxa absent from Bracken species leaves: {missing_taxa[:10]}")

    raw_depths = pd.Series(
        full_counts.to_numpy(dtype=np.int64, copy=False).sum(
            axis=0, dtype=np.int64
        ),
        index=full_counts.columns,
    )
    rarefied_full = rarefy_counts(
        full_counts,
        depth=depth,
        seed=seed,
        sample_stable=sample_stable_rarefaction,
    )
    rarefied_depths = pd.Series(
        rarefied_full.to_numpy(dtype=np.int64, copy=False).sum(
            axis=0, dtype=np.int64
        ),
        index=rarefied_full.columns,
    )
    if not np.all(rarefied_depths.to_numpy() == depth):
        raise RuntimeError("Rarefaction did not produce the requested exact depth")
    fixed_counts = rarefied_full.loc[common_taxa].copy()
    retained_depths = pd.Series(
        fixed_counts.to_numpy(dtype=np.int64, copy=False).sum(
            axis=0, dtype=np.int64
        ),
        index=fixed_counts.columns,
    )
    taxonomy = build_taxonomy_table(common_taxa, paths["taxonomy_source"])

    sample_audit = selected.copy()
    sample_audit.insert(0, "sampleID", sample_audit["ids"])
    sample_audit["source_species_depth"] = sample_audit["ids"].map(raw_depths)
    sample_audit["rarefied_full_species_depth"] = sample_audit["ids"].map(
        rarefied_depths
    )
    sample_audit["retained_259_depth"] = sample_audit["ids"].map(retained_depths)
    sample_audit["retained_259_fraction"] = (
        sample_audit["retained_259_depth"] / depth
    )
    sample_audit.to_csv(output_dir / "sample_audit.tsv", sep="\t", index=False)

    assignment_audit = (
        selected[["subject", "rank", "recovery_group"]]
        .drop_duplicates()
        .sort_values(
            ["recovery_group", "rank"],
            key=lambda series: (
                series.map(GROUP_ORDER)
                if series.name == "recovery_group"
                else series
            ),
        )
    )
    assignment_audit["source_group_alias"] = assignment_audit["recovery_group"].map(
        GROUP_ALIASES
    )
    assignment_audit.to_csv(
        output_dir / "subject_assignments.tsv", sep="\t", index=False
    )

    species_audit = pd.DataFrame(
        {
            "fixed_order": np.arange(1, len(common_taxa) + 1),
            "taxon": common_taxa,
            "source_lineage": lineage_map.reindex(common_taxa).to_numpy(),
            "total_rarefied_reads_all_selected_samples": fixed_counts.sum(axis=1)
            .astype(np.int64)
            .to_numpy(),
            "sample_prevalence_all_selected_samples": (fixed_counts > 0)
            .sum(axis=1)
            .astype(int)
            .to_numpy(),
        }
    )
    for group in GROUP_ORDER:
        group_ids = selected.loc[selected["recovery_group"] == group, "ids"]
        species_audit[f"total_rarefied_reads_{group}"] = (
            fixed_counts[group_ids].sum(axis=1).astype(np.int64).to_numpy()
        )
        species_audit[f"sample_prevalence_{group}"] = (
            (fixed_counts[group_ids] > 0).sum(axis=1).astype(int).to_numpy()
        )
    species_audit.to_csv(output_dir / "species_audit.tsv", sep="\t", index=False)

    group_manifests: dict[str, Any] = {}
    for group in GROUP_ORDER:
        group_dir = output_dir / group
        group_dir.mkdir(parents=True, exist_ok=True)
        group_meta = selected.loc[selected["recovery_group"] == group].copy()
        group_ids = group_meta["ids"].tolist()

        mdsine_metadata = pd.DataFrame(
            {
                "sampleID": group_ids,
                "subject": group_meta["subject"].to_numpy(),
                "time": group_meta["day"].astype(float).to_numpy(),
                "canon": group_meta["canon"].astype(str).to_numpy(),
                "canon_index": group_meta["canons"].to_numpy(),
                "recovery_group": group,
                "source_group_alias": GROUP_ALIASES[group],
            }
        )
        mdsine_metadata.to_csv(group_dir / "metadata.tsv", sep="\t", index=False)
        fixed_counts[group_ids].to_csv(group_dir / "reads.tsv", sep="\t")
        taxonomy.to_csv(group_dir / "taxonomy.tsv", sep="\t")

        perturbations = pd.DataFrame(
            {
                "name": "ciprofloxacin",
                "subject": group_meta["subject"].drop_duplicates().tolist(),
                "start": PERTURBATION_START,
                "end": PERTURBATION_END,
                "value": 1.0,
            }
        )
        perturbations.to_csv(group_dir / "perturbations.tsv", sep="\t", index=False)

        qpcr_placeholder = pd.DataFrame(
            {
                "sampleID": group_ids,
                "structural_placeholder_not_measured_qpcr": QPCR_PLACEHOLDER_VALUE,
            }
        )
        qpcr_placeholder.to_csv(
            group_dir / "qpcr_structural_placeholder.tsv", sep="\t", index=False
        )

        pickle_path = group_dir / f"{group}_study.pkl"
        if make_study_pickles:
            create_study_pickle(group, group_dir, pickle_path)

        group_sample_audit = sample_audit.loc[
            sample_audit["recovery_group"] == group
        ].copy()
        group_manifest = {
            "group": group,
            "source_group_alias": GROUP_ALIASES[group],
            "subjects": group_meta["subject"].drop_duplicates().tolist(),
            "n_subjects": group_meta["subject"].nunique(),
            "n_samples": len(group_meta),
            "n_taxa": len(common_taxa),
            "actual_day_used_as_mdsine2_time": True,
            "canonical_stage_is_plotting_metadata_only": True,
            "ciprofloxacin_interval_days": [PERTURBATION_START, PERTURBATION_END],
            "retained_259_depth": {
                "min": int(group_sample_audit["retained_259_depth"].min()),
                "median": float(group_sample_audit["retained_259_depth"].median()),
                "max": int(group_sample_audit["retained_259_depth"].max()),
            },
            "files": {
                "metadata": str(group_dir / "metadata.tsv"),
                "reads": str(group_dir / "reads.tsv"),
                "taxonomy": str(group_dir / "taxonomy.tsv"),
                "perturbations": str(group_dir / "perturbations.tsv"),
                "qpcr_structural_placeholder": str(
                    group_dir / "qpcr_structural_placeholder.tsv"
                ),
                "study_pickle": str(pickle_path) if make_study_pickles else None,
            },
            "qpcr_semantics": {
                "is_measured_qpcr": False,
                "placeholder_value": QPCR_PLACEHOLDER_VALUE,
                "interpretation": "structural API placeholder only; not biological data",
                "runner_override": (
                    "mdsine2_ra.py replaces values with biomass mean 1e11 and uses "
                    "a fixed 10% coefficient of variation"
                ),
            },
        }
        write_json(group_manifest, group_dir / "input_manifest.json")
        group_manifests[group] = group_manifest

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "purpose": "Figure 3 HDR/LDR longitudinal MDSINE2 no-qPCR input preparation",
        "sources": {
            label: {"path": str(path), "sha256": sha256_file(path)}
            for label, path in paths.items()
        },
        "selection": {
            "groups": ["LDR", "HDR"],
            "group_aliases": GROUP_ALIASES,
            "expected_subjects_per_group": expected_subjects_per_group,
            "observed_subjects_per_group": selected.groupby("recovery_group")[
                "subject"
            ]
            .nunique()
            .to_dict(),
            "samples_per_subject": expected_samples_per_subject,
            "observed_total_samples": len(selected),
        },
        "taxon_processing": {
            "bacterial_species_leaf_definition": (
                "lineage starts k__Bacteria and terminal pipe-delimited rank starts s__"
            ),
            "n_full_bacterial_species_leaves": len(full_counts),
            "n_fixed_common_species": len(common_taxa),
            "fixed_species_order_source": str(paths["fixed_259_table"]),
            "operation_order": [
                "extract all bacterial species leaves",
                "rarefy each full species vector",
                "subset and reorder the fixed 259 species",
            ],
        },
        "rarefaction": {
            "depth": depth,
            "seed": seed,
            "algorithm": "numpy.random.Generator.multivariate_hypergeometric",
            "sampling": "without replacement",
            "rng_stream_order": (
                "independent SHA256(global seed, sample ID)-derived stream; invariant to group definition"
                if sample_stable_rarefaction
                else "LDR then HDR; assignment rank; actual day; sample ID"
            ),
            "sample_stable_across_group_definitions": sample_stable_rarefaction,
            "numpy_version": np.__version__,
        },
        "time_and_perturbation": {
            "mdsine2_time": "raw/meta.csv day",
            "canonical_stage_usage": "plotting/annotation only",
            "perturbation": "ciprofloxacin",
            "start_day": PERTURBATION_START,
            "end_day": PERTURBATION_END,
        },
        "structural_placeholders": {
            "qpcr": (
                "one value of 1.0 per sample for Study/API structure; it is not a "
                "measurement and is overridden by the no-qPCR runner"
            ),
            "sequence": (
                "50 adenines per taxon because no sequence input was supplied; not "
                "used as biological sequence evidence"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "study_pickles_created": make_study_pickles,
        "groups": group_manifests,
        "audit_files": {
            "samples": str(output_dir / "sample_audit.tsv"),
            "subjects": str(output_dir / "subject_assignments.tsv"),
            "species": str(output_dir / "species_audit.tsv"),
        },
    }
    write_json(manifest, output_dir / "input_manifest.json")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mpa", type=Path, default=DEFAULT_MPA)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--assignments", type=Path, default=DEFAULT_ASSIGNMENTS)
    parser.add_argument("--common-table", type=Path, default=DEFAULT_COMMON_TABLE)
    parser.add_argument("--taxonomy-source", type=Path, default=DEFAULT_TAXONOMY_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--depth", type=int, default=RAREFACTION_DEPTH)
    parser.add_argument("--seed", type=int, default=RAREFACTION_SEED)
    parser.add_argument(
        "--skip-study-pickle",
        action="store_true",
        help="Export tables only (unit-test/debug option; formal pipeline should create pickles).",
    )
    parser.add_argument(
        "--sample-stable-rarefaction",
        action="store_true",
        help=(
            "Derive one deterministic RNG stream per sample ID so the same sample "
            "has identical rarefied counts across alternative group definitions."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = prepare_inputs(
        mpa_path=args.mpa,
        metadata_path=args.metadata,
        assignments_path=args.assignments,
        common_table_path=args.common_table,
        taxonomy_source_path=args.taxonomy_source,
        output_dir=args.output_dir,
        depth=args.depth,
        seed=args.seed,
        make_study_pickles=not args.skip_study_pickle,
        sample_stable_rarefaction=args.sample_stable_rarefaction,
    )
    print(
        f"Prepared {manifest['selection']['observed_total_samples']} samples, "
        f"{manifest['taxon_processing']['n_fixed_common_species']} taxa -> "
        f"{Path(args.output_dir).resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
