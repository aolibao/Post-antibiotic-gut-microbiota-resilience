#!/usr/bin/env python3
"""Load the supplied Kraken2/Bracken MPA table for Figures 1, 2 and 4.

Reverse calibration against the source raster establishes that the plotted
alpha-diversity points were calculated from unrarefied bacterial species-leaf
counts, normalized to proportions only for Shannon and composition.  This
source-faithful branch reproduces the Figure 2C statistics and Figure 2D/E
marker coordinates.  A separate sensitivity export implements the manuscript
Methods statement (rarefaction to 2,693,000 reads) so the text/figure conflict
is explicit rather than silently resolved in either direction.

The original rarefaction seed and the original HDR/LDR subject list were not
reported.  The subject labels are therefore reconstructed from the supplied
``species.csv`` endpoint-diversity ranking, which is the only supplied ranking
that reproduces the correlations printed in manuscript Figure 2C when the
Bracken Shannon values are plotted.  Both the grouping score and the plotted
Bracken values are retained in ``group_table`` for auditability.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_all as audit


RAREFACTION_DEPTH = 2_693_000
RAREFACTION_SEED = audit.SEED
MPA_NAME = "combine_brk.mpa"


def _rank_mask(index: pd.Index, prefix: str, bacteria_only: bool = True) -> np.ndarray:
    terminal = index.to_series(index=index).str.rsplit("|", n=1).str[-1]
    selected = terminal.str.startswith(prefix)
    if bacteria_only:
        selected &= index.to_series(index=index).str.startswith("k__Bacteria|")
    return selected.to_numpy()


def _terminal_names(index: pd.Index, prefix: str) -> pd.Index:
    names = (
        index.to_series(index=index)
        .str.rsplit("|", n=1)
        .str[-1]
        .str.removeprefix(prefix)
        .str.replace("_", " ", regex=False)
        .str.strip()
        .replace("", "Unclassified")
    )
    return pd.Index(names.to_numpy(), name=prefix[0])


def _taxonomy_frame(paths: pd.Index) -> pd.DataFrame:
    ranks = ["k", "p", "c", "o", "f", "g", "s"]
    rows: list[dict[str, str]] = []
    for path in paths.astype(str):
        parsed = {part[0]: part for part in path.split("|") if len(part) > 2 and part[1:3] == "__"}
        rows.append({rank: parsed.get(rank, f"{rank}__Unclassified") for rank in ranks})
    frame = pd.DataFrame(rows, index=paths)
    frame.index.name = "taxon"
    return frame


def _rarefy(counts: pd.DataFrame, depth: int, seed: int) -> pd.DataFrame:
    totals = counts.sum(axis=0).to_numpy(dtype=np.int64)
    if np.any(totals < depth):
        bad = counts.columns[totals < depth].tolist()
        raise ValueError(f"{len(bad)} samples are below rarefaction depth {depth}: {bad[:5]}")
    values = counts.to_numpy(dtype=np.int64, copy=False)
    rarefied = np.empty(values.shape, dtype=np.int32)
    rng = np.random.default_rng(seed)
    for column in range(values.shape[1]):
        rarefied[:, column] = rng.multivariate_hypergeometric(values[:, column], depth)
    return pd.DataFrame(rarefied, index=counts.index, columns=counts.columns)


def _alpha_metrics(percent: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    prop = percent / 100.0
    log_prop = np.zeros(prop.shape, dtype=float)
    values = prop.to_numpy(dtype=float)
    np.log(values, out=log_prop, where=values > 0)
    shannon = pd.Series(-(values * log_prop).sum(axis=0), index=percent.columns)
    richness = (percent > 0).sum(axis=0).astype(float)
    simpson = 1.0 - prop.pow(2).sum(axis=0)
    evenness = shannon / np.log(richness.where(richness > 1))
    return shannon, richness, simpson, evenness


def _legacy_grouping(root: Path, meta: pd.DataFrame) -> pd.DataFrame:
    legacy = pd.read_csv(root / "species.csv", index_col=0)
    legacy = legacy.loc[:, meta["ids"]]
    prop = legacy.div(legacy.sum(axis=0).replace(0, np.nan), axis=1)
    shannon = -(prop * np.log(prop.where(prop > 0))).sum(axis=0)
    frame = meta.copy()
    frame["shannon"] = frame["ids"].map(shannon)
    exposed = frame.loc[frame["abx"] == "Y"]
    baseline = exposed.loc[exposed["canon"].isin(audit.BASELINE)].groupby("subject")["shannon"].mean()
    endpoint = (
        exposed.loc[exposed["canon"].isin(audit.ENDPOINT)]
        .groupby("subject")["shannon"]
        .mean()
        .sort_values(kind="mergesort")
    )
    out = pd.DataFrame(
        {
            "grouping_baseline_shannon": baseline,
            "grouping_endpoint_shannon": endpoint,
        }
    )
    out["rank"] = out["grouping_endpoint_shannon"].rank(method="first").astype(int)
    out["recovery_group"] = "Middle"
    out.loc[endpoint.index[:15], "recovery_group"] = "LDR"
    out.loc[endpoint.index[-15:], "recovery_group"] = "HDR"
    out.index.name = "subject"
    return out


def load_raw_bracken_dataset(root: Path) -> audit.Dataset:
    root = Path(root).resolve()
    raw = root / "raw"
    meta = pd.read_csv(
        raw / "meta.csv",
        dtype={"ids": str, "library": str, "subject": str, "canon": str},
    )
    mpa = pd.read_csv(raw / MPA_NAME, sep="\t", index_col=0)
    mpa.columns = mpa.columns.str.removesuffix(".brk.report")

    if meta["ids"].duplicated().any():
        raise ValueError("raw/meta.csv contains duplicate sample ids")
    if set(meta["ids"]) != set(mpa.columns):
        raise ValueError("raw metadata and combine_brk.mpa sample ids differ")
    sample_order = meta["ids"].tolist()

    species_counts = mpa.loc[_rank_mask(mpa.index, "s__"), sample_order]
    genus_counts = mpa.loc[_rank_mask(mpa.index, "g__"), sample_order]
    if (species_counts.to_numpy() < 0).any() or (genus_counts.to_numpy() < 0).any():
        raise ValueError("Bracken table contains negative counts")

    abundance = species_counts.div(species_counts.sum(axis=0), axis=1) * 100.0
    genus_rank = genus_counts.div(genus_counts.sum(axis=0), axis=1) * 100.0

    taxa = _taxonomy_frame(abundance.index)
    species_named = abundance.groupby(_terminal_names(abundance.index, "s__"), sort=False).sum()
    genus = genus_rank.groupby(_terminal_names(genus_rank.index, "g__"), sort=False).sum()

    shannon, richness, simpson, evenness = _alpha_metrics(abundance)
    metrics = meta.copy()
    metrics["day"] = pd.to_numeric(metrics["day"], errors="coerce")
    metrics["canons"] = pd.to_numeric(metrics["canons"], errors="coerce")
    metrics["shannon"] = metrics["ids"].map(shannon)
    metrics["richness"] = metrics["ids"].map(richness)
    metrics["simpson"] = metrics["ids"].map(simpson)
    metrics["evenness"] = metrics["ids"].map(evenness)
    metrics["canon"] = pd.Categorical(metrics["canon"], audit.TIME_ORDER, ordered=True)

    presence = abundance.to_numpy() > 0
    prop = abundance / 100.0
    sample_index = {sample: i for i, sample in enumerate(abundance.columns)}
    paired_jaccard: dict[str, float] = {}
    bray_to_baseline_centroid: dict[str, float] = {}
    mean_jaccard_to_baseline: dict[str, float] = {}
    for _, subject_rows in metrics.groupby("subject", observed=True):
        a1_id = str(subject_rows.loc[subject_rows["canon"] == "A1", "ids"].iloc[0])
        a1 = presence[:, sample_index[a1_id]]
        baseline_ids = subject_rows.loc[subject_rows["canon"].isin(audit.BASELINE), "ids"].tolist()
        baseline_prop = prop[baseline_ids].mean(axis=1).to_numpy()
        baseline_binary = [presence[:, sample_index[x]] for x in baseline_ids]
        for sample in subject_rows["ids"]:
            current_presence = presence[:, sample_index[sample]]
            union = np.logical_or(a1, current_presence).sum()
            paired_jaccard[sample] = 0.0 if union == 0 else 1.0 - np.logical_and(a1, current_presence).sum() / union
            current_prop = prop[sample].to_numpy()
            denom = np.abs(current_prop + baseline_prop).sum()
            bray_to_baseline_centroid[sample] = 0.0 if denom == 0 else np.abs(current_prop - baseline_prop).sum() / denom
            distances = []
            for baseline_presence in baseline_binary:
                baseline_union = np.logical_or(baseline_presence, current_presence).sum()
                distances.append(
                    0.0
                    if baseline_union == 0
                    else 1.0 - np.logical_and(baseline_presence, current_presence).sum() / baseline_union
                )
            mean_jaccard_to_baseline[sample] = float(np.mean(distances))
    metrics["jaccard_to_A1"] = metrics["ids"].map(paired_jaccard)
    metrics["bray_to_baseline_centroid"] = metrics["ids"].map(bray_to_baseline_centroid)
    metrics["mean_jaccard_to_baseline"] = metrics["ids"].map(mean_jaccard_to_baseline)

    grouping = _legacy_grouping(root, meta)
    exposed = metrics.loc[metrics["abx"] == "Y"]
    raw_baseline = exposed.loc[exposed["canon"].isin(audit.BASELINE)].groupby("subject")["shannon"].mean()
    raw_endpoint = exposed.loc[exposed["canon"].isin(audit.ENDPOINT)].groupby("subject")["shannon"].mean()
    grouping["baseline_shannon"] = raw_baseline
    grouping["endpoint_shannon"] = raw_endpoint
    metrics = metrics.merge(grouping[["recovery_group"]], left_on="subject", right_index=True, how="left")
    metrics["recovery_group"] = metrics["recovery_group"].fillna("Control")

    return audit.Dataset(root, meta, abundance, taxa, metrics, genus, species_named, grouping)


def write_raw_audit(dataset: audit.Dataset, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    source = dataset.root / "raw" / MPA_NAME
    payload = {
        "taxonomic_source": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "figure_branch": "source-faithful unrarefied bacterial species/genus leaf counts",
        "methods_statement_rarefaction_depth": RAREFACTION_DEPTH,
        "sensitivity_rarefaction_seed": RAREFACTION_SEED,
        "sensitivity_rarefaction_algorithm": "NumPy Generator.multivariate_hypergeometric (without replacement)",
        "species_features_before_name_collapse": int(dataset.abundance.shape[0]),
        "named_species_after_collapse": int(dataset.species_named.shape[0]),
        "named_genera_after_collapse": int(dataset.genus.shape[0]),
        "samples": int(dataset.meta.shape[0]),
        "observed_source_figure_richness": {
            "min": int(dataset.metrics["richness"].min()),
            "median": float(dataset.metrics["richness"].median()),
            "max": int(dataset.metrics["richness"].max()),
        },
        "grouping_reconstruction": (
            "bottom/top 15 C4-C6 mean Shannon subjects from supplied species.csv; "
            "plotted alpha diversity and taxonomic abundance come from bacterial Bracken leaves"
        ),
        "methods_figure_conflict": (
            "Methods states rarefaction to 2,693,000 reads, but source-raster coordinate calibration "
            "and exact Figure 2C statistic recovery identify unrarefied bacterial leaf counts as the "
            "actual figure input; the rarefied sensitivity table is exported separately"
        ),
        "group_counts": dataset.group_table["recovery_group"].value_counts().to_dict(),
    }
    (output / "raw_bracken_audit.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    dataset.metrics.sort_values(["subject", "canons"]).to_csv(output / "raw_bracken_sample_metrics.csv", index=False)
    dataset.group_table.sort_values("grouping_endpoint_shannon").to_csv(
        output / "raw_bracken_hdr_ldr_assignments.csv"
    )

    # Methods-faithful sensitivity analysis.  This is not fed into the source-
    # faithful figures, because doing so measurably shifts the original points.
    mpa = pd.read_csv(source, sep="\t", index_col=0)
    mpa.columns = mpa.columns.str.removesuffix(".brk.report")
    counts = mpa.loc[_rank_mask(mpa.index, "s__"), dataset.meta["ids"]]
    rare = _rarefy(counts, RAREFACTION_DEPTH, RAREFACTION_SEED)
    percent = rare.astype(float) * (100.0 / RAREFACTION_DEPTH)
    shannon, richness, simpson, evenness = _alpha_metrics(percent)
    sensitivity = dataset.meta[["ids", "subject", "canon", "abx"]].copy()
    sensitivity["shannon_rarefied"] = sensitivity["ids"].map(shannon)
    sensitivity["richness_rarefied"] = sensitivity["ids"].map(richness)
    sensitivity["simpson_rarefied"] = sensitivity["ids"].map(simpson)
    sensitivity["evenness_rarefied"] = sensitivity["ids"].map(evenness)
    sensitivity.to_csv(output / "methods_faithful_rarefied_sensitivity.csv", index=False)
