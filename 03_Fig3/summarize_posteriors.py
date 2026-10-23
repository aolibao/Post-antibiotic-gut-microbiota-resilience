#!/usr/bin/env python3
"""Summarize MDSINE2 posterior chains for the restored manuscript Figure 3.

This module deliberately keeps two network estimands separate:

* ``module_spearman_edges_<GROUP>.csv`` is an undirected correlation network
  calculated from the 16 subject-aggregated canonical-stage trajectories.  It
  is the source-caption-consistent network used in the centre of Fig. 3E/F.
* ``directed_mdsine_edges_<GROUP>.csv`` is the directed MDSINE2 interaction
  network and contains only module pairs with posterior Bayes factor >= 10.

The qPCR file used to construct the input Study is a structural placeholder.
Consequently, growth and perturbation parameters are reported in the model's
relative-abundance scale and must not be described as absolute biomass rates.

The public functions are NumPy-first so that their numerical behaviour can be
tested without importing MDSINE2.  MDSINE2 is imported lazily only when a chain
pickle is loaded; run that operation in the ``mdsine2`` conda environment.
"""

from __future__ import annotations

import argparse
import json
import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score

from mdsine2_noqpcr_loader import load_noqpcr_mcmc


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ANALYSIS_ROOT = PROJECT_ROOT / "reproduction" / "outputs" / "fig3_mdsine2"
DEFAULT_TABLE_DIR = DEFAULT_ANALYSIS_ROOT / "tables"
CANONICAL_STAGES = [
    "A1", "A2", "A3", "A4",
    "B1", "B2", "B3", "B4", "B5", "B6",
    "C1", "C2", "C3", "C4", "C5", "C6",
]


def bh_fdr(pvalues: Sequence[float] | np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, preserving NaNs and shape."""

    values = np.asarray(pvalues, dtype=float)
    shape = values.shape
    flat = values.ravel()
    valid = np.isfinite(flat)
    adjusted = np.full(flat.shape, np.nan, dtype=float)
    if not np.any(valid):
        return adjusted.reshape(shape)
    p = np.clip(flat[valid], 0.0, 1.0)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1, dtype=float)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(q, 0.0, 1.0)
    adjusted[np.flatnonzero(valid)] = restored
    return adjusted.reshape(shape)


def posterior_bayes_factor(
    on_counts: np.ndarray,
    total_counts: int | np.ndarray,
    prior_a: float,
    prior_b: float,
) -> np.ndarray:
    """Return posterior-odds/prior-odds Bayes factors for spike indicators.

    MDSINE2 parameterizes the inclusion probability with ``Beta(a, b)``.  Its
    post-hoc utility therefore computes ``(n_on * b) / (n_off * a)``.  This
    implementation has explicit behaviour at the boundaries: zero evidence is
    returned when a parameter is never included and infinity when it is always
    included.
    """

    on = np.asarray(on_counts, dtype=float)
    total = np.asarray(total_counts, dtype=float)
    if prior_a <= 0 or prior_b <= 0:
        raise ValueError("Beta prior parameters must both be positive")
    off = total - on
    with np.errstate(divide="ignore", invalid="ignore"):
        bf = (on * float(prior_b)) / (off * float(prior_a))
    bf = np.where(on <= 0, 0.0, bf)
    bf = np.where((off <= 0) & (on > 0), np.inf, bf)
    return np.asarray(bf, dtype=float)


def _quantile_summary(samples: np.ndarray) -> dict[str, np.ndarray]:
    values = np.asarray(samples, dtype=float)
    if values.ndim < 1 or values.shape[0] == 0:
        raise ValueError("Posterior samples must have a non-empty draw axis")
    return {
        "mean": np.mean(values, axis=0),
        "median": np.median(values, axis=0),
        "q025": np.quantile(values, 0.025, axis=0),
        "q975": np.quantile(values, 0.975, axis=0),
    }


def _conditional_quantile_summary(samples_with_nan: np.ndarray) -> dict[str, np.ndarray]:
    values = np.asarray(samples_with_nan, dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return {
            "conditional_median": np.nanmedian(values, axis=0),
            "conditional_q025": np.nanquantile(values, 0.025, axis=0),
            "conditional_q975": np.nanquantile(values, 0.975, axis=0),
        }


def _equal_length_stack(chains: Sequence[np.ndarray]) -> np.ndarray:
    if not chains:
        raise ValueError("At least one chain is required")
    arrays = [np.asarray(x, dtype=float) for x in chains]
    trailing = arrays[0].shape[1:]
    if any(x.ndim < 1 or x.shape[1:] != trailing for x in arrays):
        raise ValueError("All chains must have identical parameter dimensions")
    n = min(x.shape[0] for x in arrays)
    if n < 4:
        raise ValueError("At least four posterior draws per chain are required")
    return np.stack([x[:n] for x in arrays], axis=0)


def _rank_normalize(values: np.ndarray) -> np.ndarray:
    """Rank-normalize across chains/draws independently for each parameter."""

    x = np.asarray(values, dtype=float)
    if x.ndim < 2:
        raise ValueError("Expected chain and draw dimensions")
    original_shape = x.shape
    flat = x.reshape(original_shape[0] * original_shape[1], -1)
    ranks = stats.rankdata(flat, axis=0, method="average")
    n = flat.shape[0]
    probabilities = (ranks - 3.0 / 8.0) / (n + 1.0 / 4.0)
    normalized = stats.norm.ppf(probabilities)
    return normalized.reshape(original_shape)


def _split(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    half = x.shape[1] // 2
    if half < 2:
        raise ValueError("Too few draws to split chains")
    return np.concatenate([x[:, :half], x[:, -half:]], axis=0)


def _basic_rhat(values: np.ndarray) -> np.ndarray:
    x = _split(values)
    n = x.shape[1]
    means = np.mean(x, axis=1)
    between = n * np.var(means, axis=0, ddof=1)
    within = np.mean(np.var(x, axis=1, ddof=1), axis=0)
    variance = ((n - 1.0) / n) * within + between / n
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.sqrt(variance / within)
    both_constant = (within == 0) & (between == 0)
    result = np.where(both_constant, 1.0, result)
    result = np.where((within == 0) & (between > 0), np.inf, result)
    return result


def rank_normalized_split_rhat(chains: Sequence[np.ndarray]) -> np.ndarray:
    """Vehtari-style rank-normalized split-Rhat (rank and folded maximum).

    A single chain cannot identify between-chain non-convergence; in that case
    an array of NaNs is returned rather than a misleading diagnostic.
    """

    x = _equal_length_stack(chains)
    parameter_shape = x.shape[2:]
    if x.shape[0] < 2:
        return np.full(parameter_shape or (), np.nan, dtype=float)
    rank_rhat = _basic_rhat(_rank_normalize(x))
    median = np.median(x, axis=(0, 1), keepdims=True)
    folded_rhat = _basic_rhat(_rank_normalize(np.abs(x - median)))
    return np.maximum(rank_rhat, folded_rhat)


def bulk_effective_sample_size(chains: Sequence[np.ndarray]) -> np.ndarray:
    """Estimate rank-normalized bulk ESS with Geyer's paired sequence rule."""

    x = _rank_normalize(_equal_length_stack(chains))
    m, n = x.shape[:2]
    trailing = x.shape[2:]
    flat = x.reshape(m, n, -1)
    centred = flat - flat.mean(axis=1, keepdims=True)
    fft_length = 1 << (2 * n - 1).bit_length()
    transformed = np.fft.rfft(centred, n=fft_length, axis=1)
    acov = np.fft.irfft(transformed * np.conjugate(transformed), n=fft_length, axis=1)[:, :n]
    acov /= np.arange(n, 0, -1, dtype=float)[None, :, None]

    within = np.mean(acov[:, 0, :], axis=0)
    chain_means = flat.mean(axis=1)
    between = n * np.var(chain_means, axis=0, ddof=1) if m > 1 else np.zeros_like(within)
    var_plus = ((n - 1.0) / n) * within + between / n
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = 1.0 - (within[None, :] - acov.mean(axis=0)) / var_plus[None, :]
    rho[0, :] = 1.0

    ess = np.empty(flat.shape[2], dtype=float)
    for parameter in range(flat.shape[2]):
        if not np.isfinite(var_plus[parameter]) or var_plus[parameter] <= 0:
            ess[parameter] = float(m * n)
            continue
        correlations = rho[:, parameter]
        pair_sums: list[float] = []
        for lag in range(0, n - 1, 2):
            pair = float(correlations[lag] + correlations[lag + 1])
            if not np.isfinite(pair) or pair < 0:
                break
            if pair_sums and pair > pair_sums[-1]:
                pair = pair_sums[-1]
            pair_sums.append(pair)
        tau = -1.0 + 2.0 * sum(pair_sums)
        tau = max(tau, 1.0)
        ess[parameter] = min(float(m * n), float(m * n) / tau)
    return ess.reshape(trailing or ())


def consensus_modules(
    cocluster_probability: np.ndarray,
    n_modules: int,
    linkage: str = "average",
) -> np.ndarray:
    """Cluster ``1 - P(co-cluster)`` and preserve sklearn's raw labels.

    MDSINE2's historical single-chain visualization used average linkage and
    wrote the raw ``AgglomerativeClustering.fit_predict`` labels.  Those labels
    are visually meaningful in the manuscript because they became M0..Mk.
    They must therefore not be recoded by first taxon or cluster size.
    """

    probability = np.asarray(cocluster_probability, dtype=float)
    if probability.ndim != 2 or probability.shape[0] != probability.shape[1]:
        raise ValueError("Co-clustering probability must be a square matrix")
    if not 1 <= int(n_modules) <= probability.shape[0]:
        raise ValueError("n_modules is outside the valid range")
    probability = np.clip((probability + probability.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(probability, 1.0)
    distance = 1.0 - probability
    np.fill_diagonal(distance, 0.0)
    return AgglomerativeClustering(
        n_clusters=int(n_modules), metric="precomputed", linkage=linkage
    ).fit_predict(distance)


def align_labels_to_reference(
    labels: Sequence[int] | np.ndarray,
    reference_labels: Sequence[int] | np.ndarray,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Align candidate cluster IDs to a declared reference by maximum overlap.

    The Hungarian assignment maximizes the number of taxa retaining their
    reference module ID.  When cluster counts differ, unmatched candidate
    clusters receive new integer IDs in raw-label order above the largest
    reference ID.  Cluster membership is never altered; only labels change.
    """

    candidate = np.asarray(labels, dtype=int)
    reference = np.asarray(reference_labels, dtype=int)
    if candidate.ndim != 1 or reference.ndim != 1 or candidate.shape != reference.shape:
        raise ValueError("labels and reference_labels must be equal-length vectors")
    candidate_ids = np.sort(np.unique(candidate))
    reference_ids = np.sort(np.unique(reference))
    overlap = np.zeros((len(candidate_ids), len(reference_ids)), dtype=int)
    for candidate_index, candidate_id in enumerate(candidate_ids):
        for reference_index, reference_id in enumerate(reference_ids):
            overlap[candidate_index, reference_index] = int(
                np.sum((candidate == candidate_id) & (reference == reference_id))
            )
    row_indices, column_indices = linear_sum_assignment(-overlap)
    mapping = {
        int(candidate_ids[row]): int(reference_ids[column])
        for row, column in zip(row_indices, column_indices)
    }
    next_label = int(reference_ids.max() + 1) if len(reference_ids) else 0
    for candidate_id in candidate_ids:
        if int(candidate_id) not in mapping:
            mapping[int(candidate_id)] = next_label
            next_label += 1
    aligned = np.asarray([mapping[int(value)] for value in candidate], dtype=int)
    rows = []
    for candidate_id in candidate_ids:
        aligned_id = mapping[int(candidate_id)]
        reference_overlap = int(
            np.sum((candidate == candidate_id) & (reference == aligned_id))
        )
        rows.append(
            {
                "raw_candidate_label": int(candidate_id),
                "aligned_reference_label": int(aligned_id),
                "overlap_taxa": reference_overlap,
                "candidate_cluster_size": int(np.sum(candidate == candidate_id)),
                "reference_cluster_size": int(np.sum(reference == aligned_id)),
            }
        )
    return aligned, pd.DataFrame(rows)


def _mode_integer(values: np.ndarray) -> int:
    x = np.asarray(values, dtype=int).ravel()
    if x.size == 0:
        raise ValueError("Cannot calculate a mode from no values")
    counts = np.bincount(x)
    return int(np.flatnonzero(counts == counts.max())[0])


def _taxon_display_name(name: str) -> str:
    value = str(name)
    if value.startswith("s__"):
        value = value[3:]
    return value.replace("_", " ")


@dataclass
class ChainData:
    path: Path
    taxa: list[str]
    growth: np.ndarray
    perturbation: np.ndarray
    coclusters: np.ndarray
    n_clusters: np.ndarray
    interactions: np.ndarray
    perturbation_prior_a: float
    perturbation_prior_b: float
    interaction_prior_a: float
    interaction_prior_b: float
    perturbation_name: str
    burnin: int
    n_samples: int


def load_mdsine_chain(path: Path, perturbation_name: str | None = None) -> ChainData:
    """Load all Fig. 3 posterior quantities from an MDSINE2 ``mcmc.pkl``."""

    try:
        from mdsine2.names import STRNAMES
    except ImportError as exc:  # pragma: no cover - environment-dependent
        raise RuntimeError(
            "MDSINE2 is required to read mcmc.pkl. Run with "
            "/home/qqt/.conda/envs/mdsine2/bin/python and set NUMBA_CACHE_DIR."
        ) from exc

    source = Path(path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    # The official no-qPCR runner defines two serialized classes under
    # ``__main__``.  The shared loader reinstalls those names before unpickling.
    mcmc = load_noqpcr_mcmc(source)
    graph = mcmc.graph
    perturbations = list(graph.perturbations or [])
    if not perturbations:
        raise ValueError(f"No perturbation posterior in {source}")
    if perturbation_name is None:
        perturbation = perturbations[0]
    else:
        matches = [p for p in perturbations if p.name == perturbation_name]
        if not matches:
            raise KeyError(f"Perturbation {perturbation_name!r} absent from {source}")
        perturbation = matches[0]
    clustering = graph[STRNAMES.CLUSTERING_OBJ]
    interaction_probability = graph[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB]
    taxa = [str(taxon.name) for taxon in graph.data.taxa]

    return ChainData(
        path=source,
        taxa=taxa,
        growth=np.asarray(
            graph[STRNAMES.GROWTH_VALUE].get_trace_from_disk(section="posterior"),
            dtype=float,
        ),
        perturbation=np.asarray(
            perturbation.get_trace_from_disk(section="posterior"), dtype=float
        ),
        coclusters=np.asarray(
            clustering.coclusters.get_trace_from_disk(section="posterior"),
            dtype=bool,
        ),
        n_clusters=np.asarray(
            clustering.n_clusters.get_trace_from_disk(section="posterior"),
            dtype=int,
        ),
        interactions=np.asarray(
            graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk(section="posterior"),
            dtype=float,
        ),
        perturbation_prior_a=float(perturbation.probability.prior.a.value),
        perturbation_prior_b=float(perturbation.probability.prior.b.value),
        interaction_prior_a=float(interaction_probability.prior.a.value),
        interaction_prior_b=float(interaction_probability.prior.b.value),
        perturbation_name=str(perturbation.name),
        burnin=int(mcmc.burnin),
        n_samples=int(mcmc.n_samples),
    )


def _check_chain_compatibility(chains: Sequence[ChainData], group: str) -> None:
    if not chains:
        raise ValueError(f"No {group} chains supplied")
    taxa = chains[0].taxa
    for chain in chains:
        if chain.taxa != taxa:
            raise ValueError(f"Taxon order differs between {group} chains: {chain.path}")
        n_taxa = len(taxa)
        expected_vector = (n_taxa,)
        expected_matrix = (n_taxa, n_taxa)
        if chain.growth.shape[1:] != expected_vector:
            raise ValueError(f"Unexpected growth shape in {chain.path}: {chain.growth.shape}")
        if chain.perturbation.shape[1:] != expected_vector:
            raise ValueError(f"Unexpected perturbation shape in {chain.path}: {chain.perturbation.shape}")
        if chain.coclusters.shape[1:] != expected_matrix:
            raise ValueError(f"Unexpected co-clustering shape in {chain.path}: {chain.coclusters.shape}")
        if chain.interactions.shape[1:] != expected_matrix:
            raise ValueError(f"Unexpected interaction shape in {chain.path}: {chain.interactions.shape}")


def _common_prior(chains: Sequence[ChainData], prefix: str) -> tuple[float, float]:
    a_values = np.asarray([getattr(chain, f"{prefix}_prior_a") for chain in chains])
    b_values = np.asarray([getattr(chain, f"{prefix}_prior_b") for chain in chains])
    if not (np.allclose(a_values, a_values[0]) and np.allclose(b_values, b_values[0])):
        raise ValueError(f"{prefix} indicator priors differ across chains")
    return float(a_values[0]), float(b_values[0])


def infer_evidence_tier(paths: Sequence[Path]) -> str:
    joined = " ".join(str(path).lower() for path in paths)
    if "fast" in joined:
        return "development_only"
    if "pilot" in joined:
        return "pilot_single_or_short_chain"
    if len(paths) < 2:
        return "single_chain_unverified"
    return "formal_multichain"


def _diagnostic_rows(
    group: str,
    parameter_type: str,
    parameter_names: Sequence[str],
    chains: Sequence[np.ndarray],
    evidence_tier: str,
) -> list[dict[str, Any]]:
    try:
        rhat = np.asarray(rank_normalized_split_rhat(chains)).reshape(-1)
    except ValueError:
        rhat = np.full(len(parameter_names), np.nan)
    try:
        ess = np.asarray(bulk_effective_sample_size(chains)).reshape(-1)
    except ValueError:
        ess = np.full(len(parameter_names), np.nan)
    minimum_draws = min(np.asarray(chain).shape[0] for chain in chains)
    status = "ok" if len(chains) >= 2 else "rhat_unavailable_single_chain"
    return [
        {
            "group": group,
            "parameter_type": parameter_type,
            "parameter": str(name),
            "rhat_rank_folded": float(rhat[index]),
            "ess_bulk": float(ess[index]),
            "n_chains": len(chains),
            "draws_per_chain_used": minimum_draws,
            "diagnostic_status": status,
            "evidence_tier": evidence_tier,
        }
        for index, name in enumerate(parameter_names)
    ]


def summarize_group(
    group: str,
    mcmc_paths: Sequence[Path],
    output_dir: Path,
    *,
    bf_threshold: float = 10.0,
    evidence_tier: str | None = None,
    perturbation_name: str | None = None,
) -> dict[str, Any]:
    """Summarize one HDR/LDR chain collection and write stable table names."""

    group = str(group).upper()
    if group not in {"LDR", "HDR"}:
        raise ValueError("group must be LDR or HDR")
    paths = [Path(path).resolve() for path in mcmc_paths]
    tier = evidence_tier or infer_evidence_tier(paths)
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    chains = [load_mdsine_chain(path, perturbation_name=perturbation_name) for path in paths]
    _check_chain_compatibility(chains, group)
    taxa = chains[0].taxa
    n_taxa = len(taxa)

    # Co-clustering consensus; chains contribute in proportion to retained draws.
    # The source manuscript's single-chain exports used average linkage and the
    # raw sklearn labels (those integer IDs became M0, M1, ... in the figure).
    # Keep that historical convention as the primary/source-restored result.
    # For multiple chains, a pooled partition has arbitrary label IDs, so align
    # them by maximum overlap to an explicitly declared reference: the first
    # CLI chain, clustered with the pooled modal module count and average
    # linkage.  Input CLI order therefore defines module *names*, not membership.
    cocluster_sum = np.zeros((n_taxa, n_taxa), dtype=float)
    cocluster_draws = 0
    for chain in chains:
        cocluster_sum += chain.coclusters.sum(axis=0, dtype=np.float64)
        cocluster_draws += chain.coclusters.shape[0]
    cocluster_probability = cocluster_sum / float(cocluster_draws)
    cocluster_probability = (cocluster_probability + cocluster_probability.T) / 2.0
    np.fill_diagonal(cocluster_probability, 1.0)
    n_cluster_trace = np.concatenate([chain.n_clusters for chain in chains])
    n_modules = _mode_integer(n_cluster_trace)
    raw_average_assignments = consensus_modules(
        cocluster_probability, n_modules=n_modules, linkage="average"
    )
    label_reference_chain = str(chains[0].path)
    if len(chains) == 1:
        # Preserve the exact historical raw agglomerative labels.  Recoding by
        # first taxon or cluster size would change the manuscript's M0..Mk IDs.
        assignments = raw_average_assignments.copy()
        _, average_alignment = align_labels_to_reference(
            raw_average_assignments, raw_average_assignments
        )
        label_alignment_strategy = "raw_sklearn_labels_single_chain"
    else:
        reference_probability = np.asarray(
            chains[0].coclusters.mean(axis=0), dtype=float
        )
        reference_probability = (
            reference_probability + reference_probability.T
        ) / 2.0
        np.fill_diagonal(reference_probability, 1.0)
        reference_assignments = consensus_modules(
            reference_probability, n_modules=n_modules, linkage="average"
        )
        assignments, average_alignment = align_labels_to_reference(
            raw_average_assignments, reference_assignments
        )
        label_alignment_strategy = (
            "hungarian_maximum_overlap_to_first_cli_chain_average_linkage"
        )

    average_alignment.insert(0, "group", group)
    average_alignment["candidate_linkage"] = "average"
    average_alignment["reference_chain"] = label_reference_chain
    average_alignment["label_alignment_strategy"] = label_alignment_strategy
    average_alignment.to_csv(
        output_dir / f"{group}_module_label_alignment.csv", index=False
    )

    # The current MDSINE2 multi-chain post-hoc CLI uses complete linkage.  Keep
    # that convention as a clearly separated sensitivity analysis rather than
    # silently replacing the source-restored average-linkage partition.
    raw_complete_assignments = consensus_modules(
        cocluster_probability, n_modules=n_modules, linkage="complete"
    )
    complete_assignments, complete_alignment = align_labels_to_reference(
        raw_complete_assignments, assignments
    )
    complete_alignment.insert(0, "group", group)
    complete_alignment["candidate_linkage"] = "complete"
    complete_alignment["reference_linkage"] = "average"
    complete_alignment["reference_chain"] = label_reference_chain
    complete_alignment["label_alignment_strategy"] = (
        "hungarian_maximum_overlap_to_primary_average_partition"
    )
    complete_alignment.to_csv(
        output_dir / f"{group}_complete_linkage_label_alignment.csv", index=False
    )
    linkage_ari = float(
        adjusted_rand_score(assignments, raw_complete_assignments)
    )
    complete_matches_primary = complete_assignments == assignments
    linkage_sensitivity = {
        "group": group,
        "primary_linkage": "average",
        "sensitivity_linkage": "complete",
        "n_modules": int(n_modules),
        "n_taxa": int(n_taxa),
        "adjusted_rand_index": linkage_ari,
        "n_taxa_same_after_label_alignment": int(
            np.sum(complete_matches_primary)
        ),
        "n_taxa_changed_after_label_alignment": int(
            np.sum(~complete_matches_primary)
        ),
        "fraction_taxa_same_after_label_alignment": float(
            np.mean(complete_matches_primary)
        ),
        "primary_label_alignment_strategy": label_alignment_strategy,
        "primary_label_reference_chain": label_reference_chain,
        "evidence_tier": tier,
    }
    module_names = np.asarray([f"M{value}" for value in assignments], dtype=object)

    cocluster_frame = pd.DataFrame(cocluster_probability, index=taxa, columns=taxa)
    cocluster_frame.index.name = "taxon"
    cocluster_frame.to_csv(output_dir / f"{group}_cocluster_consensus.csv")
    within_probability = np.empty(n_taxa, dtype=float)
    outside_probability = np.empty(n_taxa, dtype=float)
    for taxon_index in range(n_taxa):
        same = np.flatnonzero(assignments == assignments[taxon_index])
        same = same[same != taxon_index]
        other = np.flatnonzero(assignments != assignments[taxon_index])
        within_probability[taxon_index] = (
            float(np.mean(cocluster_probability[taxon_index, same]))
            if len(same)
            else 1.0
        )
        outside_probability[taxon_index] = (
            float(np.max(cocluster_probability[taxon_index, other]))
            if len(other)
            else 0.0
        )
    assignment_frame = pd.DataFrame(
        {
            "taxon": taxa,
            "module": module_names,
            "taxon_index": np.arange(n_taxa, dtype=int),
            "module_index": assignments,
            "display_name": [_taxon_display_name(x) for x in taxa],
            "within_module_cocluster_mean": within_probability,
            "max_outside_module_cocluster": outside_probability,
            "assignment_probability_margin": (
                within_probability - outside_probability
            ),
            "consensus_linkage": "average",
            "label_alignment_strategy": label_alignment_strategy,
            "label_reference_chain": label_reference_chain,
        }
    )
    assignment_frame.to_csv(output_dir / f"{group}_module_assignments.csv", index=False)

    complete_assignment_frame = pd.DataFrame(
        {
            "taxon": taxa,
            "taxon_index": np.arange(n_taxa, dtype=int),
            "primary_module": module_names,
            "primary_module_index": assignments,
            "complete_module": np.asarray(
                [f"M{value}" for value in complete_assignments], dtype=object
            ),
            "complete_module_index_aligned": complete_assignments,
            "complete_module_index_raw": raw_complete_assignments,
            "same_as_primary_after_label_alignment": complete_matches_primary,
            "primary_linkage": "average",
            "sensitivity_linkage": "complete",
            "adjusted_rand_index_partition_level": linkage_ari,
            "evidence_tier": tier,
        }
    )
    complete_assignment_frame.to_csv(
        output_dir / f"{group}_module_assignments_complete_linkage_sensitivity.csv",
        index=False,
    )
    clustering_chain_records = []
    for chain_index, chain in enumerate(chains):
        chain_probability = chain.coclusters.mean(axis=0)
        difference = chain_probability - cocluster_probability
        clustering_chain_records.append(
            {
                "group": group,
                "chain_index": chain_index,
                "path": str(chain.path),
                "posterior_draws": chain.coclusters.shape[0],
                "n_modules_mode": _mode_integer(chain.n_clusters),
                "n_modules_median": float(np.median(chain.n_clusters)),
                "cocluster_RMSE_vs_pooled": float(np.sqrt(np.mean(difference**2))),
                "cocluster_max_abs_difference_vs_pooled": float(
                    np.max(np.abs(difference))
                ),
                "evidence_tier": tier,
            }
        )
    pd.DataFrame(clustering_chain_records).to_csv(
        output_dir / f"clustering_chain_diagnostics_{group}.csv", index=False
    )

    # Perturbation spike-and-slab marginal summaries (excluded draws are zero).
    pert_raw = np.concatenate([chain.perturbation for chain in chains], axis=0)
    pert_on = np.sum(np.isfinite(pert_raw), axis=0)
    pert_total = pert_raw.shape[0]
    pert_a, pert_b = _common_prior(chains, "perturbation")
    pert_bf = posterior_bayes_factor(pert_on, pert_total, pert_a, pert_b)
    pert_marginal = np.nan_to_num(pert_raw, nan=0.0)
    pert_summary = _quantile_summary(pert_marginal)
    pert_summary.update(_conditional_quantile_summary(pert_raw))
    perturbation_frame = pd.DataFrame(
        {
            "taxon": taxa,
            "display_name": [_taxon_display_name(x) for x in taxa],
            "module": module_names,
            "median": pert_summary["median"],
            "mean": pert_summary["mean"],
            "q025": pert_summary["q025"],
            "q975": pert_summary["q975"],
            "conditional_median": pert_summary["conditional_median"],
            "conditional_q025": pert_summary["conditional_q025"],
            "conditional_q975": pert_summary["conditional_q975"],
            "inclusion_probability": pert_on / float(pert_total),
            "bayes_factor": pert_bf,
            "bf_ge_threshold": pert_bf >= bf_threshold,
            "posterior_draws": pert_total,
            "perturbation_name": chains[0].perturbation_name,
            "evidence_tier": tier,
            "scale_note": "relative-abundance gLV; no measured qPCR",
        }
    )
    perturbation_frame.to_csv(output_dir / f"perturbation_effects_{group}.csv", index=False)

    growth_raw = np.concatenate([chain.growth for chain in chains], axis=0)
    growth_summary = _quantile_summary(growth_raw)
    growth_frame = pd.DataFrame(
        {
            "taxon": taxa,
            "display_name": [_taxon_display_name(x) for x in taxa],
            "module": module_names,
            "median": growth_summary["median"],
            "mean": growth_summary["mean"],
            "q025": growth_summary["q025"],
            "q975": growth_summary["q975"],
            "posterior_draws": growth_raw.shape[0],
            "evidence_tier": tier,
            "scale_note": "relative-abundance gLV; no measured qPCR",
        }
    )
    growth_frame.to_csv(output_dir / f"growth_effects_{group}.csv", index=False)

    # Directed item-level evidence matrix.  Rows are targets and columns are
    # sources in MDSINE2's interaction matrix convention.
    interaction_on = np.zeros((n_taxa, n_taxa), dtype=np.int64)
    interaction_total = 0
    for chain in chains:
        interaction_on += np.isfinite(chain.interactions).sum(axis=0)
        interaction_total += chain.interactions.shape[0]
    int_a, int_b = _common_prior(chains, "interaction")
    interaction_bf = posterior_bayes_factor(
        interaction_on, interaction_total, int_a, int_b
    )
    np.fill_diagonal(interaction_bf, 0.0)

    selected_target, selected_source = np.where(interaction_bf >= bf_threshold)
    if len(selected_target):
        item_edge_chains = [
            np.nan_to_num(
                chain.interactions[:, selected_target, selected_source], nan=0.0
            )
            for chain in chains
        ]
        item_edge_all = np.concatenate(item_edge_chains, axis=0)
        item_edge_summary = _quantile_summary(item_edge_all)
        taxon_edge_frame = pd.DataFrame(
            {
                "source": np.asarray(taxa, dtype=object)[selected_source],
                "target": np.asarray(taxa, dtype=object)[selected_target],
                "source_module": module_names[selected_source],
                "target_module": module_names[selected_target],
                "effect_median": item_edge_summary["median"],
                "effect_q025": item_edge_summary["q025"],
                "effect_q975": item_edge_summary["q975"],
                "inclusion_probability": interaction_on[
                    selected_target, selected_source
                ] / float(interaction_total),
                "bayes_factor": interaction_bf[selected_target, selected_source],
                "sign": np.where(item_edge_summary["median"] >= 0, "positive", "negative"),
                "orientation": "source_column_to_target_row",
                "evidence_tier": tier,
            }
        ).sort_values("bayes_factor", ascending=False, kind="mergesort")
    else:
        item_edge_chains = []
        taxon_edge_frame = pd.DataFrame(
            columns=[
                "source", "target", "source_module", "target_module",
                "effect_median", "effect_q025", "effect_q975",
                "inclusion_probability", "bayes_factor", "sign",
                "orientation", "evidence_tier",
            ]
        )
    taxon_edge_frame.to_csv(
        output_dir / f"directed_mdsine_taxon_edges_{group}.csv", index=False
    )

    # Collapse evidence to consensus modules without treating all taxon pairs as
    # independent observations.  With learned (rather than fixed) clustering,
    # item-wise post-hoc BFs need not be constant throughout a consensus block.
    # We therefore use the single most-supported item pair as the representative
    # edge, mirroring MDSINE2's fixed-cluster condensation by a representative
    # item.  The block median BF is retained as a sensitivity field.
    module_members = {
        module: np.flatnonzero(assignments == module) for module in range(n_modules)
    }
    module_records: list[dict[str, Any]] = []
    module_edge_chain_values: list[list[np.ndarray]] = []
    module_edge_names: list[str] = []
    for target_module in range(n_modules):
        targets = module_members[target_module]
        for source_module in range(n_modules):
            if source_module == target_module:
                continue
            sources = module_members[source_module]
            block_bf = interaction_bf[np.ix_(targets, sources)]
            module_bf = float(np.nanmax(block_bf))
            if not (module_bf >= bf_threshold):
                continue
            representative_flat = int(np.nanargmax(block_bf))
            representative_target_offset, representative_source_offset = np.unravel_index(
                representative_flat, block_bf.shape
            )
            representative_target = int(targets[representative_target_offset])
            representative_source = int(sources[representative_source_offset])
            per_chain = []
            for chain in chains:
                per_chain.append(
                    np.nan_to_num(
                        chain.interactions[
                            :, representative_target, representative_source
                        ],
                        nan=0.0,
                    )
                )
            pooled = np.concatenate(per_chain)
            summary = _quantile_summary(pooled)
            finite_fraction = float(
                interaction_on[representative_target, representative_source]
                / float(interaction_total)
            )
            source_name = f"M{source_module}"
            target_name = f"M{target_module}"
            module_records.append(
                {
                    "source": source_name,
                    "target": target_name,
                    "effect_median": float(summary["median"]),
                    "effect_q025": float(summary["q025"]),
                    "effect_q975": float(summary["q975"]),
                    "inclusion_probability": finite_fraction,
                    "bayes_factor": module_bf,
                    "median_item_bayes_factor": float(np.nanmedian(block_bf)),
                    "representative_source_taxon": taxa[representative_source],
                    "representative_target_taxon": taxa[representative_target],
                    "module_evidence_aggregation": "maximum_item_bf_representative_pair",
                    "n_source_taxa": len(sources),
                    "n_target_taxa": len(targets),
                    "sign": "positive" if summary["median"] >= 0 else "negative",
                    "orientation": "directed_source_to_target",
                    "bf_threshold": bf_threshold,
                    "evidence_tier": tier,
                }
            )
            module_edge_chain_values.append(per_chain)
            module_edge_names.append(f"{source_name}->{target_name}")
    module_edge_columns = [
        "source", "target", "effect_median", "effect_q025", "effect_q975",
        "inclusion_probability", "bayes_factor",
        "median_item_bayes_factor", "representative_source_taxon",
        "representative_target_taxon", "module_evidence_aggregation",
        "n_source_taxa", "n_target_taxa", "sign", "orientation",
        "bf_threshold", "evidence_tier",
    ]
    directed_module_frame = pd.DataFrame(
        module_records, columns=module_edge_columns
    )
    if not directed_module_frame.empty:
        directed_module_frame = directed_module_frame.sort_values(
            "bayes_factor", ascending=False, kind="mergesort"
        )
    directed_module_frame.to_csv(
        output_dir / f"directed_mdsine_edges_{group}.csv", index=False
    )

    diagnostics: list[dict[str, Any]] = []
    diagnostics.extend(
        _diagnostic_rows(
            group, "growth", taxa, [chain.growth for chain in chains], tier
        )
    )
    diagnostics.extend(
        _diagnostic_rows(
            group,
            "perturbation_marginal",
            taxa,
            [np.nan_to_num(chain.perturbation, nan=0.0) for chain in chains],
            tier,
        )
    )
    diagnostics.extend(
        _diagnostic_rows(
            group,
            "n_modules",
            ["n_modules"],
            [chain.n_clusters[:, None] for chain in chains],
            tier,
        )
    )
    for edge_name, edge_chains in zip(module_edge_names, module_edge_chain_values):
        diagnostics.extend(
            _diagnostic_rows(
                group,
                "directed_module_interaction_marginal",
                [edge_name],
                [values[:, None] for values in edge_chains],
                tier,
            )
        )

    return {
        "group": group,
        "taxa": taxa,
        "assignments": assignment_frame,
        "perturbation": perturbation_frame,
        "growth": growth_frame,
        "diagnostics": pd.DataFrame(diagnostics),
        "n_modules": n_modules,
        "n_cluster_trace": n_cluster_trace,
        "cocluster_probability": cocluster_probability,
        "directed_module_edges": directed_module_frame,
        "directed_taxon_edges": taxon_edge_frame,
        "linkage_sensitivity": linkage_sensitivity,
        "label_alignment_strategy": label_alignment_strategy,
        "label_reference_chain": label_reference_chain,
        "evidence_tier": tier,
        "paths": [str(path) for path in paths],
        "chain_draws": [int(chain.growth.shape[0]) for chain in chains],
        "burnin": [chain.burnin for chain in chains],
        "n_samples": [chain.n_samples for chain in chains],
    }


def _read_table(path: Path, *, index_col: int | None = None) -> pd.DataFrame:
    source = Path(path)
    if source.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(source, sep="\t", index_col=index_col)
    return pd.read_csv(source, sep=None, engine="python", index_col=index_col)


def load_abundance_inputs(
    reads_path: Path,
    metadata_path: Path,
    group: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load counts as taxa x samples and standardized sample metadata.

    The function detects accidental transposition by comparing metadata sample
    identifiers to both axes.  Counts are converted to proportions after the
    orientation check; zero-depth samples are rejected.
    """

    reads = _read_table(Path(reads_path), index_col=0)
    metadata = _read_table(Path(metadata_path), index_col=None)
    sample_column = next(
        (
            column
            for column in ["sampleID", "sample_id", "sample", "ids", "library"]
            if column in metadata.columns
        ),
        None,
    )
    if sample_column is None:
        raise ValueError(f"No sample identifier column found in {metadata_path}")
    samples = metadata[sample_column].astype(str)
    reads.index = reads.index.astype(str)
    reads.columns = reads.columns.astype(str)
    column_overlap = len(set(samples) & set(reads.columns))
    index_overlap = len(set(samples) & set(reads.index))
    if index_overlap > column_overlap:
        reads = reads.T
        reads.index = reads.index.astype(str)
        reads.columns = reads.columns.astype(str)
        column_overlap = len(set(samples) & set(reads.columns))
    if column_overlap != len(samples):
        missing = sorted(set(samples) - set(reads.columns))
        raise ValueError(f"Reads/metadata sample mismatch; missing: {missing[:10]}")
    if samples.duplicated().any():
        raise ValueError("Metadata contains duplicate sample identifiers")
    counts = reads.loc[:, samples].apply(pd.to_numeric, errors="raise")
    if (counts.to_numpy() < 0).any():
        raise ValueError("Reads table contains negative counts")
    depths = counts.sum(axis=0)
    if (depths <= 0).any():
        raise ValueError("Reads table contains zero-depth samples")
    abundance = counts.div(depths, axis=1)

    standardized = metadata.copy()
    standardized["sampleID"] = samples.to_numpy()
    subject_column = next(
        (column for column in ["subject", "subject_id", "participant"] if column in metadata.columns),
        None,
    )
    canon_column = next(
        (column for column in ["canon", "canonical_stage", "stage"] if column in metadata.columns),
        None,
    )
    time_column = next(
        (column for column in ["time", "day", "actual_day"] if column in metadata.columns),
        None,
    )
    if subject_column is None or canon_column is None:
        raise ValueError("Metadata must contain subject and canonical-stage columns")
    standardized["subject"] = metadata[subject_column].astype(str)
    standardized["canon"] = metadata[canon_column].astype(str)
    standardized["time"] = (
        pd.to_numeric(metadata[time_column], errors="coerce")
        if time_column is not None
        else np.nan
    )
    if group is not None:
        standardized["recovery_group"] = str(group).upper()
    elif "recovery_group" not in standardized:
        raise ValueError("group must be supplied when metadata lacks recovery_group")
    standardized = standardized.set_index("sampleID", drop=False).loc[samples]
    return abundance, standardized


def module_spearman_network(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    assignments: pd.DataFrame,
    group: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate an undirected module network from 16 canonical-stage means.

    Samples are first averaged within subject and canonical stage.  Module
    abundance is then averaged across subjects, leaving exactly one observation
    per canonical stage for each module.  The Spearman unit is therefore the
    stage, not hundreds of repeated samples from the same participants.
    """

    required_assignment = {"taxon", "module"}
    if not required_assignment.issubset(assignments.columns):
        raise ValueError("assignments must contain taxon and module columns")
    mapping = assignments.set_index("taxon")["module"].astype(str)
    shared_taxa = [taxon for taxon in mapping.index if taxon in abundance.index]
    if not shared_taxa:
        raise ValueError("No assigned taxa match the abundance table")
    missing_stages = sorted(set(CANONICAL_STAGES) - set(metadata["canon"].astype(str)))
    if missing_stages:
        raise ValueError(f"Missing canonical stages for {group}: {missing_stages}")
    module_abundance = abundance.loc[shared_taxa].groupby(mapping.loc[shared_taxa], sort=False).sum().T
    module_abundance.index.name = "sampleID"
    joined = module_abundance.join(metadata[["subject", "canon"]], how="inner")
    subject_stage = joined.groupby(["subject", "canon"], observed=True)[module_abundance.columns].mean()
    stage = subject_stage.groupby("canon", observed=True).mean().reindex(CANONICAL_STAGES)

    trajectory = (
        stage.reset_index()
        .melt(id_vars="canon", var_name="module", value_name="mean_relative_abundance")
    )
    trajectory["group"] = str(group).upper()
    trajectory["stage_index"] = trajectory["canon"].map(
        {stage_name: index for index, stage_name in enumerate(CANONICAL_STAGES)}
    )
    trajectory["aggregation"] = "subject-stage mean, then subject mean"

    edge_columns = [
        "source",
        "target",
        "rho",
        "p",
        "n_stages",
        "group",
        "directed",
        "observation_unit",
    ]
    records: list[dict[str, Any]] = []
    modules = list(stage.columns.astype(str))
    for left_index in range(len(modules)):
        for right_index in range(left_index + 1, len(modules)):
            source = modules[left_index]
            target = modules[right_index]
            result = stats.spearmanr(stage[source], stage[target], nan_policy="omit")
            records.append(
                {
                    "source": source,
                    "target": target,
                    "rho": float(result.statistic),
                    "p": float(result.pvalue),
                    "n_stages": int(stage[[source, target]].dropna().shape[0]),
                    "group": str(group).upper(),
                    "directed": False,
                    "observation_unit": "canonical_stage",
                }
            )
    # Preserve a stable schema even if a short developmental chain yields a
    # single consensus module.  Without explicit columns, ``to_csv`` writes a
    # headerless one-byte file; the renderer then fails in ``pd.read_csv`` and
    # the cross-definition R audit cannot distinguish a valid empty network
    # from a corrupt export.
    edges = pd.DataFrame(records, columns=edge_columns)
    if not edges.empty:
        edges["q"] = bh_fdr(edges["p"].to_numpy())
        edges["selected_q05"] = edges["q"] < 0.05
        edges["selected_nominal_p05"] = edges["p"] < 0.05
        edges["sign"] = np.where(edges["rho"] >= 0, "positive", "negative")
    else:
        edges = edges.assign(
            q=pd.Series(dtype=float),
            selected_q05=pd.Series(dtype=bool),
            selected_nominal_p05=pd.Series(dtype=bool),
            sign=pd.Series(dtype=str),
        )
    return edges, trajectory


def _merge_group_effects(
    ldr: pd.DataFrame,
    hdr: pd.DataFrame,
    *,
    include_bf: bool,
) -> pd.DataFrame:
    keys = ["taxon"]
    ldr_columns = ["taxon", "display_name", "median", "q025", "q975"]
    hdr_columns = ["taxon", "display_name", "median", "q025", "q975"]
    if include_bf:
        ldr_columns += ["inclusion_probability", "bayes_factor"]
        hdr_columns += ["inclusion_probability", "bayes_factor"]
    left = ldr[ldr_columns].rename(
        columns={column: f"LDR_{column}" for column in ldr_columns if column != "taxon"}
    )
    right = hdr[hdr_columns].rename(
        columns={column: f"HDR_{column}" for column in hdr_columns if column != "taxon"}
    )
    combined = left.merge(right, on=keys, how="outer", validate="one_to_one")
    combined.insert(1, "display_name", combined["LDR_display_name"].fillna(combined["HDR_display_name"]))
    combined = combined.drop(columns=["LDR_display_name", "HDR_display_name"])
    if include_bf:
        # Short aliases are part of the renderer-facing stable schema; retain
        # the descriptive columns as well for human-readable audit tables.
        combined["LDR_BF"] = combined["LDR_bayes_factor"]
        combined["HDR_BF"] = combined["HDR_bayes_factor"]
    combined["HDR_minus_LDR_median"] = combined["HDR_median"] - combined["LDR_median"]
    return combined


def summarize_all(
    ldr_paths: Sequence[Path],
    hdr_paths: Sequence[Path],
    output_dir: Path = DEFAULT_TABLE_DIR,
    *,
    ldr_reads: Path | None = None,
    ldr_metadata: Path | None = None,
    hdr_reads: Path | None = None,
    hdr_metadata: Path | None = None,
    bf_threshold: float = 10.0,
    evidence_tier: str | None = None,
) -> dict[str, Any]:
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    ldr = summarize_group(
        "LDR", ldr_paths, output_dir,
        bf_threshold=bf_threshold, evidence_tier=evidence_tier,
    )
    hdr = summarize_group(
        "HDR", hdr_paths, output_dir,
        bf_threshold=bf_threshold, evidence_tier=evidence_tier,
    )
    if set(ldr["taxa"]) != set(hdr["taxa"]):
        raise ValueError("HDR and LDR fits must use the identical fixed taxon set")
    taxon_order_status = (
        "identical" if ldr["taxa"] == hdr["taxa"] else "same_set_different_order"
    )

    perturbation_combined = _merge_group_effects(
        ldr["perturbation"], hdr["perturbation"], include_bf=True
    )
    growth_combined = _merge_group_effects(
        ldr["growth"], hdr["growth"], include_bf=False
    )
    perturbation_combined.to_csv(
        output_dir / "perturbation_effects_combined.csv", index=False
    )
    growth_combined.to_csv(output_dir / "growth_effects_combined.csv", index=False)
    diagnostics = pd.concat([ldr["diagnostics"], hdr["diagnostics"]], ignore_index=True)
    diagnostics.to_csv(output_dir / "posterior_diagnostics.csv", index=False)
    linkage_sensitivity = pd.DataFrame(
        [ldr["linkage_sensitivity"], hdr["linkage_sensitivity"]]
    )
    linkage_sensitivity.to_csv(
        output_dir / "clustering_linkage_sensitivity.csv", index=False
    )

    spearman_status: dict[str, str] = {}
    input_sets = {
        "LDR": (ldr_reads, ldr_metadata, ldr["assignments"]),
        "HDR": (hdr_reads, hdr_metadata, hdr["assignments"]),
    }
    for group, (reads_path, metadata_path, assignments) in input_sets.items():
        if reads_path is None or metadata_path is None:
            spearman_status[group] = "not_requested"
            continue
        abundance, metadata = load_abundance_inputs(reads_path, metadata_path, group)
        edges, trajectory = module_spearman_network(
            abundance, metadata, assignments, group
        )
        edges.to_csv(output_dir / f"module_spearman_edges_{group}.csv", index=False)
        trajectory.to_csv(
            output_dir / f"module_stage_trajectories_{group}.csv", index=False
        )
        spearman_status[group] = "written"

    manifest = {
        "schema_version": 1,
        "groups": {
            "LDR": {
                "chains": ldr["paths"],
                "posterior_draws_by_chain": ldr["chain_draws"],
                "n_modules_mode": ldr["n_modules"],
                "primary_consensus_linkage": "average",
                "module_label_alignment_strategy": ldr["label_alignment_strategy"],
                "module_label_reference_chain": ldr["label_reference_chain"],
                "complete_linkage_sensitivity_ari": ldr["linkage_sensitivity"][
                    "adjusted_rand_index"
                ],
                "evidence_tier": ldr["evidence_tier"],
            },
            "HDR": {
                "chains": hdr["paths"],
                "posterior_draws_by_chain": hdr["chain_draws"],
                "n_modules_mode": hdr["n_modules"],
                "primary_consensus_linkage": "average",
                "module_label_alignment_strategy": hdr["label_alignment_strategy"],
                "module_label_reference_chain": hdr["label_reference_chain"],
                "complete_linkage_sensitivity_ari": hdr["linkage_sensitivity"][
                    "adjusted_rand_index"
                ],
                "evidence_tier": hdr["evidence_tier"],
            },
        },
        "bayes_factor_threshold": bf_threshold,
        "posterior_interval": [0.025, 0.975],
        "cross_group_taxon_order": taxon_order_status,
        "effect_summary": (
            "Spike-and-slab marginal: indicator-off draws are zero; conditional summaries "
            "are exported separately in group-specific tables."
        ),
        "clustering_consensus": (
            "Primary/source-restored assignments use average linkage on 1 minus the "
            "posterior co-clustering matrix. Single-chain summaries preserve raw "
            "sklearn agglomerative labels so historical M0..Mk identifiers remain "
            "unchanged. For multiple chains, pooled labels are aligned by a "
            "maximum-overlap Hungarian assignment to the first CLI chain clustered "
            "with average linkage and the pooled modal module count; CLI order affects "
            "only module IDs, not pooled partition membership."
        ),
        "clustering_sensitivity": (
            "Complete-linkage assignments (the current MDSINE2 multi-chain post-hoc "
            "CLI convention) are exported separately with adjusted Rand index and "
            "taxon-level agreement against the primary average-linkage partition."
        ),
        "interaction_orientation": "MDSINE2 matrix row=target, column=source",
        "main_figure_network": (
            "Undirected Spearman correlations among 16 subject-aggregated canonical-stage "
            "module trajectories; full p and BH q values exported."
        ),
        "separate_model_network": "Directed MDSINE2 module interactions with BF >= threshold",
        "directed_module_edge_aggregation": (
            "For learned clustering, each consensus module pair uses the maximum-BF "
            "taxon pair as its representative; the block median BF and representative "
            "taxon identifiers are also exported."
        ),
        "spearman_status": spearman_status,
        "qPCR_interpretation": (
            "No measured qPCR. Structural placeholders were overridden by the no-qPCR "
            "runner; parameters must not be interpreted as absolute biomass rates."
        ),
    }
    (output_dir / "posterior_summary_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ldr-mcmc", type=Path, action="append", required=True)
    parser.add_argument("--hdr-mcmc", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_TABLE_DIR)
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
    parser.add_argument("--bf-threshold", type=float, default=10.0)
    parser.add_argument(
        "--evidence-tier",
        choices=[
            "development_only", "pilot_single_or_short_chain",
            "single_chain_unverified", "formal_multichain",
        ],
        default=None,
    )
    parser.add_argument(
        "--without-spearman", action="store_true",
        help="Skip abundance-derived Spearman tables.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    inputs = {
        "ldr_reads": None if args.without_spearman else args.ldr_reads,
        "ldr_metadata": None if args.without_spearman else args.ldr_metadata,
        "hdr_reads": None if args.without_spearman else args.hdr_reads,
        "hdr_metadata": None if args.without_spearman else args.hdr_metadata,
    }
    manifest = summarize_all(
        args.ldr_mcmc,
        args.hdr_mcmc,
        args.output_dir,
        bf_threshold=args.bf_threshold,
        evidence_tier=args.evidence_tier,
        **inputs,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
