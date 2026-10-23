#!/usr/bin/env python3
"""Audit plausible raw-Bracken reconstructions of Supplementary Figure 3.

The source caption defines a positive hit as a within-group, within-timepoint
Spearman correlation with ``P < 0.05`` and ``rho > 0``.  What it does not say
is whether abundance entered the correlation as a count or a proportion, or
whether Shannon came from the recovered Bracken table or the supplied legacy
``species.csv`` table.  This script evaluates that complete 2 x 2 response
grid.  It also checks path-level versus terminal-name-collapsed taxa and, for
genus, a species-descendant aggregation; modes are deliberately retained even
when they prove identical, because that is an auditable result rather than an
assumed preprocessing detail.

Outputs are ordinary CSV files so the source labels, every seven-way phase
intersection, and all 16 degree bins can be inspected without the figure.
The script does not modify any manuscript figure.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.stats import t as student_t


TIME_ORDER = [
    "A1", "A2", "A3", "A4",
    "B1", "B2", "B3", "B4", "B5", "B6",
    "C1", "C2", "C3", "C4", "C5", "C6",
]
ENDPOINT = ["C4", "C5", "C6"]
PHASE_COLUMNS = {
    "A": np.arange(0, 4),
    "B": np.arange(4, 10),
    "C": np.arange(10, 16),
}
PATTERN_ORDER = ["A", "B", "C", "AB", "AC", "BC", "ABC"]

# Values printed above the source bars.  Patterns were read from the dot
# matrix, not inferred from the descending bar order.
SOURCE_UPSET = {
    ("genus", "LDR"): {
        "A": 70, "B": 309, "C": 41, "AB": 325,
        "AC": 39, "BC": 231, "ABC": 862,
    },
    ("genus", "HDR"): {
        "A": 124, "B": 107, "C": 27, "AB": 297,
        "AC": 46, "BC": 132, "ABC": 1294,
    },
    ("species", "LDR"): {
        "A": 598, "B": 1715, "C": 316, "AB": 1149,
        "AC": 248, "BC": 653, "ABC": 1125,
    },
    ("species", "HDR"): {
        "A": 424, "B": 1336, "C": 246, "AB": 1230,
        "AC": 109, "BC": 1406, "ABC": 3211,
    },
}

# Panels C/D do not print values.  These estimates were digitized from the
# embedded 1341 x 1745 JPEG and constrained only by the exact identity
# ``sum(degree bins) == sum(the seven source UpSet intersections)``.  At the
# source resolution, one vertical pixel is about 0.82 genus or 6.5 species;
# the uncertainty column in the export prevents these from being presented as
# author-supplied exact counts.  The genus solution is effectively pixel-
# resolved; species values remain approximate.
SOURCE_DEGREE_ESTIMATES = {
    ("genus", "LDR"): [261, 262, 269, 224, 207, 165, 121, 120, 71, 71, 44, 26, 25, 11, 0, 0],
    ("genus", "HDR"): [161, 165, 171, 156, 147, 135, 129, 132, 144, 127, 143, 119, 114, 103, 81, 0],
    ("species", "LDR"): [2130, 1524, 1023, 593, 345, 137, 39, 13, 0, 0, 0, 0, 0, 0, 0, 0],
    ("species", "HDR"): [1328, 1035, 852, 774, 722, 670, 683, 644, 540, 332, 195, 130, 39, 13, 5, 0],
}
SOURCE_DEGREE_UNCERTAINTY = {"genus": 1, "species": 7}


def rank_mask(index: pd.Index, prefix: str) -> np.ndarray:
    paths = index.to_series(index=index)
    terminal = paths.str.rsplit("|", n=1).str[-1]
    return (terminal.str.startswith(prefix) & paths.str.startswith("k__Bacteria|")).to_numpy()


def terminal_name(path: str, prefix: str) -> str:
    value = path.rsplit("|", 1)[-1].removeprefix(prefix).replace("_", " ").strip()
    return value or "Unclassified"


def genus_name_from_species_path(path: str) -> str:
    genera = [part for part in path.split("|") if part.startswith("g__")]
    if not genera:
        return "Unclassified"
    return terminal_name(genera[-1], "g__")


def shannon_from_counts(counts: pd.DataFrame) -> pd.Series:
    proportions = counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1)
    values = proportions.to_numpy(dtype=float)
    log_values = np.zeros_like(values)
    np.log(values, out=log_values, where=values > 0)
    return pd.Series(-(values * log_values).sum(axis=0), index=counts.columns)


def recovery_groups(meta: pd.DataFrame, legacy_shannon: pd.Series) -> dict[str, str]:
    scored = meta.copy()
    scored["shannon"] = scored["ids"].map(legacy_shannon)
    endpoint = (
        scored.loc[
            scored["abx"].eq("Y") & scored["canon"].isin(ENDPOINT)
        ]
        .groupby("subject")["shannon"]
        .mean()
        .sort_values(kind="mergesort")
    )
    result = {str(subject): "LDR" for subject in endpoint.index[:15]}
    result.update({str(subject): "HDR" for subject in endpoint.index[-15:]})
    return result


def feature_spearman(matrix: np.ndarray, response: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized scipy-style Spearman rho and two-sided t-approximation P."""
    ranked_x = rankdata(matrix, axis=1)
    ranked_y = rankdata(response)
    centered_x = ranked_x - ranked_x.mean(axis=1, keepdims=True)
    centered_y = ranked_y - ranked_y.mean()
    denominator = np.sqrt(
        np.sum(centered_x * centered_x, axis=1)
        * np.sum(centered_y * centered_y)
    )
    rho = np.divide(
        centered_x @ centered_y,
        denominator,
        out=np.full(matrix.shape[0], np.nan),
        where=denominator > 0,
    )
    df = response.size - 2
    t_value = rho * np.sqrt(
        np.divide(df, np.maximum(1e-15, 1.0 - rho * rho))
    )
    p_value = 2.0 * student_t.sf(np.abs(t_value), df)
    p_value[~np.isfinite(rho)] = np.nan
    return rho, p_value


def significant_matrix(
    feature_table: pd.DataFrame,
    response: pd.Series,
    meta: pd.DataFrame,
    groups: dict[str, str],
    group: str,
) -> np.ndarray:
    subjects = sorted(subject for subject, value in groups.items() if value == group)
    significant = np.zeros((feature_table.shape[0], len(TIME_ORDER)), dtype=bool)
    for column, canon in enumerate(TIME_ORDER):
        sample_meta = meta.loc[
            meta["subject"].isin(subjects) & meta["canon"].eq(canon)
        ].sort_values("subject")
        if sample_meta.shape[0] != 15:
            raise ValueError(f"Expected 15 {group} samples at {canon}; found {sample_meta.shape[0]}")
        matrix = feature_table[sample_meta["ids"]].to_numpy(dtype=float)
        y = sample_meta["ids"].map(response).to_numpy(dtype=float)
        rho, p_value = feature_spearman(matrix, y)
        significant[:, column] = np.isfinite(p_value) & (p_value < 0.05) & (rho > 0)
    return significant


def pattern_counts(significant: np.ndarray) -> dict[str, int]:
    phase_hits = {
        phase: significant[:, columns].any(axis=1)
        for phase, columns in PHASE_COLUMNS.items()
    }
    counts: dict[str, int] = {}
    for pattern in PATTERN_ORDER:
        selected = np.ones(significant.shape[0], dtype=bool)
        for phase in "ABC":
            selected &= phase_hits[phase] if phase in pattern else ~phase_hits[phase]
        counts[pattern] = int(selected.sum())
    return counts


def exact_timepoint_patterns(significant: np.ndarray) -> pd.DataFrame:
    degree = significant.sum(axis=1)
    retained = significant[degree > 0]
    if retained.shape[0] == 0:
        return pd.DataFrame(columns=["pattern", "degree", "n"])
    unique, counts = np.unique(retained, axis=0, return_counts=True)
    rows = []
    for flags, count in zip(unique, counts):
        pattern = "|".join(
            canon for canon, selected in zip(TIME_ORDER, flags) if selected
        )
        rows.append(
            {"pattern": pattern, "degree": int(flags.sum()), "n": int(count)}
        )
    return pd.DataFrame(rows).sort_values(
        ["degree", "n", "pattern"], ascending=[True, False, True]
    )


def prepare_feature_tables(
    direct_counts: pd.DataFrame,
    prefix: str,
) -> dict[tuple[str, str], pd.DataFrame]:
    named = direct_counts.copy()
    named.index = pd.Index(
        [terminal_name(str(path), prefix) for path in direct_counts.index],
        name="feature",
    )
    collapsed_counts = named.groupby(level=0, sort=False).sum()
    modes = {
        "direct_path": direct_counts,
        "collapsed_terminal": collapsed_counts,
    }
    result: dict[tuple[str, str], pd.DataFrame] = {}
    for taxon_mode, counts in modes.items():
        result[(taxon_mode, "counts")] = counts
        result[(taxon_mode, "relative_abundance")] = counts.div(
            counts.sum(axis=0).replace(0, np.nan), axis=1
        )
    return result


def audit(data_root: Path, output: Path) -> None:
    raw = data_root / "raw"
    meta = pd.read_csv(
        raw / "meta.csv",
        dtype={"ids": str, "subject": str, "canon": str},
    )
    mpa = pd.read_csv(raw / "combine_brk.mpa", sep="\t", index_col=0)
    mpa.columns = mpa.columns.str.removesuffix(".brk.report")
    sample_order = meta["ids"].tolist()

    direct_counts = {
        "genus": mpa.loc[rank_mask(mpa.index, "g__"), sample_order].astype(float),
        "species": mpa.loc[rank_mask(mpa.index, "s__"), sample_order].astype(float),
    }
    raw_shannon = shannon_from_counts(direct_counts["species"])
    legacy_counts = pd.read_csv(data_root / "species.csv", index_col=0).loc[:, sample_order]
    legacy_shannon = shannon_from_counts(legacy_counts)
    responses = {
        "raw_bracken_shannon": raw_shannon,
        "legacy_species_csv_shannon": legacy_shannon,
    }
    groups = recovery_groups(meta, legacy_shannon)

    summary_rows: list[dict[str, object]] = []
    upset_rows: list[dict[str, object]] = []
    degree_rows: list[dict[str, object]] = []
    exact_pattern_rows: list[dict[str, object]] = []
    inventory_rows: list[dict[str, object]] = []

    for level, prefix in [("genus", "g__"), ("species", "s__")]:
        tables = prepare_feature_tables(direct_counts[level], prefix)
        # A second biologically meaningful interpretation of "collapsed
        # genus" is to sum species leaves by their genus ancestor, rather than
        # use the MPA genus-rank rows.  Keep it as an explicit sensitivity
        # mode; it is not interchangeable with merely collapsing duplicate
        # terminal names.
        if level == "genus":
            species_to_genus = direct_counts["species"].copy()
            species_to_genus.index = pd.Index(
                [
                    genus_name_from_species_path(str(path))
                    for path in species_to_genus.index
                ],
                name="feature",
            )
            species_to_genus = species_to_genus.groupby(level=0, sort=False).sum()
            tables[("species_descendants_to_genus", "counts")] = species_to_genus
            tables[("species_descendants_to_genus", "relative_abundance")] = (
                species_to_genus.div(
                    species_to_genus.sum(axis=0).replace(0, np.nan), axis=1
                )
            )
        for (taxon_mode, abundance_scale), feature_table in tables.items():
            inventory_rows.append(
                {
                    "level": level,
                    "taxon_mode": taxon_mode,
                    "abundance_scale": abundance_scale,
                    "n_features": feature_table.shape[0],
                    "n_samples": feature_table.shape[1],
                }
            )
            for response_name, response in responses.items():
                for group in ["LDR", "HDR"]:
                    significant = significant_matrix(
                        feature_table, response, meta, groups, group
                    )
                    observed = pattern_counts(significant)
                    expected = SOURCE_UPSET[(level, group)]
                    absolute_errors = {
                        pattern: abs(observed[pattern] - expected[pattern])
                        for pattern in PATTERN_ORDER
                    }
                    expected_total = sum(expected.values())
                    observed_total = sum(observed.values())
                    source_degree = np.asarray(
                        SOURCE_DEGREE_ESTIMATES[(level, group)], dtype=int
                    )
                    observed_degree = np.asarray(
                        [
                            np.count_nonzero(significant.sum(axis=1) == value)
                            for value in range(1, 17)
                        ],
                        dtype=int,
                    )
                    summary_rows.append(
                        {
                            "level": level,
                            "group": group,
                            "taxon_mode": taxon_mode,
                            "abundance_scale": abundance_scale,
                            "shannon_response": response_name,
                            "source_total": expected_total,
                            "observed_total": observed_total,
                            "total_difference": observed_total - expected_total,
                            "pattern_L1": sum(absolute_errors.values()),
                            "pattern_relative_L1": (
                                sum(absolute_errors.values()) / expected_total
                            ),
                            "pattern_RMSE": float(
                                np.sqrt(
                                    np.mean(
                                        [
                                            (observed[p] - expected[p]) ** 2
                                            for p in PATTERN_ORDER
                                        ]
                                    )
                                )
                            ),
                            "degree_digitized_L1": int(
                                np.abs(observed_degree - source_degree).sum()
                            ),
                            "degree_digitized_relative_L1": float(
                                np.abs(observed_degree - source_degree).sum()
                                / source_degree.sum()
                            ),
                            "degree_digitized_RMSE": float(
                                np.sqrt(
                                    np.mean(
                                        (observed_degree - source_degree) ** 2
                                    )
                                )
                            ),
                        }
                    )
                    for pattern in PATTERN_ORDER:
                        upset_rows.append(
                            {
                                "level": level,
                                "group": group,
                                "taxon_mode": taxon_mode,
                                "abundance_scale": abundance_scale,
                                "shannon_response": response_name,
                                "pattern": pattern,
                                "source_n": expected[pattern],
                                "observed_n": observed[pattern],
                                "difference": observed[pattern] - expected[pattern],
                                "absolute_error": absolute_errors[pattern],
                            }
                        )
                    degree = significant.sum(axis=1)
                    for value in range(1, 17):
                        degree_rows.append(
                            {
                                "level": level,
                                "group": group,
                                "taxon_mode": taxon_mode,
                                "abundance_scale": abundance_scale,
                                "shannon_response": response_name,
                                "degree": value,
                                "n_features": int(np.count_nonzero(degree == value)),
                            }
                        )
                    exact_patterns = exact_timepoint_patterns(significant)
                    for row in exact_patterns.itertuples(index=False):
                        exact_pattern_rows.append(
                            {
                                "level": level,
                                "group": group,
                                "taxon_mode": taxon_mode,
                                "abundance_scale": abundance_scale,
                                "shannon_response": response_name,
                                "pattern": row.pattern,
                                "degree": row.degree,
                                "n": row.n,
                            }
                        )

    output.mkdir(parents=True, exist_ok=True)
    summary_frame = pd.DataFrame(summary_rows).sort_values(
        ["level", "group", "pattern_relative_L1", "taxon_mode"]
    )
    summary_frame.to_csv(output / "sfig3_variant_summary.csv", index=False)
    pd.DataFrame(upset_rows).to_csv(output / "sfig3_variant_upset_counts.csv", index=False)
    pd.DataFrame(degree_rows).to_csv(output / "sfig3_variant_degree_counts.csv", index=False)
    pd.DataFrame(exact_pattern_rows).to_csv(
        output / "sfig3_variant_exact_timepoint_patterns.csv", index=False
    )
    pd.DataFrame(inventory_rows).drop_duplicates().to_csv(
        output / "sfig3_variant_feature_inventory.csv", index=False
    )
    source_upset_rows = []
    for (level, group), patterns in SOURCE_UPSET.items():
        for pattern in PATTERN_ORDER:
            source_upset_rows.append(
                {
                    "level": level,
                    "group": group,
                    "pattern": pattern,
                    "A": "A" in pattern,
                    "B": "B" in pattern,
                    "C": "C" in pattern,
                    "source_n": patterns[pattern],
                    "provenance": "value printed above source UpSet bar",
                }
            )
    pd.DataFrame(source_upset_rows).to_csv(
        output / "sfig3_source_upset_counts.csv", index=False
    )
    source_degree_rows = []
    for (level, group), values in SOURCE_DEGREE_ESTIMATES.items():
        for degree, value in enumerate(values, start=1):
            source_degree_rows.append(
                {
                    "level": level,
                    "group": group,
                    "degree": degree,
                    "digitized_n": value,
                    "uncertainty_plus_minus_features": SOURCE_DEGREE_UNCERTAINTY[level],
                    "provenance": (
                        "raster estimate constrained to exact source UpSet total; "
                        "not an author-supplied label"
                    ),
                }
            )
    pd.DataFrame(source_degree_rows).to_csv(
        output / "sfig3_source_degree_digitization.csv", index=False
    )

    # A fair global comparison uses one common taxon construction across all
    # four source panels.  ``direct_path`` is selected because terminal-name
    # collapse is provably identical here, while species-descendant collapse
    # is defined only for genus.
    global_scores = (
        summary_frame.loc[summary_frame["taxon_mode"].eq("direct_path")]
        .groupby(["abundance_scale", "shannon_response"], as_index=False)
        .agg(
            source_total=("source_total", "sum"),
            upset_pattern_L1=("pattern_L1", "sum"),
            digitized_degree_L1=("degree_digitized_L1", "sum"),
        )
    )
    global_scores["upset_pattern_relative_L1"] = (
        global_scores["upset_pattern_L1"] / global_scores["source_total"]
    )
    global_scores["digitized_degree_relative_L1"] = (
        global_scores["digitized_degree_L1"] / global_scores["source_total"]
    )
    global_scores["mean_relative_L1"] = (
        global_scores["upset_pattern_relative_L1"]
        + global_scores["digitized_degree_relative_L1"]
    ) / 2.0
    global_scores = global_scores.sort_values("mean_relative_L1")
    global_scores.to_csv(output / "sfig3_variant_global_scores.csv", index=False)

    best = global_scores.iloc[0]
    (output / "README.md").write_text(
        f"""# Supplementary Figure 3 variant audit

The source-labelled seven-way phase intersections were transcribed from
Supplementary Figure 3 in `抗生素稿件 isme j/Supplymentary material.pdf`
(page 6; embedded 1341 x 1745 JPEG). Panels C/D were digitized from that
raster and constrained to the exact labelled intersection totals. Genus bins
have about ±1-feature resolution and species bins about ±7; they are not
author-supplied exact values.

The best common rule is **relative abundance correlated against raw-Bracken
Shannon**. Its global UpSet-pattern relative L1 error is
{best['upset_pattern_relative_L1']:.3f}, and its digitized-degree relative L1
error is {best['digitized_degree_relative_L1']:.3f}. No tested rule recovers
the source counts exactly. Raw MPA paths and terminal-name collapse are
identical because the supplied bacterial rank rows contain no duplicate
terminal genus or species labels. Species-descendant aggregation changes the
genus result only marginally and is not a clear basis for changing the figure.

Each correlation uses 15 subjects and the source rule `rho > 0` plus two-sided
`P < 0.05`. P is the usual t approximation on ranked values with 13 degrees
of freedom, matching the existing renderer. The manuscript does not identify
the original correlation implementation or its tie handling, so this remains
a reconstruction detail. L1/RMSE values are descriptive reconstruction-error
scores, not hypothesis tests, and the variant screen is exploratory.

`sfig3_variant_exact_timepoint_patterns.csv` gives exact reconstructed panel
E/F pattern counts for every tested rule. The corresponding source bars are
too dense and unlabeled to recover exact individual pattern counts from the
raster; no pseudo-exact source values are asserted.
""",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent
        / "outputs"
        / "integrated_fig124_raw"
        / "tables"
        / "sfig3_variant_audit",
    )
    args = parser.parse_args()
    audit(args.data_root.resolve(), args.output.resolve())
    print(f"sFig3 variant audit written to {args.output.resolve()}")


if __name__ == "__main__":
    main()
