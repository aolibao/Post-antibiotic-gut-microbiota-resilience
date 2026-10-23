#!/usr/bin/env python3
"""Participant-level phase contrasts for Reviewer 1, Comment 6.

The analysis uses the 268 species eligible for the original Fig. 2F display,
CLR-transforms their Bracken relative abundances, averages repeated samples
within participant and prespecified phase, and compares the original HDR/LDR
reference groups.  A single BH family covers all species and all three
contrasts.  The participant, not the sample, is the independent unit.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "reproduction"))
from load_raw_bracken import load_raw_bracken_dataset  # noqa: E402


OUT = ROOT / "reviewer_response" / "R1.6"
ELIGIBLE = (
    ROOT
    / "reproduction"
    / "outputs"
    / "integrated_fig124_raw"
    / "tables"
    / "Fig2F_LEfSe_all_eligible_species.csv"
)

PHASES = {
    "baseline_HDR_minus_LDR": ["A1", "A2", "A3", "A4"],
    "antibiotic_change_HDR_minus_LDR": ["B1", "B2", "B3", "B4", "B5"],
    "late_recovery_change_HDR_minus_LDR": ["C4", "C5", "C6"],
}


def bh_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * len(p) / np.arange(1, len(p) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1.0)
    return out


def welch_result(hdr: np.ndarray, ldr: np.ndarray) -> dict[str, float]:
    n1, n0 = len(hdr), len(ldr)
    m1, m0 = float(np.mean(hdr)), float(np.mean(ldr))
    v1, v0 = float(np.var(hdr, ddof=1)), float(np.var(ldr, ddof=1))
    se = np.sqrt(v1 / n1 + v0 / n0)
    df = (v1 / n1 + v0 / n0) ** 2 / (
        (v1 / n1) ** 2 / (n1 - 1) + (v0 / n0) ** 2 / (n0 - 1)
    )
    diff = m1 - m0
    if se == 0:
        p = 1.0 if diff == 0 else 0.0
        lo = hi = diff
    else:
        tval = diff / se
        p = float(2 * stats.t.sf(abs(tval), df))
        crit = float(stats.t.ppf(0.975, df))
        lo, hi = diff - crit * se, diff + crit * se
    pooled = np.sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
    d = 0.0 if pooled == 0 else diff / pooled
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    return {
        "n_HDR": n1,
        "n_LDR": n0,
        "mean_HDR": m1,
        "mean_LDR": m0,
        "HDR_minus_LDR": diff,
        "ci_low": lo,
        "ci_high": hi,
        "welch_df": df,
        "p_value": p,
        "hedges_g": correction * d,
    }


def main() -> None:
    dataset = load_raw_bracken_dataset(ROOT)
    eligible = pd.read_csv(ELIGIBLE)["species"].drop_duplicates().tolist()
    missing = sorted(set(eligible) - set(dataset.species_named.index))
    if missing:
        raise RuntimeError(f"{len(missing)} eligible species absent from Bracken matrix: {missing[:5]}")

    abund = dataset.species_named.loc[eligible].astype(float)
    positive = abund.to_numpy()[abund.to_numpy() > 0]
    pseudocount = float(positive.min() / 2.0)
    logged = np.log(abund + pseudocount)
    clr = logged.subtract(logged.mean(axis=0), axis=1)

    meta = dataset.metrics.loc[
        (dataset.metrics["abx"] == "Y")
        & (dataset.metrics["recovery_group"].isin(["HDR", "LDR"])),
        ["ids", "subject", "canon", "recovery_group"],
    ].copy()
    meta["canon"] = meta["canon"].astype(str)
    group_map = meta.drop_duplicates("subject").set_index("subject")["recovery_group"]

    baseline_samples = meta.loc[meta["canon"].isin(PHASES["baseline_HDR_minus_LDR"])]
    baseline = (
        clr[baseline_samples["ids"]]
        .T.assign(subject=baseline_samples.set_index("ids").loc[clr[baseline_samples["ids"]].columns, "subject"].values)
        .groupby("subject")
        .mean()
    )

    rows: list[dict[str, object]] = []
    participant_rows: list[pd.DataFrame] = []
    for contrast, labels in PHASES.items():
        phase_samples = meta.loc[meta["canon"].isin(labels)]
        phase = (
            clr[phase_samples["ids"]]
            .T.assign(subject=phase_samples.set_index("ids").loc[clr[phase_samples["ids"]].columns, "subject"].values)
            .groupby("subject")
            .mean()
        )
        value = phase if contrast == "baseline_HDR_minus_LDR" else phase - baseline.loc[phase.index]
        value["group"] = value.index.map(group_map)
        long = value.reset_index().melt(
            id_vars=["subject", "group"], var_name="species", value_name="participant_phase_value"
        )
        long["contrast"] = contrast
        participant_rows.append(long)
        for species in eligible:
            hdr = value.loc[value["group"] == "HDR", species].to_numpy()
            ldr = value.loc[value["group"] == "LDR", species].to_numpy()
            rows.append({"contrast": contrast, "species": species, **welch_result(hdr, ldr)})

    result = pd.DataFrame(rows)
    result["q_global"] = bh_adjust(result["p_value"].to_numpy())
    result["significant_global_FDR"] = result["q_global"] < 0.05
    result = result.sort_values(["contrast", "q_global", "p_value", "species"])
    result.to_csv(OUT / "tables" / "R1_6_participant_phase_CLR_global_FDR.csv", index=False)
    pd.concat(participant_rows, ignore_index=True).to_csv(
        OUT / "tables" / "R1_6_participant_phase_CLR_source_values.csv", index=False
    )

    summary = (
        result.groupby("contrast", sort=False)
        .agg(
            species_tested=("species", "size"),
            global_q_lt_0_05=("significant_global_FDR", "sum"),
            nominal_p_lt_0_05=("p_value", lambda x: int((x < 0.05).sum())),
            min_global_q=("q_global", "min"),
        )
        .reset_index()
    )
    summary["independent_unit"] = "participant"
    summary["n_HDR"] = 15
    summary["n_LDR"] = 15
    summary["global_family_size"] = len(result)
    summary["pseudocount_percent"] = pseudocount
    summary.to_csv(OUT / "tables" / "R1_6_participant_phase_CLR_summary.csv", index=False)

    top = result.groupby("contrast", sort=False, group_keys=False).head(10)
    top.to_csv(OUT / "tables" / "R1_6_participant_phase_CLR_top10.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
