#!/usr/bin/env python3
"""Reproducible evidence audit for Reviewer 1, major comment 5.

This script uses only the supplied participant metadata for numerical summaries.
External study/accession facts are recorded separately with their source URLs and
verification status; unverified manuscript assertions are never promoted to facts.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW_META = ROOT / "raw" / "meta.csv"
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
LOGS = HERE / "logs"

CANON_ORDER = ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4",
               "B5", "B6", "C1", "C2", "C3", "C4", "C5", "C6"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_and_validate() -> pd.DataFrame:
    df = pd.read_csv(RAW_META)
    required = {"ids", "library", "subject", "day", "canon", "abx", "age", "sex"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if df["ids"].duplicated().any():
        raise ValueError("Sample IDs are not unique")
    if set(df["abx"].dropna().unique()) != {"Y", "N"}:
        raise ValueError("Unexpected antibiotic-exposure coding")
    counts = df.groupby("subject", observed=True).size()
    if not counts.eq(16).all():
        raise ValueError("Not every participant has exactly 16 samples")
    observed = set(df["canon"].dropna().unique())
    if observed != set(CANON_ORDER):
        raise ValueError(f"Unexpected canonical labels: {sorted(observed)}")
    return df


def participant_frame(df: pd.DataFrame) -> pd.DataFrame:
    def first_nonmissing(x: pd.Series):
        values = x.dropna().unique()
        if len(values) > 1:
            raise ValueError(f"Participant-level field varies within participant: {values}")
        return values[0] if len(values) else np.nan

    return (
        df.groupby(["subject", "abx"], as_index=False, observed=True)
        .agg(age=("age", first_nonmissing), sex=("sex", first_nonmissing), samples_n=("ids", "size"))
    )


def write_tables(df: pd.DataFrame) -> dict[str, Path]:
    TABLES.mkdir(parents=True, exist_ok=True)
    people = participant_frame(df)

    cohort = pd.DataFrame([
        {
            "cohort": "all supplied participants",
            "participants_n": people["subject"].nunique(),
            "samples_n": len(df),
            "interpretation": "complete local longitudinal metadata; 16 samples per participant",
        },
        {
            "cohort": "ciprofloxacin-exposed",
            "participants_n": people.loc[people.abx.eq("Y"), "subject"].nunique(),
            "samples_n": int(df.abx.eq("Y").sum()),
            "interpretation": "primary human analysis cohort",
        },
        {
            "cohort": "antibiotic-unexposed",
            "participants_n": people.loc[people.abx.eq("N"), "subject"].nunique(),
            "samples_n": int(df.abx.eq("N").sum()),
            "interpretation": "descriptive reference only; source-paper n=4 versus repository/local n=5 remains unresolved",
        },
    ])
    cohort_path = TABLES / "R1_5_cohort_summary.csv"
    cohort.to_csv(cohort_path, index=False)

    exposed = df.loc[df.abx.eq("Y")].copy()
    calendar = (
        exposed.groupby("canon", observed=True)["day"]
        .agg(subjects_n="count", actual_day_min="min", actual_day_q25=lambda x: x.quantile(.25),
             actual_day_median="median", actual_day_q75=lambda x: x.quantile(.75), actual_day_max="max")
        .reindex(CANON_ORDER)
        .reset_index()
    )
    calendar.insert(1, "nominal_reference", "first dose = day 0; last treatment day = day 4")
    for column in ["actual_day_min", "actual_day_q25", "actual_day_median", "actual_day_q75", "actual_day_max"]:
        calendar[column.replace("actual_day", "day_from_last_dose")] = calendar[column] - 4
    calendar_path = TABLES / "R1_5_sample_calendar_relative_to_doses.csv"
    calendar.to_csv(calendar_path, index=False)

    exposed_people = people.loc[people.abx.eq("Y")]
    fields = [
        ("age", int(exposed_people.age.notna().sum()), "present in local metadata"),
        ("sex", int(exposed_people.sex.notna().sum()), "present in local metadata"),
        ("BMI", 0, "field absent from local metadata"),
        ("diet", 0, "field absent from local metadata"),
        ("bowel preparation", 0, "field absent from local metadata"),
        ("concomitant medications", 0, "field absent from local metadata"),
    ]
    covariates = pd.DataFrame(fields, columns=["field", "participants_nonmissing", "status"])
    covariates.insert(1, "participants_total", len(exposed_people))
    covariates["participants_missing"] = covariates.participants_total - covariates.participants_nonmissing
    covariates["completeness_percent"] = 100 * covariates.participants_nonmissing / covariates.participants_total
    covariate_path = TABLES / "R1_5_host_covariate_availability.csv"
    covariates.to_csv(covariate_path, index=False)

    accessions = pd.DataFrame([
        {
            "accession_or_record": "PRJNA974858",
            "manuscript_or_source_role": "previously published human longitudinal shotgun metagenomic cohort",
            "verification_status_2026_10_12": "public record verified",
            "verified_scope": "60 healthy adults exposed to ciprofloxacin; 960 samples; five unexposed participants also described by BioProject",
            "required_action": "cite Yaffe et al. 2025 and state explicitly that the human analysis reuses public data",
            "official_url": "https://www.ncbi.nlm.nih.gov/bioproject/974858",
        },
        {
            "accession_or_record": "PRJNA1470252",
            "manuscript_or_source_role": "manuscript asserts newly generated mouse 16S rRNA data",
            "verification_status_2026_10_12": "not independently verified in a public record",
            "verified_scope": "none",
            "required_action": "author must confirm public release, organism, assay, sample count and figure mapping",
            "official_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1470252",
        },
        {
            "accession_or_record": "PRJNA1472735",
            "manuscript_or_source_role": "manuscript asserts newly generated metagenomic data",
            "verification_status_2026_10_12": "not independently verified in a public record",
            "verified_scope": "none",
            "required_action": "author must identify cohort/specimen, assay, sample count, public release and figure mapping",
            "official_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1472735",
        },
        {
            "accession_or_record": "Yaffe et al. 2025",
            "manuscript_or_source_role": "primary publication for the reused human cohort",
            "verification_status_2026_10_12": "article metadata verified",
            "verified_scope": "Nature 641, 182-191; DOI 10.1038/s41586-025-08781-x; PMID 40269166; PMCID PMC12168174",
            "required_action": "add full citation to the reference list and cite it at first mention of the public cohort",
            "official_url": "https://www.nature.com/articles/s41586-025-08781-x",
        },
    ])
    accession_path = TABLES / "R1_5_accession_and_provenance_audit.csv"
    accessions.to_csv(accession_path, index=False)

    facts = pd.DataFrame([
        ["population", "60 healthy U.S. adults", "Yaffe et al. 2025 / PRJNA974858", "verified"],
        ["ciprofloxacin regimen", "500 mg orally twice daily for 5 days", "Yaffe et al. 2025 / PRJNA974858", "verified"],
        ["clinical indication", "protocol exposure in healthy volunteers; not treatment of a documented infection", "source study design", "verified"],
        ["prior antibiotics", "no antibiotic use for at least 6 months before sampling", "Yaffe et al. 2025", "verified"],
        ["source-study ethics", "Stanford University IRB #25268; written informed consent", "Yaffe et al. 2025 Methods", "verified source-study fact"],
        ["current secondary-analysis ethics", "approval, exemption or not-human-subjects determination", "not present in supplied package", "AUTHOR_INPUT_NEEDED"],
        ["unexposed reference count", "source Methods reportedly n=4; NCBI/local metadata n=5", "source article versus PRJNA974858/raw metadata", "AUTHOR_INPUT_NEEDED"],
        ["concomitant medications", "not present in local metadata", "raw/meta.csv schema", "AUTHOR_INPUT_NEEDED"],
    ], columns=["item", "finding", "evidence", "status"])
    facts_path = TABLES / "R1_5_verified_facts_and_gaps.csv"
    facts.to_csv(facts_path, index=False)

    return {
        "cohort": cohort_path,
        "calendar": calendar_path,
        "covariates": covariate_path,
        "accessions": accession_path,
        "facts": facts_path,
    }


def make_figure(paths: dict[str, Path]) -> list[Path]:
    calendar = pd.read_csv(paths["calendar"])
    cov = pd.read_csv(paths["covariates"])

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 8.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(8.27, 6.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[2.15, 1.1], width_ratios=[1.55, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    # A: actual-day calendar for exposed participants.
    y = np.arange(len(calendar))
    ax_a.axvspan(0, 4, color="#F4A261", alpha=.18, zorder=0)
    ax_a.axvline(0, color="#B4472D", lw=1, ls="--")
    ax_a.axvline(4, color="#B4472D", lw=1, ls=":")
    ax_a.hlines(y, calendar.actual_day_min, calendar.actual_day_max,
               color="#97A6B2", lw=1.2, zorder=1)
    ax_a.hlines(y, calendar.actual_day_q25, calendar.actual_day_q75,
               color="#2A6F97", lw=4.5, zorder=2)
    ax_a.scatter(calendar.actual_day_median, y, s=19, color="#16324F", zorder=3)
    ax_a.set_yticks(y, calendar.canon)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Recorded day relative to first ciprofloxacin dose")
    ax_a.set_title("A  Actual sampling calendar (n = 60)", loc="left", weight="bold")
    ax_a.text(.01, -.13, "Thin line: range; thick line: IQR; dot: median. Shading: treatment days 0–4.",
              transform=ax_a.transAxes, fontsize=7, color="#4F5963")
    ax_a.spines[["top", "right"]].set_visible(False)
    ax_a.grid(axis="x", color="#DCE2E6", lw=.55)

    # B: availability of requested host covariates.
    cov_plot = cov.iloc[::-1].copy()
    yb = np.arange(len(cov_plot))
    ax_b.barh(yb, cov_plot.participants_total, color="#E6EAED", height=.65)
    ax_b.barh(yb, cov_plot.participants_nonmissing, color="#2A6F97", height=.65)
    ax_b.set_yticks(yb, cov_plot.field)
    ax_b.set_xlim(0, 60)
    ax_b.set_xlabel("Participants with non-missing data")
    ax_b.set_title("B  Local covariate availability", loc="left", weight="bold")
    for yi, n in zip(yb, cov_plot.participants_nonmissing):
        ax_b.text(min(n + 1.0, 56), yi, f"{int(n)}/60", va="center", fontsize=7)
    ax_b.spines[["top", "right", "left"]].set_visible(False)
    ax_b.grid(axis="x", color="#DCE2E6", lw=.55)

    # C: provenance status kept deliberately simple and auditable.
    ax_c.axis("off")
    ax_c.set_title("C  Data provenance and reporting status", loc="left", weight="bold", pad=7)
    boxes = [
        (0.01, "PUBLIC HUMAN DATA", "Yaffe et al., Nature 2025\nPRJNA974858 · verified", "#DCEFE8", "#176B5B"),
        (0.35, "CURRENT STUDY", "Secondary human microbiome analyses\nNew mouse work claimed in manuscript", "#E5EEF7", "#2A5C85"),
        (0.69, "REPORTING LIMITS", "BMI · diet · bowel preparation · medications\nunavailable; not adjusted", "#FCE8E6", "#9C3B32"),
    ]
    for x, header, body, face, edge in boxes:
        ax_c.add_patch(plt.Rectangle((x, .14), .30, .68, transform=ax_c.transAxes,
                                     facecolor=face, edgecolor=edge, lw=1.1))
        ax_c.text(x + .015, .69, header, transform=ax_c.transAxes, fontsize=7.4,
                  weight="bold", color=edge, va="top")
        ax_c.text(x + .015, .55, body, transform=ax_c.transAxes, fontsize=7.3,
                  color="#263238", va="top", linespacing=1.45)
    ax_c.annotate("", xy=(.345, .48), xytext=(.315, .48), xycoords=ax_c.transAxes,
                  arrowprops=dict(arrowstyle="->", lw=1.1, color="#66727A"))
    ax_c.annotate("", xy=(.685, .48), xytext=(.655, .48), xycoords=ax_c.transAxes,
                  arrowprops=dict(arrowstyle="->", lw=1.1, color="#66727A"))

    fig.suptitle("Reviewer 1.5 · Human-cohort provenance, sampling and metadata audit",
                 x=.01, ha="left", fontsize=11.5, weight="bold")
    pdf = FIGURES / "R1_5_provenance_sampling_covariate_audit.pdf"
    tif = FIGURES / "R1_5_provenance_sampling_covariate_audit.tif"
    png = FIGURES / "R1_5_provenance_sampling_covariate_audit_preview.png"
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(tif, dpi=600, pil_kwargs={"compression": "tiff_lzw"}, bbox_inches="tight")
    fig.savefig(png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return [pdf, tif, png]


def write_manifest(paths: dict[str, Path], figures: list[Path], df: pd.DataFrame) -> None:
    outputs = list(paths.values()) + figures
    manifest = {
        "analysis": "Reviewer 1 major comment 5 provenance/calendar/covariate audit",
        "input": str(RAW_META.relative_to(ROOT)),
        "input_sha256": sha256(RAW_META),
        "samples_n": int(len(df)),
        "participants_n": int(df.subject.nunique()),
        "exposed_participants_n": int(df.loc[df.abx.eq("Y"), "subject"].nunique()),
        "unexposed_participants_n": int(df.loc[df.abx.eq("N"), "subject"].nunique()),
        "external_verification_date": "2026-10-12",
        "outputs": [{"path": str(p.relative_to(HERE)), "sha256": sha256(p)} for p in outputs],
    }
    (LOGS / "R1_5_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def main() -> None:
    LOGS.mkdir(parents=True, exist_ok=True)
    df = load_and_validate()
    paths = write_tables(df)
    figures = make_figure(paths)
    write_manifest(paths, figures, df)
    print("R1.5 analysis complete")
    print(f"participants={df.subject.nunique()} samples={len(df)}")


if __name__ == "__main__":
    main()
