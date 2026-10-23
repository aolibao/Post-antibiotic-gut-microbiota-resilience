#!/usr/bin/env python3
"""Create the participant/sample flow requested in Reviewer 1 comment 4."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import pandas as pd
from PIL import Image


PROJECT = Path(__file__).resolve().parents[2]
WORKSPACE = Path(__file__).resolve().parents[1]
META = PROJECT / "raw" / "meta.csv"
TABLES = WORKSPACE / "03_tables"
FIGURES = WORKSPACE / "04_figures"


def derive_flow(meta: pd.DataFrame) -> pd.DataFrame:
    required = {"ids", "subject", "abx", "canon"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Missing metadata columns: {sorted(missing)}")
    if meta["ids"].duplicated().any():
        raise ValueError("Sample identifiers are not unique")
    per_subject = meta.groupby("subject").agg(
        samples=("ids", "size"),
        exposure=("abx", lambda x: ",".join(sorted(set(x.astype(str))))),
        timepoints=("canon", "nunique"),
    )
    if not (per_subject["samples"] == 16).all() or not (per_subject["timepoints"] == 16).all():
        raise ValueError("The supplied metadata are not complete at 16 samples/time points per subject")
    exposed_subjects = int((per_subject["exposure"] == "Y").sum())
    unexposed_subjects = int((per_subject["exposure"] == "N").sum())
    rows = [
        {
            "stage": "available_longitudinal_cohort",
            "subjects_n": int(per_subject.shape[0]),
            "samples_n": int(meta.shape[0]),
            "role": "all supplied human participants",
            "reason_or_note": "No participant or sample exclusion is evident in the supplied metadata",
        },
        {
            "stage": "ciprofloxacin_exposed_primary_cohort",
            "subjects_n": exposed_subjects,
            "samples_n": int((meta["abx"] == "Y").sum()),
            "role": "primary full-cohort continuous recovery analysis",
            "reason_or_note": "Exposure-defined analytic cohort; all 16 samples retained per subject",
        },
        {
            "stage": "antibiotic_unexposed_reference",
            "subjects_n": unexposed_subjects,
            "samples_n": int((meta["abx"] == "N").sum()),
            "role": "descriptive reference only",
            "reason_or_note": (
                "n=5 in raw/meta.csv and PRJNA974858; the Nature source Methods says n=4, "
                "so the discrepancy requires author reconciliation; no confirmatory inference"
            ),
        },
        {
            "stage": "original_endpoint_shannon_LDR_reference",
            "subjects_n": 15,
            "samples_n": 15 * 16,
            "role": "reference binary grouping; secondary to full-cohort inference",
            "reason_or_note": "Original LDR: bottom 15 by mean C4-C6 Shannon; retain label with the defining metric stated",
        },
        {
            "stage": "middle_late_diversity_subjects",
            "subjects_n": 30,
            "samples_n": 30 * 16,
            "role": "retained in primary full-cohort analysis",
            "reason_or_note": "Excluded only from the original extreme-stratum illustration",
        },
        {
            "stage": "original_endpoint_shannon_HDR_reference",
            "subjects_n": 15,
            "samples_n": 15 * 16,
            "role": "reference binary grouping; secondary to full-cohort inference",
            "reason_or_note": "Original HDR: top 15 by mean C4-C6 Shannon; retain label with the defining metric stated",
        },
    ]
    return pd.DataFrame(rows)


def box(ax, x, y, w, h, title, subtitle, color):
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=0.9,
        edgecolor=color,
        facecolor="white",
    )
    ax.add_patch(patch)
    ax.text(x, y + 0.025, title, ha="center", va="center", fontsize=8.0, weight="bold", color=color)
    ax.text(x, y - 0.035, subtitle, ha="center", va="center", fontsize=6.7, color="#333333", linespacing=1.15)


def arrow(ax, start, end):
    ax.add_patch(
        FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9, linewidth=0.8, color="#707070")
    )


def render(flow: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    navy, blue, orange, grey = "#294C60", "#227C9D", "#D97706", "#6B7280"

    box(ax, 0.50, 0.88, 0.50, 0.13, "Available longitudinal cohort", "65 participants · 1,040 samples\n16 time points per participant", navy)
    box(ax, 0.30, 0.62, 0.38, 0.15, "Ciprofloxacin-exposed cohort", "60 participants · 960 samples\nContinuous + binary sensitivity analyses", blue)
    box(ax, 0.78, 0.62, 0.31, 0.15, "Unexposed reference", "5 participants · 80 samples\nDescriptive only", orange)
    arrow(ax, (0.44, 0.81), (0.34, 0.70))
    arrow(ax, (0.58, 0.81), (0.73, 0.70))

    box(ax, 0.15, 0.32, 0.27, 0.17, "LDR reference group", "Bottom 15 by C4–C6 Shannon\n240 samples", "#1677A8")
    box(ax, 0.50, 0.32, 0.27, 0.17, "Middle 30", "Retained in full-cohort analysis\n480 samples", grey)
    box(ax, 0.85, 0.32, 0.27, 0.17, "HDR reference group", "Top 15 by C4–C6 Shannon\n240 samples", "#C23B22")
    arrow(ax, (0.28, 0.54), (0.17, 0.41))
    arrow(ax, (0.31, 0.54), (0.49, 0.41))
    arrow(ax, (0.34, 0.54), (0.82, 0.41))

    ax.text(
        0.50,
        0.10,
        "The original 15 + 15 endpoint-Shannon-defined HDR/LDR groups are the reference binary comparison, not cohort exclusions.\n"
        "Local metadata/BioProject contain 5 controls, whereas the source article Methods says 4; authors must reconcile this discrepancy.",
        ha="center",
        va="center",
        fontsize=6.9,
        color="#333333",
    )
    ax.text(0.015, 0.985, "Reviewer analysis · cohort flow", va="top", fontsize=9.0, weight="bold")
    fig.tight_layout(pad=0.7)

    FIGURES.mkdir(parents=True, exist_ok=True)
    stem = FIGURES / "cohort_flow"
    fig.savefig(stem.with_suffix(".pdf"), dpi=600, bbox_inches="tight", facecolor="white")
    tmp = stem.with_name(stem.name + ".rgba.tmp.tif")
    fig.savefig(tmp, dpi=600, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(tmp) as image:
        image.convert("RGB").save(stem.with_suffix(".tif"), dpi=(600, 600), compression="tiff_lzw")
    tmp.unlink()
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(META)
    flow = derive_flow(meta)
    flow.to_csv(TABLES / "cohort_flow.csv", index=False)
    calendar = (
        meta.groupby(["abx", "canons", "canon"], dropna=False)
        .agg(
            subjects_n=("subject", "nunique"),
            samples_n=("ids", "size"),
            actual_day_min=("day", "min"),
            actual_day_q25=("day", lambda x: x.quantile(0.25)),
            actual_day_median=("day", "median"),
            actual_day_q75=("day", lambda x: x.quantile(0.75)),
            actual_day_max=("day", "max"),
        )
        .reset_index()
        .sort_values(["abx", "canons"])
    )
    calendar.to_csv(TABLES / "sample_calendar.csv", index=False)
    subject_meta = meta.sort_values(["subject", "canons"]).groupby("subject", as_index=False).first()
    coverage = []
    for abx, block in subject_meta.groupby("abx", dropna=False):
        for field in ["age", "sex"]:
            coverage.append(
                {
                    "exposure": abx,
                    "field": field,
                    "subjects_total": int(block.shape[0]),
                    "subjects_nonmissing": int(block[field].notna().sum()),
                    "subjects_missing": int(block[field].isna().sum()),
                }
            )
    for absent in ["BMI", "diet", "bowel_preparation", "clinical_indication", "concomitant_medications"]:
        coverage.append(
            {
                "exposure": "all",
                "field": absent,
                "subjects_total": int(subject_meta.shape[0]),
                "subjects_nonmissing": 0,
                "subjects_missing": int(subject_meta.shape[0]),
            }
        )
    pd.DataFrame(coverage).to_csv(TABLES / "host_covariate_availability.csv", index=False)
    render(flow)
    manifest = {
        "source": str(META),
        "source_sha256": hashlib.sha256(META.read_bytes()).hexdigest(),
        "subjects": int(meta["subject"].nunique()),
        "samples": int(meta.shape[0]),
        "confirmation_needed": "Verify against the source cohort that no exclusions occurred before these supplied tables.",
        "outputs": [
            str(TABLES / "cohort_flow.csv"),
            str(TABLES / "sample_calendar.csv"),
            str(TABLES / "host_covariate_availability.csv"),
            str(FIGURES / "cohort_flow.pdf"),
            str(FIGURES / "cohort_flow.tif"),
        ],
    }
    (TABLES / "cohort_flow_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

