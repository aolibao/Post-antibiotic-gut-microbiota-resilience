#!/usr/bin/env python3
"""Validate Figure 3 inputs, exports, traceability and visual integrity."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"

EXPECTED_PIXELS = (7805, 11162)
EXPECTED_DPI = 1000


@dataclass
class Check:
    name: str
    status: str
    observed: str
    expected: str
    severity: str = "error"


def _check(condition: bool, name: str, observed: object, expected: object, severity: str = "error") -> Check:
    return Check(name, "PASS" if condition else ("WARN" if severity == "warning" else "FAIL"), str(observed), str(expected), severity)


def validate_inputs(output: Path) -> list[Check]:
    checks: list[Check] = []
    manifest_path = output / "inputs" / "input_manifest.json"
    checks.append(_check(manifest_path.exists(), "input manifest exists", manifest_path.exists(), True))
    if not manifest_path.exists():
        return checks
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    text = json.dumps(manifest).lower()
    checks.append(_check("structural" in text and "qpcr" in text, "qPCR placeholder is disclosed", "structural" in text and "qpcr" in text, True))
    nonmeasured = any(
        marker in text
        for marker in (
            "not_measured",
            "not measured",
            "not biological",
            "non-biological",
            '"is_measured_qpcr": false',
        )
    )
    checks.append(_check(nonmeasured, "qPCR placeholder is non-measured", "non-measured marker present" if nonmeasured else "missing", "explicit non-measured marker"))
    for group in ("LDR", "HDR"):
        group_dir = output / "inputs" / group
        reads_path = group_dir / "reads.tsv"
        metadata_path = group_dir / "metadata.tsv"
        pickle_path = group_dir / f"{group}_study.pkl"
        for name, path in [("reads", reads_path), ("metadata", metadata_path), ("Study pickle", pickle_path)]:
            checks.append(_check(path.exists(), f"{group} {name} exists", path.exists(), True))
        if reads_path.exists() and metadata_path.exists():
            reads = pd.read_csv(reads_path, sep="\t", index_col=0)
            meta = pd.read_csv(metadata_path, sep="\t")
            checks.extend(
                [
                    _check(reads.shape[0] == 259, f"{group} fixed taxon count", reads.shape[0], 259),
                    _check(reads.shape[1] == 240, f"{group} sample count", reads.shape[1], 240),
                    _check(meta["subject"].nunique() == 15, f"{group} subject count", meta["subject"].nunique(), 15),
                    _check(meta.groupby("subject").size().eq(16).all(), f"{group} complete longitudinal series", meta.groupby("subject").size().value_counts().to_dict(), "16 samples for every subject"),
                    _check(set(meta["canon"]) == {f"A{i}" for i in range(1, 5)} | {f"B{i}" for i in range(1, 7)} | {f"C{i}" for i in range(1, 7)}, f"{group} canonical stages", sorted(meta["canon"].unique()), "A1-A4, B1-B6, C1-C6"),
                ]
            )
    return checks


def validate_posterior(output: Path, allow_development: bool) -> list[Check]:
    checks: list[Check] = []
    manifest_path = output / "tables" / "posterior_summary_manifest.json"
    checks.append(_check(manifest_path.exists(), "posterior summary manifest exists", manifest_path.exists(), True))
    if not manifest_path.exists():
        return checks
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    for group in ("LDR", "HDR"):
        info = payload.get("groups", {}).get(group, {})
        tier = str(info.get("evidence_tier", "missing"))
        formal = tier.startswith("formal")
        severity = "warning" if allow_development else "error"
        checks.append(_check(formal, f"{group} posterior is formal multichain", tier, "formal_multichain", severity=severity))
        chains = info.get("chains", [])
        draws = info.get("posterior_draws_by_chain", [])
        checks.append(_check(len(chains) == 3 if formal else allow_development, f"{group} chain count", len(chains), 3, severity=severity))
        enough_draws = len(draws) == 3 and all(int(value) >= 1500 for value in draws)
        checks.append(_check(enough_draws if formal else allow_development, f"{group} retained draws", draws, "three chains, >=1500 draws each", severity=severity))
    diagnostics_path = output / "tables" / "posterior_diagnostics.csv"
    checks.append(_check(diagnostics_path.exists(), "posterior diagnostics exist", diagnostics_path.exists(), True))
    if diagnostics_path.exists():
        diagnostics = pd.read_csv(diagnostics_path)
        for group in ("LDR", "HDR"):
            subset = diagnostics.loc[diagnostics["group"] == group]
            tier_values = subset.get("evidence_tier", pd.Series(dtype=str)).astype(str)
            formal = bool(len(tier_values) and tier_values.str.startswith("formal").all())
            if formal:
                rhat = pd.to_numeric(subset["rhat_rank_folded"], errors="coerce")
                ess = pd.to_numeric(subset["ess_bulk"], errors="coerce")
                finite_fraction = float(rhat.notna().mean())
                pass_fraction = float((rhat.dropna() <= 1.10).mean()) if rhat.notna().any() else 0.0
                checks.append(_check(finite_fraction >= 0.95, f"{group} finite R-hat coverage", f"{finite_fraction:.3f}", ">=0.95"))
                checks.append(_check(pass_fraction >= 0.90, f"{group} R-hat <=1.10", f"{pass_fraction:.3f}", ">=0.90"))
                checks.append(_check(float(ess.median()) >= 100, f"{group} median bulk ESS", f"{float(ess.median()):.1f}", ">=100"))
    return checks


def _pdf_raster_ppi(path: Path) -> list[tuple[float, float]]:
    if shutil.which("pdfimages") is None:
        return []
    run = subprocess.run(["pdfimages", "-list", str(path)], check=False, text=True, capture_output=True)
    values: list[tuple[float, float]] = []
    for line in run.stdout.splitlines():
        fields = line.split()
        if len(fields) < 14 or not fields[0].isdigit():
            continue
        try:
            values.append((float(fields[-3]), float(fields[-2])))
        except ValueError:
            continue
    return values


def validate_exports(output: Path, allow_development: bool) -> list[Check]:
    checks: list[Check] = []
    figure_dir = output / "figures"
    paths = {suffix: figure_dir / f"Fig3_restored.{suffix}" for suffix in ("svg", "pdf", "tif", "png")}
    for suffix, path in paths.items():
        checks.append(_check(path.exists() and path.stat().st_size > 0, f"{suffix.upper()} export exists", path.stat().st_size if path.exists() else 0, ">0 bytes"))
    if paths["tif"].exists():
        with Image.open(paths["tif"]) as image:
            dpi = tuple(round(float(x)) for x in image.info.get("dpi", (0, 0)))
            compression = str(image.info.get("compression", ""))
            checks.extend(
                [
                    _check(image.size == EXPECTED_PIXELS, "TIFF pixel dimensions", image.size, EXPECTED_PIXELS),
                    _check(dpi == (EXPECTED_DPI, EXPECTED_DPI), "TIFF resolution", dpi, (EXPECTED_DPI, EXPECTED_DPI)),
                    _check(image.mode == "RGB", "TIFF color mode", image.mode, "opaque RGB"),
                    _check("tiff_lzw" in compression.lower() or "lzw" in compression.lower(), "TIFF lossless compression", compression, "LZW", severity="warning"),
                ]
            )
    if paths["svg"].exists():
        svg = paths["svg"].read_text(encoding="utf-8", errors="replace")
        text_nodes = len(re.findall(r"<text\b", svg))
        checks.append(_check(text_nodes > 0, "SVG text remains editable", text_nodes, ">0 <text> nodes"))
    if paths["pdf"].exists():
        ppi = _pdf_raster_ppi(paths["pdf"])
        minimum = min((min(x, y) for x, y in ppi), default=float("nan"))
        checks.append(
            _check(
                not ppi or minimum >= 595,
                "PDF rasterized layer resolution",
                "no raster layers" if not ppi else f"minimum {minimum:.1f} ppi",
                ">=595 ppi",
                severity="warning" if not ppi else "error",
            )
        )
    if paths["png"].exists():
        with Image.open(paths["png"]) as preview:
            rgb = np.asarray(preview.convert("RGB"))
        border = np.concatenate([rgb[:, :2].reshape(-1, 3), rgb[:, -2:].reshape(-1, 3)], axis=0)
        nonwhite = int((border.min(axis=1) < 245).sum())
        checks.append(_check(nonwhite == 0, "left/right artist clipping guard", nonwhite, "0 non-white border pixels"))
    render_manifest = output / "render_manifest.json"
    checks.append(_check(render_manifest.exists(), "render manifest exists", render_manifest.exists(), True))
    if render_manifest.exists():
        payload = json.loads(render_manifest.read_text(encoding="utf-8"))
        correction = str(payload.get("caption_correction", ""))
        checks.append(_check("E=HDR" in correction and "F=LDR" in correction, "caption E/F correction recorded", correction, "E=HDR; F=LDR"))
        semantics = json.dumps(payload.get("network_semantics", {})).lower()
        checks.append(_check("spearman" in semantics and "bayes factor >10" in semantics, "network semantics separated", semantics, "main Spearman; separate directed BF>10"))
        development = bool(payload.get("development_only", True))
        checks.append(
            _check(
                not development,
                "render uses formal posterior",
                development,
                False,
                severity="warning" if allow_development else "error",
            )
        )
    return checks


def make_comparison(project_root: Path, output: Path) -> Path | None:
    source = project_root / "抗生素稿件 isme j" / "Figures" / "Figure3.tif"
    restored = output / "figures" / "Fig3_restored.png"
    if not source.exists() or not restored.exists():
        return None
    with Image.open(source) as im:
        source_rgb = im.convert("RGB")
        ratio = 860 / source_rgb.width
        source_preview = source_rgb.resize((860, round(source_rgb.height * ratio)), Image.Resampling.LANCZOS)
    with Image.open(restored) as im:
        restored_rgb = im.convert("RGB")
        ratio = 860 / restored_rgb.width
        restored_preview = restored_rgb.resize((860, round(restored_rgb.height * ratio)), Image.Resampling.LANCZOS)
    fig, axes = plt.subplots(1, 2, figsize=(12, 8.7), facecolor="white")
    for ax, image, title in zip(axes, [source_preview, restored_preview], ["Source manuscript", "Restored / current inference state"]):
        ax.imshow(image)
        ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
        ax.set_axis_off()
    fig.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.01, wspace=0.025)
    path = output / "comparisons" / "Fig3_source_vs_restored.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)
    return path


def validate(project_root: Path, output: Path, strict: bool = False, allow_development: bool = False) -> dict[str, object]:
    qa = output / "qa"
    qa.mkdir(parents=True, exist_ok=True)
    checks = validate_inputs(output) + validate_posterior(output, allow_development) + validate_exports(output, allow_development)
    comparison = make_comparison(project_root, output)
    if comparison is not None:
        checks.append(_check(comparison.exists(), "source/restored comparison generated", comparison, "existing PNG"))
    table = pd.DataFrame([asdict(x) for x in checks])
    table.to_csv(qa / "validation_checks.tsv", sep="\t", index=False)
    counts = table["status"].value_counts().to_dict()
    report = {
        "status_counts": counts,
        "strict": strict,
        "allow_development": allow_development,
        "checks": [asdict(x) for x in checks],
        "comparison": str(comparison) if comparison else None,
    }
    (qa / "validation_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    failed = int((table["status"] == "FAIL").sum())
    if strict and failed:
        raise SystemExit(f"Figure 3 validation failed: {failed} error-level check(s)")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "outputs" / "fig3_mdsine2")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument(
        "--allow-development",
        action="store_true",
        help="Downgrade non-formal posterior checks to warnings for layout-only QA.",
    )
    args = parser.parse_args()
    report = validate(args.project_root.resolve(), args.output.resolve(), args.strict, args.allow_development)
    print(json.dumps(report["status_counts"], indent=2))


if __name__ == "__main__":
    main()
