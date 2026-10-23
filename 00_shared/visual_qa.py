#!/usr/bin/env python3
"""Quantitative visual diagnostics for the high-fidelity figure layer.

These metrics assess layout/density only.  They must never be used to alter data
values or to claim analytical agreement.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageFilter
from scipy.fft import dctn
from scipy.spatial.distance import jensenshannon
from skimage.metrics import structural_similarity


def fit_rgb(path: Path, width: int = 1200) -> np.ndarray:
    with Image.open(path) as opened:
        image = opened.convert("RGB")
        height = max(1, round(width * image.height / image.width))
        image = image.resize((width, height), Image.Resampling.LANCZOS)
        return np.asarray(image)


def ink_metrics(rgb: np.ndarray) -> dict[str, float]:
    ink = np.any(rgb < 245, axis=2)
    gray = np.asarray(Image.fromarray(rgb).convert("L"))
    dark = gray < 100
    blank_rows = ink.mean(axis=1) < .001
    longest = run_length(blank_rows)
    return {
        "ink_fraction": float(ink.mean()),
        "dark_fraction": float(dark.mean()),
        "blank_row_fraction": float(blank_rows.mean()),
        "longest_blank_row_run_fraction": float(longest / len(blank_rows)),
    }


def run_length(values: np.ndarray) -> int:
    best = current = 0
    for value in values:
        current = current + 1 if value else 0
        best = max(best, current)
    return best


def color_hist(rgb: np.ndarray) -> np.ndarray:
    mask = np.any(rgb < 245, axis=2)
    pixels = rgb[mask]
    if not len(pixels):
        return np.ones(4096) / 4096
    bins = (pixels // 16).astype(int)
    index = bins[:, 0] * 256 + bins[:, 1] * 16 + bins[:, 2]
    hist = np.bincount(index, minlength=4096).astype(float) + 1e-12
    return hist / hist.sum()


def phash(rgb: np.ndarray) -> np.ndarray:
    image = Image.fromarray(rgb).convert("L").resize((32, 32), Image.Resampling.LANCZOS)
    coeff = dctn(np.asarray(image, dtype=float), norm="ortho")[:8, :8]
    threshold = np.median(coeff[1:, :])
    return coeff.ravel() > threshold


def compare(original: Path, recreated: Path) -> dict[str, object]:
    a = fit_rgb(original)
    b = fit_rgb(recreated)
    target_h = min(a.shape[0], b.shape[0])
    a2 = np.asarray(Image.fromarray(a).resize((1200, target_h), Image.Resampling.LANCZOS))
    b2 = np.asarray(Image.fromarray(b).resize((1200, target_h), Image.Resampling.LANCZOS))
    ga = np.asarray(Image.fromarray(a2).convert("L"))
    gb = np.asarray(Image.fromarray(b2).convert("L"))
    blur_a = np.asarray(Image.fromarray(ga).filter(ImageFilter.GaussianBlur(radius=8)))
    blur_b = np.asarray(Image.fromarray(gb).filter(ImageFilter.GaussianBlur(radius=8)))
    with Image.open(original) as io, Image.open(recreated) as ir:
        aspect_o = io.width / io.height
        aspect_r = ir.width / ir.height
        pixel_o = [io.width, io.height]
        pixel_r = [ir.width, ir.height]
    out: dict[str, object] = {
        "original": str(original),
        "recreated": str(recreated),
        "original_pixels": pixel_o,
        "recreated_pixels": pixel_r,
        "original_aspect": aspect_o,
        "recreated_aspect": aspect_r,
        "aspect_relative_error": abs(aspect_r / aspect_o - 1),
        "grayscale_ssim_diagnostic": float(structural_similarity(ga, gb, data_range=255)),
        "blurred_ssim_diagnostic": float(structural_similarity(blur_a, blur_b, data_range=255)),
        "phash_hamming_64_diagnostic": int(np.count_nonzero(phash(a2) != phash(b2))),
        "nonwhite_color_hist_jsd_diagnostic": float(jensenshannon(color_hist(a2), color_hist(b2)) ** 2),
    }
    out.update({f"original_{k}": v for k, v in ink_metrics(a2).items()})
    out.update({f"recreated_{k}": v for k, v in ink_metrics(b2).items()})
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--high-fidelity", type=Path, default=Path(__file__).resolve().parent / "outputs" / "high_fidelity")
    args = parser.parse_args()
    source = args.data_root / "抗生素稿件 isme j" / "Figures"
    rows = []
    for number in [1, 2, 4]:
        row = compare(source / f"Figure{number}.tif", args.high_fidelity / "figures" / f"Figure{number}_high_fidelity.png")
        row["figure"] = f"Figure{number}"
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(args.high_fidelity / "visual_qa.csv", index=False)
    (args.high_fidelity / "visual_qa.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8")
    print(table[["figure", "aspect_relative_error", "blurred_ssim_diagnostic", "original_ink_fraction", "recreated_ink_fraction"]].to_string(index=False))


if __name__ == "__main__":
    main()
