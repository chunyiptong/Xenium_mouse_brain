#!/usr/bin/env python3
"""Predict a grayscale 2D TIFF and export points.tsv + details.json.gz for R.

Example:
    python spotiflow_to_r.py --img /path/Per2.tif --out /path/results

R:
    points <- read.delim('/path/results/points.tsv')
    details <- jsonlite::fromJSON(gzfile('/path/results/details.json.gz'),
                                 simplifyMatrix = TRUE)

Spotiflow API: https://weigertlab.org/spotiflow/api.html
"""
import argparse
import gzip
import json
import math
from pathlib import Path

import numpy as np
import tifffile


def json_value(value):
    """Preserve nested array axes; encode unavailable/nonfinite values as null."""
    if isinstance(value, np.ndarray):
        if np.issubdtype(value.dtype, np.inexact) and not np.isfinite(value).all():
            array = value.astype(object)
            array[~np.isfinite(value)] = None
            return array.tolist()
        return value.tolist()
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, dict):
        return {key: json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def export_results(points, details, out):
    """Export actual model outputs, including all attributes of details."""
    points = np.asarray(points)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError(f"Expected points shaped (N, 2) in y,x order; got {points.shape}")
    n = len(points)

    def spot_values(name):
        value = getattr(details, name, None)
        if value is None:
            return np.full(n, np.nan)
        array = np.asarray(value).reshape(-1)
        if len(array) != n:
            raise ValueError(f"details.{name} does not match the number of spots")
        return array

    table = np.column_stack((points, spot_values("prob"), spot_values("intens")))
    out = Path(out)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"Output directory must be new or empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    # Coordinates stay 0-based. 'score' matches analyze_per2_2d.R --points-tsv.
    np.savetxt(out / "points.tsv", table, delimiter="\t", fmt="%.17g",
               header="y_px\tx_px\tscore\tintensity", comments="")
    # Heatmap/flow/subpix can be large, so compress the complete JSON losslessly.
    with gzip.open(out / "details.json.gz", "wt", encoding="utf-8") as f:
        json.dump(json_value(vars(details)), f, ensure_ascii=False,
                  allow_nan=False, separators=(",", ":"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--img", required=True, type=Path, help="Input grayscale 2D TIFF path")
    parser.add_argument("--out", type=Path, help="Output folder; default: <image_stem>_spotiflow")
    args = parser.parse_args()
    img_path = args.img.expanduser().resolve()
    out = args.out.expanduser().resolve() if args.out else img_path.with_name(img_path.stem + "_spotiflow")
    if not img_path.is_file():
        parser.error(f"Image does not exist: {img_path}")
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error(f"Output directory must be new or empty: {out}")

    img = tifffile.imread(img_path)
    if img.ndim != 2 or not img.size or not np.isfinite(img).all():
        parser.error(f"Expected one finite grayscale 2D plane; got shape {img.shape}")

    from spotiflow.model import Spotiflow

    model = Spotiflow.from_pretrained("general")
    points, details = model.predict(img)
    export_results(points, details, out)
    print(f"Detected {len(points)} spots. Saved to: {out}")
    print("points.tsv: y_px, x_px, score (= details.prob), intensity (= details.intens)")
    print("details.json.gz: all details fields, including heatmap, flow and subpix")


if __name__ == "__main__":
    main()
