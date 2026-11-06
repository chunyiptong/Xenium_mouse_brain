#!/usr/bin/env python3
"""Detect spots in a grayscale 2D TIFF and export results.RData + points.tsv.

Install in your Python environment:
    python -m pip install spotiflow numpy tifffile rdata

Example:
    python spotiflow_to_rdata.py --img /path/Per2.tif --out /path/results
    python spotiflow_to_rdata.py --img /path/Per2.tif --prob-thresh 0.2

R:
    load('/path/results/results.RData')
    colnames(points) <- c('y_px', 'x_px')
    head(points)
    details$prob
    dim(details$flow)

points is the original N-by-2 coordinate matrix (y, x), in 0-based pixels.
details is a named R list containing all fields returned by Spotiflow.
The TSV additionally combines coordinates, score and intensity in one table.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import tifffile


def r_compatible(value):
    """Convert NumPy floats to R doubles without flattening or swapping axes."""
    if isinstance(value, np.ndarray):
        if value.dtype.kind == "f":
            return value.astype(np.float64, copy=False)
        return value
    if isinstance(value, np.generic):
        return r_compatible(value.item())
    if isinstance(value, SimpleNamespace):
        return r_compatible(vars(value))
    if isinstance(value, dict):
        return {key: r_compatible(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [r_compatible(item) for item in value]
    # None becomes R NULL; NaN/Inf remain NaN/Inf, rather than JSON null.
    return value


def check_output(out):
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"Output directory must be new or empty: {out}")


def export_results(points, details, out):
    """Write native R objects, plus a table compatible with --points-tsv."""
    import rdata

    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("Expected finite points shaped (N, 2), in y,x order")
    fields = details if isinstance(details, dict) else vars(details)

    def spot_values(name):
        value = fields.get(name)
        if value is None:
            return np.full(len(points), np.nan)
        values = np.asarray(value, dtype=np.float64).reshape(-1)
        if len(values) != len(points):
            raise ValueError(f"details.{name} does not match the number of spots")
        return values

    table = np.column_stack((points, spot_values("prob"), spot_values("intens")))
    payload = {"points": points, "details": r_compatible(fields)}
    out = Path(out).expanduser()
    check_output(out)
    out.mkdir(parents=True, exist_ok=True)
    # The temporary file prevents a failed serialization from looking complete.
    temporary = out / ".spotiflow.results.RData.tmp"
    try:
        rdata.write_rda(temporary, payload, compression="gzip")
        temporary.replace(out / "spotiflow.results.RData")
    finally:
        temporary.unlink(missing_ok=True)
    np.savetxt(out / "points.tsv", table, delimiter="\t", fmt="%.17g",
               header="y_px\tx_px\tscore\tintensity", comments="")
    return out / "spotiflow.results.RData"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--img", required=True, type=Path, help="Input grayscale 2D TIFF path")
    parser.add_argument("--out", type=Path,
                        help="Output folder; default: <image_stem>_spotiflow_rdata beside the image")
    parser.add_argument("--prob-thresh", type=float, default=None,
                        help="Detection threshold in [0, 1]; omitted uses the model default")
    args = parser.parse_args(argv)
    if args.prob_thresh is not None and not (
        np.isfinite(args.prob_thresh) and 0 <= args.prob_thresh <= 1
    ):
        parser.error("--prob-thresh must be a finite number between 0 and 1")
    img_path = args.img.expanduser().resolve()
    out = (args.out.expanduser().resolve() if args.out else
#           img_path.with_name(img_path.stem + "_spotiflow_rdata"))
           img_path.with_name("spotiflow_rdata"))
    if not img_path.is_file():
        parser.error(f"Image does not exist: {img_path}")
    try:
        check_output(out)
    except ValueError as error:
        parser.error(str(error))

    img = tifffile.imread(img_path)
    if img.ndim != 2 or not img.size or not np.isfinite(img).all():
        parser.error(f"Expected one finite grayscale 2D plane; got shape {img.shape}")
    try:
        import rdata  # Check the export dependency before running inference.
        from spotiflow.model import Spotiflow
    except ImportError as error:
        parser.error(f"Missing dependency: {error}. Run: python -m pip install spotiflow rdata")

    model = Spotiflow.from_pretrained("general")
    print("Detection threshold:", "model default" if args.prob_thresh is None else args.prob_thresh)
    points, details = model.predict(img, prob_thresh=args.prob_thresh)
    result_path = export_results(points, details, out)
    print(f"Detected {len(points)} spots. Saved: {result_path}")
    print("results.RData: points (y,x matrix) and details (named list); read with R load()")
    print("points.tsv: y_px, x_px, score (= details.prob), intensity (= details.intens)")


if __name__ == "__main__":
    main()
