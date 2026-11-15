#!/usr/bin/env python3
"""Compare the R port with the original Python workflow on synthetic fixtures.

Requires the original script's Python dependencies (not Spotiflow), Rscript,
and the R packages ijtiff, png, jsonlite. No real image results are generated.
"""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image
import roifile
from scipy import ndimage as ndi
import tifffile

WORKFLOW = Path(__file__).resolve().parents[1]


def run(command):
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise AssertionError(f"Command failed: {command}\n{result.stdout}\n{result.stderr}")
    return result


def compare_table(left, right):
    with left.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        a, fields = list(reader), reader.fieldnames
    with right.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        b, other_fields = list(reader), reader.fieldnames
    assert fields == other_fields, (left.name, fields, other_fields)
    assert len(a) == len(b), left.name
    for i, (ar, br) in enumerate(zip(a, b)):
        for field in fields:
            av, bv = ar[field], br[field]
            try:
                af, bf = float(av), float(bv)
            except ValueError:
                assert av == bv, (left.name, i, field, av, bv)
            else:
                np.testing.assert_allclose(af, bf, rtol=1e-11, atol=1e-10,
                                           equal_nan=True, err_msg=f"{left.name}:{i}:{field}")


def main():
    with tempfile.TemporaryDirectory(prefix="per2-r-parity-") as temporary:
        root = Path(temporary)
        h, w = 39, 53
        raw = np.full((h, w), 100, np.uint16)
        raw[10, 10] += 500
        raw[20, 33] += 900
        raw[21, 35] += 600
        raw[0, 0] = 4095
        labels = np.zeros((h, w), np.uint16)
        labels[5:16, 5:17] = 1
        labels[17:28, 28:40] = 300
        labels[31:36, 45:50] = 60000  # Retain this zero-spot object.
        tifffile.imwrite(root / "raw.tif", raw)
        tifffile.imwrite(root / "mask.tif", labels)
        Image.fromarray(raw).save(root / "raw.png")
        Image.fromarray(labels).save(root / "mask.png")
        high_labels = labels.astype(np.uint32)
        high_labels[high_labels == 60000] = 70000
        tifffile.imwrite(root / "mask32.tif", high_labels)
        rois = [
            roifile.ImagejRoi.frompoints(np.array([[3, 3], [3, 30], [34, 30], [34, 3]]), name="SCN_left"),
            roifile.ImagejRoi.frompoints(np.array([[28.25, 15.5], [48.5, 15.5], [48.5, 35.25], [28.25, 35.25]], np.float32), name="SCN_右"),
        ]
        for roi in rois:
            roi.roitype = roifile.ROI_TYPE.POLYGON
        roifile.roiwrite(root / "rois.zip", rois)
        roifile.roiwrite(root / "single.roi", rois[0])
        points = np.array([[10, 10], [20, 33], [21.25, 35.5], [0, 0],
                           [38.75, 52.75], [3, 3], [16.5, 12], [30, 34]], float)
        scores = np.linspace(.5, .99, len(points))
        np.savetxt(root / "points.tsv", np.column_stack((points, scores)), delimiter="\t",
                   header="y_px\tx_px\tscore", comments="")
        (root / "empty.tsv").write_text("y_px\tx_px\tscore\n")
        cases = [
            ("nuclei", "tif", "mask.tif", "points.tsv", 0, "rois.zip"),
            ("nuclei", "png", "mask.png", "points.tsv", 2.1, "rois.zip"),
            ("cells", "tif", "mask32.tif", "points.tsv", 0, "rois.zip"),
            ("nuclei", "tif", "mask.tif", "empty.tsv", 0, "rois.zip"),
            ("nuclei", "tif", "mask.tif", "empty.tsv", 1.5, "rois.zip"),
            ("cells", "png", "mask.png", "points.tsv", 0, "single.roi"),
        ]
        for case, (kind, extension, mask, coords, expansion, roi_file) in enumerate(cases):
            args = ["--per2", str(root / f"raw.{extension}"), "--mask", str(root / mask),
                    "--mask-kind", kind, "--expand-um", str(expansion),
                    "--rois", str(root / roi_file), "--image-id", "synthetic",
                    "--image-mode", "single_plane", "--source-description", "Synthetic parity test",
                    "--pixel-size-x", "0.4", "--pixel-size-y", "0.9", "--points-tsv", str(root / coords)]
            py_out, r_out = root / f"py_{case}", root / f"r_{case}"
            run([sys.executable, str(WORKFLOW / "analyze_per2_2d.py"), *args, "--out", str(py_out)])
            run(["Rscript", "--vanilla", str(WORKFLOW / "analyze_per2_2d.R"), *args, "--out", str(r_out)])
            for table in ("per2_spots.txt", "cells.txt", "contours.txt"):
                compare_table(py_out / table, r_out / table)
            np.testing.assert_array_equal(tifffile.imread(py_out / "assignment_labels.tif"),
                                          tifffile.imread(r_out / "assignment_labels.tif"))
            with tifffile.TiffFile(r_out / "assignment_labels.tif") as tif:
                assert tif.pages[0].dtype == np.dtype("uint32")
            with Image.open(r_out / "qc_overlay.png") as image:
                assert image.size == (2400, 1600)
                image.verify()
            info = json.loads((r_out / "run_info.json").read_text())
            assert info["shape_yx"] == [h, w]
            assert info["parameters"]["expand_um"] == expansion
            assert info["model"] == "external_coordinates"
            if case == 0:
                with (r_out / "per2_spots.txt").open() as f:
                    first = next(csv.DictReader(f, delimiter="\t"))
                assert float(first["peak_raw"]) == 600
                assert float(first["background_local"]) == 100
                assert float(first["integrated_bgsub"]) == 500
            protected = subprocess.run(["Rscript", "--vanilla", str(WORKFLOW / "analyze_per2_2d.R"),
                                        *args, "--out", str(r_out)], capture_output=True, text=True)
            assert protected.returncode != 0 and "nonempty" in protected.stderr
            print(f"PASS full workflow {case + 1}/{len(cases)}", flush=True)

        # Random anisotropic distance transforms, including exact-distance ties.
        rng = np.random.default_rng(712)
        transforms = []
        for i in range(24):
            shape = [(17, 23), (1, 23), (17, 1)][i % 3]
            mask = np.where(rng.random(shape) < .15, rng.integers(1, 9, shape), 0)
            mask.flat[0] = 1
            sx, sy = [(1., 1.), (.4, .9), (2., .5)][i % 3]
            distance, index = ndi.distance_transform_edt(mask == 0, sampling=(sy, sx), return_indices=True)
            expanded = mask.copy()
            hits = (mask == 0) & (distance <= 3.1)
            expanded[hits] = mask[tuple(index[:, hits])]
            transforms.append(dict(mask=mask.tolist(), expected=expanded.tolist(), sx=sx, sy=sy))
        fixtures = root / "transforms.json"
        fixtures.write_text(json.dumps(transforms))
        r_test = root / "check_transforms.R"
        r_test.write_text('''args <- commandArgs(TRUE)
source(args[1])
tests <- jsonlite::fromJSON(args[2], simplifyVector = FALSE)
for (case in tests) {
  mask <- do.call(rbind, case$mask)
  expected <- do.call(rbind, case$expected)
  storage.mode(mask) <- storage.mode(expected) <- "integer"
  actual <- expand_labels(mask, case$sx, case$sy, 3.1)
  stopifnot(identical(actual, expected))
}
cat("PASS 24 anisotropic and tie distance transforms\\n")
''')
        print(run(["Rscript", "--vanilla", str(r_test), str(WORKFLOW / "analyze_per2_2d.R"),
                   str(fixtures)]).stdout, end="")
    print("All R/Python comparisons passed.")


if __name__ == "__main__":
    main()
