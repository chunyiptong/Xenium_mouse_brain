"""Validate native RData exports with R's load(); no model download required."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace

import numpy as np
import tifffile

workflow = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("exporter", workflow / "spotiflow_to_rdata.py")
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)

with tempfile.TemporaryDirectory(prefix="spotiflow-rdata-") as temp:
    root = Path(temp)
    points = np.array([[0.25, 2.5], [1.75, 3.]], dtype=np.float32)
    details = SimpleNamespace(
        prob=np.array([.8, .9], dtype=np.float32), intens=np.array([123., 456.]),
        heatmap=np.arange(12, dtype=np.float32).reshape(3, 4),
        flow=np.arange(36, dtype=np.float32).reshape(3, 4, 3),
        subpix=np.arange(24, dtype=np.float32).reshape(3, 4, 2), fit_params=None)
    exporter.export_results(points, details, root / "normal")
    empty = SimpleNamespace(prob=np.array([], dtype=np.float32), intens=np.array([]),
        heatmap=details.heatmap, flow=None, subpix=None, fit_params=False)
    exporter.export_results(np.empty((0, 2)), empty, root / "empty")
    special = SimpleNamespace(prob=np.array([.8, np.nan]), intens=None,
        heatmap=np.array([[np.nan, np.inf], [-np.inf, 1.]], dtype=np.float32),
        flow=None, subpix=None, fit_params={"amplitude": np.array([1., 2.], dtype=np.float32)})
    exporter.export_results(points, special, root / "special")
    try:
        exporter.export_results(points, details, root / "normal")
    except ValueError:
        pass
    else:
        raise AssertionError("Existing results were not protected")

    r_script = root / "validate.R"
    r_script.write_text('''root <- commandArgs(TRUE)[1]
read_result <- function(name) {
  env <- new.env()
  loaded <- load(file.path(root, name, "results.RData"), envir=env)
  stopifnot(setequal(loaded, c("points", "details")))
  as.list(env)
}
x <- read_result("normal")
stopifnot(identical(dim(x$points),c(2L,2L)), x$points[1,1] == .25,
  x$points[1,2] == 2.5, abs(x$details$prob[2]-.9) < 1e-6,
  x$details$intens[1] == 123, is.null(x$details$fit_params),
  identical(dim(x$details$heatmap),c(3L,4L)),
  identical(dim(x$details$flow),c(3L,4L,3L)),
  identical(dim(x$details$subpix),c(3L,4L,2L)),
  x$details$heatmap[2,3] == 6, x$details$flow[2,3,2] == 19,
  x$details$subpix[2,3,2] == 13)
tab <- read.delim(file.path(root,"normal","points.tsv"))
stopifnot(identical(names(tab),c("y_px","x_px","score","intensity")),
  all(tab$y_px == x$points[,1]), all(tab$x_px == x$points[,2]))
e <- read_result("empty")
stopifnot(identical(dim(e$points),c(0L,2L)), length(e$details$prob) == 0L,
  is.null(e$details$flow), identical(e$details$fit_params,FALSE),
  nrow(read.delim(file.path(root,"empty","points.tsv"))) == 0L)
s <- read_result("special")
stopifnot(is.nan(s$details$prob[2]), is.null(s$details$intens),
  is.nan(s$details$heatmap[1,1]), s$details$heatmap[1,2] == Inf,
  s$details$heatmap[2,1] == -Inf, s$details$fit_params$amplitude[2] == 2)
cat("PASS: R load(), coordinate/array axes, nested details, empty spots, NaN/Inf/NULL, TSV\\n")
''')
    subprocess.run(["Rscript", "--vanilla", str(r_script), str(root)], check=True)

    # Exercise rejection paths before the optional Spotiflow import.
    tifffile.imwrite(root / "raw.tif", np.zeros((3, 4), dtype=np.uint16))
    tifffile.imwrite(root / "rgb.tif", np.zeros((3, 4, 3), dtype=np.uint16), photometric="rgb")
    for args, message in [
        (["--img", str(root / "raw.tif"), "--prob-thresh", "nan"], "finite number"),
        (["--img", str(root / "raw.tif"), "--prob-thresh", "1.5"], "between 0 and 1"),
        (["--img", str(root / "raw.tif"), "--out", str(root / "normal")], "new or empty"),
        (["--img", str(root / "rgb.tif")], "grayscale 2D"),
    ]:
        result = subprocess.run([sys.executable, str(workflow / "spotiflow_to_rdata.py"), *args],
                                capture_output=True, text=True)
        assert result.returncode != 0 and message in result.stderr, result.stderr
    print("PASS: CLI validation and overwrite protection")
