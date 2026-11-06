"""Validate Python exports in R without running or downloading the model."""
import importlib.util
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace

import numpy as np

workflow = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("exporter", workflow / "spotiflow_to_r.py")
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)

with tempfile.TemporaryDirectory(prefix="spotiflow-r-export-") as temp:
    root = Path(temp)
    details = SimpleNamespace(
        prob=np.array([0.8, 0.9]), intens=np.array([123., 456.]),
        heatmap=np.arange(12).reshape(3, 4).astype(np.float32),
        flow=np.arange(36).reshape(3, 4, 3).astype(np.float32),
        subpix=np.arange(24).reshape(3, 4, 2).astype(np.float32), fit_params=None)
    points = np.array([[0.25, 2.5], [1.75, 3.]])
    exporter.export_results(points, details, root / "normal")
    empty = SimpleNamespace(prob=np.array([]), intens=np.array([]),
        heatmap=details.heatmap, flow=None, subpix=None, fit_params=False)
    exporter.export_results(np.empty((0, 2)), empty, root / "empty")
    missing = SimpleNamespace(prob=np.array([0.8, np.nan]), intens=None,
        heatmap=np.array([[np.nan, np.inf], [-np.inf, 1.]]), flow=None, subpix=None)
    exporter.export_results(points, missing, root / "missing")
    try:
        exporter.export_results(points, details, root / "normal")
    except ValueError:
        pass
    else:
        raise AssertionError("Did not protect existing results")
    script = root / "validate.R"
    script.write_text('''source(commandArgs(TRUE)[1])
root <- commandArgs(TRUE)[2]
x <- read_spotiflow(file.path(root, "normal"))
stopifnot(nrow(x$points) == 2L, x$points$y_px[1] == 0.25,
          x$points$x_px[1] == 2.5, x$details$prob[2] == 0.9,
          x$details$intens[1] == 123,
          identical(dim(x$details$heatmap), c(3L,4L)),
          identical(dim(x$details$flow), c(3L,4L,3L)),
          identical(dim(x$details$subpix), c(3L,4L,2L)),
          x$details$heatmap[2,3] == 6,
          x$details$flow[2,3,2] == 19,
          x$details$subpix[2,3,2] == 13)
e <- read_spotiflow(file.path(root, "empty"))
stopifnot(nrow(e$points) == 0L, length(e$details$prob) == 0L, is.null(e$details$flow))
m <- read_spotiflow(file.path(root, "missing"))
stopifnot(is.na(m$points$score[2]), all(is.na(m$points$intensity)),
          sum(is.na(m$details$heatmap)) == 3L)
cat("PASS: R round-trip, array axes/values, empty spots, missing values, overwrite protection\\n")
''')
    subprocess.run(["Rscript", "--vanilla", str(script),
                    str(workflow / "read_spotiflow.R"), str(root)], check=True)
