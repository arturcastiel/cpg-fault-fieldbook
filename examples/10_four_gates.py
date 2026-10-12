"""Four different questions: validity, corridor, chain and physical rim."""

import numpy as np
from cpg2unstructured import morph_pillars, read_grdecl
from fault_projections import targets_from_polyline
from fieldbook.diagnostics import evaluate
from fieldbook.fixtures import quickstart

grid = read_grdecl(quickstart())
trace = np.array([[1., 1.], [3., 3.]])
targets = targets_from_polyline(grid, trace)
result = morph_pillars(
    grid, targets, boundary="pinned", smoothness="elastic",
    w0=0.05, polish="barrier", on_invalid="raise",
)
report = evaluate(grid, result.grid, trace, targets, min_sj=0.3)
assert report.ok, report
print(report)
