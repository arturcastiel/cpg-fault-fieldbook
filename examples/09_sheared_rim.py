"""Show that 'slide' can move a sheared physical outline despite ok=True."""

from dataclasses import replace

import numpy as np
from cpg2unstructured import pillar_lattice, read_grdecl
from fault_projections import adapt_grid
from fieldbook.fixtures import quickstart

base = read_grdecl(quickstart())
coord = base.coord.copy()
coord[:, 0] += 0.5 * coord[:, 1]
coord[:, 3] += 0.5 * coord[:, 4]
grid = replace(base, coord=coord)
trace = np.array([[1.1, 2.2], [6.7, 3.4]])
errors = {}
for boundary in ("pinned", "slide", "free"):
    result = adapt_grid(
        grid, trace, boundary=boundary, w0=0.05,
        regenerate_faults=False, on_invalid="raise",
    )
    xy = pillar_lattice(result.grid).top
    error = max(
        float(np.max(np.abs(xy[0, :, 0] - .5 * xy[0, :, 1])) / np.sqrt(1.25)),
        float(np.max(np.abs(xy[-1, :, 0] - .5 * xy[-1, :, 1] - 5)) / np.sqrt(1.25)),
        float(np.max(np.abs(xy[:, 0, 1]))),
        float(np.max(np.abs(xy[:, -1, 1] - 5))),
    )
    assert result.ok
    errors[boundary] = error
    print(boundary, "physical outline error:", round(error, 4))
assert errors["pinned"] < 1e-9
assert errors["slide"] > 0.1  # defect in the pinned upstream revision
