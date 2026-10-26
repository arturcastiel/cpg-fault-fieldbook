"""Stage 2 with two planar facets meeting along one depth-varying seam."""

import numpy as np
from cpg2unstructured import pillar_lattice, read_grdecl
from fault_projections import adapt_grid
from fault_projections.targets import targets_from_facets
from fieldbook.fixtures import quickstart

grid = read_grdecl(quickstart())
facets = [
    np.array([[2.2, 0., 0.], [2.2, 2.5, 0.],
              [2.8, 2.5, 3.], [2.8, 0., 3.]]),
    np.array([[2.2, 2.5, 0.], [2.4, 5., 0.],
              [3.0, 5., 3.], [2.8, 2.5, 3.]]),
]
targets = targets_from_facets(grid, facets)
assert targets.n > 0 and not targets.is_rigid
result = adapt_grid(
    grid, facets, boundary="pinned", w0=0.05,
    regenerate_faults=False, on_invalid="raise",
)
assert result.ok and result.faults_after is None
original, adapted = pillar_lattice(grid), pillar_lattice(result.grid)
tilt = np.max(np.abs((adapted.top - original.top) - (adapted.bot - original.bot)))
assert tilt > 0
print("facets:", len(facets), "target rows:", targets.n,
      "mean residual:", round(float(result.distance_after.mean()), 4))
print("3D adaptation did not export FAULTS.")
