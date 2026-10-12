"""Crossing faults: duplicate pillar claims survive in a simultaneous solve."""

import numpy as np
from cpg2unstructured import read_grdecl
from fault_projections import adapt_grid, targets_from_polylines
from fieldbook.fixtures import quickstart

grid = read_grdecl(quickstart())
traces = [
    np.array([[1., 1.], [3., 3.]]),
    np.array([[1., 3.], [3., 1.]]),
]
targets = targets_from_polylines(grid, traces)
assert targets.n == 10
assert len(set(map(tuple, targets.ij))) == 8
result = adapt_grid(
    grid, traces, boundary="pinned", w0=0.05,
    regenerate_faults=False, on_invalid="raise",
)
assert result.ok and result.distance_after.shape == (targets.n,)
print("claim rows:", targets.n, "unique pillars:", len(set(map(tuple, targets.ij))))
print(result.summary())
