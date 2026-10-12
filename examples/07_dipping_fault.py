"""Stage 2: a planar 3D polygon moves top and bottom XY independently."""

import numpy as np
from cpg2unstructured import pillar_lattice, read_grdecl
from fault_projections import adapt_grid
from fieldbook.fixtures import quickstart

grid = read_grdecl(quickstart())
polygon = np.array([
    [2.2, 0., 0.], [2.2, 5., 0.],
    [2.8, 5., 3.], [2.8, 0., 3.],
])
result = adapt_grid(
    grid, polygon, boundary="pinned", w0=0.05,
    regenerate_faults=False, on_invalid="raise",
)
before, after = pillar_lattice(grid), pillar_lattice(result.grid)
top_shift, bottom_shift = after.top - before.top, after.bot - before.bot
tilt = float(np.max(np.abs(top_shift - bottom_shift)))
assert result.ok and tilt > 0
assert result.faults_after is None  # 3D adaptation does not regenerate 2D FAULTS
print("largest differential XY displacement:", round(tilt, 4))
