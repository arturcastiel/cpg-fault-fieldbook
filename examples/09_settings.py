"""Compare four solver choices on one identical fault and grid."""

import numpy as np
from cpg2unstructured import active_columns, cell_quality, check_validity, pillar_lattice, read_grdecl
from fault_projections import adapt_grid
from fieldbook.fixtures import quickstart

grid = read_grdecl(quickstart())
trace = np.array([[1., 1.], [3., 3.]])
choices = {
    "elastic": dict(w0=0.05, smoothness="elastic"),
    "laplacian": dict(w0=0.05, smoothness="laplacian"),
    "stiffer": dict(w0=2., smoothness="elastic"),
    "barrier": dict(w0=0.05, smoothness="elastic", polish="barrier"),
}
original = pillar_lattice(grid)
for name, settings in choices.items():
    result = adapt_grid(
        grid, trace, boundary="pinned", regenerate_faults=False,
        on_invalid="raise", **settings,
    )
    adapted = pillar_lattice(result.grid)
    top = cell_quality(adapted, plane="top")
    bottom = cell_quality(adapted, plane="bot")
    gate = check_validity(
        adapted, original, active=active_columns(grid), min_scaled_jacobian=0.3,
    )
    assert result.ok
    print(name, "residual", round(float(result.distance_after.mean()), 4),
          "SJ top/bot", round(top.sj_min, 3), round(bottom.sj_min, 3),
          "new cells below 0.3", gate.n_degenerate, "quality gate", gate.ok)
