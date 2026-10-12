"""Use the same targets for morphing and independent post-morph checks."""

import numpy as np
from cpg2unstructured import morph_pillars, read_grdecl
from fault_projections import targets_from_polyline
from fieldbook.diagnostics import evaluate
from fieldbook.fixtures import quickstart

original = read_grdecl(quickstart())
trace = np.array([[1., 1.], [3., 3.]])
targets = targets_from_polyline(original, trace)
assert targets.n > 0 and targets.is_rigid
result = morph_pillars(
    original, targets, boundary="pinned", smoothness="elastic",
    w0=0.05, polish="barrier", on_invalid="raise",
)
gates = evaluate(original, result.grid, trace, targets, min_sj=0.3)
assert result.ok and gates.ok, (result.summary(), gates)
print("targets:", targets.n, result.summary())
print("gates:", gates)
