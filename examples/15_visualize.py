"""Generate identical-scale before/after SVGs; inspect the numerical gates too."""

import argparse
from pathlib import Path

import numpy as np
from cpg2unstructured import morph_pillars, read_grdecl
from fault_projections import targets_from_polyline
from fieldbook.diagnostics import evaluate
from fieldbook.fixtures import quickstart
from fieldbook.visualize import render_pair

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", required=True, type=Path)
args = parser.parse_args()
grid = read_grdecl(quickstart())
trace = np.array([[1., 1.], [3., 3.]])
targets = targets_from_polyline(grid, trace)
result = morph_pillars(
    grid, targets, boundary="pinned", w0=0.05,
    polish="barrier", on_invalid="raise",
)
gates = evaluate(grid, result.grid, trace, targets)
assert gates.ok, gates
before, after = render_pair(grid, result.grid, trace, args.out)
print("before:", before, "after:", after)
print("gates:", gates)
