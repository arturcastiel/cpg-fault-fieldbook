"""Create a traceable synthetic case report, explicitly NOT simulator approval."""

import argparse
import json
from pathlib import Path

import numpy as np
from cpg2unstructured import morph_pillars, read_grdecl
from fault_projections import targets_from_polyline
from fieldbook.diagnostics import evaluate
from fieldbook.fixtures import quickstart
from fieldbook.provenance import case_record

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
source = quickstart()
grid = read_grdecl(source)
trace = np.array([[1., 1.], [3., 3.]])
settings = {
    "boundary": "pinned", "smoothness": "elastic", "w0": 0.05,
    "polish": "barrier", "on_invalid": "raise",
}
targets = targets_from_polyline(grid, trace)
result = morph_pillars(grid, targets, **settings)
gates = evaluate(grid, result.grid, trace, targets, min_sj=0.3)
record = case_record(source, trace, settings, result, gates)
args.out.mkdir(parents=True, exist_ok=True)
path = args.out / "one-trace.json"
path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print("case:", path, "verdict:", record["verdict"])
assert gates.ok, gates
