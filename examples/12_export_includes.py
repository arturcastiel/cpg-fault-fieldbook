"""Write a grid INCLUDE only to an explicitly requested output directory."""

import argparse
from pathlib import Path

import numpy as np
from fault_projections import adapt_deck
from fieldbook.fixtures import quickstart

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
result = adapt_deck(
    quickstart(), np.array([[1., 1.], [3., 3.]]), args.out,
    prefix="fault_adapted", boundary="pinned", polish="barrier",
    regenerate_faults=False, on_invalid="raise",
)
assert result.ok and result.grid_path is not None
assert result.grid_path.is_file() and result.faults_path is None
print("grid INCLUDE:", result.grid_path)
print("FAULTS INCLUDE:", result.faults_path)
print("Original deck was not rewritten; add the new INCLUDE to a reviewed copy.")
