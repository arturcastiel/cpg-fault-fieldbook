"""Read-only intake of a user deck before considering any adaptation."""

import argparse
from pathlib import Path

from cpg2unstructured import active_columns, cell_quality, pillar_lattice, read_egrid, read_grdecl

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--deck", required=True, type=Path)
args = parser.parse_args()
if not args.deck.is_file():
    parser.error(f"deck does not exist: {args.deck}")
grid = (read_egrid(args.deck) if args.deck.suffix.upper() in {".EGRID", ".FEGRID"}
        else read_grdecl(args.deck))
lat = pillar_lattice(grid)
active = active_columns(grid)
print("source:", grid.source, "logical cells:", grid.cart_dims)
print("active cells:", int(grid.actnum.sum()), "/", grid.n_cells,
      "active columns:", int(active.sum()), "/", active.size)
print("coordinate frame: MAPAXES present:", grid.mapaxes is not None,
      "already applied:", grid.mapaxes_applied)
print("top quality (all columns):", cell_quality(lat, plane="top").line())
print("bottom quality (all columns):", cell_quality(lat, plane="bot").line())
print("parsed properties:", sorted(grid.props))
print("unrecognized keywords:", grid.unrecognized)
print("NO topology, adaptation, output deck, or simulator check was performed.")
