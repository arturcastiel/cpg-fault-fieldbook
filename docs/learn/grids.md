# 1. Read and navigate a corner-point grid

**What to expect:** the bundled GRDECL has 5 × 5 × 3 logical cells, 70
active cells, 36 pillars and a `PORO` property. Its XY coordinates span
0–5 in arbitrary units. Run `python -m fieldbook.cli run 01_read_grid`.

{{ example:01_read_grid.py }}

`read_grdecl` returns a `CornerPointGrid`. `from_grdecl` turns it into a
corner cube; only then can `build_topology` build explicit mesh faces
and cell adjacency. `Grdecl` is a type alias for `CornerPointGrid`,
so this pivot also accepts `read_egrid` results. `read_egrid` applies
MAPAXES by default; when writing transformed coordinates, use the
writer's explicit MAPAXES override **only** if the simulator coordinate
frame was checked.

Run `python -m fieldbook.cli run 02_egrid` to compare the bundled binary
MAPAXES fixture in raw-local and transformed-world frames:

{{ example:02_egrid.py }}

The writer refuses to silently export MAPAXES-transformed world
coordinates into a deck expecting local coordinates. Run
`python -m fieldbook.cli run 17_mapaxes_guard`; it asserts the
specific error and leaves no output file:

{{ example:17_mapaxes_guard.py }}

Run `python -m fieldbook.cli run 03_mesh_graph`:

{{ example:03_mesh_graph.py }}

The `Faces` attribute is spelled `neighbours`. Topology construction
may fail on complex field decks, including the recorded Norne case;
pillar adaptation with `regenerate_faults=False` does not require
full mesh topology.
