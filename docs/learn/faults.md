# 2. Project fault geometry without moving the grid

**2.5D:** build explicit topology, project it into XY, snap an authored
polyline, stitch projected nodes along graph edges and emit FAULTS
through all layers. Run `python -m fieldbook.cli run 03_project_trace`.

{{ example:03_project_trace.py }}

The 2D result is a staircase on the **original** grid. `FAULTS` records
are one-based Eclipse face references; the NumPy inputs are physical
XY coordinates in the grid's current frame. Projection does not
deform cells or compute transmissibility multipliers.

**Planar 3D:** project a coplanar polygon onto the explicit cell-face
mesh, with X/Y/Z and real K ranges. Run
`python -m fieldbook.cli run 04_project_surface`.

{{ example:04_project_surface.py }}

**Coherent quad surface:** the upstream `quad_mesh_fault.py` builds
two quads that share a ridge and uses one common projected skeleton.
Run `python -m fieldbook.cli run 12_quad_surface --out output/quad`.
This needs a sibling source checkout; it writes FAULTS and VTK only
to the requested directory. Independently projected quads can tear
at the seam. This is a demonstration, not a guarantee for every
faulted or pinched field grid.

{{ example:12_quad_surface.py }}
