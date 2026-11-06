# Workflow 2 · Full 3D projection

Put a 3D fault surface on the grid's cell faces, in X, Y **and** Z,
with real K ranges, without moving the grid. A dipping fault becomes
a staircase in depth as well as in plan.

![Workflow 2: full 3D projection, from a fault surface to 3D FAULTS](../assets/workflows/flow-2-full-3d.png)

## Steps

```python
import numpy as np
from cpg2unstructured import build_topology, from_grdecl, read_grdecl
from fault_projections import QuadMesh, faults_text_3d, project_fault_surface, write_faults_3d

mesh = build_topology(from_grdecl(read_grdecl("FIELD.GRDECL")))

# a planar polygon: (N, 3) coplanar XYZ; here a 45 degree plane
polygon = np.array([[310.0, -10.0, 1995.0], [310.0, 410.0, 1995.0],
                    [410.0, 410.0, 2095.0], [410.0, -10.0, 2095.0]])
surface = project_fault_surface(polygon, mesh)
# planar facets:  surface = project_faceted_surface(quads, mesh)   # list of (4, 3) quads
# a quad mesh:    surface = project_quad_mesh_surface(QuadMesh(verts, quads), mesh)

surface.face_ids                           # the selected cell faces (any input)
surface.cell_ijk, surface.direction        # polygon/facets: cell pairs, 0 = X, 1 = Y, 2 = Z
                                           # (empty for a quad mesh)
write_faults_3d("faults3d.grdecl", {"F1": surface}, mesh)   # faults_text_3d returns the text
```

1. **Build the cell-face mesh** of the grid
    ([lesson 03](../learn/grids.md)).
2. **Choose the call for your surface:** one planar polygon
    (lesson 05), planar facets (lesson 25 C) or a quad mesh with shared
    seams (lessons 13, 25 D). Sticks or a point cloud: fit a plane
    (`fit_plane`, lesson 25 B) or build facets (`facet_surface`,
    lesson 25 C) first.
3. **Write 3D `FAULTS`.** Z faces are a layer's top or base.

## Before and after

![Full 3D projection: a 45 degree plane, its 2D trace projection, and its 3D projection](../assets/workflows/ill-2-full-3d.png)

The same 45° plane, two ways. Projecting only its top trace
([workflow 1](stair-step.md)) gives 8 vertical faces through every
layer, at x = 300 m all the way down. Projecting the surface gives
48 X faces and 16 Z faces (0 Y faces; the plane strikes north), and
the staircase follows the dip. The grid is the same in all three
panels.

## Run the lessons

| Lesson | What it shows | Run | Read |
| --- | --- | --- | --- |
| 05 | a planar polygon onto X/Y/Z faces | `python -m fieldbook.cli run 05_project_surface` | [Fault projection](../learn/faults.md) |
| 13 | a quad mesh with a shared ridge (needs a sibling source checkout) | `python -m fieldbook.cli run 13_quad_surface --out output/quad` | [Fault projection](../learn/faults.md) |
| 25 B, C, D | sticks → `fit_plane`; a point cloud → facets; a kinked quad mesh | `python -m fieldbook.cli run 25_gallery_geometry_inputs` | [Gallery 3](../learn/gallery.md#3-from-fault-data-to-a-slid-grid-lesson-25) |

Source: [05_project_surface.py](https://github.com/arturcastiel/cpg-fault-fieldbook/blob/main/examples/05_project_surface.py),
[13_quad_surface.py](https://github.com/arturcastiel/cpg-fault-fieldbook/blob/main/examples/13_quad_surface.py),
[25_gallery_geometry_inputs.py](https://github.com/arturcastiel/cpg-fault-fieldbook/blob/main/examples/25_gallery_geometry_inputs.py).
Method page: [Surface projection](../functions/surface-projection.md).

![From fault data to a slid grid: arc, sticks, point cloud, quad mesh, dome, dipping layers](../assets/gallery/geometry-inputs.png)

## Checks and limits

- The grid is unchanged: a dipping fault is a staircase in depth.
  To make the cells follow it, tilt the pillars: stage 2 of
  [workflow 3](deformation.md), or [workflow 6](dip-angle.md).
- Independently projected quads can tear at a seam; a `QuadMesh`
  with shared vertices avoids it
  ([Surface projection](../functions/surface-projection.md#project_quad_mesh_surface)).
- See [Limits](../reference/limits.md).
