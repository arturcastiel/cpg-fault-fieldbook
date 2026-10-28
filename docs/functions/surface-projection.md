# Three-dimensional fault projection

These routines produce cell-face surfaces on an existing
`cpg2unstructured.UnstructuredGrid`. They differ from plan-view trace targets:
they select discrete 3-D cell faces and do not morph the CPG pillars.

## `project_fault_surface`

### Syntax

```python
project_fault_surface(polygon, grid, graph=None) -> FaultSurface
```

`polygon` is either a `Polygon3D` or raw `(N,3)` XYZ vertices describing a
closed planar fault polygon. The grid is an `UnstructuredGrid`. `graph` is an
optional legacy argument; the current implementation does not need it.

The function fits the polygon plane, maps its footprint into plane coordinates,
and selects a minimal cell-face surface. It also derives the returned boundary
from the selected surface itself.

### Output

`FaultSurface` contains:

- `face_ids`: `(M,)` face IDs in `grid.faces`;
- `cell_ijk`: `(M,2,3)` structured indices of the two cells adjacent to each
  selected face;
- `direction`: `(M,)`, 0=i/X, 1=j/Y, 2=k/Z;
- `area`: total selected face area;
- `boundary_loop`: closed grid-node loop where a loop exists;
- `p0`, `normal`, `u`, `v`: plane origin, normal, and in-plane basis.

### Example

```python
import numpy as np
from fault_projections.surface import project_fault_surface

polygon = np.array([
    [1, 0, 0], [1, 3, 0], [1, 3, 2], [1, 0, 2],
], dtype=float)
surface = project_fault_surface(polygon, grid)
print(surface.face_ids.shape, surface.area, surface.boundary_loop.shape)
```

The values follow the axis-aligned rectangular fixture in
`tests/test_surface.py`; adapt polygon coordinates to the actual grid.

### Limits and maturity

The input is planar; this is not a general curved-surface projection. A
polygon/plane that misses the grid can yield an empty selection. Tested cases
cover axis-aligned, tilted, partial, and triangular planar inputs and verify
that returned boundary edges lie on the selected surface. Real-grid topology
and pinch-out structure still constrain the result.

## `project_quad_mesh_surface`

### Syntax

```python
project_quad_mesh_surface(
    mesh: QuadMesh,
    grid,
    graph=None,
    *,
    margin: float = 3.0,
    validate: bool = True,
    polish: bool = False,
) -> FaultSurface
```

`mesh` is a `QuadMesh` whose `vertices` are `(V,3)` and whose cyclic
`quads` are `(Q,4)` indices. `grid` is the rock-matrix `UnstructuredGrid`.
An optional `graph` avoids rebuilding the mesh graph. `margin` is passed to
the per-patch collar fill. `validate=True` runs strict manifold/orientation/
planarity checks; `polish=False` leaves the near-minimal collar fill as
assembled. `polish=True` runs a boundary-preserving area reducer.

The algorithm projects the shared mesh skeleton once, seeds each patch with
`project_fault_surface`, fills to the prescribed projected boundary, then
combines faces by symmetric difference. Shared patch boundaries therefore
cancel rather than creating independent seams.

The returned `FaultSurface.face_ids` are the combined surface and
`boundary_loop` is the projected mesh boundary. `area` is the total selected
area. For this combined mesh result, `cell_ijk` is empty `(0,2,3)`,
`direction` is empty `(0,)`, and plane fields `p0`, `normal`, `u`, `v` are
zero vectors because the surface has no single plane.

### Source-based example

```python
from fault_projections.quad_surface import project_quad_mesh_surface
from fault_projections.skeleton import QuadMesh

# `mesh` is a validated QuadMesh, `grid` an UnstructuredGrid.
surface = project_quad_mesh_surface(
    mesh, grid,
    validate=True,
    polish=False,
)
print(surface.face_ids.size, surface.area)
```

The full fixture construction and ridge geometry are in
[`tests/test_quad_surface.py`](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/tests/test_quad_surface.py);
that test suite demonstrates a sharp two-patch ridge with no interior holes.
This excerpt is source-based and has not been run as a standalone snippet.

### Failure modes and maturity

- With validation enabled, invalid meshes raise mesh validation exceptions
  (for example, non-planar, non-manifold, disconnected, or non-orientable
  inputs).
- If a patch boundary does not bound a cell-face surface (for example, an
  inactive-region/pinch-out obstruction), the routine raises
  `BoundaryDoesNotBound` with the patch index.
- `validate=False` bypasses input validation; it does not make malformed
  meshes safe.
- The implementation is backed by synthetic Cartesian tests for planar
  equivalence, ridge watertightness, closed boundary, and validation. Broader
  real-grid and pinch-out coverage remains limited.

## Source and tests

- [`project_fault_surface`](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/fault_projections/surface.py)
  · [Planar surface tests](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/tests/test_surface.py)
- [`project_quad_mesh_surface`](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/fault_projections/quad_surface.py)
- [Quad mesh input and validation](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/fault_projections/skeleton.py)
  · [Validation implementation](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/fault_projections/validate.py)
- [Watertight quad-surface tests](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/tests/test_quad_surface.py)
- [Related `cpg2python` surface primitives](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/cpg2unstructured/minimal_surface.py)
- [Fieldbook limitations](../legacy-field-guide.html)

