# Topology and geometric checks

This page follows the grid-input pipeline:

```text
CornerPointGrid --from_grdecl--> CornerCube --build_topology--> UnstructuredGrid
       \                                                         /
        +---------------- pillar_lattice -----------------------+
                                  |
                        check_validity / cell_quality
```

## `from_grdecl`

### Syntax

```python
from_grdecl(g: Grdecl, pinch_tol: float = 0.0) -> CornerCube
```

`g` is the `CornerPointGrid` returned by **either** `read_grdecl` or
`read_egrid`. The legacy `Grdecl` name aliases the same type.

### Output

`CornerCube.nodes` is `(Nnodes, 3)` XYZ; `P` is `(2*NX, 2*NY, 2*NZ)` node
indices, with `-1` for slots belonging to inactive/pinched cells;
`cart_dims` is `(NX, NY, NZ)`; and `actnum_3d` is `(NX, NY, NZ)` Boolean.
`pinched_3d`, also `(NX, NY, NZ)`, marks collapsed cells detected by this
conversion. The corner coordinates are interpolated along each `COORD`
pillar at `ZCORN` depths.

`pinch_tol=0.0` drops only cells whose thickness is zero at all four pillar
corners. A positive tolerance also treats sufficiently thin collapsed cells as
pinches; those cells are removed from the active mask and carried forward so
the topology builder can form NNC bridges where applicable.

### Example

```python
from cpg2unstructured.cpg import from_grdecl
from cpg2unstructured.egrid_reader import read_egrid

parsed = read_egrid("FIELD.EGRID")
cube = from_grdecl(parsed, pinch_tol=0.0)
print(cube.nodes.shape, cube.P.shape, cube.actnum_3d.shape)
```

The EGRID-to-cube handoff is also demonstrated in
[the sibling EGRID example](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/examples/read_egrid_demo.py).

### Limits and maturity

Node deduplication quantizes XYZ to approximately `1e-9` coordinate units;
nearly coincident pillars can therefore be welded into shared nodes. Pinch
detection is geometric thickness detection and is not a full simulator
interpretation of all pinch-out policies. The conversion is implemented and
exercised by GRDECL, EGRID, pinch, and component tests.

## `build_topology`

### Syntax

```python
build_topology(cube: CornerCube) -> UnstructuredGrid
```

Builds unique i-, j-, and k-faces, cell-to-face incidence and active cell
indices. It removes inactive and zero-volume cells and remaps the remaining
cell IDs contiguously. `UnstructuredGrid` exposes:

| Field | Shape / meaning |
|---|---|
| `nodes` | `(Nnodes, 3)` XYZ |
| `faces.nodes` | `(Nfaces, 4)` node IDs |
| `faces.neighbours` | `(Nfaces, 2)` cell IDs; `-1` marks a boundary |
| `faces.direction` | `(Nfaces,)`, 0=i, 1=j, 2=k |
| `cell_faces` | `(Nactive, 6)` face IDs per cell |
| `cell_corners` | `(Nactive, 8)` node IDs per cell |
| `original_ijk` | `(Nactive, 3)` zero-based original `(i,j,k)` |
| `nnc` | `(M,2)` cell pairs for pinch bridges, otherwise empty |
| `component` | `(Nactive,)` connected-component labels |

This is a topology container, not automatically populated with centroids,
areas, or volumes. Those are separate geometry APIs.

### Example

```python
from cpg2unstructured.cpg import from_grdecl
from cpg2unstructured.grdecl_reader import read_grdecl
from cpg2unstructured.unstructured import build_topology

parsed = read_grdecl("tests/fixtures/simple_2x2x2.grdecl")
mesh = build_topology(from_grdecl(parsed))
print(mesh.n_cells, mesh.faces.nodes.shape, mesh.cell_faces.shape)
```

### Failure modes and maturity

- Non-matching-Z i/j faults are supported in general position (0.4.0):
  inactive cells and pinch-outs beside a fault get boundary faces, and a
  fault face where layer edges cross gets extra nodes (`faces.nodes` then
  has more than 4 columns; `mesh3d.face_edges` raises for such faces).
- Non-matching-Z k-direction faults are not implemented and raise
  `NonMatchingZError`.
- Face deduplication raises `RuntimeError` if it finds a face with more than
  two neighboring cells. The sibling wiki records this as a real-deck
  limitation; topology reconstruction must not be presumed to work on every
  field model.
- Synthetic Cartesian, pinch, and selected i/j-fault cases are covered by
  tests. The method is implemented and test-backed within those boundaries,
  not a universal Eclipse topology converter.

## `check_validity`

### Syntax

```python
check_validity(
    candidate: PillarLattice,
    original: PillarLattice,
    *,
    min_fraction: float = 0.10,
    min_scaled_jacobian: float | None = None,
    active: numpy.ndarray | None = None,
) -> ValidityReport
```

`candidate` and `original` are `PillarLattice` instances with matching
`(NX+1, NY+1, 2)` `top` and `bot` plan coordinates. The candidate is checked
against the original from which it was deformed.

Checks include endpoint-plane folds at **both** top and bottom, pillar welds,
adjacent-pillar separation, depth crossings, and non-adjacent proximity.
`min_fraction` is the minimum separation relative to original spacing for
quality checks. `active`, when supplied, is `(NX, NY)` Boolean and limits
findings to live columns and their pillars. Use `active_columns(grid)` to
build it. Without a mask the full lattice, including inactive regions, is
checked.

`min_scaled_jacobian` adds an optional degeneracy floor. It is `None` by
default, so a geometrically non-folded but nearly collapsed cell may still
pass. The check is a plan-view deformation gate, not a complete 3-D mesh,
reservoir simulator, or engineering qualification.

### Source-based use

```python
from cpg2unstructured.pillars import active_columns, pillar_lattice
from cpg2unstructured.validity import check_validity

before = pillar_lattice(grid)
after = pillar_lattice(candidate_grid)
report = check_validity(
    after, before,
    active=active_columns(grid),
    min_scaled_jacobian=0.3,
)
if not report.ok:
    print(report.summary())
```

`ValidityReport` includes `ok`, Boolean masks `folded_top` and `folded_bot`
(`(NX,NY)`), `welded` and `too_close` (`(NX+1,NY+1)`), counts, and an optional
`degenerate` mask. The active mask changes only which findings are counted; it
does not alter the underlying candidate geometry.

### Maturity

Implemented and directly exercised by validity-scope and morph tests. The
default checks crossings/welds/separation; the optional Jacobian threshold is
an explicit stricter policy that callers must select.

## `cell_quality`

### Syntax

```python
cell_quality(lattice: PillarLattice, *, plane: str = "top") -> QualityReport
```

`plane` is intended to be `"top"` or `"bot"`. The function evaluates the
selected endpoint plane's plan-view quads, each `(NX, NY)` in count. The
`QualityReport` summarizes:

- scaled Jacobian: min/mean/count below `0.3`;
- aspect ratio: maximum and 95th percentile;
- skew: maximum and mean;
- face non-orthogonality in degrees: maximum, 95th percentile, and mean;
- `n_cells`, the number of plan-view column quads.

It returns summary statistics, not the per-cell arrays. For array-level
metrics, call `scaled_jacobian`, `aspect_ratio`, `skew`, or
`non_orthogonality` directly from `cpg2unstructured.quality`.

### Example

```python
from cpg2unstructured.pillars import pillar_lattice
from cpg2unstructured.quality import cell_quality

lattice = pillar_lattice(grid)
top_report = cell_quality(lattice, plane="top")
bottom_report = cell_quality(lattice, plane="bot")
print(top_report.line())
```

Scaled Jacobian and skew measure plan-view shape; non-orthogonality is the
angle between neighboring cell-center direction and shared-face normal.
These are useful diagnostics, not a single pass/fail threshold for all
simulators. The function is implemented and backed by analytic-shape tests
(squares, stretched rectangles, shear, and inversion).

## Source and tests

- [`from_grdecl`](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/cpg2unstructured/cpg.py) ·
  [`build_topology`](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/cpg2unstructured/unstructured.py)
- [Topology/component tests](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/tests/test_components.py) ·
  [Faulted topology tests](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/tests/test_build_topology_faulted.py)
- [`check_validity`](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/cpg2unstructured/validity.py) ·
  [Validity-scope tests](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/tests/test_validity_scope.py)
- [`cell_quality`](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/cpg2unstructured/quality.py) ·
  [Quality tests](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/tests/test_quality.py)
- [Context and field limitations](../legacy-field-guide.html)

