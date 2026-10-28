# Trace targets and cell corridors

The target functions in `fault_projections.targets` select pillars and compute
XY targets for the `cpg2unstructured` morph solver. `cell_corridor` provides a
separate cell-based containment/connectivity check; it does not generate
targets or write Eclipse `FAULTS`.

## `targets_from_polyline`

### Syntax

```python
targets_from_polyline(
    grid,
    polyline,
    *,
    radius: float | None = None,
    plane: str = "top",
    dip_targets=None,
    consolidate: bool = True,
    tol: float = 1e-9,
) -> PillarTargets
```

`grid` is a `CornerPointGrid`; `polyline` is an `(N,2)` array-like of model XY
vertices (`N >= 2`). Coordinates must use the grid's coordinate frame.
`plane="top"` (default) selects from the top endpoint lattice; `"bot"` selects
from the bottom.

By default, each segment claims pillars at crossed lattice edges and
near-incident vertices. The chosen endpoint of each crossed edge is the nearer
pillar, avoiding a zero-length collapsed edge. The result consists of ordered
pillar/segment incidences connected into a chain; this allows fault-face
construction to rely on both endpoints of each chain segment being targeted.

### Options and output

- `radius=None`: topological edge-crossing selection. A numeric `radius`
  instead selects pillars within that distance of each segment, widening the
  target band (units are grid-coordinate units).
- `consolidate=True`: collapse a pillar's competing segment targets to a
  nearest target. `False` retains one row per incidence; kink pillars can then
  appear more than once and retain multiple residuals.
- `dip_targets`: optional callable `(ij, xy_top_target) -> xy_bot_target`.
  When present, supplies per-row bottom targets and enables independent
  endpoint movement. With `None`, targets are rigid translations.
- `tol=1e-9`: geometric tolerance used for segment/lattice-edge crossing.

The returned `PillarTargets` has `(M,2)` zero-based `ij`, `(M,2)` `xy_top`,
optional `(M,2)` `xy_bot`, and optional `(M,)` weights. If the trace misses the
grid, the result can be empty (`n == 0`) rather than an error.

### Example

```python
import numpy as np
from fault_projections.targets import targets_from_polyline

trace = np.array([[140.0, 0.0], [140.0, 800.0]])
targets = targets_from_polyline(grid, trace)
print(targets.n, targets.ij.shape, targets.xy_top.shape, targets.is_rigid)

raw_incidences = targets_from_polyline(grid, trace, consolidate=False)
wider_band = targets_from_polyline(grid, trace, radius=150.0)
```

These call patterns follow `tests/test_targets.py`; grid extents and units must
match the supplied trace. The numbers above are from the source test's
synthetic 8-by-8 grid, not a result claimed for an arbitrary model.

## `targets_from_polylines`

### Syntax

```python
targets_from_polylines(grid, polylines, **kw) -> PillarTargets
```

Accepts one `(N,2)` polyline or an iterable of `(Ni,2)` polylines. Keyword
arguments are forwarded to `targets_from_polyline`. Targets from all traces
are merged into a single system so shared pillars carry all of their target
incidences rather than applying faults sequentially. Empty trace sets or
traces that claim no pillars produce an empty `PillarTargets`.

```python
from fault_projections.targets import targets_from_polylines

traces = [
    np.array([[140.0, 0.0], [140.0, 800.0]]),
    np.array([[0.0, 260.0], [800.0, 260.0]]),
]
combined = targets_from_polylines(grid, traces)
```

This input form follows the tests. For lists of coordinate arrays, use a list
whose each item is a 2-D `(N,2)` array; a single trace can be passed directly.

## `cell_corridor`

### Syntax

```python
cell_corridor(
    grid,
    line: numpy.ndarray,
    *,
    plane: str = "top",
    per_cell: float = 8.0,
) -> CellCorridor
```

`line` is an `(N,2)` or `(N,3)` polyline; Z is ignored and the XY coordinates
are tested against the grid plan view. `grid` may be an original or morphed
`CornerPointGrid`. `plane` chooses the lattice endpoint used for locating
columns; `per_cell` controls densification, with 8 samples per local cell-size
as the default. Finer values add more samples.

The result contains:

| Member | Shape / meaning |
|---|---|
| `cells` | `(M,2)` zero-based `(i,j)` columns in trace order |
| `interfaces` | `M-1` `FaultFace` records between adjacent corridor cells; Eclipse face indices are 1-based |
| `samples` | densified XY samples |
| `located` | `(S,2)` cell indices per sample; `(-1,-1)` means off-grid |
| `inserted`, `jumps` | repair-inserted columns and suspected sampling skips |

`ok` requires face-connected cells, trace containment, and a resolved sample
walk (`jumps == 0`). Use `n_outside` to detect trace samples outside the grid;
it is possible for a corridor to be connected while part of the trace leaves
the grid. `.mask(cart_dims)` builds an `(NX,NY)` Boolean selection.

### Post-morph check

```python
from fault_projections.corridor import cell_corridor

before = cell_corridor(grid, trace)
after = cell_corridor(adapted_grid, trace)  # check the deformed grid
if not after.ok:
    print(after.summary())
```

`AdaptResult.corridor(trace)` is a convenience for checking the result's grid.
This is a plan-view fault strip check, not a 3-D fault-surface or simulator
validation.

## Limits and maturity

Target selection and corridor construction are implemented with focused
synthetic tests, including misses, kinks, multiple traces, edge-boundary
vertices, undersampling, and checks after adaptation. However, target-chain
continuity does not itself prove cell-face corridor containment: use both
target/result inspection and `cell_corridor` as separate checks. The corridor
classifies plan-view columns and interfaces; it does not account for layerwise
fault-surface topology.

## Source and tests

- [`targets_from_polyline` / `targets_from_polylines`](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/fault_projections/targets.py)
- [Polyline input normalization](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/fault_projections/polylines.py)
- [Target tests](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/tests/test_targets.py)
- [`cell_corridor`](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/fault_projections/corridor.py)
  · [Corridor tests](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/tests/test_corridor.py)
- [Adaptation example](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/examples/adapt_corridor.py)
- [Fieldbook limitations and recommendations](../legacy-field-guide.html)

