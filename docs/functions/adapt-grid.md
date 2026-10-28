# End-to-end adaptation: `adapt_grid` and `adapt_deck`

`fault_projections.adapt` connects fault geometry to `targets_from_*` and the
`cpg2unstructured` morph solver. `adapt_grid` works in memory; `adapt_deck`
reads a GRDECL and writes separate include files without rewriting the master
deck.

## `adapt_grid`

### Syntax

```python
adapt_grid(
    grid: CornerPointGrid,
    fault,
    *,
    regenerate_faults: bool = True,
    taper_cells: float = 3.0,
    slide_rim_ends: bool = False,
    top_trace=None,
    **morph_kw,
) -> AdaptResult
```

Accepted `fault` forms are dispatched by geometry:

| Form | Interpretation | Target behavior |
|---|---|---|
| `(N,2)` XY vertices, or a list of `(Ni,2)` lines | Plan-view polyline(s) | Rigid pillar translation from `targets_from_polylines` |
| `(N,3)` vertices or `Polygon3D` | Closed planar dipping polygon | Independent top/bottom targets from `targets_from_polygon` |
| List of `(4,3)` planar quads | Curved, faceted surface | Targets from individual facets |

A polygon or facets alone claim the pillars near the fault as a staircase
band, not an edge-connected chain, so they give no FAULTS. Pass
`top_trace=` to run both stages in one call. Give either the fault's `(N,2)`
trace on the top surface, or `"auto"` to derive it with
`surface_top_trace(grid, surface)`. Stage 1 snaps a pillar chain to the trace
and reads the FAULTS off it; stage 2 then tilts the pillars into the fault.
`result.stage1` is the trace stage alone, and `result.top_trace` the trace
used. Lesson 25 B–D and lesson 30 C use it.

`morph_kw` is forwarded to `morph_pillars` (for example, `boundary`,
`smoothness`, `support_rings`, and validity options). `taper_cells` applies to
polygon/facet target generation. `regenerate_faults` applies only to the
2-D-polyline path; 3-D surface input uses a separate 3-D surface/fault export
workflow.

### Output: `AdaptResult`

`grid` is the morphed grid; `morph` is the underlying `MorphResult`;
`distance_before` and `distance_after` are target residual arrays.
For 2-D traces with fault regeneration enabled, `faults_before` and
`faults_after` carry Eclipse face records and `faces_invariant` reports whether
their face sets agree. Only faces with active cells on both sides are
compared (`live_faces`). A trace through an inactive region therefore does not
count as a change: there, re-projection walks round the hole while the chain
runs through it. `grid_path` and `faults_path` are unset for in-memory
adaptation.

`ok` forwards the morph validity result. `summary()` reports distance changes,
morph details, and fault-face invariance when available. It does not check
cell-corridor containment or certify the deck for a simulator.

### Source-based usage

```python
import numpy as np
from fault_projections.adapt import adapt_grid

trace = np.array([[140.0, 0.0], [140.0, 800.0]])
result = adapt_grid(
    grid, trace,
    boundary="pinned",
    regenerate_faults=False,
    on_invalid="raise",
)
print(result.summary())
```

The in-memory adaptation tests use synthetic 8-by-8 grids with this trace.
Set `regenerate_faults=False` when only the grid adaptation is needed or when
full topology-based face regeneration is unsuitable. For the deformed
trace-cell strip, separately call
[`cell_corridor`](targets-corridors.md#cell_corridor).

## `adapt_deck`

### Syntax

```python
adapt_deck(
    deck: str | pathlib.Path,
    polylines,
    out_dir: str | pathlib.Path,
    *,
    prefix: str = "adapted",
    regenerate_faults: bool = True,
    fault_names=None,
    slide=None,
    slide_options=None,
    top_trace=None,
    **morph_kw,
) -> AdaptResult
```

This file workflow is intended for GRDECL text decks. `polylines` takes 2-D
polyline traces, or one dipping fault (a polygon or facets) with `top_trace=`
(see above). The dipping fault's FAULTS are then named by `fault_names`
(default `FAULT_1`).
It reads `deck`, calls `adapt_grid`, creates `out_dir` if needed, and writes:

- `<prefix>_grid.grdecl`, a grid include with the modified `COORD`;
- `<prefix>_faults.grdecl` when 2-D fault regeneration is enabled;
- `<prefix>_slide_report.txt` when `slide` is given.

`slide` (a `Throw` or a list, naming the faults in `fault_names` or the
mapping keys) slides the regenerated faults right after the morph, by ZCORN
only, before anything is written: the grid include then holds the adapted and
slid grid, `result.slide` the [`SlideResult`](fault-slide.md) and
`result.slide.before` the adapted grid alone. `slide_options` go to
`slide_grid`. It needs regenerated faults (2-D traces, or a surface with
`top_trace=`) and `regenerate_faults=True`; a failed slide raises
`SlideError` and writes nothing. Lesson 30 runs it and checks that it equals
`adapt_deck` followed by `slide_grid`. Its case C slides a 65° polygon in one
call, `adapt_deck(deck, polygon, out, top_trace="auto", fault_names="DIP",
slide=Throw("DIP", 15.0))`: the hanging wall loses within 1.2 % of the wedge
volume.

The master deck is not rewritten. Add an appropriate Eclipse `INCLUDE` to the
master deck to use the generated grid include; inspect `faults_path` and
`grid_path` on the returned `AdaptResult`.

```python
from fault_projections.adapt import adapt_deck

result = adapt_deck(
    "FIELD.GRDECL",
    trace,
    "adapted-output",
    prefix="case_a",
    regenerate_faults=False,
    boundary="pinned",
    on_invalid="raise",
)
print(result.grid_path)
```

This call pattern follows `tests/test_adapt.py`; the snippet has not been
executed for this documentation.

## Failure modes and limitations

- `adapt_grid` treats a trace that claims no pillars as a no-op with an empty
  distance array. Check `.summary()` and the corridor, rather than assuming a
  nonempty morph.
- When 2-D FAULTS regeneration is requested, the implementation rebuilds
  projected topology/graph state from the grid. The sibling wiki documents
  real-deck `build_topology` failures such as faces with more than two
  neighbors, as well as unsupported k-direction and complex fault-face cases.
  Turning off regeneration avoids that part; it does not repair grid topology.
- Adaptation always moves XY pillar positions; it does not perform a
  geomechanical solve, calculate transmissibility, or simulate flow.
- The 3-D planar and faceted adaptation forms do not use the 2-D FAULTS
  regeneration branch. Use the relevant 3-D surface exporter instead.
- Generated include files are not a complete rewritten simulator deck.
  Validate coordinate frames, original and adapted cell quality, trace
  containment, and simulator parsing independently.

## Maturity

The interface and 2-D/3-D dispatch have focused synthetic tests, including
unchanged master files, include round-tripping, face-set invariance, dipping
polygons, and faceted faults. Broad field-deck regeneration remains limited by
upstream topology support. This is an implemented workflow with scoped
evidence, not a claim of production readiness.

## Source and tests

- [`adapt_grid`, `adapt_deck`, and `AdaptResult`](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/fault_projections/adapt.py)
- [Adaptation tests](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/tests/test_adapt.py)
- [2-D corridor example](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/examples/adapt_corridor.py)
  · [Dipping-fault example](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/examples/adapt_dipping.py)
- [Morph solver](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/cpg2unstructured/morph.py)
- [Known field limitations](../legacy-field-guide.html)

