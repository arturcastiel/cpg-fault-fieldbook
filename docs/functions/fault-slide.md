# Fault slides in Z: `slide_grid`, `slide_deck` and `zcorn_ops`

A fault slide gives existing `FAULTS` a throw by editing ZCORN only. Every
corner moves along its own pillar, so COORD never changes and only faces listed
in `FAULTS` open. `fault_projections.slide` decides the shifts;
`cpg2unstructured.zcorn_ops` applies and checks them.

## `slide_grid`

### Syntax

```python
slide_grid(
    grid: CornerPointGrid,
    throws: Throw | Sequence[Throw],
    *,
    faults=None,             # every face that may open: {name: faces}; default the throws' faces
    method: str = "kernel",  # far field of faults with tips: "kernel" or "laplacian"
    on_invalid: str = "raise",   # "raise", "bisect" or "ignore"
    max_bisect: int = 20,
    report: bool = True,     # re-process the grid and fill the report
    tol: float = 0.0,
) -> SlideResult
```

`Throw(fault, throw, hanging_side=None, hanging=1.0, profile="cosine",
taper=None, decay="auto", growth=None, snap=False, name="")` describes one fault:

| Field | Meaning |
|---|---|
| `fault` | A name in `faults`, or the fault's faces (one branch-free chain; a ring is allowed) |
| `throw` | Positive = hanging wall down. A peak value that tapers to zero at each tip, an `(n, 2)` table of `(arc length, throw)`, or one value per chain pillar |
| `hanging_side` | A plan point on the hanging wall, `"left"`/`"right"`, or `None` to read it from the pillar dip |
| `hanging` | Share of the throw taken by the hanging wall; the footwall moves up by the rest |
| `decay` | `"auto"`: a fault that separates the grid slides a block rigidly; one that ends inside the grid fades with distance |
| `growth` | `NZ + 1` factors, one per horizon: a growth fault (layers change thickness) |
| `snap` | Lower the throw at pillars so layer boundaries never cross inside a face |

Several throws add up: sliding two faults at once gives the same grid as
sliding them one after the other.

### Output: `SlideResult`

`grid` is the slid grid and `before` the input. `ok` is True when no cell
inverted, no layers overlap, no face outside `FAULTS` opened (`unexpected`) and
the grid re-processes. `modes` gives `"rigid"`, `"decay"` or `"block"` per fault;
`realized_throw(name)` the throw along the chain; `scale` the fraction of the
throws applied after bisection. With `report=True`, `juxtaposition`,
`connections` (cell pairs across faults, with areas) and `volume_change` are
filled. `summary()` gives one line, `report()` a full text report.

### Example

```python
import numpy as np
from cpg2unstructured import CornerPointGrid
from fault_projections import Throw, slide_grid
from fault_projections.faults import FaultFace

nx, ny, nz, dx, dz = 6, 4, 3, 100.0, 10.0
coord = [(i * dx, j * dx, 0.0, i * dx, j * dx, 100.0)
         for j in range(ny + 1) for i in range(nx + 1)]
depth = np.repeat(np.arange(nz) * dz, 2) + np.tile([0.0, dz], nz)
zcorn = np.broadcast_to(depth[:, None, None], (2 * nz, 2 * ny, 2 * nx)).ravel()
grid = CornerPointGrid((nx, ny, nz), np.array(coord), zcorn.astype(float),
                       np.ones(nx * ny * nz, np.int8))

f1 = [FaultFace(3, j, "X") for j in range(1, ny + 1)]   # rim to rim
res = slide_grid(grid, Throw(f1, 15.0, hanging_side=(500.0, 200.0), name="F1"))
print(res.summary())
# OK: 1 fault(s) slid (F1 rigid) · 4 face(s) opened · max |volume change| 0 %
```

This snippet was run against the pinned revisions.

## `slide_deck`

```python
slide_deck(deck, throws, out_dir, *, prefix="slid", faults=None, vtk=False,
           png=False, **slide_kw) -> SlideResult
```

Reads a deck, slides it and writes `<prefix>_grid.grdecl` (COORD, ZCORN,
ACTNUM) and `<prefix>_report.txt`. `faults` defaults to the deck's own `FAULTS`
(INCLUDE is followed), so a throw can name a fault: `Throw("F1", 20.0)`. The
`FAULTS` set does not change, so replace the grid in the deck with
`INCLUDE '<prefix>_grid.grdecl' /`. `png=True` also writes `<prefix>_slide.png`:
a plan map of the top shift with the fault chains, and the realized throw
along each fault (`viz.plot_slide(res)`, `write_slide_png(res, path)`; needs
matplotlib).

Lesson 30 slides Norne's E_01 (5 m, fading to its tips) and BC (3 m, rigid)
through a deck: 8 faces open, none outside `FAULTS`, a round trip returns to
5e-13 m and cell volumes change by up to 2.4 % on the tilted pillars. To adapt
and slide in one call, use [`adapt_deck(slide=...)`](adapt-grid.md#adapt_deck);
case C does it for a dipping polygon with `top_trace="auto"`. With matplotlib
installed, cases B and C also write a PNG.

{{ example:30_deck_slide_workflows.py }}

## `zcorn_ops` (cpg2unstructured)

| Function | What it does |
|---|---|
| `zcorn_cube(grid)` / `with_zcorn(grid, cube)` | ZCORN as a `(2NX, 2NY, 2NZ)` cube of corner slots, and back |
| `apply_zcorn_displacement(grid, dz)` | Shift corners along their pillars: a scalar, per column corner `(2NX, 2NY)`, per horizon `(2NX, 2NY, NZ+1)` or per slot |
| `unwelded_faces(grid, tol=0.0)` | `(open_i, open_j)` masks of lateral faces whose corners no longer match |
| `zcorn_validity(grid, tol=0.0)` | `ZcornValidity`: inverted cells and overlapping layers |

## Failure modes and limitations

- Faults must lie on pillar lines (a deck's `FAULTS`, or the faces
  `adapt_grid` regenerates). k-direction faults are not supported.
- A throw that fails the gate raises `SlideError`, unless
  `on_invalid="bisect"` or `"ignore"`.
- On tilted pillars, cells next to a moved pillar change volume by the fault
  wedge; `volume_change` reports it.
- A fault that ends inside the grid fades by a distance kernel or a harmonic
  solve, not by a geomechanical model.

## Source and tests

- [`slide_grid`, `slide_deck`, `Throw`, `SlideResult`](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/fault_projections/slide.py)
  · [Slide tests](https://github.com/arturcastiel/fault-projections/blob/247507e70be854bbc18d6b02c27efeb59f998149/tests/test_slide.py)
- [`zcorn_ops`](https://github.com/arturcastiel/cpg2python/blob/7d98aebbd64b03ee4623cf0c97e01dd19cda9f5b/cpg2unstructured/zcorn_ops.py)
- Lessons 19–30 ([Z slide](../learn/z-slide.md), [gallery](../learn/gallery.md))
