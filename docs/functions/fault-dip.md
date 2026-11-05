# Fault dip: `dip_grid`, `dip_deck` and `Dip`

A fault dip sets the angle of existing `FAULTS` by tilting pillars, so only
COORD changes. Every corner keeps its depth and slides sideways with its
pillar; ZCORN, ACTNUM and the `FAULTS` face set never change. Apply a dip
before a slide: `dip_deck(..., slide=...)` does both.

The lesson with every case worked through is
[7. Set a fault's dip angle](../learn/fault-dip.md).

## `dip_grid`

### Syntax

```python
dip_grid(
    grid: CornerPointGrid,
    dips: Dip | Sequence[Dip],
    *,
    faults=None,               # every fault face: {name: faces}; default the dips' faces
    on_invalid: str = "raise", # "raise", "bisect" or "ignore"
    widen: int = 4,            # how many times the fades may be doubled (0: never)
    max_bisect: int = 20,
    min_jacobian: float = 0.05,  # gate floor (applied only where a cell got worse)
    min_fraction: float = 0.10,  # gate floor on edge length, fraction of the original
    passes: int = 6,           # rounds for several dips to hold their angles together
    report: bool = True,       # re-process the grid (part of the gate), fill volumes
) -> DipResult
```

`faults` decides which faults separate the grid, which in turn decides each
dip's mode. Give all of a deck's faults when the dipped fault abuts others.

`Dip(fault, angle=60.0, hanging_side=None, pivot="middle", decay="auto",
footwall="auto", smoothing=None, name="")` describes one fault:

| Field | Meaning |
|---|---|
| `fault` | A name in `faults`, or the fault's faces (one branch-free chain; a ring is allowed) |
| `angle` | Degrees below horizontal, in (0, 90]: one value, an `(n, 2)` table of `(arc length, angle)` along the oriented chain, or one value per chain pillar |
| `hanging_side` | A plan point on the hanging wall, `"left"`/`"right"` of the chain, or `None` to read it from the current lean (a vertical fault needs it) |
| `pivot` | The depth that keeps its plan position: `"top"`, `"middle"` (default) or `"bottom"` of the cells on each fault pillar, or a depth |
| `decay` | `"auto"`: rigid for a fault that separates the grid, else a fade over `max(L/2, π·reach, 2h)`. `None`: rigid (an error for a fault with a tip). A distance: decay mode with that fade, even for a separating fault |
| `footwall` | Fade of the footwall. `"auto"`: `max(π·reach, 2h)` in rigid mode, the `decay` distance in decay mode. `None`/`0`: no fade (rigid only) |
| `smoothing` | Gaussian width over which the strike is averaged along the chain. Default: twice the median face length |
| `name` | Name for a fault given as faces |

Here `L` is the fault length, `h` the median face length, and `reach` the
largest corner shift on the fault at the full dip.

### What moves

At each pillar, the change is a horizontal shift that grows linearly with
depth away from the pivot depth `z0`. A corner at depth `z` moves by

```text
Δxy = w · (z − z0) · k
```

- **On the fault**, `w = 1` and `k` sets the angle exactly. With `s` the
  pillar's present lean per unit depth, `t` the smoothed strike and `n` the
  dip direction (toward the hanging wall), the new lean is
  `(s·t) t + cot(angle) n`. The along-strike component of the lean is kept;
  the component along the dip direction is set to `cot(angle)`. So
  `k = (s·t) t + cot(angle) n − s`.
- **Off the fault**, `k` and `z0` are spread from the chain pillars by the
  weights `ℓ_m / (r² + σ²)`, with `ℓ_m` the arc length each chain pillar
  stands for and `σ` the smoothing. Curved and stair-stepped faults thus give
  a smooth field.
- **The weight `w`** fades with a cosine ramp:
  - in **rigid** mode, `w = 1` on every pillar of the hanging wall (the flood
    fill of the side the fault separates, *including the rim and any blocks
    beyond other faults*), fading to 0 over `footwall` into the footwall;
  - in **decay** mode, `w` fades with distance from the fault on both sides:
    over `decay` on the hanging side and over `footwall` on the footwall.
- **With several dips**, each dip's `w` is pinned to 0 at the other dipped
  faults, locally (within the other fault's footwall fade). The dips are
  applied in turn, and repeated up to `passes` times until every fault holds
  its angle to 1e-9.

The `"auto"` distances keep the strain in the fade zone at or below one half:
`π` times the largest corner shift, and at least two cells.

**Volumes.** In rigid mode on a straight fault, the hanging wall is sheared
parallel to the fault, so its cells keep their volumes (lesson 31 measures
1e-14). With the middle pivot, the top of a column moves opposite to its base
by the same amount, so every column keeps its **total** volume. Volume still
moves between the layers of a column (case A: −21.65 % to +17.32 % in one
footwall column). The other pivots change column volumes.

### The gate

The dipped grid is checked where cells are. At the top, middle and base of
every active cell, the plan-view corner angles (scaled Jacobian, signed
against the original orientation) and edge lengths are compared with the
original. A cell fails if:

- a corner folds;
- its scaled Jacobian falls below `min_jacobian` **and** below what it was
  before (a cell that was already poor is not blamed for staying poor);
- an edge falls below `min_fraction` of its original length;
- or, with `report=True`, the grid no longer re-processes (`build_topology`).

On failure the fades are doubled, up to `widen` times (×2, ×4, …; the angle
on the fault stays exact). From the first doubling on, each fade is at least
its auto distance: `max(d·f, auto·f/2)` at factor `f`. A successful widening
warns. Widening never changes a rigid hanging wall; it only lengthens the
fades. If the grid still fails, `on_invalid` decides:

- `"raise"`: a `DipError` with the failing counts;
- `"bisect"`: the largest valid fraction of the change is applied, with a
  warning. **The requested angle is then not reached**; `scale` reports the
  fraction;
- `"ignore"`: the result is returned with `ok` False.

### Output: `DipResult`

| Field | Meaning |
|---|---|
| `grid`, `before` | The dipped grid and the input |
| `ok` | The gate passed and the grid re-processes |
| `realized_dip(name, grid=None)` | `(M, 2)`: arc length and dip in degrees at each chain pillar, measured from COORD along the dip direction; above 90 means overturned. Pass `grid` to measure another grid, e.g. one read back from disk |
| `angle_before`, `requested` | Per fault: the dip before, and the requested dip, per chain pillar |
| `modes` | `"rigid"` or `"decay"` per fault |
| `fade` | `(decay, footwall)` distances in use, after widening |
| `widen`, `scale`, `bisection_steps`, `passes` | Widening factor, fraction of the change applied, bisection steps, dip rounds |
| `d_top`, `d_bot` | `(NX+1, NY+1, 2)` applied shifts of the pillar tops and bases |
| `check` | `DipCheck`: `n_folded`, `n_degenerate`, `n_squeezed`, `min_jacobian_before`, `min_jacobian`, `bad` (an `(NX, NY, NZ)` mask), `ok`, `summary()` |
| `pillar_validity` | `check_validity` of the COORD planes. Informative only: COORD ends may lie off the cells |
| `volume_change` | `(NX, NY, NZ)` relative change; NaN where inactive |
| `normals`, `chains`, `hanging_sign`, `faults` | The geometry used |
| `slide`, `grid_path`, `report_path` | Filled by `dip_deck` |

`summary()` gives one line and `report()` a full text report: the before,
requested and realized angle at every chain pillar.

### Example

```python
import numpy as np
from cpg2unstructured import CornerPointGrid
from fault_projections import Dip, dip_grid
from fault_projections.faults import FaultFace

nx, ny, nz, dx, dz = 6, 4, 3, 100.0, 10.0
coord = [(i * dx, j * dx, 0.0, i * dx, j * dx, nz * dz)       # vertical pillars, 0 to 30 m
         for j in range(ny + 1) for i in range(nx + 1)]
depth = np.repeat(np.arange(nz) * dz, 2) + np.tile([0.0, dz], nz)
zcorn = np.broadcast_to(depth[:, None, None], (2 * nz, 2 * ny, 2 * nx)).ravel()
grid = CornerPointGrid((nx, ny, nz), np.array(coord), zcorn.astype(float),
                       np.ones(nx * ny * nz, np.int8))

f1 = [FaultFace(3, j, "X") for j in range(1, ny + 1)]   # rim to rim, x = 300 m
res = dip_grid(grid, Dip(f1, 70.0, hanging_side=(500.0, 200.0), name="F1"))
print(res.summary())
# OK: 1 fault(s) dipped (F1 rigid) · F1 70..70 deg · cells valid (min corner Jacobian 1 -> 1) · max |volume change| 1.82 %
print(res.realized_dip("F1"))           # (arc length, degrees) per fault pillar
# [[  0.  70.]
#  [100.  70.]
#  [200.  70.]
#  [300.  70.]
#  [400.  70.]]
print(np.array(res.grid.coord).reshape(ny + 1, nx + 1, 6)[0, 3])
# [294.54044649   0.           0.         305.45955351   0.
#   30.        ]
```

Check it by hand: the pillar runs `305.460 − 294.540 = 10.919` m over 30 m,
and `atan(30 / 10.919) = 70.000°`. The middle pivot (15 m) moved the top
5.46 m west and the base 5.46 m east. This snippet was run against the pinned
revisions.

## `dip_deck`

```python
dip_deck(deck, dips, out_dir, *, prefix="dipped", faults=None, slide=None,
         slide_options=None, **dip_kw) -> DipResult
```

Reads a deck, dips it and writes `<prefix>_grid.grdecl` (COORD, ZCORN, ACTNUM)
and `<prefix>_report.txt`.

- `faults` defaults to the deck's own `FAULTS` (INCLUDE is followed), so a dip
  can name a fault: `Dip("F1", 70.0)`.
- `slide=Throw(...)` slides the dipped grid with
  [`slide_grid`](fault-slide.md#slide_grid) (`slide_options` go to it). The
  written grid is then the slid one, and `DipResult.slide` holds the
  `SlideResult`. A throw that names a deck fault gets the fault's `FAULTS` K
  range as `Throw.layers`, as in `slide_deck`.
- A fault whose `FAULTS` faces are all `Z` raises `DipError`.
- Other keywords go to `dip_grid`.
- The `FAULTS` set does not change, so replace the grid in the deck with
  `INCLUDE '<prefix>_grid.grdecl' /`. The input deck is never written.

Lesson 31 runs every case below and re-checks each result from COORD, a
re-processed topology or the bytes on disk:

| Case | What it shows |
|---|---|
| A | A rigid 60° dip on a rim-to-rim fault; `decay=300` keeps the rim vertical |
| B | 40° needs more than a cell of lean: the naive tilt folds 16 columns, `dip_grid` stays valid |
| C | The four pivots and what they do to column volumes |
| D | A fault with a tip: decay mode, fading around the tip |
| E | An angle table, 55° to 80° along strike |
| F | A stair-stepped diagonal fault dips along its smoothed strike |
| G | A graben: two faults at 60° and 70° dipping toward each other |
| H | The gate: widening, `DipError`, and bisection that misses the angle |
| I | `dip_deck` with `slide=`: angle and throw read back from disk |
| J | Norne (when present): E_01 → 60° passes; BC → 70° fails with the middle pivot and passes with `pivot="top"` |

The full source and the numbers each case prints are on the
[lesson page](../learn/fault-dip.md).

## Failure modes and limitations

- **Faults must lie on pillar lines** (a deck's `FAULTS`). A dip cannot move a
  fault in plan or change which faces it uses.
- **Angles are limited to (0°, 90°].** Overturned faults can be measured
  (`realized_dip` above 90) but not requested.
- **A rigid hanging wall leans entirely**, including the grid rim and any
  fault blocks beyond it. Pass a `decay` distance to keep rims vertical.
- **Cell volumes are not conserved one by one.** The rigid hanging wall
  keeps them; footwall and fade-zone cells do not. With the middle pivot only
  column totals are kept.
- **Faults whose present dip varies sharply along strike can fail the gate**
  whatever the fade. Norne's BC jumps from 59° to 84° within four pillars, and
  at 70° with the middle pivot one top-layer cell folds even when widened ×16;
  `pivot="top"` passes. The neighbouring pillars receive a smoothed blend of
  very different tilts, and widening only lengthens the fades.
- **Large changes on short fades** raise `DipError` unless `widen` or
  `on_invalid` allow otherwise. Bisection keeps the grid valid, but then the
  angle is not reached.
- **The result is geometric only:** no restoration, no geomechanics.

## Source and tests

- [`dip_grid`, `dip_deck`, `Dip`, `DipResult`, `DipCheck`](https://github.com/arturcastiel/fault-projections/blob/75f6dd28e8d6fd3df0d866c71ca2849505e0c5f5/fault_projections/dip.py)
  · [Dip tests](https://github.com/arturcastiel/fault-projections/blob/75f6dd28e8d6fd3df0d866c71ca2849505e0c5f5/tests/test_dip.py)
- Lesson 31 ([fault dip](../learn/fault-dip.md))
