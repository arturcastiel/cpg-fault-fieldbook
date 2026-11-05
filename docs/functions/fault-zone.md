# Fault zones: `fault_zone`, `FaultZone`, `write_zone`, and `zone=`

A fault zone is the set of cells within a distance `D` of a fault, on either
side. There are two separate features:

- **Selection:** `fault_zone` finds the cells in `−D … +D` of each fault. It
  never changes the grid.
- **Confinement (optional):** `zone=` on [`dip_grid`](fault-dip.md) and
  [`slide_grid`](fault-slide.md) lets only cells of the zone deform. **The
  default is no zone:** the whole domain may deform, and the grid's outer
  boundary holds (`dip_grid(boundary="fixed")`).

The lesson with every case worked through is
[9. Fault zones](../learn/fault-zone.md).

## `fault_zone`

### Syntax

```python
fault_zone(
    grid: CornerPointGrid,
    faults,                      # {name: faces}, a list of face lists, or one face list
    distance,                    # D, or (D_footwall, D_hanging_wall)
    *,
    hanging_side=None,           # point, "left"/"right", {name: value}; None reads the dip
    slices: int | None = None,   # depth slabs for the candidate search (not the accuracy)
) -> FaultZone
```

- **`faults`:** a deck's FAULTS (from `faults_from_deck`) work as they are.
  A fault may have several chains. Z faces (no plan trace) are skipped.
- **`distance`:** in model units, from each cell's centre along the normal of
  the fault surface.
- **`hanging_side`:** which side is `+`. It works as `Throw.hanging_side`
  does. With `None`, it is read from the fault's present dip; a vertical
  fault then takes `+` on the left of its chain, with a warning, and is
  listed in `FaultZone.assumed`.
- **`slices`:** the default is about square slabs (4 to 64). It only decides
  which faces are searched near a cell; the distance itself is exact.

### How the distance is measured

The fault surface is the ruled surface between the fault's pillars, extended
past the top and base of the grid by `0.01·H + D`. Between two pillars that
are not parallel, the face is *twisted* (curved). For each cell centre:

1. Candidate faces are found with a KD-tree on face patches (depth slabs),
   pruned by distance bounds.
2. On each candidate patch, the nearest point on each ruling (the straight
   line between the two pillars at one depth) has a closed form. The best
   position along the face is found by 33 samples followed by golden-section
   refinement. The result is never worse than the best sample.
3. The sign comes from the side of the surface's normal, flipped so `+` is
   the hanging wall.

Past a fault's tip or the grid's rim, the distance is to the surface's edge,
so the zone wraps round a tip. On Norne the result agrees with a brute force
over 2.6 million surface points to 0.011 m (the brute force's sampling
error), and is never farther.

### Output: `FaultZone`

| Field | Meaning |
|---|---|
| `distance` | `(D_footwall, D_hanging_wall)` |
| `signed[name]` | `(NX, NY, NZ)` signed distance per fault; NaN beyond the zone (not measured) |
| `hanging_sign[name]` | `+1` when the hanging wall is left of the chain |
| `assumed` | Vertical faults whose `+` side was assumed (left) |
| `names` | The fault names |
| `nearest` | `(NX, NY, NZ)` index of the nearest fault, −1 if none. Faults within 1e-6 of each other tie, and the first listed wins |
| `side(name=None)` | int8: −1 footwall band, +1 hanging-wall band, 0 outside. Without `name`, the nearest fault decides |
| `minus`, `plus`, `cells` | Bool masks: footwall band, hanging-wall band, either |
| `counts()` | `{name: (footwall cells, hanging-wall cells)}`, each fault on its own |
| `property()` | `NX*NY*NZ` ints in deck order (I fastest): 0 outside, 1 footwall band, 2 hanging-wall band |
| `summary()` | One line |

## `write_zone`

```python
write_zone(zone: FaultZone, path, keyword: str = "FLTZONE") -> Path
```

Writes `property()` as a GRDECL include, run-length encoded (`n*value`).
Rename `keyword` (for example to `MULTNUM` or `FIPNUM`) to use it in a deck.

## `zone=` on `dip_grid` and `slide_grid`

```python
dip_grid(grid, dips, ..., zone=None)     # also dip_deck(..., zone=...)
slide_grid(grid, throws, ..., zone=None) # also slide_deck(..., zone=...)
```

`zone` is `D`, `(D_footwall, D_hanging_wall)` or a `FaultZone`. A distance
is measured around the dipped or slid faults on the input grid.

- **Dip:** a pillar moves only if every cell on it (all layers, all four
  columns) is in the zone.
- **Slide:** a ZCORN corner moves only if every cell sharing it is in the
  zone.
- **Inside the zone,** the change fades with a cosine in
  `d_fault / (d_fault + d_frozen)`: full on the fault, zero at the nearest
  frozen pillar or corner.
- **Outside the zone,** every corner of every cell is bit-identical.
- **`dip_grid(slide=..., zone=...)`** confines the slide to the same zone.
- **`ZoneError`** (a `ValueError`) if the zone freezes a pillar or corner on
  the fault: the zone is narrower than the cells around the fault.
- **The usual gate still applies.** A zone too narrow for the change folds
  cells, and the dip or slide is refused (`DipError` or `SlideError`).
- `DipResult.zone`, `SlideResult.zone`: the zone used. `DipResult.moved`
  marks the cells with a corner that moved.

### Example

```python
import numpy as np
from cpg2unstructured import CornerPointGrid
from fault_projections import Dip, dip_grid, fault_zone, write_zone
from fault_projections.faults import FaultFace

nx, ny, nz, dx, dz = 8, 4, 3, 50.0, 10.0
coord = [(i * dx, j * dx, 0.0, i * dx, j * dx, nz * dz)       # vertical pillars, 0 to 30 m
         for j in range(ny + 1) for i in range(nx + 1)]
depth = np.repeat(np.arange(nz) * dz, 2) + np.tile([0.0, dz], nz)
zcorn = np.broadcast_to(depth[:, None, None], (2 * nz, 2 * ny, 2 * nx)).ravel()
grid = CornerPointGrid((nx, ny, nz), np.array(coord), zcorn.astype(float),
                       np.ones(nx * ny * nz, np.int8))

faults = {"F1": [FaultFace(4, j, "X") for j in range(1, ny + 1)]}   # x = 200 m, rim to rim
z = fault_zone(grid, faults, 80.0, hanging_side="right")
print(z.summary())
# fault zone -80 / +80: 24 footwall + 24 hanging-wall cells; F1: 24 footwall + 24 hanging-wall cells
print(z.signed["F1"][:, 0, 0])          # signed distance along row 0, top layer
# [ nan  nan -75. -25.  25.  75.  nan  nan]
print(z.side()[:, 0, 0])                # -1 footwall band, +1 hanging-wall band, 0 outside
# [ 0  0 -1 -1  1  1  0  0]
write_zone(z, "zone.grdecl")            # FLTZONE: 0 outside, 1 footwall, 2 hanging wall
# -- fault zone -80 / +80: 0 outside, 1 footwall band, 2 hanging-wall band
# FLTZONE
# 2*0 2*1 2*2 4*0 2*1 2*2 4*0 2*1 2*2 4*0 2*1 2*2 ...

res = dip_grid(grid, Dip("F1", 70.0, hanging_side="right"), faults=faults, zone=z)
print(res.summary())
# OK: 1 fault(s) dipped (F1 rigid) · F1 70..70 deg · kept inside the zone -80/+80 (48 of 48 zone cells moved) · cells valid (min corner Jacobian 1 -> 1, 3-D 1 -> 0.94) · max |volume change| 5.46 %
print(int(res.moved.sum()), int((res.moved & ~z.cells).sum()))   # moved; moved outside
# 48 0
```

Check it by hand: the cell centres on row 0 lie at x = 25, 75, …, 375 m, so
their distances from x = 200 m are −175, −125, −75, −25, 25, 75, … m. Only
|d| ≤ 80 is in the zone: two columns on each side, 4 rows × 3 layers each,
24 + 24 cells. This snippet was run against the pinned revisions.

## Failure modes and limitations

- **Membership is by cell centre.** A cell whose centre lies within `D` is
  in the zone, however far its corners reach.
- **The zone is measured on the grid given.** After a dip the fault has
  moved; measure a new zone on the result if needed.
- **A confined edit is harder on the cells it moves.** The whole change is
  taken up between the fault and the zone's edge. Too narrow a zone raises
  `ZoneError` (fault pillars frozen) or fails the gate (cells fold); lesson
  33 J shows both.
- **Faults must lie on pillar lines** (a deck's `FAULTS`). Z faces are
  skipped.
- **Run time** grows with the number of cells near faults. All 61 Norne
  faults at D = 100 m take about 9 s.

## Source and tests

- [`fault_zone`, `FaultZone`, `write_zone`, `ZoneError`](https://github.com/arturcastiel/fault-projections/blob/48dd1ece3614b875c524e9a6a6bbefac99ddf5d7/fault_projections/zone.py)
  · [Zone tests](https://github.com/arturcastiel/fault-projections/blob/48dd1ece3614b875c524e9a6a6bbefac99ddf5d7/tests/test_zone.py)
- Lesson 33 ([fault zones](../learn/fault-zone.md))
