"""Control the throw of a fault that ends inside the grid: profile, far field and share.

Around the tip of such a fault the two sides are still joined, so a uniform throw would tear
welded faces past the tip: the throw must fall to zero there. How it falls, and how far the
motion spreads from the fault, are choices. This lesson builds one throw field on the pillars
from three controls and applies it by editing ZCORN only:

- profile: the relative throw rises from 0 at each tip, linearly or as a cosine (smooth at
  both ends), over a taper length; past the taper the full throw applies;
- far field: either the whole side moves with its nearest fault point, or the motion fades
  to zero at a decay distance from the fault;
- share: the hanging wall (east) goes down by its share of the throw, the footwall (west)
  goes up by the rest.

Six cases on a 16x12x3 grid with a vertical fault on pillar line i = 8. Each keeps COORD and
FAULTS, opens exactly the FAULTS faces, keeps every cell volume, and re-processes to contact
areas that match the analytic overlaps. Pass --out DIR to keep the decks.
"""

import argparse
from dataclasses import dataclass, replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import (CornerPointGrid, build_topology, from_grdecl, read_grdecl,
                              write_grid_include)
from cpg2unstructured.geometry import cell_geometry, face_geometry
from fault_projections.faults import faults_text
from fault_projections.from_faults import edge_face, faults_from_deck

NX, NY, NZ, DX, DZ = 16, 12, 3, 50.0, 20.0
LINE = 8                # the fault lies on pillar line i = 8 (x = 400 m)
TMAX = 16.0             # full relative throw; below DZ, so layer tops never cross inside a face


@dataclass(frozen=True)
class Case:
    name: str
    rows: tuple                   # first and last fault pillar row; row 0 or NY is on the rim
    profile: str = "cosine"       # rise from a tip: "linear" or "cosine"
    taper: float | None = None    # rise length from each tip [m]; None: to the far end or middle
    decay: float | None = None    # None: the whole side moves; else 0 at this distance [m]
    hanging: float = 1.0          # share of the throw taken by the hanging wall (east, down)

    def describe(self):
        far = "whole side" if self.decay is None else f"decay {self.decay:g} m"
        taper = "full length" if self.taper is None else f"taper {self.taper:g} m"
        return f"{self.profile}, {taper}, {far}, hanging wall {self.hanging:.0%}"


CASES = [
    Case("A_linear_whole_side", (0, 8), profile="linear"),
    Case("B_short_taper", (0, 8), taper=150.0),
    Case("C_cosine_whole_side", (0, 8)),
    Case("D_decay", (0, 8), decay=300.0),
    Case("E_decay_split", (0, 8), decay=300.0, hanging=0.5),
    Case("F_two_tips", (2, 10), decay=300.0, hanging=0.7),
]


def layer_cake():
    """Flat DZ m layers on vertical pillars."""
    coord = [(i * DX, j * DX, 0.0, i * DX, j * DX, NZ * DZ)
             for j in range(NY + 1) for i in range(NX + 1)]
    depth = np.repeat(np.arange(NZ) * DZ, 2) + np.tile([0.0, DZ], NZ)
    zcorn = np.broadcast_to(depth[:, None, None], (2 * NZ, 2 * NY, 2 * NX)).ravel()
    return CornerPointGrid((NX, NY, NZ), np.array(coord), zcorn.astype(float),
                           np.ones(NX * NY * NZ, np.int8))


def corner_cube(grid):
    """ZCORN as z[k-slot, j-slot, i-slot]; slot 2c / 2c+1 = near / far corner of cell c."""
    return grid.zcorn.reshape(2 * NZ, 2 * NY, 2 * NX).copy()


def relative_throw(case, y):
    """Throw at along-strike position y, clamped onto the fault (so 0 beyond a tip)."""
    y0, y1 = case.rows[0] * DX, case.rows[1] * DX
    tips = [t for t, row in ((y0, case.rows[0]), (y1, case.rows[1])) if 0 < row < NY]
    rise = case.taper or (y1 - y0) / len(tips)
    u = np.clip(min(abs(np.clip(y, y0, y1) - t) for t in tips) / rise, 0.0, 1.0)
    return TMAX * (u if case.profile == "linear" else 0.5 * (1.0 - np.cos(np.pi * u)))


def pillar_shifts(case):
    """Down-shift (m) of each pillar's corners: (east cells, west cells), each (NY+1, NX+1)."""
    east, west = np.zeros((2, NY + 1, NX + 1))
    y0, y1 = case.rows[0] * DX, case.rows[1] * DX
    for j in range(NY + 1):
        for i in range(NX + 1):
            x, y = i * DX, j * DX
            d = np.hypot(x - LINE * DX, y - np.clip(y, y0, y1))      # distance to the fault
            fade = 1.0 if case.decay is None else \
                0.5 * (1.0 + np.cos(np.pi * min(d / case.decay, 1.0)))
            throw = relative_throw(case, y) * fade
            down, up = case.hanging * throw, (case.hanging - 1.0) * throw
            east[j, i] = down if i >= LINE else up
            west[j, i] = up if i <= LINE else down
    return east, west


def apply(grid, case):
    """Add the pillar shifts to every corner, east-cell or west-cell value by the cell's side."""
    east, west = pillar_shifts(case)
    row = (np.arange(2 * NY) + 1) // 2                     # j-slot -> pillar row
    col = (np.arange(2 * NX) + 1) // 2                     # i-slot -> pillar column
    is_east = (np.arange(2 * NX) // 2 >= LINE)[None, :]
    shift = np.where(is_east, east[row][:, col], west[row][:, col])
    return replace(grid, zcorn=(corner_cube(grid) + shift[None]).ravel())


def unwelded_faces(grid):
    """Cell faces whose two sides no longer share corners: {('X', i, j) | ('Y', i, j)}."""
    z, out = corner_cube(grid), set()
    for i in range(1, NX):
        for j in np.flatnonzero((z[:, :, 2 * i] != z[:, :, 2 * i - 1]).any(0)) // 2:
            out.add(("X", i, int(j)))
    for j in range(1, NY):
        for i in np.flatnonzero((z[:, 2 * j, :] != z[:, 2 * j - 1, :]).any(0)) // 2:
            out.add(("Y", int(i), j))
    return out


def write_deck(grid, faults, path, header):
    write_grid_include(grid, path, header=header)
    with open(path, "a") as deck:
        deck.write("\n" + faults_text([faults], NZ, names=["P"]))
    return path


def reprocess(path):
    """{row j: {(k west | None, k east | None): area}} on the fault line, and cell volumes."""
    topology = build_topology(from_grdecl(read_grdecl(path)))
    areas = face_geometry(topology.nodes, topology.faces, topology.cell_corners)
    volumes = cell_geometry(topology.nodes, topology.faces, topology.cell_faces, areas).volumes
    ijk, out = topology.original_ijk, {}
    for face in np.flatnonzero(topology.faces.direction == 0):
        a, b = topology.faces.neighbours[face]
        cells = sorted([ijk[a]] + ([ijk[b]] if b >= 0 else []), key=lambda c: c[0])
        if b >= 0 and (cells[0][0], cells[1][0]) == (LINE - 1, LINE):
            key = (int(cells[0][2]), int(cells[1][2]))
        elif b < 0 and cells[0][0] in (LINE - 1, LINE):
            k = int(cells[0][2])
            key = (k, None) if cells[0][0] == LINE - 1 else (None, k)
        else:
            continue
        row = out.setdefault(int(cells[0][1]), {})
        row[key] = row.get(key, 0.0) + float(areas.areas[face])
    return out, volumes


def expected_face(throw_a, throw_b):
    """Analytic contacts on one fault face whose throw runs linearly from throw_a to throw_b."""
    mean = (throw_a + throw_b) / 2
    want = {(k, k): DZ - mean for k in range(NZ)}
    if mean > 0:
        want.update({(k + 1, k): mean for k in range(NZ - 1)})
        want.update({(0, None): mean, (None, NZ - 1): mean})
    return {key: DX * length for key, length in want.items()}


def run(case, grid, out):
    faults = [edge_face((LINE, j), (LINE, j + 1)) for j in range(*case.rows)]
    path0 = write_deck(grid, faults, out / f"{case.name}_0_input.grdecl", "lesson 22: input")
    on_deck = faults_from_deck(path0)["P"]
    slid = apply(read_grdecl(path0), case)
    path1 = write_deck(slid, on_deck, out / f"{case.name}_1_slid.grdecl",
                       f"lesson 22 {case.name}: {case.describe()}")
    after = read_grdecl(path1)
    assert np.array_equal(after.coord, grid.coord) and faults_from_deck(path1)["P"] == on_deck
    assert unwelded_faces(after) == {("X", LINE, j) for j in range(*case.rows)}, \
        "the slide opened a face that is not in FAULTS, or left a fault face welded"

    east, west = pillar_shifts(case)
    throw = east[:, LINE] - west[:, LINE]
    tips = [row for row in case.rows if 0 < row < NY]
    assert all(throw[row] == 0.0 for row in tips), "the throw must be 0 at a tip"

    found, volumes = reprocess(path1)
    _, volumes0 = reprocess(path0)
    assert np.allclose(volumes, volumes0, rtol=1e-12, atol=0.0), "a cell changed volume"
    for j in range(*case.rows):
        want = {key: area for key, area in expected_face(throw[j], throw[j + 1]).items() if area}
        got = {key: area for key, area in found[j].items() if area > 1e-9}
        assert got.keys() == want.keys(), (case.name, j, sorted(got, key=str))
        assert all(np.isclose(got[key], want[key]) for key in want), (case.name, j)

    peak = int(np.argmax(throw))
    print(f"{case.name:20s} {case.describe()}")
    print(f"  {'':20s}throw on the fault, rows {case.rows[0]}..{case.rows[1]}: "
          + " ".join(f"{t:.1f}" for t in throw[case.rows[0]:case.rows[1] + 1]))
    print(f"  {'':20s}row {peak}: hanging wall down {east[peak, LINE]:.1f} m at the fault, "
          f"{east[peak, NX]:.1f} m at the east rim; footwall up {abs(west[peak, LINE]):.1f} m at "
          f"the fault, {abs(west[peak, 0]):.1f} m at the west rim")


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    grid = layer_cake()

    # A rigid slide of the whole east side tears the welded faces past the tip.
    z = corner_cube(grid)
    z[:, :, 2 * LINE:] += TMAX
    torn = unwelded_faces(replace(grid, zcorn=z.ravel())) - {("X", LINE, j) for j in range(8)}
    assert torn
    print(f"rigid {TMAX:g} m slide of the east side: opens {len(torn)} faces past the tip that "
          "are not in FAULTS, so it is refused; taper the throw instead")

    for case in CASES:
        run(case, grid, out)
    print("every case: COORD and FAULTS unchanged, exactly the FAULTS faces open, no cell "
          "changed volume, contact areas match the analytic overlaps")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
