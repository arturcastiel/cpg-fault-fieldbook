"""A fault that ends inside the grid: taper the throw to zero at the tip.

FAULTS P cuts pillar line i=3 only for rows j=0..3 of a 6x6x3 layer cake; its tip is pillar
(3, 4). Around the tip the east side is still joined to the west side, so no rigid block can
slide. Instead the throw is a field on pillars: 8 m at the grid edge, falling linearly to 0 at
the tip, applied to every corner east of the fault line. Off the fault every corner stays
welded, so the only new discontinuity is on the FAULTS faces. Pass --out DIR to keep the decks.
"""

import argparse
from collections import deque
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import (CornerPointGrid, NonMatchingZError, build_topology, from_grdecl,
                              read_grdecl, write_grid_include)
from cpg2unstructured.geometry import face_geometry
from fault_projections.faults import faults_text
from fault_projections.from_faults import edge_face, faults_from_deck

NX, NY, NZ, DX, DZ = 6, 6, 3, 100.0, 10.0
LINE, TIP = 3, 4                       # fault on pillar line i=3 for rows j < 4; tip pillar (3, 4)
SLANT = np.hypot(1.0, 40.0 / 80.0)     # fault pillars dip 40 m east over 80 m depth


def layer_cake():
    coord = []
    for j in range(NY + 1):
        for i in range(NX + 1):
            h = 20.0 if i == LINE else 0.0
            coord.append((i * DX - h, j * DX, 0.0, i * DX + h, j * DX, 80.0))
    depth = np.repeat(np.arange(NZ) * DZ, 2) + np.tile([0.0, DZ], NZ)
    zcorn = np.broadcast_to(depth[:, None, None], (2 * NZ, 2 * NY, 2 * NX)).ravel()
    return CornerPointGrid((NX, NY, NZ), np.array(coord), zcorn.astype(float),
                           np.ones(NX * NY * NZ, np.int8))


def corner_cube(grid):
    """ZCORN as z[k-slot, j-slot, i-slot]; slot 2c / 2c+1 = near / far corner of cell c."""
    return grid.zcorn.reshape(2 * NZ, 2 * NY, 2 * NX).copy()


def closes_a_block(faces):
    """Flood fill from the east side of `faces`; True if it cannot reach the west side."""
    sides = [((f.j - 1, f.i - 1), (f.j - 1, f.i)) for f in faces]
    cut = {frozenset(pair) for pair in sides}
    footwall, seed = sides[0]
    seen, queue = {seed}, deque([seed])
    while queue:
        j, i = queue.popleft()
        for n in ((j, i + 1), (j, i - 1), (j + 1, i), (j - 1, i)):
            if 0 <= n[0] < NY and 0 <= n[1] < NX and n not in seen \
                    and frozenset({(j, i), n}) not in cut:
                seen.add(n)
                queue.append(n)
    return footwall not in seen


def tapered_slide(grid, max_throw):
    """Move corners east of LINE down their pillars by a throw that depends on the pillar row."""
    throw = np.clip(max_throw * (TIP - np.arange(NY + 1)) / TIP, 0.0, None)   # per pillar row
    pillar_row = (np.arange(2 * NY) + 1) // 2                                # j-slot -> row
    east = np.arange(2 * NX) // 2 >= LINE                                    # i-slot -> side
    z = corner_cube(grid) + np.where(east[None, :], throw[pillar_row][:, None], 0.0)
    return replace(grid, zcorn=z.ravel())


def unwelded_faces(grid):
    """Cell faces whose two sides no longer share corners: {('X'|'Y', line, row)}."""
    z, out = corner_cube(grid), set()
    for i in range(1, NX):
        for j in np.flatnonzero((z[:, :, 2 * i] != z[:, :, 2 * i - 1]).any(0)) // 2:
            out.add(("X", i, int(j)))
    for j in range(1, NY):
        for i in np.flatnonzero((z[:, 2 * j, :] != z[:, 2 * j - 1, :]).any(0)) // 2:
            out.add(("Y", j, int(i)))
    return out


def write_deck(grid, faults, path, header):
    write_grid_include(grid, path, header=header)
    with open(path, "a") as deck:
        deck.write("\n" + faults_text([faults], NZ, names=["P"]))
    return path


def contacts_by_row(path):
    """{row j: {(k west, k east): area}} on the fault; None marks an exposed piece."""
    topology = build_topology(from_grdecl(read_grdecl(path)))
    areas = face_geometry(topology.nodes, topology.faces, topology.cell_corners).areas
    ijk, out = topology.original_ijk, {}
    for face in np.flatnonzero(topology.faces.direction == 0):
        a, b = topology.faces.neighbours[face]
        cells = sorted([ijk[a]] + ([ijk[b]] if b >= 0 else []), key=lambda c: c[0])
        j = int(cells[0][1])
        if j >= TIP:
            continue
        if b >= 0 and (cells[0][0], cells[1][0]) == (LINE - 1, LINE):
            key = (int(cells[0][2]), int(cells[1][2]))
        elif b < 0 and cells[0][0] in (LINE - 1, LINE):
            k = int(cells[0][2])
            key = (k, None) if cells[0][0] == LINE - 1 else (None, k)
        else:
            continue
        out.setdefault(j, {})[key] = out.setdefault(j, {}).get(key, 0.0) + float(areas[face])
    return out


def expected_row(mean_throw):
    """Throw < DZ everywhere: each overlap is linear in the throw, so use its row mean."""
    want = {(k, k): DZ - mean_throw for k in range(NZ)}
    want.update({(k + 1, k): mean_throw for k in range(NZ - 1)})
    want.update({(0, None): mean_throw, (None, NZ - 1): mean_throw})
    return {key: DX * SLANT * length for key, length in want.items()}


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    fault = [edge_face((LINE, j), (LINE, j + 1)) for j in range(TIP)]
    source = write_deck(layer_cake(), fault, out / "input.grdecl", "lesson 20: input")
    grid, on_deck = read_grdecl(source), faults_from_deck(source)["P"]
    assert not closes_a_block(on_deck)
    print(f"rigid slide impossible: east of P reaches the west side around the tip "
          f"pillar ({LINE}, {TIP})")

    slid = tapered_slide(grid, 8.0)
    path = write_deck(slid, on_deck, out / "tapered_8m.grdecl", "lesson 20: tapered 8 m")
    z = corner_cube(slid)
    profile = [z[0, min(2 * row, 2 * NY - 1), 2 * LINE] - z[0, min(2 * row, 2 * NY - 1), 2 * LINE - 1]
               for row in range(NY + 1)]
    print("throw on fault-line pillars, rows 0..6 [m]:", " ".join(f"{t:g}" for t in profile))
    assert np.array_equal(slid.coord, grid.coord) and faults_from_deck(path)["P"] == on_deck
    assert unwelded_faces(slid) == {("X", LINE, j) for j in range(TIP)}
    print(f"only the {TIP} FAULTS faces opened; every other face, including past the tip, "
          f"is still welded")

    found = contacts_by_row(path)
    for j in range(TIP):
        want = expected_row((profile[j] + profile[j + 1]) / 2)
        assert found[j].keys() == want.keys(), (j, sorted(found[j], key=str))
        assert all(np.isclose(found[j][key], want[key]) for key in want)
        listed = ", ".join(f"{'-' if w is None else f'k{w}'}|{'-' if e is None else f'k{e}'} "
                           f"{a:,.0f}" for (w, e), a in sorted(found[j].items(), key=str))
        print(f"  row {j} throw {profile[j]:g}->{profile[j + 1]:g} m  contacts W|E m2: {listed}")
    print("re-processed contact areas match the analytic overlaps; '-' = exposed fault piece")

    deep = write_deck(tapered_slide(grid, 15.0), on_deck, out / "tapered_15m.grdecl",
                      "lesson 20: tapered 15 m (layer boundaries cross along strike)")
    try:
        contacts_by_row(deep)
    except NonMatchingZError as error:
        print("Expected current limit for a 15 m taper:", error)
    else:
        raise AssertionError("crossing layer boundaries along the fault were accepted")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
