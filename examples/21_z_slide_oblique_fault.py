"""Slide a fault whose trace is oblique to the grid axes.

A vertical fault trace runs from (300, 0) to (600, 600), 27 degrees off the J axis of an 8x6x3
layer cake, so no single pillar line carries it. fault_faces_for turns the trace into a
staircase of X and Y FAULTS faces. Two decks carry the same FAULTS: (A) the regular grid, whose
staircase is 900 m long, and (B) the grid after adapt_grid, whose fault pillars sit on the
trace. In both, the south-east block slides 7 m down its pillars by editing ZCORN only; the
re-processed contact areas follow the fault length of each deck. Pass --out DIR to keep the decks.
"""

import argparse
from collections import deque
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import (CornerPointGrid, build_topology, from_grdecl, read_grdecl,
                              write_grid_include)
from cpg2unstructured.geometry import cell_geometry, face_geometry
from fault_projections import adapt_grid
from fault_projections.adapt import fault_faces_for
from fault_projections.faults import faults_text
from fault_projections.from_faults import faults_from_deck

NX, NY, NZ, DX, DZ = 8, 6, 3, 100.0, 10.0
TRACE = np.array([[300.0, 0.0], [600.0, 600.0]])
THROW = 7.0
SEED, FAR = (0, NX - 1), (NY - 1, 0)       # (j, i) columns on the moving / fixed side


def layer_cake():
    """Flat 10 m layers on vertical pillars."""
    coord = [(i * DX, j * DX, 0.0, i * DX, j * DX, 80.0)
             for j in range(NY + 1) for i in range(NX + 1)]
    depth = np.repeat(np.arange(NZ) * DZ, 2) + np.tile([0.0, DZ], NZ)   # 0,10,10,20,20,30
    zcorn = np.broadcast_to(depth[:, None, None], (2 * NZ, 2 * NY, 2 * NX)).ravel()
    return CornerPointGrid((NX, NY, NZ), np.array(coord), zcorn.astype(float),
                           np.ones(NX * NY * NZ, np.int8))


def corner_cube(grid):
    """ZCORN as z[k-slot, j-slot, i-slot]; slot 2c / 2c+1 = near / far corner of cell c."""
    return grid.zcorn.reshape(2 * NZ, 2 * NY, 2 * NX).copy()


def face_sides(f):
    """The (j, i) columns on either side of a 1-based FAULTS face."""
    return (f.j - 1, f.i - 1), (f.j - 1, f.i) if f.face == "X" else (f.j, f.i - 1)


def face_pillars(f):
    """The two (i, j) pillars that bound a 1-based FAULTS face."""
    return ((f.i, f.j - 1), (f.i, f.j)) if f.face == "X" else ((f.i - 1, f.j), (f.i, f.j))


def moving_block(faces):
    """Columns (j, i) reachable from SEED without crossing `faces`; FAR must stay outside."""
    cut = {frozenset(face_sides(f)) for f in faces}
    block = np.zeros((NY, NX), bool)
    block[SEED] = True
    queue = deque([SEED])
    while queue:
        j, i = queue.popleft()
        for n in ((j, i + 1), (j, i - 1), (j + 1, i), (j - 1, i)):
            if 0 <= n[0] < NY and 0 <= n[1] < NX and not block[n] \
                    and frozenset({(j, i), n}) not in cut:
                block[n] = True
                queue.append(n)
    assert not block[FAR], "the FAULTS do not cut the grid in two"
    return block


def slide(grid, faces, throw):
    """Move every corner of the moving block `throw` m down its own pillar."""
    z = corner_cube(grid)
    z[:, np.repeat(np.repeat(moving_block(faces), 2, 0), 2, 1)] += throw
    return replace(grid, zcorn=z.ravel())


def fault_length(grid, faces):
    """Plan-view length of the FAULTS faces, from the pillar tops in COORD."""
    top = grid.coord.reshape(NY + 1, NX + 1, 6)[:, :, :2]
    return sum(float(np.linalg.norm(top[b[1], b[0]] - top[a[1], a[0]]))
               for a, b in map(face_pillars, faces))


def unwelded_faces(grid):
    """FAULTS-style keys (i, j, face) of cell faces whose two sides no longer share corners."""
    z, out = corner_cube(grid), set()
    for i in range(1, NX):
        for j in np.flatnonzero((z[:, :, 2 * i] != z[:, :, 2 * i - 1]).any(0)) // 2:
            out.add((i, int(j) + 1, "X"))
    for j in range(1, NY):
        for i in np.flatnonzero((z[:, 2 * j, :] != z[:, 2 * j - 1, :]).any(0)) // 2:
            out.add((int(i) + 1, j, "Y"))
    return out


def write_deck(grid, faces, path, header):
    write_grid_include(grid, path, header=header)
    with open(path, "a") as deck:
        deck.write("\n" + faults_text([faces], NZ, names=["OBLIQUE"]))
    return path


def reprocess(path, block):
    """Re-read the deck; return {(k fixed, k moving): area} across the fault and cell volumes."""
    topology = build_topology(from_grdecl(read_grdecl(path)))
    areas = face_geometry(topology.nodes, topology.faces, topology.cell_corners)
    volumes = cell_geometry(topology.nodes, topology.faces, topology.cell_faces, areas).volumes
    ijk, found = topology.original_ijk, {}
    for face in np.flatnonzero(topology.faces.direction <= 1):          # X and Y faces
        a, b = topology.faces.neighbours[face]
        if b < 0:
            continue
        ca, cb = ijk[a], ijk[b]
        if block[ca[1], ca[0]] == block[cb[1], cb[0]]:
            continue
        fixed, moving = (ca, cb) if block[cb[1], cb[0]] else (cb, ca)
        key = (int(fixed[2]), int(moving[2]))
        found[key] = found.get(key, 0.0) + float(areas.areas[face])
    return found, volumes


def overlap(k_fixed, k_moving):
    """Depth overlap (m) of two layers when the moving side sits THROW m deeper."""
    top = max(k_fixed * DZ, k_moving * DZ + THROW)
    return max(0.0, min((k_fixed + 1) * DZ, (k_moving + 1) * DZ + THROW) - top)


def run(name, grid, faces, out):
    """Write, slide and re-process one deck; return its fault length."""
    path0 = write_deck(grid, faces, out / f"{name}_0_input.grdecl", f"lesson 21 {name}: input")
    on_deck = faults_from_deck(path0)["OBLIQUE"]
    assert on_deck == faces, "FAULTS did not survive the round trip"
    before = read_grdecl(path0)
    block = moving_block(on_deck)
    path1 = write_deck(slide(before, on_deck, THROW), on_deck, out / f"{name}_1_slid.grdecl",
                       f"lesson 21 {name}: slid {THROW:g} m")
    after = read_grdecl(path1)
    assert np.array_equal(after.coord, before.coord), "COORD changed"
    assert faults_from_deck(path1)["OBLIQUE"] == faces, "FAULTS changed"
    assert unwelded_faces(after) == {(f.i, f.j, f.face) for f in faces}, \
        "the slide opened a face that is not in FAULTS"

    length = fault_length(after, faces)
    found, volumes = reprocess(path1, block)
    expected = {(kf, km): length * overlap(kf, km)
                for kf in range(NZ) for km in range(NZ) if overlap(kf, km) > 0}
    assert found.keys() == expected.keys(), (sorted(found), sorted(expected))
    for key, area in found.items():
        assert np.isclose(area, expected[key], rtol=1e-9), (key, area, expected[key])
    _, volumes0 = reprocess(path0, block)
    assert np.allclose(volumes, volumes0, rtol=1e-12), "vertical pillars keep cell volumes"

    n_x = sum(f.face == "X" for f in faces)
    pairs = ", ".join(f"k{kf}|k{km} {a:,.0f}" for (kf, km), a in sorted(found.items()))
    print(f"{name}: {len(faces)} faces ({n_x} X, {len(faces) - n_x} Y), fault length "
          f"{length:6.1f} m, {block.sum()} columns slid {THROW:g} m")
    print(f"  contacts fixed|moving m2: {pairs}")
    return length


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    regular = layer_cake()
    staircase = fault_faces_for(regular, TRACE)
    adapted = adapt_grid(regular, TRACE, boundary="pinned", w0=0.05, polish="barrier",
                         on_invalid="raise")
    assert adapted.ok and adapted.faces_invariant, adapted.summary()
    assert adapted.faults_after == staircase, "adaptation should keep the FAULTS faces"
    top = adapted.grid.coord.reshape(-1, 6)
    assert np.allclose(top[:, :2], top[:, 3:5]), "adapted pillars stay vertical"
    print(f"trace length {np.linalg.norm(TRACE[1] - TRACE[0]):.1f} m; adapted fault pillars "
          f"within {adapted.distance_after.max():.1f} m of it "
          f"(regular grid: {adapted.distance_before.max():.1f} m)")

    length_a = run("A_regular", regular, staircase, out)
    length_b = run("B_adapted", adapted.grid, adapted.faults_after, out)
    print(f"same FAULTS, same slide; the regular staircase has {length_a / length_b - 1:.0%} "
          "more fault area than the adapted grid")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the four decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
