"""Slide a fault that is oblique to the grid axes and dips.

The fault trace at the top (z = 0) runs from (300, 0) to (600, 600), 27 degrees off the J axis
of an 8x6x3 layer cake (40 m layers), so no single pillar line carries it. fault_faces_for turns
the trace into a staircase of X and Y FAULTS faces. The staircase pillars, and every other
pillar of the south-east hanging wall, keep their tops and lean 75 m east at the base (120 m),
so the fault dips about 61 degrees to the south-east; the footwall pillars stay vertical.
Two decks carry the same FAULTS: (A) the regular grid and (B) the grid after adapt_grid, which
moves the fault pillars rigidly onto the trace and keeps their lean. In both, the hanging wall
slides 25 m down its pillars by editing ZCORN only. Its pillars are parallel, so it moves as one
rigid piece along the fault and no cell changes volume; the re-processed contact areas follow
each deck's fault surface. Pass --out DIR to keep the four decks.
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

NX, NY, NZ, DX, DZ = 8, 6, 3, 100.0, 40.0
TRACE = np.array([[300.0, 0.0], [600.0, 600.0]])   # fault trace at the top, z = 0
LEAN = 75.0                             # hanging-wall pillars move 75 m east by the base
THROW = 25.0
SEED, FAR = (0, NX - 1), (NY - 1, 0)       # (j, i) columns on the moving / fixed side


def layer_cake(leaning=()):
    """Flat DZ m layers; pillars in `leaning` keep their top and move LEAN m east at the base."""
    coord = [(i * DX, j * DX, 0.0, i * DX + (LEAN if (i, j) in leaning else 0.0), j * DX, NZ * DZ)
             for j in range(NY + 1) for i in range(NX + 1)]
    depth = np.repeat(np.arange(NZ) * DZ, 2) + np.tile([0.0, DZ], NZ)   # 0,40,40,80,80,120
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


def area_per_metre(grid, faces):
    """Fault area per metre of depth, from COORD.

    Both pillars of a FAULTS face lean the same way, so the face is a parallelogram: its area
    between two depths is |edge x lean| times the depth difference, with lean = d(xyz)/dz.
    """
    coord = grid.coord.reshape(NY + 1, NX + 1, 6)
    lean = (coord[..., 3:] - coord[..., :3]) / (coord[..., 5:] - coord[..., 2:3])
    total = 0.0
    for a, b in map(face_pillars, faces):
        assert np.allclose(lean[a[1], a[0]], lean[b[1], b[0]]), "fault pillars must be parallel"
        edge = coord[b[1], b[0], :3] - coord[a[1], a[0], :3]
        total += float(np.linalg.norm(np.cross(edge, lean[a[1], a[0]])))
    return total


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
    """Re-read the deck; return {(k fixed, k moving): area} across the fault, cell volumes."""
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


def dip_degrees(grid, faces):
    """Dip of the plane through the top trace and the fault pillars' lean."""
    coord = grid.coord.reshape(NY + 1, NX + 1, 6)
    i, j = face_pillars(faces[0])[0]
    lean = coord[j, i, 3:] - coord[j, i, :3]
    normal = np.cross(np.append(TRACE[1] - TRACE[0], 0.0), lean)
    return float(np.degrees(np.arccos(abs(normal[2]) / np.linalg.norm(normal))))


def run(name, grid, faces, out):
    """Write, slide and re-process one deck; return its fault area per metre of depth."""
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

    per_metre = area_per_metre(after, faces)
    found, volumes = reprocess(path1, block)
    expected = {(kf, km): per_metre * overlap(kf, km)
                for kf in range(NZ) for km in range(NZ) if overlap(kf, km) > 0}
    assert found.keys() == expected.keys(), (sorted(found), sorted(expected))
    for key, area in found.items():
        assert np.isclose(area, expected[key], rtol=1e-9), (key, area, expected[key])

    # Every hanging-wall pillar leans the same way, so every moving corner slides by the same
    # vector: the block translates rigidly along the fault and no cell changes volume.
    _, volumes0 = reprocess(path0, block)
    assert np.allclose(volumes, volumes0, rtol=1e-9, atol=0.0), "a cell changed volume"

    n_x = sum(f.face == "X" for f in faces)
    pairs = ", ".join(f"k{kf}|k{km} {a:,.0f}" for (kf, km), a in sorted(found.items()))
    print(f"{name}: {len(faces)} faces ({n_x} X, {len(faces) - n_x} Y), fault area "
          f"{per_metre:6.1f} m2 per m of depth, {block.sum()} columns slid {THROW:g} m")
    print(f"  contacts fixed|moving m2: {pairs}")
    print("  no cell changed volume: the hanging wall moved as one rigid piece")
    return per_metre


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    staircase = fault_faces_for(layer_cake(), TRACE)
    hanging = {(i + di, j + dj) for j, i in zip(*np.nonzero(moving_block(staircase)))
               for di in (0, 1) for dj in (0, 1)}            # every pillar of a moving column
    regular = layer_cake(hanging)
    assert fault_faces_for(regular, TRACE) == staircase, "the lean should not move the tops"
    adapted = adapt_grid(regular, TRACE, boundary="pinned", w0=0.05, polish="barrier",
                         on_invalid="raise")
    assert adapted.ok and adapted.faces_invariant, adapted.summary()
    assert adapted.faults_after == staircase, "adaptation should keep the FAULTS faces"
    lean = lambda g: g.coord.reshape(-1, 6)[:, 3:5] - g.coord.reshape(-1, 6)[:, :2]  # noqa: E731
    assert np.allclose(lean(adapted.grid), lean(regular)), "stage 1 keeps every pillar's lean"
    print(f"trace length {np.linalg.norm(TRACE[1] - TRACE[0]):.1f} m at z = 0; fault dips "
          f"{dip_degrees(regular, staircase):.1f} deg to the south-east")
    print(f"adapted fault pillars within {adapted.distance_after.max():.1f} m of the trace "
          f"(regular grid: {adapted.distance_before.max():.1f} m)")

    area_a = run("A_regular", regular, staircase, out)
    area_b = run("B_adapted", adapted.grid, adapted.faults_after, out)
    print(f"same FAULTS, same slide; the regular staircase has {area_a / area_b - 1:.0%} "
          "more fault area than the adapted grid")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the four decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
