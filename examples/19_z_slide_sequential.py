"""Slide two faults in Z, one after the other, by editing ZCORN only.

Two tilted normal faults cut a 9x4x3 layer cake into blocks W | M | E. Sliding F1 moves its
whole hanging wall (M and E) 7 m down along the pillars; sliding F2 then moves E another 12 m.
COORD, ACTNUM and FAULTS never change. After each step the deck is written and re-processed,
because contacts, exposed fault pieces and areas are derived from ZCORN. Pass --out DIR to keep
the three decks.
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
from fault_projections.faults import faults_text
from fault_projections.from_faults import edge_face, faults_from_deck

NX, NY, NZ, DX, DZ = 9, 4, 3, 100.0, 10.0
FAULT_LINE = {"F1": 3, "F2": 6}          # pillar line i carrying each fault
THROW = {"F1": 7.0, "F2": 12.0}


def layer_cake():
    """Flat 10 m layers; both fault pillar lines dip east (x +-20 m over 80 m depth)."""
    coord = []
    for j in range(NY + 1):
        for i in range(NX + 1):
            h = 20.0 if i in FAULT_LINE.values() else 0.0
            coord.append((i * DX - h, j * DX, 0.0, i * DX + h, j * DX, 80.0))
    depth = np.repeat(np.arange(NZ) * DZ, 2) + np.tile([0.0, DZ], NZ)   # 0,10,10,20,20,30
    zcorn = np.broadcast_to(depth[:, None, None], (2 * NZ, 2 * NY, 2 * NX)).ravel()
    return CornerPointGrid((NX, NY, NZ), np.array(coord), zcorn.astype(float),
                           np.ones(NX * NY * NZ, np.int8))


def corner_cube(grid):
    """ZCORN as z[k-slot, j-slot, i-slot]; slot 2c / 2c+1 = near / far corner of cell c."""
    return grid.zcorn.reshape(2 * NZ, 2 * NY, 2 * NX).copy()


def hanging_wall(faces):
    """Columns (j, i) reachable from the east side of `faces` without crossing them."""
    sides = [((f.j - 1, f.i - 1), (f.j - 1, f.i) if f.face == "X" else (f.j, f.i - 1))
             for f in faces]
    cut = {frozenset(pair) for pair in sides}
    footwall, seed = sides[0]
    block = np.zeros((NY, NX), bool)
    block[seed] = True
    queue = deque([seed])
    while queue:
        j, i = queue.popleft()
        for n in ((j, i + 1), (j, i - 1), (j + 1, i), (j - 1, i)):
            if 0 <= n[0] < NY and 0 <= n[1] < NX and not block[n] \
                    and frozenset({(j, i), n}) not in cut:
                block[n] = True
                queue.append(n)
    assert not block[footwall], "fault does not close a block; use a tapered slide (lesson 20)"
    return block


def slide(grid, faces, throw):
    """Move every corner of the hanging-wall columns `throw` m down its own pillar."""
    z = corner_cube(grid)
    z[:, np.repeat(np.repeat(hanging_wall(faces), 2, 0), 2, 1)] += throw
    return replace(grid, zcorn=z.ravel())


def write_deck(grid, faults, path, header):
    write_grid_include(grid, path, header=header)
    with open(path, "a") as deck:
        deck.write("\n" + faults_text(list(faults.values()), NZ, names=list(faults)))
    return path


def reprocess(path):
    """Read the deck back and derive the fault contacts from its ZCORN."""
    topology = build_topology(from_grdecl(read_grdecl(path)))
    areas = face_geometry(topology.nodes, topology.faces, topology.cell_corners)
    volumes = cell_geometry(topology.nodes, topology.faces, topology.cell_faces, areas).volumes
    return topology, areas.areas, volumes


def contacts(topology, areas, line):
    """{(k west, k east): area} on pillar line `line`; None marks an exposed piece."""
    ijk, out = topology.original_ijk, {}
    for face in np.flatnonzero(topology.faces.direction == 0):
        a, b = topology.faces.neighbours[face]
        cells = sorted([ijk[a]] + ([ijk[b]] if b >= 0 else []), key=lambda c: c[0])
        if b >= 0 and (cells[0][0], cells[1][0]) == (line - 1, line):
            key = (int(cells[0][2]), int(cells[1][2]))
        elif b < 0 and cells[0][0] in (line - 1, line):
            k = int(cells[0][2])
            key = (k, None) if cells[0][0] == line - 1 else (None, k)
        else:
            continue
        out[key] = out.get(key, 0.0) + float(areas[face])
    return out


def overlapping_pairs(throw):
    """Layer pairs whose depth intervals overlap when the east side sits `throw` m deeper."""
    return {(kw, ke) for kw in range(NZ) for ke in range(NZ)
            if min(kw + 1, ke + 1 + throw / DZ) - max(kw, ke + throw / DZ) > 1e-9}


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    faults = {name: [edge_face((i, j), (i, j + 1)) for j in range(NY)]
              for name, i in FAULT_LINE.items()}
    paths = [write_deck(layer_cake(), faults, out / "0_input.grdecl", "lesson 19: input")]
    grid = read_grdecl(paths[0])
    for step, name in enumerate(FAULT_LINE, start=1):
        on_deck = faults_from_deck(paths[-1])     # the slide is driven by the deck's FAULTS
        grid = slide(grid, on_deck[name], THROW[name])
        paths.append(write_deck(grid, on_deck, out / f"{step}_after_{name}.grdecl",
                                f"lesson 19: after sliding {name} by {THROW[name]:g} m"))

    grids = [read_grdecl(path) for path in paths]
    assert all(np.array_equal(g.coord, grids[0].coord) for g in grids), "COORD changed"
    assert all(faults_from_deck(path) == faults for path in paths), "FAULTS changed"
    reverse = slide(slide(grids[0], faults["F2"], THROW["F2"]), faults["F1"], THROW["F1"])
    assert np.array_equal(reverse.zcorn, grids[2].zcorn), "slides should commute"

    volume0 = reprocess(paths[0])[2]
    for step, (path, g) in enumerate(zip(paths, grids)):
        z = corner_cube(g)
        topology, areas, volumes = reprocess(path)
        print(f"stage {step}: {path.name}")
        for name, i in FAULT_LINE.items():
            throw = np.unique(z[:, :, 2 * i] - z[:, :, 2 * i - 1])
            assert throw.size == 1, "a rigid slide keeps the throw constant along the fault"
            found = contacts(topology, areas, i)
            assert {k for k in found if None not in k} == overlapping_pairs(throw[0])
            listed = ", ".join(f"{'-' if w is None else f'k{w}'}|{'-' if e is None else f'k{e}'}"
                               f" {a:,.0f}" for (w, e), a in sorted(found.items(), key=str))
            print(f"  {name} throw {throw[0]:4.1f} m  contacts W|E m2: {listed}")
        changed = sorted({int(c[0]) for c in topology.original_ijk[
            ~np.isclose(volumes, volume0, rtol=1e-12)]})
        print(f"  columns i whose cell volume changed: {changed or 'none'}")
    print("COORD and FAULTS unchanged; F1-then-F2 == F2-then-F1; '-' = exposed fault piece")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the three decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
