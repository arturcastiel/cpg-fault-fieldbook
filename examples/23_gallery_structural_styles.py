"""Gallery 1: structural styles built by adapting a grid to fault traces, then sliding blocks.

Each case starts from a layer cake (24x16x4, 50 m columns, layers 15/25/20/30 m, top dipping
1 m per 50 m to the east), adapts it to one or more fault traces with adapt_grid, and moves
the blocks the FAULTS close along their pillars by editing ZCORN only:

- A graben: two curved traces; the middle block drops 25 m;
- B asymmetric horst: the same faults; the flanks drop 15 m (west) and 20 m (east);
- C domino: three oblique traces; each block steps 12 m below its western neighbour;
- D ring collapse: a closed elliptic trace; the block inside drops 30 m (a caldera);
- E reverse fault: the pillars of the east side lean 40 m east by the base, so the fault on
  pillar line i = 12 dips about 66 degrees east; the hanging wall goes UP 20 m along it.

Each case asserts that COORD and FAULTS stay put, exactly the FAULTS faces open, no cell
changes volume, and the re-processed contact areas match the analytic overlaps (or, on the
leaning pillars, the face area per metre times the remaining overlap). Pass --out DIR to keep
the decks (<case>_0_input.grdecl, <case>_1_slid.grdecl).
"""

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import pillar_lattice
from fault_projections import adapt_grid
from fault_projections.from_faults import edge_face

from fieldbook.zslide import (analytic_contacts, blocks, close_ring, column_of_point, compare,
                              interior, layer_cake, named_chains, process, slide_blocks,
                              unwelded_faces, write_deck)

NX, NY, DX = 24, 16, 50.0
DZ = [15.0, 25.0, 20.0, 30.0]
MORPH = dict(boundary="pinned", w0=0.05, polish="barrier", on_invalid="raise")


def base_grid():
    return layer_cake(NX, NY, len(DZ), DX, DZ, top=lambda x, y: 1000.0 + 0.02 * x)


def adapted(grid, traces, ring=False):
    """Adapt to the traces; return the grid and the FAULTS faces named per chain."""
    res = adapt_grid(grid, traces, **MORPH)
    faces = interior(res.grid, res.faults_after)
    if ring:
        faces = close_ring(res.grid, faces)
    shift = np.linalg.norm(pillar_lattice(res.grid).top - pillar_lattice(grid).top, axis=-1)
    print(f"  adapted: {len(faces)} FAULTS faces, pillars moved up to {shift.max():.1f} m, "
          f"mean distance to the traces {res.distance_before.mean():.1f} -> "
          f"{res.distance_after.mean():.1f} m")
    return res.grid, named_chains(res.grid, faces)


def check(name, grid, faults, throws, out, contacts="analytic", per_metre=None):
    """Slide, write both decks, and assert the invariants."""
    faces = [f for chain in faults.values() for f in chain]
    label = blocks(grid, faces)
    slid = slide_blocks(grid, faces, throws)
    assert np.array_equal(slid.coord, grid.coord)
    assert unwelded_faces(slid) == set(faces), "exactly the FAULTS faces must open"
    write_deck(out / f"{name}_0_input.grdecl", grid, faults, header=f"{name}: input")
    write_deck(out / f"{name}_1_slid.grdecl", slid, faults, header=f"{name}: slid")
    before, after = process(grid), process(slid)
    dims = grid.cart_dims
    assert np.allclose(before.volume_grid(dims), after.volume_grid(dims)), "volumes changed"
    found = after.contacts(faces)
    if contacts == "analytic":
        total = compare(found, analytic_contacts(slid, faces))
    else:
        total = sum(found.values())
        assert np.isclose(total, per_metre), (total, per_metre)
    moved = ", ".join(f"block {label[seed]} {t:+g} m" for seed, t in throws.items())
    print(f"{name:18s} {label.max() + 1} blocks; {moved}; fault contact {total:,.0f} m2")


def leaning_cake():
    """Flat layers; pillars from x = 600 m eastward lean 40 m east by the base."""
    lean = 40.0 / sum(DZ)
    return layer_cake(NX, NY, len(DZ), DX, DZ,
                      lean=lambda x, y: (lean * sum(DZ), 0.0) if x >= 600.0 - 1e-6 else (0, 0))


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    grid = base_grid()

    print("A/B: two curved traces")
    west = np.array([[380.0, 0.0], [420.0, 400.0], [400.0, 800.0]])
    east = np.array([[780.0, 0.0], [760.0, 400.0], [800.0, 800.0]])
    g2, faults = adapted(grid, [west, east])
    mid = column_of_point(g2, (600.0, 400.0))
    check("A_graben", g2, faults, {mid: 25.0}, out)
    check("B_horst", g2, faults,
          {column_of_point(g2, (100.0, 400.0)): 15.0, column_of_point(g2, (1100.0, 400.0)): 20.0},
          out)

    print("C: three oblique traces")
    traces = [np.array([[x, 0.0], [x + 300.0, 800.0]]) for x in (150.0, 500.0, 850.0)]
    g3, faults = adapted(grid, traces)
    seeds = [(80.0 + 350.0 * n, 400.0) for n in range(4)]
    check("C_domino", g3, faults,
          {column_of_point(g3, p): 12.0 * n for n, p in enumerate(seeds) if n}, out)

    print("D: closed elliptic trace")
    th = np.linspace(0.0, 2.0 * np.pi, 25)
    ring = np.c_[600.0 + 230.0 * np.cos(th), 400.0 + 190.0 * np.sin(th)]
    g4, faults = adapted(grid, ring, ring=True)
    assert blocks(g4, faults["F1"]).max() == 1, "the ring must close one block"
    check("D_ring_collapse", g4, faults, {column_of_point(g4, (600.0, 400.0)): 30.0}, out)

    print("E: fault on leaning pillars")
    g5 = leaning_cake()
    faces = [edge_face((12, j), (12, j + 1)) for j in range(NY)]
    h, t = sum(DZ), -20.0
    coord = g5.coord.reshape(NY + 1, NX + 1, 6)
    lean = (coord[0, 12, 3:] - coord[0, 12, :3]) / h                # d(xyz)/dz on the fault
    per_metre = sum(float(np.linalg.norm(np.cross([0.0, DX, 0.0], lean))) for _ in faces)
    dip = np.degrees(np.arctan2(h, 40.0))
    print(f"  fault dips {dip:.0f} degrees east; footwall west, hanging wall east goes up")
    check("E_reverse_leaning", g5, {"F1": faces}, {(NX - 1, 0): t}, out,
          contacts="per_metre", per_metre=per_metre * (h - abs(t)))

    print("every case: COORD and FAULTS unchanged, exactly the FAULTS faces open, no cell "
          "changed volume, contact areas as predicted")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
