"""Gallery 2: fault networks — crossing, abutting, splaying and overlapping faults slide together.

Each case starts from the layer cake of gallery 1 (24x16x4, 50 m columns, layers 15/25/20/30 m,
top dipping 1 m per 50 m to the east), adapts it to all the traces of a fault network at once
with adapt_grid (a pillar claimed by two traces is solved for both), reads each trace's FAULTS
off its own pillar chain, and slides the network by editing ZCORN only:

- A X crossing: two traces cross at one pillar and cut four blocks, which drop 0, 10, 25 and
  15 m; the crossing pillar then carries four depths per horizon, one per welded group (a
  structural framework with intersecting faults);
- B T junction: a fault abuts a through-going fault; three blocks, and the main fault's throw
  steps at the junction by exactly the abutting fault's throw (displacement transfer);
- C Y splay: a splay branches off a through-going fault and dies inside its hanging wall; the
  main fault slides its block rigidly, the splay gets a tapered throw field, and the two edits
  commute (a splay that compartmentalises only part of a block);
- D relay ramp: two en-echelon faults, each with two tips, overlap in plan; both throw fields
  in one call keep the grid one block and tilt the ramp between them (the ramp is the
  connection, or leak point, between the two fault blocks);
- E compartments: five traces cut eight compartments, each with its own throw; a table of
  columns, bulk volume and throw per compartment (also E_compartments.csv), and a FIPNUM
  region per compartment in both decks for in-place volume reports.

Each case asserts that COORD and FAULTS (read back from both decks) stay put, exactly the
FAULTS faces open, no cell changes volume, and the re-processed decks give the contact areas
of the analytic overlaps. Pass --out DIR to keep the decks (<case>_0_input.grdecl,
<case>_1_slid.grdecl).
"""

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import format_keyword, pillar_lattice
from fault_projections import adapt_grid, targets_from_polylines
from fault_projections.adapt import fault_faces_from_chain
from fault_projections.from_faults import faults_from_deck, trace_faults

from fieldbook.zslide import (Throw, analytic_contacts, blocks, column_of_point, compare,
                              corner_cube, faces_of, layer_cake, process,
                              realized_throw, slide_blocks, slide_field, unwelded_faces,
                              write_deck)

NX, NY, DX = 24, 16, 50.0
DZ = [15.0, 25.0, 20.0, 30.0]
MORPH = dict(boundary="pinned", w0=0.05, polish="barrier", on_invalid="raise")


def base_grid():
    return layer_cake(NX, NY, len(DZ), DX, DZ, top=lambda x, y: 1000.0 + 0.02 * x)


def adapted(grid, traces):
    """Adapt to all traces at once; return the grid and the FAULTS faces of each trace."""
    res = adapt_grid(grid, list(traces.values()), **MORPH)
    # faults_after is read off all chains joined end to end; each trace's own chain names it.
    faults = {name: fault_faces_from_chain(targets_from_polylines(grid, [t]))
              for name, t in traces.items()}
    faces = faces_of(faults)
    assert len(set(faces)) == len(faces), "two traces share a face"
    assert set(faces) == set(res.faults_after)
    assert all(trace_faults(res.grid, f).is_simple for f in faults.values())
    shift = np.linalg.norm(pillar_lattice(res.grid).top - pillar_lattice(grid).top, axis=-1)
    print(f"  adapted: {len(faces)} FAULTS faces on {len(faults)} traces, pillars moved up to "
          f"{shift.max():.1f} m, mean distance to the traces {res.distance_before.mean():.1f} -> "
          f"{res.distance_after.mean():.1f} m")
    return res.grid, faults


def junction(grid, faults):
    """The one pillar where three or more FAULTS faces meet."""
    branch = trace_faults(grid, faces_of(faults)).branch_pillars
    assert len(branch) == 1, branch
    return tuple(int(v) for v in branch[0])


def depths_at(grid, pillar, k=0):
    """The distinct depths of the top of layer k on one pillar: one per welded group."""
    z, (pi, pj) = corner_cube(grid), pillar
    return sorted({round(float(z[2 * k, 2 * j + pj - j, 2 * i + pi - i]), 6)
                   for i in (pi - 1, pi) for j in (pj - 1, pj) if 0 <= i < NX and 0 <= j < NY})


def check(name, grid, faults, slid, out, extra=""):
    """Write both decks, re-process them and assert the invariants.

    Returns the cell volumes before and after, and the total contact area across the FAULTS.
    """
    faces = faces_of(faults)
    assert np.array_equal(slid.coord, grid.coord)
    assert unwelded_faces(slid) == set(faces), "exactly the FAULTS faces must open"
    decks = []
    for n, (stage, g) in enumerate((("input", grid), ("slid", slid))):
        deck = write_deck(out / f"{name}_{n}_{stage}.grdecl", g, faults,
                          header=f"{name}: {stage}", extra=extra)
        on_deck = {key: set(v) for key, v in faults_from_deck(deck).items()}
        assert on_deck == {key: set(v) for key, v in faults.items()}, "FAULTS changed"
        decks.append(deck)
    before, after = process(decks[0]), process(decks[1])
    volumes = before.volume_grid(grid.cart_dims), after.volume_grid(grid.cart_dims)
    assert np.allclose(*volumes), "volumes changed"
    return volumes, compare(after.contacts(faces), analytic_contacts(slid, faces))


def seeded(grid, faults, points):
    """``{seed column: throw}`` from ``{plan point: throw}``, one seed per block."""
    label = blocks(grid, faces_of(faults))
    seeds = {column_of_point(grid, p): t for p, t in points.items()}
    assert len({label[s] for s in seeds}) == len(seeds) == label.max() + 1
    return seeds


def case_a(grid, out):
    print("A: two traces crossing in an X")
    g, faults = adapted(grid, {"F1": np.array([[400.0, 0.0], [800.0, 800.0]]),
                               "F2": np.array([[0.0, 600.0], [1200.0, 200.0]])})
    throws = {(100.0, 400.0): 0.0, (600.0, 50.0): 10.0, (1100.0, 400.0): 25.0,
              (600.0, 750.0): 15.0}                                   # W, S, E, N
    slid = slide_blocks(g, faces_of(faults), seeded(g, faults, throws))
    cross = junction(g, faults)
    depths = depths_at(slid, cross)
    assert np.allclose(np.diff(depths), [10.0, 5.0, 10.0]), depths    # 0, 10, 15, 25 m
    _, total = check("A_x_crossing", g, faults, slid, out)
    print(f"A_x_crossing       4 blocks W 0, S +10, E +25, N +15 m; fault contact {total:,.0f} m2")
    print(f"  crossing pillar {cross} carries {len(depths)} depths of the top: "
          f"{', '.join(f'{d:.1f}' for d in depths)} m, one per welded group")


def case_b(grid, out):
    print("B: a fault abutting a through-going fault (T junction)")
    g, faults = adapted(grid, {"MAIN": np.array([[0.0, 600.0], [1200.0, 200.0]]),
                               "ABUT": np.array([[600.0, 400.0], [800.0, 800.0]])})
    tee = junction(g, faults)
    abut = trace_faults(g, faults["ABUT"]).chains[0]
    assert tee in {tuple(abut.ij[0]), tuple(abut.ij[-1])}, "ABUT must end on MAIN"
    throws = {(600.0, 50.0): 0.0, (300.0, 650.0): 20.0, (900.0, 650.0): 8.0}   # S, NW, NE
    slid = slide_blocks(g, faces_of(faults), seeded(g, faults, throws))
    assert len(depths_at(slid, tee)) == 3
    main = realized_throw(g, slid, faults["MAIN"], hanging_side=(600.0, 750.0))[:, 1]
    across = realized_throw(g, slid, faults["ABUT"], hanging_side=(300.0, 650.0))[:, 1]
    assert np.allclose(across, 12.0) and np.isclose(main.max() - main.min(), 12.0)
    _, total = check("B_t_junction", g, faults, slid, out)
    print(f"B_t_junction       3 blocks S 0, NW +20, NE +8 m; fault contact {total:,.0f} m2")
    print(f"  MAIN throw {main.max():.0f} m west of the junction {tee}, {main.min():.0f} m east: "
          f"the step is ABUT's {np.median(across):.0f} m")


def case_c(grid, out):
    print("C: a splay branching off a through-going fault (Y)")
    g, faults = adapted(grid, {"MAIN": np.array([[0.0, 600.0], [1200.0, 200.0]]),
                               "SPLAY": np.array([[600.0, 400.0], [750.0, 700.0]])})
    faces = faces_of(faults)
    seeds = seeded(g, faults, {(600.0, 50.0): 0.0, (600.0, 750.0): 30.0})    # S, N
    # The splay ends inside block N, so it needs a field. The field also adds to MAIN's throw
    # east of the junction (30 -> 38 m), kept clear of the 40 m between horizons 0 and 2:
    # slide_field snaps only the splay's own faces, not MAIN's.
    splay = Throw(faces=tuple(faults["SPLAY"]), throw=8.0, hanging_side=(800.0, 500.0))
    rigid = slide_blocks(g, faces, seeds)
    slid = slide_field(rigid, splay, faults=faces)
    # Both edits add a depth shift per corner, so the order does not matter.
    other = slide_blocks(slide_field(g, splay, faults=faces), faces, seeds)
    assert np.allclose(other.zcorn, slid.zcorn, rtol=0.0, atol=1e-9)
    along = realized_throw(rigid, slid, faults["SPLAY"], hanging_side=splay.hanging_side)
    assert np.isclose(along[0, 1], 8.0) and along[-1, 1] == 0.0       # junction, tip
    wye = junction(g, faults)
    depths = depths_at(slid, wye)                    # S, N west and N east of the splay
    assert np.allclose(np.diff(depths), [30.0, 8.0]), depths
    _, total = check("C_y_splay", g, faults, slid, out)
    print(f"C_y_splay          2 blocks S 0, N +30 m, then the SPLAY field; fault contact "
          f"{total:,.0f} m2")
    print(f"  SPLAY {along[0, 1]:.1f} m at the junction {wye} -> {along[-1, 1]:.1f} m at its tip "
          f"{along[-1, 0]:.0f} m away; MAIN {depths[1] - depths[0]:.0f} m, "
          f"{depths[2] - depths[0]:.0f} m at the junction east of SPLAY")


def case_d(grid, out):
    print("D: two overlapping en-echelon faults (relay ramp)")
    g, faults = adapted(grid, {"F1": np.array([[390.0, 100.0], [410.0, 300.0], [395.0, 500.0]]),
                               "F2": np.array([[540.0, 300.0], [560.0, 500.0], [545.0, 700.0]])})
    faces = faces_of(faults)
    assert blocks(g, faces).max() == 0, "two faults with tips cut no block"
    specs = [Throw(faces=tuple(faults["F1"]), throw=12.0, hanging_side=(500.0, 300.0)),
             Throw(faces=tuple(faults["F2"]), throw=12.0, hanging_side=(650.0, 500.0))]
    slid = slide_field(g, specs, faults=faces)
    peaks = [realized_throw(g, slid, faults[name], spec.hanging_side)[:, 1]
             for name, spec in zip(faults, specs)]
    # zero at both tips; the cosine peak sits at mid-length, between two pillars
    assert all(11.99 < p.max() <= 12.0 and np.allclose(p[[0, -1]], 0.0) for p in peaks)
    shift = corner_cube(slid)[0] - corner_cube(g)[0]
    ramp = []
    for y in (325.0, 375.0, 425.0, 475.0):
        i, j = column_of_point(g, (475.0, y))
        ramp.append(shift[2 * j:2 * j + 2, 2 * i:2 * i + 2].mean())
    assert np.all(np.diff(ramp) < 0.0), "the ramp rises from F1's hanging wall to F2's footwall"
    _, total = check("D_relay_ramp", g, faults, slid, out)
    tilt = (ramp[0] - ramp[-1]) / 150.0
    print(f"D_relay_ramp       1 block; F1, F2 peak {max(peaks[0]):.1f}, {max(peaks[1]):.1f} m; "
          f"fault contact {total:,.0f} m2")
    print(f"  ramp between the faults: top down {ramp[0]:.1f} m at y = 325, {ramp[-1]:.1f} m at "
          f"y = 475; it dips {np.degrees(np.arctan(tilt)):.1f} deg ({1e3 * tilt:.0f} m/km) south")


def case_e(grid, out):
    print("E: five traces cut the grid into compartments")
    # Three vertices per trace: adapt_grid classifies its input with np.asarray(list of traces),
    # which fails on traces with different vertex counts.
    g, faults = adapted(grid, {
        "F1": np.array([[300.0, 0.0], [330.0, 400.0], [280.0, 800.0]]),
        "F2": np.array([[850.0, 0.0], [880.0, 400.0], [830.0, 800.0]]),
        "F3": np.array([[0.0, 450.0], [600.0, 420.0], [1200.0, 470.0]]),
        "F4": np.array([[580.0, 421.0], [630.0, 200.0], [650.0, 0.0]]),
        "F5": np.array([[849.0, 650.0], [1000.0, 690.0], [1200.0, 700.0]])})
    faces = faces_of(faults)
    throws = {(150.0, 200.0): 0.0, (150.0, 700.0): 10.0, (450.0, 200.0): 35.0,
              (600.0, 700.0): 20.0, (750.0, 200.0): 28.0, (1050.0, 200.0): 5.0,
              (1050.0, 550.0): 15.0, (1050.0, 750.0): 25.0}
    seeds = seeded(g, faults, throws)
    label = blocks(g, faces)
    assert label.max() + 1 == 8
    slid = slide_blocks(g, faces, seeds)
    fipnum = format_keyword("FIPNUM", np.tile(label.T.ravel() + 1, len(DZ)), integer=True)
    (before, after), total = check("E_compartments", g, faults, slid, out, extra=fipnum)
    rows = []
    for seed, throw in seeds.items():
        mask = label == label[seed]
        bulk0, bulk1 = before[mask].sum(), after[mask].sum()
        assert np.isclose(bulk0, bulk1, rtol=1e-12), "a compartment changed volume"
        rows.append((label[seed] + 1, int(mask.sum()), bulk1, throw))
    rows.sort()
    assert np.isclose(sum(r[2] for r in rows), NX * NY * DX * DX * sum(DZ))
    with open(out / "E_compartments.csv", "w") as csv:
        csv.write("fipnum,columns,bulk_volume_m3,throw_m\n")
        csv.writelines(f"{n},{c},{v:.3f},{t:g}\n" for n, c, v, t in rows)
    print("  FIPNUM  columns  bulk volume [1e6 m3]  throw [m]")
    for n, c, v, t in rows:
        print(f"  {n:6d}  {c:7d}  {v / 1e6:20.3f}  {t:+9.0f}")
    print(f"E_compartments     {len(rows)} compartments; bulk volume {after.sum() / 1e6:.3f}e6 m3 "
          f"before and after; fault contact {total:,.0f} m2")


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    grid = base_grid()
    case_a(grid, out)
    case_b(grid, out)
    case_c(grid, out)
    case_d(grid, out)
    case_e(grid, out)
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
