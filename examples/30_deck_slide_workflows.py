"""Deck workflows: adapt and slide in one call, then slide a field deck's own FAULTS.

- A adapt_deck(slide=...): a 24x16x4 layer cake on disk with staircase FAULTS is adapted to
  two traces and slid in the same call: F1 runs rim to rim, so its east block drops 20 m
  rigidly (carrying F2); F2 ends inside the grid, so its 12 m throw fades to its tip. The
  call writes the grid (adapted and slid), the regenerated FAULTS and the slide report as
  INCLUDE files; the lesson checks that it equals adapt_deck followed by slide_grid, that
  COORD is the adapted one, and that exactly the FAULTS faces open;
- B slide_deck on a field deck: a deck's FAULTS are read through INCLUDE and given a throw by
  name; only ZCORN is written. With no --deck, Norne is used when it is on this machine
  (env CPG_NORNE_GRDECL, the IRAP_1005.GRDECL grid; ACTNUM_0704.prop beside it and
  ../FAULT/FAULT_JUN_05.INC) with E_01 5 m and BC 3 m; otherwise case B is skipped;
- C a dipping fault in one call: a 65 degree planar polygon on a deck with no FAULTS.
  adapt_deck(polygon, top_trace="auto", slide=...) derives the fault's trace on the top
  surface, snaps a pillar chain to it (the FAULTS), tilts the pillars into the plane and
  drops the hanging wall 15 m along them: a dip-slip. The block loses the wedge
  throw x thickness x cot(dip) per metre of fault.

With matplotlib installed, B and C also write a PNG of the top shift and the throw profile.

Pass --out DIR to keep the files, --deck DECK and --throw NAME=METRES (repeat) for your deck.
"""

import argparse
import importlib.util
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import pillar_lattice, read_grdecl
from fault_projections import Throw, adapt_deck, slide_deck, slide_grid, write_slide_png
from fault_projections.adapt import fault_faces_for
from fault_projections.from_faults import faults_from_deck

from fieldbook.zslide import faces_of, layer_cake, process, unwelded_faces, write_deck

NX, NY, DX = 24, 16, 50.0
DZ = [15.0, 25.0, 20.0, 30.0]
MORPH = dict(boundary="pinned", w0=0.05, polish="barrier", on_invalid="raise")
TRACES = {"F1": np.array([[430.0, 0.0], [520.0, 300.0], [560.0, 560.0], [540.0, 800.0]]),
          "F2": np.array([[900.0, 0.0], [880.0, 250.0], [860.0, 450.0]])}
THROWS = [Throw("F1", 20.0, hanging_side=(1100.0, 400.0)),
          Throw("F2", 12.0, hanging_side=(700.0, 200.0))]
NORNE = "/mnt/c/projects/opm-machine/opm-tests/norne/INCLUDE/GRID/IRAP_1005.GRDECL"
PNG = importlib.util.find_spec("matplotlib") is not None


def case_a(out):
    name = "A_adapt_and_slide"
    work = out / name
    work.mkdir(exist_ok=True)
    raw = layer_cake(NX, NY, len(DZ), DX, DZ, top=lambda x, y: 1000.0 + 0.02 * x)
    stair = {n: fault_faces_for(raw, t) for n, t in TRACES.items()}
    deck = write_deck(work / "raw.grdecl", raw, stair, header=f"{name}: staircase FAULTS")

    res = adapt_deck(deck, TRACES, work, slide=THROWS, **MORPH)
    assert res.ok and res.faces_invariant and not res.slide.unexpected
    assert res.slide.modes == {"F1": "rigid", "F2": "decay"}

    # the one call equals the two steps: adapt_deck, then slide_grid on what it wrote
    plain = adapt_deck(deck, TRACES, work / "two_steps", **MORPH)
    faults = faults_from_deck(plain.faults_path)
    two = slide_grid(plain.grid, THROWS, faults=faults, report=False)
    assert np.array_equal(res.grid.zcorn, two.grid.zcorn)
    assert np.array_equal(res.grid.coord, plain.grid.coord), "the slide keeps adapted COORD"

    # read back through a master deck; exactly the FAULTS faces opened, no volume changed
    includes = (res.grid_path.name, res.faults_path.name)
    (work / "case.grdecl").write_text("".join(f"INCLUDE\n  '{p}' /\n" for p in includes))
    grid = read_grdecl(work / "case.grdecl")
    assert np.array_equal(grid.coord, res.grid.coord) and np.array_equal(grid.zcorn,
                                                                         res.grid.zcorn)
    assert faults_from_deck(work / "case.grdecl") == faults
    assert unwelded_faces(grid) == set(res.slide.opened) <= set(faces_of(faults))
    dims = grid.cart_dims
    before, after = process(res.slide.before), process(grid)
    assert np.allclose(after.volume_grid(dims), before.volume_grid(dims))

    moved = int(np.any(res.slide.before.coord != raw.coord, axis=1).sum())
    peak = {n: float(np.abs(res.slide.realized_throw(n)[:, 1]).max()) for n in TRACES}
    print(f"{name:20s} adapt_deck(slide=...) moved {moved} of {(NX + 1) * (NY + 1)} pillars, slid "
          + ", ".join(f"{n} {res.slide.modes[n]} {peak[n]:g} m" for n in TRACES))
    print(f"{'':20s} FAULTS " + ", ".join(f"{n} {len(stair[n])} -> {len(faults[n])}"
                                          for n in TRACES)
          + f" faces; {len(res.slide.opened)} faces opened, all in FAULTS; equals "
          "adapt_deck + slide_grid bit for bit")
    written = (res.grid_path, res.faults_path, res.slide_report_path)
    print(f"{'':20s} wrote {', '.join(p.name for p in written)}")


def norne_deck(out):
    grid = Path(os.environ.get("CPG_NORNE_GRDECL", NORNE))
    files = [grid, grid.parent / "ACTNUM_0704.prop",
             grid.parent.parent / "FAULT" / "FAULT_JUN_05.INC"]
    if not all(f.is_file() for f in files):
        return None
    deck = out / "NORNE_GRID.GRDECL"
    deck.write_text("".join(f"INCLUDE\n  '{f.as_posix()}' /\n\n" for f in files))
    return deck


def case_b(out, deck, throws):
    name = "B_field_deck"
    if deck is None:
        deck = norne_deck(out)
        throws = throws or [Throw("E_01", 5.0), Throw("BC", 3.0)]
        if deck is None:
            print(f"{name:20s} skipped: Norne is not on this machine "
                  "(set CPG_NORNE_GRDECL or pass --deck and --throw)")
            return
    if not throws:
        raise SystemExit("--deck needs at least one --throw NAME=METRES")
    work = out / name
    res = slide_deck(deck, throws, work, prefix="slid", png=PNG)
    assert res.ok and not res.unexpected and res.topology_error is None
    grid = read_grdecl(res.grid_path)
    assert np.array_equal(grid.coord, res.before.coord), "COORD must not change"
    back = slide_grid(res.grid, [Throw(t.fault, -t.throw, hanging_side=t.hanging_side)
                                 for t in throws], faults=res.faults, report=False)
    trip = float(np.abs(back.grid.zcorn - res.before.zcorn).max())
    assert trip < 1e-9, trip
    nx, ny, nz = res.grid.cart_dims
    dv = np.abs(res.volume_change[np.isfinite(res.volume_change)])
    print(f"{name:20s} {deck.name}: {nx}x{ny}x{nz}, {int(res.grid.actnum.sum()):,} active "
          f"cells, {len(res.faults)} named faults")
    print(f"{'':20s} slid " + ", ".join(
        f"{t.fault} {res.modes[t.fault]} {np.abs(res.realized_throw(t.fault)[:, 1]).max():g} m"
        for t in throws)
          + f"; {len(res.opened)} faces opened, none outside FAULTS; round trip {trip:.1e} m")
    print(f"{'':20s} COORD unchanged; cell volume change max {100 * dv.max():.2f} % (tilted "
          f"pillars); {len(res.connections):,} connections in {res.report_path.name}")
    if res.png_path:
        print(f"{'':20s} wrote {res.png_path.name}")


def case_c(out):
    name = "C_dipping_fault"
    work = out / name
    work.mkdir(exist_ok=True)
    raw = layer_cake(NX, NY, len(DZ), DX, DZ, top=lambda x, y: 1000.0)
    deck = write_deck(work / "raw.grdecl", raw, {}, header=f"{name}: no FAULTS")
    dip, strike, p0 = np.radians(65.0), np.radians(20.0), np.array([560.0, 400.0, 1000.0])
    along = np.array([np.sin(strike), np.cos(strike), 0.0])
    down = np.array([np.cos(strike) / np.tan(dip), -np.sin(strike) / np.tan(dip), 1.0])
    normal = np.cross(along, down) / np.linalg.norm(np.cross(along, down))

    def on_plane(depth, s):
        return p0 + (depth - 1000.0) * down + s * along

    polygon = np.array([on_plane(990.0, -2000.0), on_plane(990.0, 2000.0),
                        on_plane(1100.0, 2000.0), on_plane(1100.0, -2000.0)])
    t = 15.0
    res = adapt_deck(deck, polygon, work, top_trace="auto", fault_names="DIP",
                     slide=Throw("DIP", t, hanging_side=(1100.0, 400.0)),
                     **{**MORPH, "boundary": "slide"})
    assert res.ok and res.faces_invariant and res.slide.modes == {"DIP": "rigid"}
    assert unwelded_faces(res.grid) == set(res.slide.opened) == set(res.faults_after)
    trace = np.c_[res.top_trace, np.full(len(res.top_trace), 1000.0)]
    assert np.abs((trace - p0) @ normal).max() < 1e-6, "the auto trace lies on the plane"

    def gap(grid):
        ij = set()
        for f in res.faults_after:
            ij |= {(f.i, f.j - 1), (f.i, f.j)} if f.face == "X" else {(f.i - 1, f.j), (f.i, f.j)}
        i, j = np.array(sorted(ij)).T
        lat = pillar_lattice(grid)
        ends = np.r_[np.c_[lat.top[i, j], lat.z_top[i, j]], np.c_[lat.bot[i, j], lat.z_bot[i, j]]]
        return np.abs((ends - p0) @ normal)

    g1, g2 = gap(res.stage1.grid), gap(res.slide.before)
    assert g1.max() > 30.0 and g2.max() < 0.1 * DX, (g1.max(), g2.max())
    dims = raw.cart_dims
    dv = float(process(res.grid).volume_grid(dims).sum()
               - process(res.slide.before).volume_grid(dims).sum())
    length = res.slide.chains["DIP"].length
    wedge = -t * sum(DZ) / np.tan(dip) * length
    assert abs(dv - wedge) < 0.03 * abs(wedge), (dv, wedge)
    print(f"{name:20s} 65 deg polygon, no FAULTS in the deck: auto top trace of "
          f"{len(res.top_trace)} points on the plane; {len(res.faults_after)} FAULTS faces")
    print(f"{'':20s} chain pillar ends to the plane {g1.max():.0f} m after stage 1 -> max "
          f"{g2.max():.1f}, mean {g2.mean():.1f} m tilted; hanging wall -{t:g} m (rigid)")
    print(f"{'':20s} volume {dv:,.0f} m3 vs the wedge {wedge:,.0f} m3 "
          f"({100 * (dv / wedge - 1):+.1f} %); {len(res.slide.opened)} faces opened, all in FAULTS")
    if PNG:
        png = write_slide_png(res.slide, work / "C_dipping_fault_slide.png")
        print(f"{'':20s} wrote {png.name}")


def parse_throw(text):
    fault, _, metres = text.partition("=")
    return Throw(fault, float(metres))


def main(out, deck, throws):
    out.mkdir(parents=True, exist_ok=True)
    case_a(out)
    case_b(out, deck, throws)
    case_c(out)
    print("every case: COORD stays the adapted one, only FAULTS faces open, files read back")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the files in this directory")
parser.add_argument("--deck", type=Path, help="a GRDECL deck whose FAULTS name the faults")
parser.add_argument("--throw", type=parse_throw, action="append", default=[],
                    metavar="NAME=METRES", help="throw of a named fault (repeat)")
args = parser.parse_args()
if args.deck is not None and not args.deck.is_file():
    parser.error(f"deck does not exist: {args.deck}")
if args.out is not None:
    main(args.out, args.deck, args.throw)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary), args.deck, args.throw)
