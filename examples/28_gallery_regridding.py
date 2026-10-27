"""Gallery 6: regridding and grid adaptation around a sliding fault, from transfer to invariance.

A slide edits ZCORN only, so the fault it moves has exactly the plan shape the grid gives it: a
staircase of cell faces on an unadapted grid, the trace itself once adapt_grid has moved the
pillars onto it. Every case uses a layer cake (layers 15/25/20/30 m, top dipping 1 m per 50 m to
the east, vertical pillars), regrids or adapts it, slides the result and re-processes it:

- A fault transfer: a coarse 12x8 deck with a staircase fault is regridded 3x finer; copying its
  FAULTS onto the fine grid keeps the staircase (fault plane 34% too large), adapting the fine
  grid to the smoothed staircase removes most of the excess (moving a model to a finer grid);
- B resolution study: one curved trace on 8x6, 16x12 and 32x24 grids: the staircase area stays
  |cos|+|sin| times too large at every resolution, the adapted area converges to the exact one
  (choosing a grid resolution around faults);
- C solver settings: five adapt_grid settings (pinned or sliding rim, elastic or Laplacian
  smoothness, barrier polish on or off, starting weight w0 0.05 or a stiff 1) on one partial
  fault, each slid by a throw field: fidelity to the trace against cell quality (tuning the morph);
- D UTM invariance: case C at the origin and at (452000, 6780000): same FAULTS, same slid ZCORN,
  same contacts and volumes up to round-off (decks in map coordinates);
- E active outline: an elliptic ACTNUM outline with an inactive patch; the FAULTS keep only faces
  between active columns and re-processing drops the inactive cells (field outlines, pinch-outs);
- F damage zone: the cell corridor of the slid fault grows into a two-band damage zone whose
  PERMX is cut 10x and 3x and written into the deck (fault-zone permeability).

Pass --out DIR to keep the decks (<case>_0_input.grdecl, <case>_1_slid.grdecl).
"""

import argparse
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import (NonMatchingZError, boundary_mask, cell_quality, format_keyword,
                              pillar_lattice, read_grdecl)
from fault_projections import adapt_grid, targets_from_polyline
from fault_projections.adapt import fault_faces_from_chain
from fault_projections.corridor import cell_corridor
from fault_projections.from_faults import fault_polyline, faults_from_deck, trace_faults
from scipy.ndimage import distance_transform_cdt

from fieldbook.zslide import (Throw, analytic_contacts, blocks, column_of_point, compare,
                              face_columns, interior, layer_cake, process, slide_blocks,
                              slide_field, unwelded_faces, write_deck)

DZ = [15.0, 25.0, 20.0, 30.0]
H = sum(DZ)
THROW = 20.0
MORPH = dict(boundary="pinned", w0=0.05, polish="barrier", on_invalid="raise")
S = np.linspace(0.0, 1.0, 41)
PARTIAL = np.c_[400.0 + 300.0 * S ** 2, 560.0 * S]            # south rim to a tip at (700, 560)
CROSSING = np.c_[350.0 + 225.0 * (1.0 - np.cos(np.pi * S)), 800.0 * S]   # south rim to north rim
UTM = (452000.0, 6780000.0)


def cake(nx, ny, lx=1200.0, ly=800.0, origin=(0.0, 0.0), actnum=None):
    """An nx*ny*4 layer cake over lx*ly m; the top dips 1 m per 50 m to the east."""
    x0 = origin[0]
    return layer_cake(nx, ny, len(DZ), lx / nx, DZ, dy=ly / ny, origin=origin, actnum=actnum,
                      top=lambda x, y: 1000.0 + 0.02 * (x - x0))


def length(line):
    return float(np.linalg.norm(np.diff(line, axis=0), axis=1).sum())


def plan_length(grid, faces):
    """Plan length of a fault: the pillar-to-pillar widths of its faces."""
    return sum(chain.length for chain in trace_faults(grid, faces).chains)


def staircase(grid, trace):
    """FAULTS of a trace on an unadapted grid: the faces between the pillars it claims."""
    return interior(grid, fault_faces_from_chain(targets_from_polyline(grid, trace)))


def adapted(grid, trace, **morph):
    res = adapt_grid(grid, trace, **(morph or MORPH))
    return res, interior(res.grid, res.faults_after)


def active(grid):
    nx, ny, nz = grid.cart_dims
    return grid.actnum.reshape(nz, ny, nx).transpose(2, 1, 0).astype(bool)


def expected_contacts(grid, faces):
    """Analytic contacts between active cells (re-processing drops the inactive ones)."""
    act = active(grid)
    return {(f, ka, kb): area for (f, ka, kb), area in analytic_contacts(grid, faces).items()
            if act[face_columns(f)[0] + (ka,)] and act[face_columns(f)[1] + (kb,)]}


def finish(name, grid, slid, faces, out, *, opened=None, extra=""):
    """Write both decks (if named), assert the slide invariants, return the contact area."""
    assert np.array_equal(slid.coord, grid.coord), "a slide must not move COORD"
    assert unwelded_faces(slid) == set(opened or faces), "exactly the cut faces must open"
    if name:
        for tag, g in (("0_input", grid), ("1_slid", slid)):
            write_deck(out / f"{name}_{tag}.grdecl", g, {"F1": faces}, extra=extra,
                       header=f"lesson 28 {name}: {tag[2:]}")
    before, after = process(grid), process(slid)
    dims = grid.cart_dims
    assert np.allclose(before.volume_grid(dims), after.volume_grid(dims), equal_nan=True), \
        "a cell changed volume"
    return after, compare(after.contacts(faces), expected_contacts(slid, faces))


def smooth(xy, width, step):
    """Moving average of a polyline over ``width`` m of arc length, every ``step`` m; ends kept."""
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))])
    n = int(np.ceil(s[-1] / step))
    t = np.linspace(0.0, s[-1], 8 * n + 1)
    dense = np.c_[np.interp(t, s, xy[:, 0]), np.interp(t, s, xy[:, 1])]
    out = []
    for u in np.linspace(0.0, s[-1], n + 1):
        half = min(width / 2.0, u, s[-1] - u)               # the window shrinks at the ends
        out.append(dense[np.abs(t - u) <= half + 1e-9].mean(axis=0))
    return np.array(out)


def case_a(out):
    name = "A_fault_transfer"
    trace = np.array([[300.0, 0.0], [700.0, 800.0]])
    coarse = cake(12, 8)
    deck = write_deck(out / f"{name}_coarse.grdecl", coarse, {"F1": staircase(coarse, trace)},
                      header=f"lesson 28 {name}: the coarse 12x8 deck")
    xy = fault_polyline(coarse, faults_from_deck(deck)["F1"])      # all the coarse deck knows
    fine = cake(36, 24)
    stair = staircase(fine, xy)
    res, faces = adapted(fine, smooth(xy, 300.0, 100.0))
    true = length(trace)
    lengths = {}
    for label, grid, f in (("staircase", fine, stair), ("adapted", res.grid, faces)):
        slid = slide_blocks(grid, f, {column_of_point(grid, (1100.0, 400.0)): THROW})
        _, contact = finish(name if label == "adapted" else None, grid, slid, f, out)
        lengths[label] = plan_length(grid, f)
        assert np.isclose(contact, lengths[label] * (H - THROW)), (label, contact)
        if label == "staircase":
            write_deck(out / f"{name}_staircase_slid.grdecl", slid, {"F1": f},
                       header=f"lesson 28 {name}: re-staircased and slid")
    assert np.isclose(lengths["staircase"], length(xy)), "re-staircasing copies the coarse steps"
    assert abs(lengths["adapted"] - true) < 0.1 * abs(lengths["staircase"] - true)
    err = {k: 100.0 * (v / true - 1.0) for k, v in lengths.items()}
    print(f"{name:20s} 12x8 staircase ({len(xy) - 1} faces) -> 36x24 grid; true trace "
          f"{true:.0f} m")
    print(f"{'':20s} fault plan {lengths['staircase']:.0f} m re-staircased "
          f"({err['staircase']:+.1f}%), {lengths['adapted']:.0f} m adapted to smoothed steps "
          f"({err['adapted']:+.1f}%)")
    print(f"{'':20s} contact after a {THROW:g} m slide {lengths['staircase'] * (H - THROW):,.0f}"
          f" -> {lengths['adapted'] * (H - THROW):,.0f} m2")


def case_b(out):
    name = "B_resolution_study"
    trace = np.c_[300.0 + 300.0 * (1.0 - np.cos(np.pi * S)), 900.0 * S]
    exact = length(trace) * H
    factor = np.abs(np.diff(trace, axis=0)).sum() / length(trace)  # mean of |cos| + |sin|
    rows = []
    for n in (8, 16, 32):
        grid = cake(n, 3 * n // 4, ly=900.0)
        stair = staircase(grid, trace)
        res, faces = adapted(grid, trace)
        assert set(faces) == set(stair), "adapting moves pillars, not FAULTS faces"
        area_stair = sum(process(grid).contacts(stair).values())
        area_adapt = sum(process(res.grid).contacts(faces).values())
        rows.append((n, 3 * n // 4, 1200.0 / n, area_stair, area_adapt))
    slid = slide_blocks(res.grid, faces, {column_of_point(res.grid, (1100.0, 450.0)): THROW})
    finish(name, res.grid, slid, faces, out)
    with open(out / f"{name}.csv", "w") as csv:
        csv.write("nx,ny,cell_m,staircase_m2,adapted_m2,exact_m2\n")
        csv.writelines(f"{nx},{ny},{h:g},{a:.3f},{b:.3f},{exact:.3f}\n"
                       for nx, ny, h, a, b in rows)
    err = np.array([abs(b / exact - 1.0) for *_, b in rows])
    for *_, a, b in rows:
        assert abs(a / exact - factor) < 0.02 * factor, (a / exact, factor)
    assert np.all(np.diff(err) < 0.0) and err[-1] < 0.03, err
    print(f"{name:20s} exact {exact:,.0f} m2; staircase / adapted fault area per grid "
          f"(|cos|+|sin| = {factor:.3f}):")
    print(f"{'':20s} " + "; ".join(f"{nx}x{ny} {a / exact:.3f} / {b / exact:.3f}"
                                  for nx, ny, _, a, b in rows))


def partial_case(origin=(0.0, 0.0), **morph):
    """Case C's grid adapted to the partial fault and slid by a throw field (18 m peak)."""
    o = np.asarray(origin)
    grid = cake(24, 16, origin=origin)
    res, faces = adapted(grid, PARTIAL + o, **morph)
    spec = Throw(faces=tuple(faces), throw=18.0, hanging_side=tuple(o + (900.0, 300.0)))
    return grid, res, faces, slide_field(res.grid, spec, faults=faces)


SETTINGS = {
    "default": MORPH,
    "slide": dict(MORPH, boundary="slide"),
    "laplacian": dict(MORPH, smoothness="laplacian"),
    "stiff_no_polish": dict(MORPH, polish="none", w0=1.0),
    "no_polish": dict(MORPH, polish="none"),
}


def case_c(out):
    name = "C_solver_settings"
    stats, moves = {}, []
    for n, (tag, morph) in enumerate(SETTINGS.items()):
        grid, res, faces, slid = partial_case(**morph)
        assert res.ok and trace_faults(res.grid, faces).is_simple
        quality = cell_quality(pillar_lattice(res.grid), plane="top")
        move = np.linalg.norm(pillar_lattice(res.grid).top - pillar_lattice(grid).top, axis=-1)
        finish(None if n else name, res.grid, slid, faces, out)
        if n:
            write_deck(out / f"{name}_{tag}_slid.grdecl", slid, {"F1": faces},
                       header=f"lesson 28 {name}: {tag}, slid")
        else:
            reference = (res, faces, slid)
        assert quality.sj_min > 0.0, f"{tag}: a cell degenerated"
        d = res.distance_after
        stats[tag] = (d.mean(), d.max(), quality.sj_min, quality.sj_below_0_3)
        moves.append(move.max())
    dist = {tag: s[0] for tag, s in stats.items()}
    sj = {tag: s[2] for tag, s in stats.items()}
    assert dist["no_polish"] < dist["stiff_no_polish"] < dist["default"], "fidelity order"
    assert sj["no_polish"] < sj["stiff_no_polish"] < 0.3 <= sj["default"], "quality order"
    assert all(stats[tag][3] == 0 for tag in ("default", "slide", "laplacian"))
    text = [f"{tag.replace('_', ' ')} {m:.2f}/{x:.2f} m {s:.3g}" + (f" ({k})" if k else "")
            for tag, (m, x, s, k) in stats.items()]
    print(f"{name:20s} distance mean/max, sj_min (cells < 0.3); pillars move <= {max(moves):.0f} "
          "m:")
    print(f"{'':20s} " + ", ".join(text[:3]))
    print(f"{'':20s} " + ", ".join(text[3:]))
    return reference


def case_d(out, reference):
    name = "D_utm_invariance"
    res0, faces0, slid0 = reference
    grid, res, faces, slid = partial_case(UTM)
    after, contact = finish(name, res.grid, slid, faces, out)
    after0 = process(slid0)
    assert faces == faces0, "the FAULTS must not depend on the origin"
    # float64 resolves 6.78e6 m to ~1e-9 m; the morph solve and the volume sums lose a few
    # more digits at that offset, so coordinates get 1 micron and the rest 1e-8 relative.
    dxy = np.abs(pillar_lattice(res.grid).top - UTM - pillar_lattice(res0.grid).top).max()
    dz = np.abs(slid.zcorn - slid0.zcorn).max()
    v, v0 = after.volume_grid(grid.cart_dims), after0.volume_grid(grid.cart_dims)
    c, c0 = after.contacts(faces), after0.contacts(faces0)
    dc = max(abs(c[k] - c0[k]) for k in c0)
    assert dxy < 1e-6 and dz <= 1e-9 * np.abs(slid0.zcorn).max(), (dxy, dz)
    assert np.allclose(v, v0, rtol=1e-8, atol=0.0) and c.keys() == c0.keys()
    assert all(np.isclose(c[k], c0[k], rtol=1e-8, atol=1e-6) for k in c0)
    print(f"{name:20s} origin ({UTM[0]:.0f}, {UTM[1]:.0f}) vs (0, 0): same {len(faces)} FAULTS; "
          "largest differences:")
    print(f"{'':20s} pillar xy {dxy:.1e} m, slid ZCORN {dz:.1e} m, cell volume "
          f"{np.nanmax(np.abs(v - v0) / v0):.1e} relative,")
    print(f"{'':20s} contact {dc:.1e} m2 of {contact:,.0f} m2")


def case_e(out):
    name = "E_actnum_outline"
    x, y = np.meshgrid((np.arange(24) + 0.5) * 50.0, (np.arange(16) + 0.5) * 50.0, indexing="ij")
    inside = ((x - 600.0) / 560.0) ** 2 + ((y - 400.0) / 370.0) ** 2 <= 1.0
    act = np.repeat(inside[:, :, None], len(DZ), axis=2)
    act[17:21, 6:10, 1] = False                       # a pinch-out patch in the hanging wall
    grid = cake(24, 16, actnum=act)
    res, cut = adapted(grid, CROSSING)
    rim = boundary_mask(grid, scope="active")
    rim_move = np.linalg.norm(pillar_lattice(res.grid).top - pillar_lattice(grid).top,
                              axis=-1)[rim].max()
    column = act.any(axis=2)
    faces = [f for f in cut if all(column[c] for c in face_columns(f))]
    slid = slide_blocks(res.grid, cut, {column_of_point(res.grid, (1000.0, 400.0)): THROW})
    after, contact = finish(name, res.grid, slid, faces, out, opened=cut)
    assert len(after.volumes) == act.sum(), "re-processing must drop the inactive cells"
    assert np.array_equal(np.isnan(after.volume_grid(grid.cart_dims)), ~act)
    assert all(not all(column[c] for c in face_columns(f)) for f in set(cut) - set(faces))
    # A pinch-out against the fault: cpg2unstructured needs all-active faulted stacks.
    act2 = act.copy()
    for c in face_columns(faces[len(faces) // 2]):
        act2[c + (1,)] = False
    try:
        process(replace(slid, actnum=act2.transpose(2, 1, 0).ravel().astype(np.int8)))
        raise AssertionError("an inactive cell in a faulted stack should be rejected")
    except NonMatchingZError:
        pass
    print(f"{name:20s} {column.sum()}/{column.size} columns, {len(after.volumes)} cells active, "
          f"all re-processed; outline pillars moved {rim_move:.1f} m")
    print(f"{'':20s} {len(cut)} faces cut, the {len(faces)} between active columns in FAULTS; "
          f"contact {contact:,.0f} m2")
    print(f"{'':20s} an inactive cell beside the slid fault: process raises NonMatchingZError")


def case_f(out):
    name = "F_damage_zone"
    grid = cake(24, 16)
    res, faces = adapted(grid, CROSSING)
    slid = slide_blocks(res.grid, faces, {column_of_point(res.grid, (1000.0, 400.0)): THROW})
    corridor = cell_corridor(slid, CROSSING)
    assert corridor.ok, corridor.summary()
    assert np.array_equal(corridor.cells, cell_corridor(res.grid, CROSSING).cells)
    nx, ny = grid.cart_dims[:2]
    steps = distance_transform_cdt(~corridor.mask((nx, ny)), metric="taxicab")
    label = blocks(slid, faces)
    walls = {int(label[tuple(c)]) for c in corridor.cells}
    assert walls == {0, 1}, "the corridor must straddle the fault"
    assert all(steps[c] <= 1 for f in faces for c in face_columns(f)), "both walls are damaged"
    mult = np.where(steps <= 1, 0.1, np.where(steps == 2, 1.0 / 3.0, 1.0))
    permx = mult[:, :, None] * np.array([250.0, 80.0, 150.0, 40.0])[None, None, :]
    values = permx.transpose(2, 1, 0).ravel()
    finish(name, res.grid, slid, faces, out, extra=format_keyword("PERMX", values))
    back = read_grdecl(out / f"{name}_1_slid.grdecl").props["PERMX"]
    assert np.array_equal(back, values), "PERMX must read back unchanged"
    counts = [int((steps == d).sum()) for d in (0, 1, 2)]
    print(f"{name:20s} corridor {counts[0]} columns ({corridor.summary().split(' · ')[1]}) on "
          "both walls of the slid fault")
    print(f"{'':20s} PERMX x0.1 within 1 cell ({counts[0] + counts[1]} columns), x1/3 at 2 "
          f"({counts[2]}); {back.size} values read back")


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    case_a(out)
    case_b(out)
    reference = case_c(out)
    case_d(out, reference)
    case_e(out)
    case_f(out)
    print("every case: COORD and FAULTS unchanged, exactly the cut faces open, no cell changed "
          "volume, contacts as predicted")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
