"""Gallery 5: juxtaposition, fault seal, connectivity, uncertainty and wells on a slid grid.

Every case starts from one sand/shale layer cake (24x16x6, 50 m columns, top at 2000 m dipping
1 m per 50 m to the east): Cap shale 20 m, Sand A 30 m (the main reservoir), Shale 1 10 m,
Sand B 15 m, Shale 2 25 m, Sand C 20 m, with VSH 0.7 in the shales and 0.1 in the sands.
adapt_grid moves the pillars onto a fault trace; the slide then edits ZCORN only, and
re-processing the slid grid gives the cell-to-cell contact areas across the fault, from which
the reservoir questions are answered:

- A Allan diagram: an oblique fault crossing the grid, east side down 40 m (more than any layer
  is thick); the footwall-to-hanging-wall juxtaposition table and its sand-sand windows, for
  juxtaposition seal and leak-point analysis;
- B fault seal: the throw grows northward from 20 to 50 m (a throw table); the shale gouge ratio
  of every cell face in contact and a Manzocchi transmissibility multiplier written as
  MULTX / MULTY, for fault-rock seal in a simulation deck;
- C connectivity sweep: throws 0-80 m; the sand-sand contact area and the throw at which the
  footwall Sand A loses its last sand window, for cross-fault connectivity screening;
- D throw ensemble: 40 lognormal throws; P90 / P50 / P10 sand-sand area and the share of
  realizations in which Sand A stays connected, for throw uncertainty (the realization closest
  to P50 is re-processed);
- E well tops: a growth fault (20 m throw at the top, 44 m at the base); formation tops in a
  hanging-wall and a footwall well before and after, for well correlation across the fault;
- F partial fault: a fault ending inside the grid; hanging-wall wells further and further from
  it see the shift decay to zero, for well ties near a fault tip.

Shale gouge ratio at a point of the fault (Yielding et al. 1997): SGR = sum(VSH_i * dz_i) / D
over the interval of thickness D (the local throw) that slipped past the point; a cell face
gets its mean over the part in contact. Multiplier (Manzocchi et al. 1999), D in m, k in mD:
k_f = 10^(-4 SGR - 0.25 log10(D) (1 - SGR)^5), t_f = D / 66,
T = 1 / (1 + 2 t_f / L (k_m / k_f - 1)) with L = 50 m and k_m = 100 mD; T = 1 if k_f >= k_m.
MULTX(i, j, k) scales the connections of the +X face of cell (i, j, k), MULTY the +Y face.

Pass --out DIR to keep the decks (<case>_0_input.grdecl, <case>_1_slid.grdecl).
"""

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import format_keyword, pillar_lattice, read_grdecl
from fault_projections import adapt_grid
from fault_projections.from_faults import face_edge, faults_from_deck, trace_faults

from fieldbook.zslide import (Throw, analytic_contacts, blocks, column_of_point, compare,
                              contacts_by_layer, corner_cube, face_columns, interior, layer_cake,
                              named_chains, process, slide_blocks, slide_field, unwelded_faces,
                              write_deck)

NX, NY, DX = 24, 16, 50.0
NAMES = ["Cap shale", "Sand A", "Shale 1", "Sand B", "Shale 2", "Sand C"]
DZ = np.array([20.0, 30.0, 10.0, 15.0, 25.0, 20.0])
VSH = np.array([0.7, 0.1, 0.7, 0.1, 0.7, 0.1])
NZ, SAND, MAIN = len(DZ), VSH < 0.5, 1
Z = np.concatenate([[0.0], np.cumsum(DZ)])          # horizons below the top of the cake
H = Z[-1]
TRACE = np.array([[480.0, 0.0], [600.0, 400.0], [700.0, 800.0]])
HANGING = (1000.0, 400.0)                           # a plan point east of the fault
K_MATRIX, CELL_L = 100.0, DX
MORPH = dict(boundary="pinned", w0=0.05, polish="barrier", on_invalid="raise")


def base_grid():
    return layer_cake(NX, NY, NZ, DX, DZ, top=lambda x, y: 2000.0 + 0.02 * x)


def adapted(grid, trace):
    res = adapt_grid(grid, trace, **MORPH)
    faces = interior(res.grid, res.faults_after)
    print(f"  adapted: {len(faces)} FAULTS faces, mean distance to the trace "
          f"{res.distance_before.mean():.1f} -> {res.distance_after.mean():.1f} m")
    return res.grid, faces


class Fault:
    """An adapted grid, its FAULTS, the hanging-wall block and the fault length in plan."""

    def __init__(self, grid, faces):
        self.grid, self.faces, self.faults = grid, faces, named_chains(grid, faces)
        label = blocks(grid, faces)
        assert label.max() == 1, "the fault must cut the grid in two"
        self.seed = column_of_point(grid, HANGING)
        self.hw = label == label[self.seed]
        self.flip = [f for f in faces if self.hw[face_columns(f)[0]]]  # read FW -> HW
        top = pillar_lattice(grid).top
        self.width = {f: float(np.linalg.norm(top[p[0], p[1]] - top[q[0], q[1]]))
                      for f, (p, q) in zip(faces, map(face_edge, faces))}
        self.length = sum(self.width.values())
        self.volumes = process(grid).volume_grid(grid.cart_dims)

    def slide(self, throw):
        """Move the hanging-wall block down ``throw`` (one value, or one per horizon)."""
        slid = slide_blocks(self.grid, self.faces, {self.seed: throw})
        assert np.array_equal(slid.coord, self.grid.coord), "COORD changed"
        assert unwelded_faces(slid) == (set(self.faces) if np.any(throw) else set())
        return slid

    def table(self, contacts):
        """Juxtaposition table ``[k footwall, k hanging wall]`` of contact areas (m2)."""
        out = np.zeros((NZ, NZ))
        for (kf, kh), area in contacts_by_layer(contacts, self.flip).items():
            out[kf, kh] = area
        return out

    def decks(self, name, slid, out, note, extra=""):
        write_deck(out / f"{name}_0_input.grdecl", self.grid, self.faults, header=f"{name}: input")
        return write_deck(out / f"{name}_1_slid.grdecl", slid, self.faults,
                          header=f"{name}: {note}", extra=extra)


def overlap(throw):
    """Thickness (m) of footwall layer k facing hanging-wall layer m, the HW ``throw`` down."""
    zh = Z + throw
    return np.clip(np.minimum.outer(Z[1:], zh[1:]) - np.maximum.outer(Z[:-1], zh[:-1]), 0.0, None)


def sand_sand(tab):
    return float(tab[np.ix_(SAND, SAND)].sum())


def print_allan(tab):
    """Rows: footwall layers; columns: hanging-wall layers; * marks a sand-sand window."""
    print(f"  {'FW | HW, 1e3 m2':17s}" + "".join(f"{n:>9s} " for n in NAMES).rstrip())
    for kf in range(NZ):
        cells = [(f"{tab[kf, kh] / 1e3:9.1f}" + ("*" if SAND[kf] and SAND[kh] else " "))
                 if tab[kf, kh] > 1e-6 else "        . " for kh in range(NZ)]
        print((f"  {NAMES[kf]:17s}" + "".join(cells)).rstrip())


def case_a(fault, out):
    name, throw = "A_allan_diagram", 40.0
    slid = fault.slide(throw)
    fault.decks(name, slid, out, f"hanging wall down {throw:g} m")
    after = process(slid)
    found = after.contacts(fault.faces)
    total = compare(found, analytic_contacts(slid, fault.faces))
    tab = fault.table(found)
    assert np.allclose(tab, fault.length * overlap(throw)), "table != length x layer overlaps"
    assert np.isclose(total, fault.length * (H - throw)), (total, fault.length * (H - throw))
    assert np.allclose(after.volume_grid(slid.cart_dims), fault.volumes), "volumes changed"
    print_allan(tab)
    assert tab[MAIN, SAND].sum() < 1e-6, "the footwall Sand A should face shale only"
    windows = [f"{NAMES[kf]}|{NAMES[kh]}" for kf, kh in zip(*np.nonzero(tab > 1e-6))
               if SAND[kf] and SAND[kh]]
    print(f"{name:22s} throw {throw:g} m: contact {total:,.0f} m2 = {fault.length:.0f} m of "
          f"fault x {H - throw:g} m; sand-sand {sand_sand(tab):,.0f} m2")
    print(f"  windows (FW|HW) {', '.join(windows)}; the footwall Sand A faces shale only")


def cum_vsh(s):
    """VSH times thickness summed from the top of the cake down to stratigraphic depth s."""
    return np.interp(s, Z, np.concatenate([[0.0], np.cumsum(VSH * DZ)]))


def mean_sgr(lo, hi, throw):
    """Mean SGR over footwall-frame depths lo..hi; trapezoids are exact between the kinks."""
    s = np.unique(np.clip(np.concatenate([[lo, hi], Z, Z + throw]), lo, hi))
    g = (cum_vsh(s) - cum_vsh(s - throw)) / throw
    return float(np.sum(0.5 * (g[1:] + g[:-1]) * np.diff(s)) / (hi - lo))


def multiplier(sgr, throw):
    k_fault = 10.0 ** (-4.0 * sgr - 0.25 * np.log10(throw) * (1.0 - sgr) ** 5)
    if k_fault >= K_MATRIX:
        return 1.0
    return 1.0 / (1.0 + 2.0 * throw / 66.0 / CELL_L * (K_MATRIX / k_fault - 1.0))


def face_throws(fault, before, after):
    """Throw of every FAULTS face: hanging-wall minus footwall shift, mean of its two pillars."""
    shift = corner_cube(after)[0] - corner_cube(before)[0]
    corner = lambda c, p: shift[2 * c[1] + p[1] - c[1], 2 * c[0] + p[0] - c[0]]  # noqa: E731
    out = {}
    for f in fault.faces:
        a, b = face_columns(f)
        hw, fw = (a, b) if fault.hw[a] else (b, a)
        out[f] = float(np.mean([corner(hw, p) - corner(fw, p) for p in face_edge(f)]))
    return out


def case_b(fault, out):
    name = "B_sgr_multipliers"
    chain = trace_faults(fault.grid, fault.faces).chains[0]
    ends = (20.0, 50.0) if chain.xy[0, 1] < chain.xy[-1, 1] else (50.0, 20.0)
    spec = Throw(faces=tuple(fault.faces), throw=((0.0, ends[0]), (chain.length, ends[1])),
                 hanging_side=HANGING)
    slid = slide_field(fault.grid, spec, faults=fault.faces)
    assert np.array_equal(slid.coord, fault.grid.coord) and unwelded_faces(slid) == set(fault.faces)
    after = process(slid)
    found = after.contacts(fault.faces)
    compare(found, analytic_contacts(slid, fault.faces))
    assert np.allclose(after.volume_grid(slid.cart_dims), fault.volumes), "volumes changed"
    per_cell = {}                                   # contact area of each first-side cell face
    for (f, k, _), area in found.items():
        per_cell[(f, k)] = per_cell.get((f, k), 0.0) + area

    throws = face_throws(fault, fault.grid, slid)
    mult = {"X": np.ones((NX, NY, NZ)), "Y": np.ones((NX, NY, NZ))}
    rows = []
    for f in fault.faces:
        t, (i, j) = throws[f], face_columns(f)[0]
        on_hw = bool(fault.hw[i, j])
        lift = t if on_hw else 0.0                  # the cell's layers in the footwall frame
        per_face = []
        for k in range(NZ):
            lo, hi = max(Z[k] + lift, t), min(Z[k + 1] + lift, H)
            # The contact interval of the SGR is the one the re-processed grid sees.
            assert np.isclose(per_cell.get((f, k), 0.0), fault.width[f] * max(hi - lo, 0.0))
            if hi - lo > 1e-9:
                sgr = mean_sgr(lo, hi, t)
                m = mult[f.face][i, j, k] = multiplier(sgr, t)
                per_face.append((sgr, m))
                rows.append((i + 1, j + 1, k + 1, f.face, "HW" if on_hw else "FW", t, sgr, m))
        sgr_m = np.array(sorted(per_face))
        assert np.all(np.diff(sgr_m[:, 1]) <= 0.0), f"multiplier must fall as SGR rises on {f}"
    sgr_all, m_all = np.array([r[6] for r in rows]), np.array([r[7] for r in rows])
    assert VSH.min() <= sgr_all.min() and sgr_all.max() <= VSH.max(), "SGR is a VSH average"
    assert np.all((m_all > 0.0) & (m_all <= 1.0)), "multipliers must lie in (0, 1]"

    extra = "".join(format_keyword(f"MULT{d}", mult[d].transpose(2, 1, 0)) for d in "XY")
    path = fault.decks(name, slid, out, f"throw {ends[0]:g}-{ends[1]:g} m along strike",
                       extra=extra)
    (out / f"{name}_faces.csv").write_text("i,j,k,face,cell_side,throw_m,sgr,mult\n" + "".join(
        f"{i},{j},{k},{d},{side},{t:.3f},{sgr:.4f},{m:.6g}\n"
        for i, j, k, d, side, t, sgr, m in rows))
    back = read_grdecl(path).props
    for d in "XY":
        assert np.array_equal(back[f"MULT{d}"].reshape(NZ, NY, NX).transpose(2, 1, 0), mult[d])
    assert faults_from_deck(path) == fault.faults, "FAULTS did not survive the round trip"
    tv = np.array(list(throws.values()))
    print(f"{name:22s} throws {tv.min():.1f}-{tv.max():.1f} m S to N; SGR "
          f"{sgr_all.min():.2f}-{sgr_all.max():.2f} on {len(rows)} cell faces in contact")
    print(f"  multipliers {m_all.min():.1e}-{m_all.max():.1e}; MULTX on {np.sum(mult['X'] < 1)} "
          f"cells and MULTY on {np.sum(mult['Y'] < 1)} read back unchanged")


def case_c(fault, out):
    name = "C_connectivity_sweep"
    throws = np.arange(0.0, 81.0, 5.0)
    rows = []
    for t in throws:
        tab = fault.table(analytic_contacts(fault.slide(t), fault.faces))
        assert np.allclose(tab, fault.length * overlap(t)), t
        assert np.isclose(tab.sum(), fault.length * (H - t)), t
        rows.append((t, sand_sand(tab), tab[MAIN, SAND].sum(), tab[MAIN, MAIN]))
    rows = np.array(rows)
    self_area = rows[:, 3]
    assert np.allclose(self_area, fault.length * np.clip(DZ[MAIN] - throws, 0.0, None))
    assert np.all(np.diff(self_area) <= 0.0), "Sand A against itself can only shrink"
    cut = float(throws[np.argmax(rows[:, 2] < 1e-6)])
    assert cut == DZ[MAIN] and np.all(rows[throws >= cut, 2] < 1e-6), cut
    slid = fault.slide(cut)
    found = process(slid).contacts(fault.faces)
    compare(found, analytic_contacts(slid, fault.faces))
    assert fault.table(found)[MAIN, SAND].sum() < 1e-6, "re-processing still sees a window"
    fault.decks(name, slid, out, f"hanging wall down {cut:g} m: footwall Sand A sealed")
    np.savetxt(out / f"{name}.csv", rows, delimiter=",", fmt="%.6g", comments="",
               header="throw_m,sand_sand_m2,fw_sand_a_to_sand_m2,sand_a_to_sand_a_m2")
    print("  throw m     " + "".join(f"{t:5.0f}" for t in throws))
    print("  sand-sand   " + "".join(f"{a / 1e3:5.1f}" for a in rows[:, 1]))
    print("  FW A - sand " + "".join(f"{a / 1e3:5.1f}" for a in rows[:, 2]))
    low = int(np.argmin(rows[:, 1]))
    print(f"{name:22s} areas 1e3 m2; footwall Sand A's last sand window closes at {cut:g} m, "
          "its thickness")
    print(f"  sand-sand area is lowest at {throws[low]:g} m ({rows[low, 1]:,.0f} m2), then "
          f"deeper sands reopen windows ({rows[-1, 1]:,.0f} m2 at {throws[-1]:g} m)")


def case_d(fault, out):
    name, n = "D_throw_ensemble", 40
    throws = np.random.default_rng(27).lognormal(np.log(30.0), 0.5, n)
    areas, connected = np.empty(n), np.empty(n, bool)
    for r, t in enumerate(throws):
        tab = fault.table(analytic_contacts(fault.slide(t), fault.faces))
        assert np.allclose(tab, fault.length * overlap(t)), t
        areas[r], connected[r] = sand_sand(tab), tab[MAIN, SAND].sum() > 1e-6
    assert np.array_equal(connected, throws < DZ[MAIN]), "Sand A seals once throw > thickness"
    # Exceedance convention: P90 is the value exceeded by 90 % of the realizations.
    p90, p50, p10 = np.percentile(areas, [10.0, 50.0, 90.0])
    assert p90 <= p50 <= p10 and p90 < p10, (p90, p50, p10)
    r = int(np.argmin(np.abs(areas - p50)))
    slid = fault.slide(throws[r])
    found = process(slid).contacts(fault.faces)
    compare(found, analytic_contacts(slid, fault.faces))
    assert np.isclose(sand_sand(fault.table(found)), areas[r])
    fault.decks(name, slid, out, f"realization {r}: hanging wall down {throws[r]:.2f} m")
    np.savetxt(out / f"{name}.csv", np.c_[throws, areas, connected], delimiter=",",
               fmt="%.6g", comments="", header="throw_m,sand_sand_m2,fw_sand_a_connected")
    t90, t50, t10 = np.percentile(throws, [10.0, 50.0, 90.0])
    print(f"{name:22s} {n} lognormal throws, P90/P50/P10 {t90:.1f}/{t50:.1f}/{t10:.1f} m "
          "(Pxx: exceeded by xx %)")
    print(f"  sand-sand P90/P50/P10 {p90:,.0f}/{p50:,.0f}/{p10:,.0f} m2; FW Sand A connected in "
          f"{connected.sum()}/{n} realizations")


def well_tops(grid, column):
    """Horizon depths at the centre of a column: top of every layer, then the base."""
    i, j = column
    z = corner_cube(grid)[:, 2 * j:2 * j + 2, 2 * i:2 * i + 2].mean(axis=(1, 2))
    return np.append(z[0::2], z[-1])


def case_e(fault, out):
    name = "E_well_tops"
    throws = np.linspace(20.0, 44.0, NZ + 1)        # older horizons carry more throw
    slid = fault.slide(throws)
    fault.decks(name, slid, out, "growth fault, throw 20 m at the top to 44 m at the base")
    after = process(slid)
    compare(after.contacts(fault.faces), analytic_contacts(slid, fault.faces))
    fw, hw = (column_of_point(fault.grid, xy) for xy in ((300.0, 400.0), (900.0, 400.0)))
    assert fault.hw[hw] and not fault.hw[fw]
    hw0, hw1 = well_tops(fault.grid, hw), well_tops(slid, hw)
    fw0, fw1 = well_tops(fault.grid, fw), well_tops(slid, fw)
    assert np.allclose(hw1 - hw0, throws, rtol=0.0, atol=1e-9), "HW tops moved by the throw"
    assert np.array_equal(fw1, fw0), "footwall tops must not move"
    vol0, vol1 = fault.volumes[fault.hw].sum(), after.volume_grid(slid.cart_dims)[fault.hw].sum()
    grow = (throws[-1] - throws[0]) / H
    assert np.isclose(vol1 / vol0 - 1.0, grow), "HW thickens by the throw growth"
    print(f"  {'tops, m':16s}{f'FW {fw}':>10s}{f'HW {hw}':>11s}{'HW slid':>10s}{'shift':>8s}")
    for n, label in enumerate([f"top {x}" for x in NAMES] + [f"base {NAMES[-1]}"]):
        print(f"  {label:16s}{fw1[n]:10.1f}{hw0[n]:11.1f}{hw1[n]:10.1f}{hw1[n] - hw0[n]:+8.1f}")
    print(f"{name:22s} HW tops moved by their horizon's throw, FW tops unchanged, HW volume "
          f"+{grow:.0%}")


def case_f(grid, out):
    name = "F_partial_fault_wells"
    g, faces = adapted(grid, np.array([[300.0, 250.0], [900.0, 550.0]]))
    faults = named_chains(g, faces)
    assert len(faults) == 1 and blocks(g, faces).max() == 0, "one fault that closes no block"
    peak, side = 30.0, np.array([300.0, -600.0]) / np.hypot(300.0, 600.0)  # unit normal, SE
    slid = slide_field(g, Throw(faces=tuple(faces), throw=peak, hanging_side=(700.0, 200.0)))
    assert np.array_equal(slid.coord, g.coord) and unwelded_faces(slid) == set(faces)
    after = process(slid)
    compare(after.contacts(faces), analytic_contacts(slid, faces))
    assert np.allclose(after.volume_grid(g.cart_dims), process(g).volume_grid(g.cart_dims))
    write_deck(out / f"{name}_0_input.grdecl", g, faults, header=f"{name}: input")
    write_deck(out / f"{name}_1_slid.grdecl", slid, faults, header=f"{name}: peak {peak:g} m")
    centre = np.array([600.0, 400.0])
    offsets = [50.0, 150.0, 250.0, 400.0]
    shift = []
    for d in offsets + [-150.0]:
        col = column_of_point(g, centre + d * side)
        moved = well_tops(slid, col) - well_tops(g, col)
        assert np.allclose(moved, moved[0], rtol=0.0, atol=1e-9), "a well moves rigidly"
        shift.append(moved[0])
    hw, fw = np.array(shift[:-1]), shift[-1]
    decay = 0.5 * trace_faults(g, faces).chains[0].length
    assert 0.0 < hw[0] < peak and np.all(np.diff(hw) < 0.0), hw
    assert hw[-1] == 0.0 and offsets[-1] > decay + DX and fw == 0.0, (hw, fw)
    print(f"{name:22s} {len(faces)} faces, peak throw {peak:g} m, motion fading to 0 at "
          f"{decay:.0f} m from the fault")
    print(f"  HW wells {'/'.join(f'{d:g}' for d in offsets)} m off the fault: every top down "
          f"{'/'.join(f'{s:.1f}' for s in hw)} m; FW well {fw:.1f} m")


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    grid = base_grid()
    print("A-E: an oblique fault crossing the grid, hanging wall to the east")
    fault = Fault(*adapted(grid, TRACE))
    case_a(fault, out)
    case_b(fault, out)
    case_c(fault, out)
    case_d(fault, out)
    case_e(fault, out)
    print("F: a fault that ends inside the grid")
    case_f(grid, out)
    print("every case: COORD and FAULTS unchanged, only the FAULTS faces open, re-processed "
          "contacts = analytic")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
