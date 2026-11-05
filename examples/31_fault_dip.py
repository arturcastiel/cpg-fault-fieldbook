"""Give a fault a dip angle by tilting pillars (COORD only), and check every claim by hand.

A fault in a corner-point grid is a ruled surface between its pillars, so its dip *is* the
lean of those pillars. dip_grid rotates each fault pillar, in the vertical plane normal to the
fault's strike, until it dips at the requested angle toward the hanging wall, and moves the
neighbouring pillars smoothly so no cell folds. ZCORN, ACTNUM and FAULTS never change.

Every case re-derives its numbers from the grid itself rather than trusting the result object:
angles from the six COORD numbers of each pillar, volumes from a re-processed topology,
untouched keywords by byte comparison.

- A rigid 60 deg: a rim-to-rim fault; the hanging wall shears as one block, so its cell
  volumes do not change; the footwall absorbs the change over a cosine fade. With
  decay=<distance> the tilt fades out on both sides instead and the grid rim stays vertical;
- B why neighbours move: 40 deg by tilting only the fault pillars folds cells; dip_grid
  moves the neighbours too and every cell stays valid;
- C pivots: the depth where pillars keep their plan position (top, middle, bottom, a depth);
  only the middle pivot keeps every column's total volume;
- D a fault with a tip: the tilt fades out around the tip on both sides; volumes change there;
- E an angle table along strike: 55 deg at the south rim to 80 deg at the north rim;
- F a stair-stepped diagonal fault: the dip is set along its smoothed strike;
- G a graben: two facing faults, 60 and 70 deg, both exact, rigid and with decay;
- H the gate: an impossible footwall fade is widened (warning); with widening off it raises,
  or bisects and then does NOT reach the angle;
- I on a deck: dip_deck dips and slides 12 m in one call; the written file is read back and the
  angle, the throw, the unchanged FAULTS and the hanging-wall volumes are checked from disk;
- J a field deck (Norne, when it is on this machine; env CPG_NORNE_GRDECL as in lesson 30):
  E_01 to 60 deg passes. BC to 70 deg FAILS the gate with the default pivot: its present dip
  jumps from 59 to 84 deg within four neighbouring pillars, and one top-layer cell folds whatever
  the fade. pivot="top" passes; so does a lower angle.

With matplotlib installed and --out DIR, cases A, B, D and G also write PNG figures.
"""

import argparse
import importlib.util
import os
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import CornerPointGrid, check_validity, pillar_lattice, read_grdecl
from fault_projections import Dip, DipError, Throw, dip_deck, dip_grid
from fault_projections.adapt import fault_faces_for
from fault_projections.faults import FaultFace

from fieldbook.zslide import layer_cake, process, realized_throw, write_deck

NX, NY, DX = 16, 8, 50.0
DZ = [15.0, 25.0, 20.0, 30.0]
TOP = 2000.0
H = sum(DZ)                                  # 90 m: the pillars run from 2000 to 2090 m
EAST, WEST = (700.0, 200.0), (100.0, 200.0)  # plan points on either side of x = 400
PNG = importlib.util.find_spec("matplotlib") is not None
NORNE = "/mnt/c/projects/opm-machine/opm-tests/norne/INCLUDE/GRID/IRAP_1005.GRDECL"


def cake():
    return layer_cake(NX, NY, len(DZ), DX, DZ, top=lambda x, y: TOP)


def x_fault(i, j0, j1):
    """FAULTS faces on pillar line i (x = i * DX), rows j0..j1 (1-based, Eclipse style)."""
    return [FaultFace(i, j, "X") for j in range(j0, j1 + 1)]


def coord_of(grid, i, j):
    nx, ny, _ = grid.cart_dims
    return np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)[j, i]


def coord_dip(grid, i, j, normal):
    """Dip of pillar (i, j) in degrees, straight from its six COORD numbers: the drop over
    the horizontal run measured along ``normal`` (the dip direction)."""
    x0, y0, z0, x1, y1, z1 = coord_of(grid, i, j)
    run = (x1 - x0) * normal[0] + (y1 - y0) * normal[1]
    return float(np.degrees(np.arctan2(z1 - z0, run)))


def volumes(grid):
    return process(grid).volume_grid(grid.cart_dims)


def rel_change(before, after):
    return volumes(after) / volumes(before) - 1.0


def cells_valid(before, after):
    return check_validity(pillar_lattice(after), pillar_lattice(before))


def section(ax, grid, row, title, bad=None, faults=(), xlim=None):
    """Draw the south faces of cell row ``row`` in the x-z plane: the real cell faces."""
    from matplotlib.collections import PolyCollection

    nx, _, nz = grid.cart_dims
    coord = np.asarray(grid.coord, float).reshape(-1, nx + 1, 6)[row]
    z = np.asarray(grid.zcorn, float).reshape(2 * nz, -1, 2 * nx)[:, 2 * row]

    def x_at(i, depth):
        x0, _, z0, x1, _, z1 = coord[i]
        return x0 + (depth - z0) / (z1 - z0) * (x1 - x0)

    polys, colours = [], []
    for k in range(nz):
        for i in range(nx):
            zt0, zt1, zb0, zb1 = z[2 * k, 2 * i], z[2 * k, 2 * i + 1], \
                z[2 * k + 1, 2 * i], z[2 * k + 1, 2 * i + 1]
            polys.append([(x_at(i, zt0), zt0), (x_at(i + 1, zt1), zt1),
                          (x_at(i + 1, zb1), zb1), (x_at(i, zb0), zb0)])
            hit = bad is not None and bad[i, row]
            colours.append((0.89, 0.29, 0.2, 0.45) if hit else
                           ("#c6dbef", "#a1d99b", "#fdd0a2", "#dadaeb")[k])
    ax.add_collection(PolyCollection(polys, facecolors=colours, edgecolors="#333333",
                                     linewidths=0.6))
    for i in faults:
        x0, _, z0, x1, _, z1 = coord[i]
        ax.plot([x0, x1], [z0, z1], color="#b2182b", lw=2.2)
    ax.set_xlim(xlim or (coord[:, [0, 3]].min() - 10, coord[:, [0, 3]].max() + 10))
    ax.set_ylim(TOP + H + 5, TOP - 5)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x (m)", fontsize=8, labelpad=1)
    ax.set_ylabel("depth (m)", fontsize=8)
    ax.tick_params(labelsize=7)


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    import matplotlib.pyplot as plt

    plt.close(fig)
    print(f"{'':22s} wrote {path.name}")


# ── A rigid 60 deg ───────────────────────────────────────────────────────────────────────────


def case_a(out):
    name = "A_rigid_60"
    grid = cake()
    faults = {"F": x_fault(8, 1, NY)}                       # x = 400 m, south rim to north rim
    res = dip_grid(grid, Dip("F", 60.0, hanging_side=EAST), faults=faults)
    assert res.ok and res.modes == {"F": "rigid"} and res.scale == 1.0

    got = [coord_dip(res.grid, 8, j, (1.0, 0.0)) for j in range(NY + 1)]
    assert np.allclose(got, 60.0, atol=1e-9), got
    assert np.array_equal(res.grid.zcorn, grid.zcorn), "ZCORN must not change"
    assert np.array_equal(res.grid.actnum, grid.actnum)
    v0, v1 = volumes(grid), volumes(res.grid)
    dv, col = v1 / v0 - 1.0, v1.sum(2) / v0.sum(2) - 1.0
    assert np.abs(dv[8:]).max() < 1e-12, "a rigid hanging wall shears parallel to the fault"
    assert np.abs(col).max() < 1e-12, "a middle pivot keeps every column's total volume"
    lean = (pillar_lattice(res.grid).bot - pillar_lattice(res.grid).top)[:, 0, 0]
    assert np.allclose(lean[8:], H / np.tan(np.radians(60.0)))
    assert cells_valid(grid, res.grid).ok

    rim = dip_grid(grid, Dip("F", 60.0, hanging_side=EAST, decay=300.0), faults=faults)
    rim_lean = (pillar_lattice(rim.grid).bot - pillar_lattice(rim.grid).top)[:, 0, 0]
    rim_dv = rel_change(grid, rim.grid)
    assert rim.ok and rim.modes == {"F": "decay"} and rim_lean[-1] == 0.0
    assert np.allclose([coord_dip(rim.grid, 8, j, (1.0, 0.0)) for j in range(NY + 1)], 60.0)

    x0, y0, z0, x1, y1, z1 = coord_of(res.grid, 8, 4)
    print(f"{name:22s} F dips {min(got):.9f}..{max(got):.9f} deg at all {len(got)} pillars")
    print(f"{'':22s} pillar (8,4) COORD {x0:.3f} {y0:.0f} {z0:.0f}  {x1:.3f} {y1:.0f} "
          f"{z1:.0f}: atan({z1 - z0:.0f} / {x1 - x0:.3f}) = "
          f"{np.degrees(np.arctan2(z1 - z0, x1 - x0)):.6f} deg")
    print(f"{'':22s} ZCORN and ACTNUM identical; pillar lean (m) on lines 5..16: "
          + " ".join(f"{v:.1f}" for v in lean[5:]))
    print(f"{'':22s} cell volumes: hanging wall {np.abs(dv[8:]).max():.1e}, footwall up to "
          f"{100 * np.abs(dv[:8]).max():.2f} % (column 7, top to base: "
          + ", ".join(f"{100 * v:+.2f}" for v in dv[7, 0]) + " %); every column total "
          f"{np.abs(col).max():.0e}")
    print(f"{'':22s} rigid mode leans the whole hanging wall, east rim included. "
          f"decay=300 instead: rim lean {rim_lean[-1]:.1f} m, lean on lines 3..13 "
          + " ".join(f"{v:.1f}" for v in rim_lean[3:14])
          + f", cells change up to {100 * np.nanmax(np.abs(rim_dv)):.2f} % on both sides")
    if out and PNG:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 1, figsize=(9, 5.0))
        section(axes[0], grid, 0, "before: vertical pillars, F at x = 400 m", faults=[8])
        section(axes[1], res.grid, 0, "rigid (default): F at 60 deg, the hanging wall leans as "
                f"one block, footwall fade {res.fade['F'][1]:.0f} m", faults=[8])
        section(axes[2], rim.grid, 0, "decay=300: F at 60 deg, the tilt fades out on both "
                "sides, the grid rim stays vertical", faults=[8])
        save(fig, out / "31_A_rigid_60.png")
    return res


# ── B why neighbours move ────────────────────────────────────────────────────────────────────


def tilt_only_fault_pillars(grid, i, angle):
    """The naive edit: lean the fault's own pillars and nothing else (middle pivot)."""
    lean = H / np.tan(np.radians(angle))
    coord = np.asarray(grid.coord, float).reshape(NY + 1, NX + 1, 6).copy()
    coord[:, i, 0] -= lean / 2
    coord[:, i, 3] += lean / 2
    return CornerPointGrid(grid.cart_dims, coord.reshape(-1, 6), grid.zcorn, grid.actnum)


def case_b(out):
    name = "B_why_neighbours_move"
    grid = cake()
    faults = {"F": x_fault(8, 1, NY)}
    naive = tilt_only_fault_pillars(grid, 8, 40.0)
    bad = cells_valid(grid, naive)
    assert not bad.ok and bad.n_folded > 0
    assert np.allclose([coord_dip(naive, 8, j, (1.0, 0.0)) for j in range(NY + 1)], 40.0)

    res = dip_grid(grid, Dip("F", 40.0, hanging_side=EAST), faults=faults)
    good = cells_valid(grid, res.grid)
    assert res.ok and good.ok
    assert np.allclose([coord_dip(res.grid, 8, j, (1.0, 0.0)) for j in range(NY + 1)], 40.0)
    shift = np.abs(res.d_bot[..., 0]).max(axis=1)
    print(f"{name:22s} 40 deg needs a lean of {H / np.tan(np.radians(40.0)):.1f} m over "
          f"{H:.0f} m, more than one cell ({DX:.0f} m)")
    print(f"{'':22s} tilt only the fault pillars: {bad.summary()} "
          f"({bad.n_folded} of {NX * NY} columns)")
    print(f"{'':22s} dip_grid: {res.check.summary()}; footwall fade "
          f"{res.fade['F'][1]:.0f} m, base shift (m) on lines 4..9: "
          + " ".join(f"{s:.1f}" for s in shift[4:10]) + " (lines 8..16 move as one block)")
    if out and PNG:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(9, 3.5))
        folded = bad.folded_top | bad.folded_bot
        section(axes[0], naive, 0, f"tilt only the fault pillars: {bad.n_folded} folded "
                "columns (red)", bad=folded, faults=[8], xlim=(150, 850))
        section(axes[1], res.grid, 0, "dip_grid: neighbours move smoothly, every cell valid",
                faults=[8], xlim=(150, 850))
        save(fig, out / "31_B_why_neighbours_move.png")


# ── C pivots ─────────────────────────────────────────────────────────────────────────────────


def case_c():
    name = "C_pivots"
    grid = cake()
    faults = {"F": x_fault(8, 1, NY)}
    lean = H / np.tan(np.radians(70.0))
    v0 = volumes(grid)
    rows = []
    for pivot, top, base in (("top", 0.0, 1.0), ("middle", -0.5, 0.5), ("bottom", -1.0, 0.0),
                             (2030.0, -30.0 / H, 60.0 / H)):
        res = dip_grid(grid, Dip("F", 70.0, hanging_side=EAST, pivot=pivot), faults=faults)
        assert res.ok
        assert np.isclose(coord_dip(res.grid, 8, 4, (1.0, 0.0)), 70.0)
        dtop, dbot = res.d_top[8, 4, 0], res.d_bot[8, 4, 0]
        assert np.isclose(dtop, top * lean) and np.isclose(dbot, base * lean)
        stays = TOP - dtop / (dbot - dtop) * H
        col = volumes(res.grid).sum(2) / v0.sum(2) - 1.0
        rows.append(f"{str(pivot):>6s}: top {dtop:+7.2f} m, base {dbot:+7.2f} m, fixed at "
                    f"{stays:.0f} m; footwall column 7 volume {100 * col[7, 0]:+.2f} %")
    print(f"{name:22s} 70 deg (lean {lean:.2f} m); pillar (8,4) moves east by:")
    for row in rows:
        print(f"{'':22s}   {row}")


# ── D a fault with a tip ─────────────────────────────────────────────────────────────────────


def case_d(out):
    name = "D_tip_55"
    grid = cake()
    faults = {"T": x_fault(8, 1, 5)}            # rows 1..5 of 8: the tip is pillar (8, 5)
    res = dip_grid(grid, Dip("T", 55.0, hanging_side=EAST), faults=faults)
    assert res.ok and res.modes == {"T": "decay"}
    got = [coord_dip(res.grid, 8, j, (1.0, 0.0)) for j in range(6)]
    assert np.allclose(got, 55.0, atol=1e-9)
    lean = np.linalg.norm(pillar_lattice(res.grid).bot - pillar_lattice(res.grid).top, axis=-1)
    north = lean[8, 5:]                          # from the tip to the north rim
    east = lean[8:, 2]                           # from the fault into the hanging wall
    assert north[0] > north[-1] and east[0] > east[-1]
    dv = rel_change(grid, res.grid)
    assert cells_valid(grid, res.grid).ok
    print(f"{name:22s} T (rows 1..5) dips {min(got):.6f}..{max(got):.6f} deg, mode decay, "
          f"fade {res.fade['T'][0]:.0f} m")
    print(f"{'':22s} pillar lean north of the tip (j=5..8): "
          + ", ".join(f"{v:.1f}" for v in north) + " m")
    print(f"{'':22s} pillar lean east across the hanging wall (i=8..16): "
          + ", ".join(f"{v:.1f}" for v in east) + " m")
    print(f"{'':22s} cells valid; volumes change by up to {100 * np.nanmax(np.abs(dv)):.1f} % "
          "(a fading tilt stretches cells; only rigid mode keeps them)")
    if out and PNG:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        top = pillar_lattice(res.grid).top
        base = pillar_lattice(res.grid).bot
        shift = np.linalg.norm(res.d_bot - res.d_top, axis=-1)
        pcm = ax.pcolormesh(top[..., 0], top[..., 1], shift, shading="gouraud", cmap="viridis")
        for j in range(NY + 1):
            ax.plot(base[:, j, 0], base[:, j, 1], color="w", lw=0.4)
        for i in range(NX + 1):
            ax.plot(base[i, :, 0], base[i, :, 1], color="w", lw=0.4)
        ax.plot([400, 400], [0, 250], color="#b2182b", lw=2.5, label="T (FAULTS rows 1..5)")
        ax.set_aspect("equal")
        ax.set_title("D: tilt change per pillar (top to base, m); white = base of the cells",
                     fontsize=9)
        ax.legend(loc="upper left", fontsize=8)
        fig.colorbar(pcm, ax=ax, shrink=0.8)
        save(fig, out / "31_D_tip_55.png")


# ── E angle table ────────────────────────────────────────────────────────────────────────────


def case_e():
    name = "E_angle_table"
    grid = cake()
    faults = {"F": x_fault(8, 1, NY)}
    table = np.array([[0.0, 55.0], [NY * DX, 80.0]])            # (arc length, degrees)
    res = dip_grid(grid, Dip("F", table, hanging_side=EAST), faults=faults)
    assert res.ok
    arc = res.realized_dip("F")[:, 0]
    want = np.interp(arc, table[:, 0], table[:, 1])
    j = np.rint(np.asarray(res.chains["F"].ij)[:, 1]).astype(int)
    got = np.array([coord_dip(res.grid, 8, jj, (1.0, 0.0)) for jj in j])
    assert np.allclose(got, want, atol=1e-9)
    print(f"{name:22s} requested -> from COORD, along the fault (arc m: deg):")
    print(f"{'':22s}   " + "  ".join(f"{a:.0f}: {w:.3f}->{g:.3f}" for a, w, g in
                                     zip(arc, want, got)))


# ── F stair-stepped diagonal ─────────────────────────────────────────────────────────────────


def case_f():
    name = "F_staircase_65"
    grid = cake()
    trace = np.array([[200.0, 0.0], [600.0, 400.0]])           # 45 deg in plan
    faces = fault_faces_for(grid, trace)
    res = dip_grid(grid, Dip("D", 65.0, hanging_side=(600.0, 100.0)), faults={"D": faces})
    assert res.ok
    ij = np.asarray(res.chains["D"].ij)
    normals = res.normals["D"]
    got = np.array([coord_dip(res.grid, int(i), int(j), n) for (i, j), n in zip(ij, normals)])
    assert np.allclose(got, 65.0, atol=1e-9)
    ideal = np.array([np.sqrt(0.5), -np.sqrt(0.5)])           # normal of the straight trace
    along = np.array([coord_dip(res.grid, int(i), int(j), ideal) for i, j in ij])
    mid = normals[len(normals) // 2]
    print(f"{name:22s} {len(faces)} stair faces ({sum(f.face == 'X' for f in faces)} X, "
          f"{sum(f.face == 'Y' for f in faces)} Y), {len(ij)} chain pillars")
    print(f"{'':22s} along the smoothed strike: {got.min():.6f}..{got.max():.6f} deg "
          f"(mid-chain dip direction {mid[0]:+.3f}, {mid[1]:+.3f})")
    print(f"{'':22s} along the straight-trace normal (+0.707, -0.707): "
          f"{along.min():.2f}..{along.max():.2f} deg (the ends differ: the strike there is "
          "averaged over fewer steps)")


# ── G graben ─────────────────────────────────────────────────────────────────────────────────


def case_g(out):
    name = "G_graben_60_70"
    grid = cake()
    faults = {"A": x_fault(5, 1, NY), "B": x_fault(11, 1, NY)}  # x = 250 and x = 550
    centre = (400.0, 200.0)                                     # the graben: hanging wall of both

    def dips(**kw):
        return [Dip("A", 60.0, hanging_side=centre, **kw), Dip("B", 70.0, hanging_side=centre, **kw)]

    v0 = volumes(grid)
    rows = []
    for label, res in (("rigid", dip_grid(grid, dips(), faults=faults)),
                       ("decay=150", dip_grid(grid, dips(decay=150.0), faults=faults))):
        assert res.ok and cells_valid(grid, res.grid).ok
        a = [coord_dip(res.grid, 5, j, (1.0, 0.0)) for j in range(NY + 1)]
        b = [coord_dip(res.grid, 11, j, (-1.0, 0.0)) for j in range(NY + 1)]
        assert np.allclose(a, 60.0, atol=1e-9) and np.allclose(b, 70.0, atol=1e-9)
        v1 = volumes(res.grid)
        dv, col = v1 / v0 - 1.0, v1.sum(2) / v0.sum(2) - 1.0
        assert np.abs(col).max() < 1e-12
        lean = (pillar_lattice(res.grid).bot - pillar_lattice(res.grid).top)[:, 0, 0]
        rows.append((label, res, a, b, dv, lean))
    assert rows[0][5][0] != 0.0 and rows[1][5][0] == 0.0 == rows[1][5][-1]
    for label, res, a, b, dv, lean in rows:
        print(f"{name if label == 'rigid' else '':22s} {label:9s}: A (x=250, dips east) "
              f"{min(a):.6f}..{max(a):.6f} deg, B (x=550, dips west) {min(b):.6f}..{max(b):.6f}"
              f" deg; cells up to {100 * np.nanmax(np.abs(dv)):.2f} %, rim leans "
              f"{lean[0]:+.1f} / {lean[-1]:+.1f} m")
    print(f"{'':22s} rigid: each outer block is on the hanging side of the far fault, so it "
          "leans with it; decay keeps the rims vertical")
    if out and PNG:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(9, 3.5))
        section(axes[0], rows[0][1].grid, 0, "G rigid: A 60 deg and B 70 deg dip toward each "
                "other; the outer blocks lean with the far fault", faults=[5, 11])
        section(axes[1], rows[1][1].grid, 0, "G decay=150: same angles, the tilt fades out "
                "and the rims stay vertical", faults=[5, 11])
        save(fig, out / "31_G_graben.png")


# ── H the gate ───────────────────────────────────────────────────────────────────────────────


def case_h():
    name = "H_gate"
    grid = cake()
    faults = {"F": x_fault(8, 1, NY)}
    spec = Dip("F", 35.0, hanging_side=EAST, footwall=10.0)     # a 10 m fade cannot take it
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        wide = dip_grid(grid, spec, faults=faults)
    assert wide.ok and wide.widen > 1.0
    assert any("widened" in str(w.message) for w in caught)
    assert np.allclose([coord_dip(wide.grid, 8, j, (1.0, 0.0)) for j in range(NY + 1)], 35.0)
    try:
        dip_grid(grid, spec, faults=faults, widen=0)
        raise AssertionError("the gate should have refused")
    except DipError as exc:
        refused = str(exc).splitlines()[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        half = dip_grid(grid, spec, faults=faults, widen=0, on_invalid="bisect")
    reached = [coord_dip(half.grid, 8, j, (1.0, 0.0)) for j in range(NY + 1)]
    assert half.ok and half.scale < 1.0 and min(reached) > 35.0
    print(f"{name:22s} 35 deg with footwall=10 m: widened x{wide.widen:g} to "
          f"{wide.fade['F'][1]:.0f} m (warning; the first widening jumps to at least the "
          "auto distance), angle exact")
    print(f"{'':22s} widen=0: DipError: {refused[:90]}...")
    print(f"{'':22s} widen=0, on_invalid='bisect': {100 * half.scale:.1f} % of the change, "
          f"F reaches only {min(reached):.2f} deg, not 35")


# ── I on a deck ──────────────────────────────────────────────────────────────────────────────


def case_i(out):
    name = "I_deck_dip_and_slide"
    work = out / name
    work.mkdir(exist_ok=True)
    grid = cake()
    faults = {"F": x_fault(8, 1, NY)}
    deck = write_deck(work / "case.grdecl", grid, faults, header=f"{name}: vertical F")
    before = deck.read_bytes()
    res = dip_deck(deck, Dip("F", 65.0, hanging_side=EAST), work,
                   slide=Throw("F", 12.0, hanging_side=EAST))
    assert res.ok and res.slide is not None and res.slide.ok and not res.slide.unexpected
    assert deck.read_bytes() == before, "the input deck is never written"

    back = read_grdecl(res.grid_path)                         # from disk, not from memory
    got = [coord_dip(back, 8, j, (1.0, 0.0)) for j in range(NY + 1)]
    assert np.allclose(got, 65.0, atol=1e-6)
    throw = realized_throw(res.grid, back, faults["F"], EAST)[:, 1]
    assert np.allclose(throw, 12.0)
    assert "FAULTS" not in res.grid_path.read_text()
    dv = rel_change(grid, back)
    assert np.abs(dv[8:]).max() < 1e-9, "dip (rigid) then throw keeps hanging-wall volumes"
    print(f"{name:22s} dip_deck(Dip('F', 65), slide=Throw('F', 12)) wrote "
          f"{res.grid_path.name}, {res.report_path.name}")
    print(f"{'':22s} read back: F dips {min(got):.6f}..{max(got):.6f} deg, throw "
          f"{throw.min():.6f}..{throw.max():.6f} m; input deck unchanged; no FAULTS in the "
          "output (the deck's FAULTS stay valid)")
    print(f"{'':22s} hanging-wall volume change {np.abs(dv[8:]).max():.1e}; footwall up to "
          f"{100 * np.abs(dv[:8]).max():.2f} %")


# ── J a field deck ───────────────────────────────────────────────────────────────────────────


def norne_deck(out):
    grid = Path(os.environ.get("CPG_NORNE_GRDECL", NORNE))
    files = [grid, grid.parent / "ACTNUM_0704.prop",
             grid.parent.parent / "FAULT" / "FAULT_JUN_05.INC"]
    if not all(f.is_file() for f in files):
        return None
    deck = out / "NORNE_GRID.GRDECL"
    deck.write_text("".join(f"INCLUDE\n  '{f.as_posix()}' /\n\n" for f in files))
    return deck


def case_j(out):
    name = "J_field_deck_norne"
    work = out / name
    work.mkdir(exist_ok=True)
    deck = norne_deck(work)
    if deck is None:
        print(f"{name:22s} skipped: Norne is not on this machine (set CPG_NORNE_GRDECL)")
        return
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        e01 = dip_deck(deck, Dip("E_01", 60.0), work, prefix="E_01_60")
        back = read_grdecl(e01.grid_path)
        got = e01.realized_dip("E_01", back)[:, 1]
        assert e01.ok and np.allclose(got, 60.0, atol=1e-6)
        assert np.array_equal(back.zcorn, e01.before.zcorn)
        was = e01.angle_before["E_01"]
        nx, ny, nz = e01.grid.cart_dims
        dv = np.abs(e01.volume_change[np.isfinite(e01.volume_change)])
        print(f"{name:22s} {nx}x{ny}x{nz}, {int(e01.grid.actnum.sum()):,} active cells")
        print(f"{'':22s} E_01 -> 60 deg: was {was.min():.1f}..{was.max():.1f} deg on "
              f"{len(was)} pillars, now {got.min():.6f}..{got.max():.6f} (read back from "
              f"{e01.grid_path.name}); {e01.modes['E_01']} mode, fade "
              f"{e01.fade['E_01'][0]:.0f} m; {e01.check.summary()}; cell volumes up to "
              f"{100 * dv.max():.1f} %")
        try:
            dip_deck(deck, Dip("BC", 70.0), work, prefix="BC_70")
            raise AssertionError("BC at 70 deg with the middle pivot was expected to fail")
        except DipError as exc:
            print(f"{'':22s} BC -> 70 deg, pivot middle: DipError: {str(exc).split(';')[0]}")
        top = dip_deck(deck, Dip("BC", 70.0, pivot="top"), work, prefix="BC_70_top")
        got = top.realized_dip("BC")[:, 1]
        was = top.angle_before["BC"]
        assert top.ok and np.allclose(got, 70.0, atol=1e-6)
        print(f"{'':22s} BC -> 70 deg, pivot top: was {was.min():.1f}..{was.max():.1f} deg, "
              f"now {got.min():.6f}..{got.max():.6f}; fades widened x{top.widen:g}; "
              f"{top.check.summary()}")


def main(out):
    case_a(out)
    case_b(out)
    case_c()
    case_d(out)
    case_e()
    case_f()
    case_g(out)
    case_h()
    case_i(out)
    case_j(out)
    print("every case: angles read from COORD, ZCORN untouched by the dip, volumes re-processed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, help="keep the decks and figures here")
    args = parser.parse_args()
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        main(args.out)
    else:
        with TemporaryDirectory() as temporary:
            main(Path(temporary))
