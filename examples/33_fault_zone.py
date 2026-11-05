"""Find the cells within -D/+D of a fault, and (if asked) keep a dip or a slide inside them.

fault_zone measures, for every cell, the distance from the cell's centre to the fault surface
(the ruled surface between the fault's pillars, extended past the top and base of the grid)
along the surface's normal. The distance is signed: + on the hanging wall, - on the footwall.
A cell is in the zone when it lies within D of the fault, or within (D_footwall,
D_hanging_wall) for different distances on each side.

By default dip_grid and slide_grid deform the whole grid; the grid's outer sides hold still
(dip_grid(boundary="fixed")). Passing zone= keeps the edit inside the zone: every pillar and
corner of a cell outside it stays exactly where it was.

The lesson does not trust the library's numbers. Each selection is re-measured: from the
equation of the fault plane where the fault is planar, and by brute force against a dense
sampling of the pillar surface where it is curved. Each edit is judged by lesson 32's
independent cell checker, and "nothing outside moved" is a comparison of raw COORD and ZCORN.

- A vertical fault, D = 60 m: the band is every cell whose centre lies within 60 m of x = 300;
- B dipping fault (60 deg): the band slants with the fault; its distances are checked against
  the plane fitted through the fault's own COORD, and + is read from the dip;
- C asymmetric zone: 30 m on the footwall, 90 m on the hanging wall;
- D curved fault (an arc adapted into the pillars, dipped 65 deg): distances checked by brute
  force against 150,000 points of the pillar surface;
- E two crossing faults, D = 60 m: each cell goes to its nearest fault (a tie to the first);
- F GRDECL export: FLTZONE (0 outside, 1 footwall band, 2 hanging-wall band) written, read
  back and checked cell by cell in deck order;
- G a confined dip: 60 deg inside a 100 m zone, against the default whole-domain dip (rim
  fixed) and boundary="free"; not one corner outside the zone moves;
- H a confined slide: 20 m inside 75 m, against the default slide of the whole hanging wall;
- I dip and slide confined together: dip_grid(..., slide=..., zone=...);
- J refusals: a zone narrower than a cell, and a zone too tight for the change, are refused;
- K a field deck (Norne, when it is on this machine; env CPG_NORNE_GRDECL as in lesson 30):
  the zone of every fault at D = 100 m, and E_01 dipped 60 deg inside 400 m of it.

With matplotlib installed and --out DIR, every case writes one PNG: the plan and a section
BEFORE and AFTER, the footwall band in blue, the hanging-wall band in orange and, after an
edit, every cell that moved outlined in black. --case B,G runs only those cases.
"""

import argparse
import importlib.util
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import read_grdecl
from fault_projections import (Dip, DipError, Throw, ZoneError, dip_grid, fault_zone,
                               slide_grid, write_zone)
from fault_projections.dip import _deck_faults
from scipy.spatial import cKDTree

from fieldbook.zslide import realized_throw

# lesson 32's independent cell checker: raw COORD + ZCORN -> hexahedra, no library code
_spec = importlib.util.spec_from_file_location(
    "fault_dip_stress", Path(__file__).with_name("32_fault_dip_stress.py"))
l32 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(l32)
corners, audit, cake, x_fault, y_fault = l32.corners, l32.audit, l32.cake, l32.x_fault, l32.y_fault

NX, NY, DX = l32.NX, l32.NY, l32.DX        # 24 x 12 columns of 25 m, 5 layers, 100 m thick
EAST = (575.0, 150.0)
PNG = importlib.util.find_spec("matplotlib") is not None
OUT = None                       # set by main: where figures go (None: no figures)


# ── independent measurements ────────────────────────────────────────────────────────────────


def centres(grid):
    """(NX, NY, NZ, 3) cell centres: the mean of the eight corners from raw COORD + ZCORN."""
    return corners(grid).reshape(grid.cart_dims + (8, 3)).mean(-2)


def moved(before, after):
    """(NX, NY, NZ) cells with any corner not bit-identical to before."""
    c0, c1 = corners(before), corners(after)
    same = (c0 == c1) | (np.isnan(c0) & np.isnan(c1))
    return ~same.all(axis=(3, 4, 5, 6))


def band(d, lo, hi):
    """-1 / +1 / 0 from signed distances: the rule a zone applies."""
    out = np.zeros(d.shape, np.int8)
    out[(d < 0.0) & (d >= -lo)] = -1
    out[(d >= 0.0) & (d <= hi)] = 1
    return out


def fault_plane(grid, faces):
    """Least-squares plane through the top and base points of the fault's pillars (COORD):
    (point, unit normal, worst residual)."""
    nx, ny, _ = grid.cart_dims
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)
    ij = {p for f in faces for p in l32.face_pillars(f)}
    pts = np.array([coord[j, i, s:s + 3] for i, j in ij for s in (0, 3)])
    mid = pts.mean(0)
    normal = np.linalg.svd(pts - mid)[2][-1]
    return mid, normal, float(np.abs((pts - mid) @ normal).max())


def surface_samples(grid, faces, pad, n_depth=200, n_along=40):
    """Points of the ruled pillar surface: every face's two pillars, extended ``pad`` past the
    top and base of the grid, sampled at the same depths and blended linearly across."""
    nx, ny, _ = grid.cart_dims
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)
    z = corners(grid)[..., 2]
    depths = np.linspace(np.nanmin(z) - pad, np.nanmax(z) + pad, n_depth)
    s = np.linspace(0.0, 1.0, n_along)
    out = []
    for f in faces:
        line = []
        for i, j in l32.face_pillars(f):
            x0, y0, z0, x1, y1, z1 = coord[j, i]
            t = (depths - z0) / (z1 - z0)
            line.append(np.c_[x0 + t * (x1 - x0), y0 + t * (y1 - y0), depths])
        p, q = line
        out.append(((1 - s)[None, :, None] * p[:, None] + s[None, :, None] * q[:, None])
                   .reshape(-1, 3))
    return np.concatenate(out)


# ── figures ─────────────────────────────────────────────────────────────────────────────────

FOOT, HANG, OUTSIDE = "#4393c3", "#f4a582", "#ececec"
LAYERS = ("#e0e0e0", "#cfcfcf", "#e8e8e8", "#d6d6d6", "#dddddd")


def colours(side, nz_k):
    if side is None:
        return LAYERS[nz_k % len(LAYERS)]
    return FOOT if side < 0 else HANG if side > 0 else OUTSIDE


def zone_plan(ax, grid, side, moved_, faults, title, k=0):
    """The top of layer ``k`` in plan: each cell coloured by its zone side, moved cells
    outlined black, the faults at the top (solid) and the base (dashed) of their pillars."""
    from matplotlib.collections import PolyCollection

    nx, ny, nz = grid.cart_dims
    c = corners(grid)
    polys = [c[i, j, k, [0, 1, 1, 0], [0, 0, 1, 1], 0, :2] for i in range(nx) for j in range(ny)]
    fill = [colours(None if side is None else side[i, j, k], k)
            for i in range(nx) for j in range(ny)]
    ax.add_collection(PolyCollection(polys, facecolors=fill, edgecolors="#888888",
                                     linewidths=0.25))
    many = len(faults) > 4                                   # a field: thin lines
    if moved_ is not None and moved_[..., k].any():
        ax.add_collection(PolyCollection([polys[i * ny + j] for i, j in
                                          zip(*np.nonzero(moved_[..., k]))],
                                         facecolors="none", edgecolors="#111111",
                                         linewidths=0.4 if many else 0.9))
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)
    for faces in faults:
        for f in faces:
            p, q = l32.face_pillars(f)
            for off, style in ((0, "-"), (3, "--")):
                ax.plot([coord[p[1], p[0], off], coord[q[1], q[0], off]],
                        [coord[p[1], p[0], off + 1], coord[q[1], q[0], off + 1]],
                        style, color="black",
                        lw=(0.6 if off == 0 else 0.3) if many else (1.6 if off == 0 else 0.9))
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)


def zone_section(ax, grid, cut, side, moved_, faults, title, window=None):
    """True-scale section: ("row", j) shows the south faces of row j in x-z, ("col", i) the
    west faces of column i in y-z; cells coloured by zone side, moved cells outlined black,
    fault pillars in black."""
    from matplotlib.collections import PolyCollection

    nx, ny, nz = grid.cart_dims
    kind, n = cut
    c = corners(grid)
    act = np.asarray(grid.actnum).reshape(nz, ny, nx).transpose(2, 1, 0).astype(bool)
    axis = 0 if kind == "row" else 1
    cells = [(i, n) for i in range(nx)] if kind == "row" else [(n, j) for j in range(ny)]
    polys, fill, bold = [], [], []
    for i, j in cells:
        for k in range(nz):
            q = c[i, j, k, :, 0, :, :] if kind == "row" else c[i, j, k, 0, :, :, :]
            poly = [(p[axis], p[2]) for p in (q[0, 0], q[1, 0], q[1, 1], q[0, 1])]
            polys.append(poly)
            fill.append("#ffffff" if not act[i, j, k] else
                        colours(None if side is None else side[i, j, k], k))
            if moved_ is not None and moved_[i, j, k]:
                bold.append(poly)
    ax.add_collection(PolyCollection(polys, facecolors=fill, edgecolors="#777777",
                                     linewidths=0.4))
    if bold:
        ax.add_collection(PolyCollection(bold, facecolors="none", edgecolors="#111111",
                                         linewidths=0.9))
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)
    seen = set()
    for faces in faults:
        for f in faces:
            for p in l32.face_pillars(f):
                if (kind == "row" and p[1] == n or kind == "col" and p[0] == n) and p not in seen:
                    seen.add(p)
                    x0, y0, z0, x1, y1, z1 = coord[p[1], p[0]]
                    ax.plot([(x0, y0)[axis], (x1, y1)[axis]], [z0, z1], color="black", lw=2.0)
    pts = np.array([p for poly in polys for p in poly])
    pts = pts[np.isfinite(pts).all(1)]
    if window:
        inside = (pts[:, 0] >= window[0]) & (pts[:, 0] <= window[1])
        pts = pts[inside] if inside.any() else pts
    ax.set_xlim(window or (pts[:, 0].min() - 5, pts[:, 0].max() + 5))
    ax.set_ylim(pts[:, 1].max() + 3, pts[:, 1].min() - 3)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(("x" if kind == "row" else "y") + " (m)", fontsize=8, labelpad=1)
    ax.set_ylabel("depth (m)", fontsize=8)
    ax.tick_params(labelsize=7)


def figure(name, title, panels, faults, cut, k=0, window=None, plan_window=None,
           after_window=None):
    """One PNG. ``panels``: (label, grid, side, moved) - the first is BEFORE, the last AFTER,
    any between are shown as sections for comparison. Top: the plan of BEFORE and AFTER
    (``after_window`` zooms the AFTER plan only); below: a true-scale section of every
    panel."""
    if OUT is None or not PNG:
        return
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    fig = plt.figure(figsize=(13.0, 5.4 + 2.2 * len(panels)), layout="constrained")
    spec = fig.add_gridspec(1 + len(panels), 2, height_ratios=[5.0] + [2.2] * len(panels))
    for col, (label, grid, side, mv) in enumerate((panels[0], panels[-1])):
        ax = fig.add_subplot(spec[0, col])
        zone_plan(ax, grid, side, mv, faults, f"{label}, plan: top of layer {k + 1}\n"
                  "fault at the top (solid) and base (dashed) of its pillars", k)
        zoom = after_window if (col == 1 and after_window) else plan_window
        if zoom:
            ax.set_xlim(zoom[0])
            ax.set_ylim(zoom[1])
    fig.legend(handles=[Patch(facecolor=FOOT, label="footwall band (-D)"),
                        Patch(facecolor=HANG, label="hanging-wall band (+D)"),
                        Patch(facecolor=OUTSIDE, edgecolor="#888888", label="outside the zone"),
                        Patch(facecolor="none", edgecolor="#111111", label="cell that moved")],
               loc="outside right upper", fontsize=8)
    for m, (label, grid, side, mv) in enumerate(panels):
        zone_section(fig.add_subplot(spec[1 + m, :]), grid, cut, side, mv, faults,
                     f"{label}: {cut[0]} {cut[1]}", window)
    fig.suptitle(title, fontsize=11)
    fig.savefig(OUT / f"33_{name}.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote 33_{name}.png")


# ── A vertical fault ─────────────────────────────────────────────────────────────────────────


def case_a():
    print("A vertical fault on x = 300 m, D = 60 m")
    grid, faces = cake(), x_fault(12)
    z = fault_zone(grid, {"F": faces}, 60.0, hanging_side=EAST)
    d = centres(grid)[..., 0] - 300.0                      # the plane x = 300, + east
    np.testing.assert_allclose(z.signed["F"][np.isfinite(z.signed["F"])],
                               d[np.isfinite(z.signed["F"])], atol=1e-9)
    assert np.array_equal(z.side("F"), band(d, 60.0, 60.0))
    fw, hw = z.counts()["F"]
    print(f"  {fw} footwall + {hw} hanging-wall cells: every cell with |x_centre - 300| <= 60 "
          f"(columns {int(np.nonzero(z.cells[:, 0, 0])[0].min())}.."
          f"{int(np.nonzero(z.cells[:, 0, 0])[0].max())}); distances equal x - 300 to 1e-9")
    print(f"  {z.summary()}")
    figure("A_vertical", "A: vertical fault, the cells within 60 m of it",
           [("BEFORE (the grid)", grid, None, None), ("AFTER (the zone selected)", grid,
                                                      z.side(), None)], [faces], ("row", 6))


# ── B dipping fault ──────────────────────────────────────────────────────────────────────────


def dipped(angle=60.0):
    grid, faces = cake(), x_fault(12)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = dip_grid(grid, Dip(faces, angle, hanging_side=EAST, name="F"))
    return grid, res.grid, faces


def case_b():
    print("B fault dipping 60 deg east, D = 50 m: the band slants with the fault")
    _, grid, faces = dipped()
    z = fault_zone(grid, {"F": faces}, 50.0)               # + read from the dip: east
    mid, normal, resid = fault_plane(grid, faces)
    normal = normal if normal[0] > 0 else -normal          # + toward the hanging wall (east)
    d = (centres(grid) - mid) @ normal
    seen = np.isfinite(z.signed["F"])
    err = float(np.abs(z.signed["F"][seen] - d[seen]).max())
    assert resid < 1e-6 and err < 1e-6, (resid, err)
    assert np.array_equal(z.side("F"), band(d, 50.0, 50.0))
    dip = np.degrees(np.arccos(abs(normal[2])))
    top, base = (np.nonzero(z.side("F")[:, 6, k])[0] for k in (0, 4))
    print(f"  the fault's pillars lie on one plane (residual {resid:.0e} m) dipping {dip:.6f} "
          f"deg; zone distances match the plane's to {err:.0e} m")
    print(f"  + side read from the dip: east. Row 6: columns {top.min()}..{top.max()} in the top "
          f"layer, {base.min()}..{base.max()} in the bottom one")
    figure("B_dipping", "B: fault dipping 60 deg, the cells within 50 m measured along its "
           "normal", [("BEFORE (the dipped grid)", grid, None, None),
                      ("AFTER (the zone selected)", grid, z.side(), None)], [faces], ("row", 6))


# ── C asymmetric ─────────────────────────────────────────────────────────────────────────────


def case_c():
    print("C asymmetric zone on the 60 deg fault: 30 m footwall, 90 m hanging wall")
    _, grid, faces = dipped()
    z = fault_zone(grid, {"F": faces}, (30.0, 90.0))
    mid, normal, _ = fault_plane(grid, faces)
    d = (centres(grid) - mid) @ (normal if normal[0] > 0 else -normal)
    assert np.array_equal(z.side("F"), band(d, 30.0, 90.0))
    fw, hw = z.counts()["F"]
    print(f"  {fw} footwall + {hw} hanging-wall cells; every cell with -30 <= d <= 90 from the "
          f"plane's equation, none other")
    figure("C_asymmetric", "C: 30 m on the footwall, 90 m on the hanging wall",
           [("BEFORE (the dipped grid)", grid, None, None),
            ("AFTER (the zone selected)", grid, z.side(), None)], [faces], ("row", 6))


# ── D curved fault ───────────────────────────────────────────────────────────────────────────


def case_d():
    print("D curved fault: an arc adapted into the pillars, dipped 65 deg, D = 50 m")
    grid = cake()
    th = np.radians(np.linspace(-45.0, 45.0, 41))
    grid, faces = l32.adapted(grid, np.c_[40.0 + 260.0 * np.cos(th), 150.0 + 260.0 * np.sin(th)],
                              quiet=True)
    res, _ = l32.dip("D", grid, Dip(faces, 65.0, hanging_side=(560.0, 150.0), name="ARC"),
                     quiet=True)
    grid = res.grid
    z = fault_zone(grid, {"ARC": faces}, 50.0)
    pts = surface_samples(grid, faces, pad=60.0)
    brute = cKDTree(pts).query(centres(grid).reshape(-1, 3))[0].reshape(grid.cart_dims)
    seen = np.isfinite(z.signed["ARC"])
    err = float(np.abs(np.abs(z.signed["ARC"][seen]) - brute[seen]).max())
    tol = 1.0                                              # the sampling step
    sure = np.abs(brute - 50.0) > tol
    agree = (z.cells == (brute <= 50.0))[sure].all()
    assert err < tol and agree, (err, agree)
    fw, hw = z.counts()["ARC"]
    print(f"  {fw} footwall + {hw} hanging-wall cells; brute force over {len(pts):,} surface "
          f"points: distances agree to {err:.2f} m (sampling step 1 m), and so does every cell "
          f"farther than 1 m from the zone's edge ({int((~sure).sum())} nearer, not judged)")
    figure("D_curved", "D: curved fault dipping 65 deg, the cells within 50 m",
           [("BEFORE (adapted and dipped)", grid, None, None),
            ("AFTER (the zone selected)", grid, z.side(), None)], [faces], ("row", 6), k=4)


# ── E two faults ─────────────────────────────────────────────────────────────────────────────


def case_e():
    print("E two crossing faults: F1 on x = 200 m (+ east), F2 on y = 150 m (+ north), D = 60 m")
    grid = cake()
    faults = {"F1": x_fault(8), "F2": y_fault(6)}
    z = fault_zone(grid, faults, 60.0, hanging_side={"F1": (575.0, 10.0), "F2": (10.0, 290.0)})
    c = centres(grid)
    d1, d2 = c[..., 0] - 200.0, c[..., 1] - 150.0
    near = np.where(np.abs(d2) < np.abs(d1), d2, d1)       # a tie goes to the first (F1)
    assert np.array_equal(z.side(), band(near, 60.0, 60.0))
    assert np.array_equal(z.side("F1"), band(d1, 60.0, 60.0))
    assert np.array_equal(z.side("F2"), band(d2, 60.0, 60.0))
    both = (np.abs(d1) <= 60.0) & (np.abs(d2) <= 60.0)
    ties = both & (np.abs(d1) == np.abs(d2))
    counts = z.counts()
    print(f"  F1 {counts['F1'][0]} + {counts['F1'][1]}, F2 {counts['F2'][0]} + {counts['F2'][1]}"
          f" cells (footwall + hanging wall, each on its own); {int(both.sum())} cells are within"
          f" 60 m of both: the nearer fault decides their side ({int(ties.sum())} ties go to F1)")
    print(f"  union: {int(z.minus.sum())} footwall + {int(z.plus.sum())} hanging-wall cells")
    figure("E_two_faults", "E: two crossing faults, each cell assigned to the nearer one",
           [("BEFORE (the grid)", grid, None, None),
            ("AFTER (the zone selected)", grid, z.side(), None)], list(faults.values()),
           ("row", 2))


# ── F GRDECL export ──────────────────────────────────────────────────────────────────────────


def read_keyword(path, keyword):
    """Values of ``keyword`` in a GRDECL file, run-length ``n*v`` expanded (no library code)."""
    words = path.read_text().split()
    out = []
    for w in words[words.index(keyword) + 1:]:
        if w == "/":
            break
        n, _, v = w.rpartition("*")
        out.extend([int(v)] * (int(n) if n else 1))
    return np.array(out)


def case_f():
    print("F GRDECL export of the B zone (60 deg fault, D = 50 m)")
    _, grid, faces = dipped()
    z = fault_zone(grid, {"F": faces}, 50.0)
    with TemporaryDirectory() as work:
        path = write_zone(z, Path(work) / "zone.grdecl")
        text = path.read_text()
        vals = read_keyword(path, "FLTZONE")
    nx, ny, nz = grid.cart_dims
    side = z.side()
    assert len(vals) == nx * ny * nz
    for i, j, k in ((0, 0, 0), (11, 6, 0), (12, 6, 4), (23, 11, 4)):
        want = {-1: 1, 0: 0, 1: 2}[int(side[i, j, k])]
        assert vals[i + nx * (j + ny * k)] == want, (i, j, k)
    assert np.array_equal(vals.reshape(nz, ny, nx).transpose(2, 1, 0),
                          np.select([side == -1, side == 1], [1, 2], 0))
    print(f"  {len(text.splitlines())} lines, {len(text)} bytes; read back by hand: "
          f"{len(vals)} values, I fastest; {int((vals == 1).sum())} x 1 (footwall band), "
          f"{int((vals == 2).sum())} x 2 (hanging-wall band), {int((vals == 0).sum())} x 0")
    print("  first lines: " + " | ".join(text.splitlines()[:3]))
    back = np.select([vals == 1, vals == 2], [-1, 1], 0).reshape(nz, ny, nx).transpose(2, 1, 0)
    figure("F_export", "F: the B zone written to GRDECL (FLTZONE 0/1/2) and read back by hand",
           [("BEFORE: the zone in memory", grid, side, None),
            ("AFTER: FLTZONE read back from the file", grid, back, None)], [faces], ("row", 6))


# ── G confined dip ───────────────────────────────────────────────────────────────────────────


def case_g():
    print("G a dip of 60 deg: the whole domain (default), boundary='free', and inside 100 m")
    grid, faces = cake(), x_fault(12)
    spec = Dip(faces, 60.0, hanging_side=EAST, name="F")
    z = fault_zone(grid, {"F": faces}, 100.0, hanging_side=EAST)
    runs = {}
    for label, kw in (("whole domain, rim fixed (default)", {}),
                      ("whole domain, boundary='free'", {"boundary": "free"}),
                      ("inside the 100 m zone", {"zone": z})):
        res, a = l32.dip("G", grid, spec, quiet=True, **kw)
        mv = moved(grid, res.grid)
        got = l32.angles(res, "F")
        assert np.abs(got - 60.0).max() < 1e-9
        runs[label] = (res, mv)
        print(f"  {label:34s}: F {got.min():.9f}..{got.max():.9f} deg; {int(mv.sum()):4d} of "
              f"{mv.size} cells moved, {int((mv & ~z.cells).sum()):4d} of them outside the zone;"
              f" cells change volume up to {100 * np.nanmax(np.abs(a.dv)):.1f} %; {a.n_bad} broken")
    res, mv = runs["inside the 100 m zone"]
    assert not (mv & ~z.cells).any(), "a confined dip moved a cell outside its zone"
    print(f"  confined: the {int((~z.cells).sum())} cells outside the zone are bit-identical "
          f"(raw COORD and ZCORN); {int((z.cells & ~mv).sum())} zone cells next to its edge "
          "stay put too (a pillar moves only when every cell around it is in the zone)")
    figure("G_confined_dip", "G: 60 deg dip inside a 100 m zone, against the whole-domain "
           "dips", [("BEFORE (the zone)", grid, z.side(), None)]
           + [(label, r.grid, z.side(), m) for label, (r, m) in runs.items()],
           [faces], ("row", 6))


# ── H confined slide ─────────────────────────────────────────────────────────────────────────


def case_h():
    print("H a slide of 20 m: the whole hanging wall (default) and inside 75 m")
    grid, faces = cake(), x_fault(12)
    throw = Throw(faces, 20.0, hanging_side=EAST, name="F")
    z = fault_zone(grid, {"F": faces}, 75.0, hanging_side=EAST)
    runs = {}
    for label, kw in (("whole domain (default)", {}), ("inside the 75 m zone", {"zone": z})):
        res, a = l32.slide("H", grid, throw, quiet=True, **kw)
        mv = moved(grid, res.grid)
        t = realized_throw(grid, res.grid, faces, EAST)[:, 1]
        assert np.allclose(t, 20.0, atol=1e-9), t
        runs[label] = (res, mv)
        print(f"  {label:24s}: throw {t.min():.9f}..{t.max():.9f} m; {int(mv.sum()):4d} cells "
              f"moved, {int((mv & ~z.cells).sum()):4d} outside the zone; {a.n_bad} broken")
    res, mv = runs["inside the 75 m zone"]
    assert not (mv & ~z.cells).any(), "a confined slide moved a cell outside its zone"
    print("  confined: the hanging wall drops 20 m at the fault and the drop fades to 0 at the "
          "zone's edge (a cosine); outside it every corner is bit-identical")
    figure("H_confined_slide", "H: 20 m slide inside a 75 m zone, against the whole-domain "
           "slide", [("BEFORE (the zone)", grid, z.side(), None)]
           + [(label, r.grid, z.side(), m) for label, (r, m) in runs.items()],
           [faces], ("row", 6))


# ── I dip and slide confined ─────────────────────────────────────────────────────────────────


def case_i():
    print("I dip 60 deg and slide 15 m in one call, both inside 120 m")
    grid, faces = cake(), x_fault(12)
    z = fault_zone(grid, {"F": faces}, 120.0, hanging_side=EAST)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = dip_grid(grid, Dip(faces, 60.0, hanging_side=EAST, name="F"), zone=z,
                       slide=Throw(faces, 15.0, hanging_side=EAST, name="F"))
    assert res.slide is not None and res.slide.zone is not None
    out = res.slide.grid
    a = audit(grid, out)
    assert a.n_bad == 0, a.line()
    mv = moved(grid, out)
    assert not (mv & ~z.cells).any()
    got = l32.angles(res, "F")
    t = realized_throw(res.grid, out, faces, EAST)[:, 1]
    assert np.abs(got - 60.0).max() < 1e-9 and np.allclose(t, 15.0, atol=1e-9)
    print(f"  F {got.min():.6f} deg, throw {t.min():.6f}..{t.max():.6f} m; {int(mv.sum())} cells "
          f"moved, 0 outside the zone; {a.line()}")
    figure("I_dip_and_slide", "I: dip 60 deg and slide 15 m, both inside a 120 m zone",
           [("BEFORE (the zone)", grid, z.side(), None),
            ("dipped", res.grid, z.side(), moved(grid, res.grid)),
            ("AFTER (dipped and slid)", out, z.side(), mv)], [faces], ("row", 6))


# ── J refusals ───────────────────────────────────────────────────────────────────────────────


def case_j():
    print("J refusals")
    grid, faces = cake(), x_fault(12)
    try:
        dip_grid(grid, Dip(faces, 60.0, hanging_side=EAST), zone=10.0)
        raise AssertionError("a 10 m zone around 25 m cells should be refused")
    except ZoneError as exc:
        print(f"  dip inside 10 m (cells are 25 m wide): ZoneError: {str(exc)[:110]}...")
    try:
        slide_grid(grid, Throw(faces, 20.0, hanging_side=EAST), zone=5.0)
        raise AssertionError("a 5 m zone should be refused")
    except ZoneError as exc:
        print(f"  slide inside 5 m: ZoneError: {str(exc)[:110]}...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            dip_grid(grid, Dip(faces, 30.0, hanging_side=EAST), zone=40.0)
            raise AssertionError("30 deg inside 40 m should be refused")
        except DipError as exc:
            print(f"  30 deg inside 40 m: DipError: {str(exc)[:110]}...")
        loose = dip_grid(grid, Dip(faces, 30.0, hanging_side=EAST), zone=40.0,
                         on_invalid="ignore")
    a = audit(grid, loose.grid)
    assert a.n_bad > 0
    print(f"  the grid it refused, judged independently: {a.line()}")
    z40 = fault_zone(grid, {"F": faces}, 40.0, hanging_side=EAST)
    figure("J_refused", "J: 30 deg inside a 40 m zone is refused; the grid it would have "
           "returned (forced with on_invalid='ignore') folds",
           [("BEFORE: the 40 m zone", grid, z40.side(), None),
            ("AFTER (REFUSED; forced only to show why)", loose.grid, z40.side(),
             moved(grid, loose.grid))], [faces], ("row", 6))
    z = corners(grid)[..., 2]
    lean = float(np.nanmax(z) - np.nanmin(z)) / np.tan(np.radians(30.0))
    print(f"  a zone only removes room: at 30 deg the fault's base sits {lean:.0f} m from its top "
          "in plan, and the pillars within 40 m must follow it, so the cells fold; the gate "
          "says so rather than returning them")


# ── K Norne ──────────────────────────────────────────────────────────────────────────────────


def case_k():
    print("K Norne: the zone of every fault at D = 100 m; E_01 dipped 60 deg inside 400 m")
    with TemporaryDirectory() as work:
        deck = l32.norne_deck(Path(work))
        if deck is None:
            print("  skipped: Norne is not on this machine (set CPG_NORNE_GRDECL)")
            return
        grid = read_grdecl(deck)
        faults = _deck_faults(deck, grid.cart_dims[2])[0]
    nx, ny, nz = grid.cart_dims
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        z = fault_zone(grid, faults, 100.0)
    act = np.asarray(grid.actnum).reshape(nz, ny, nx).transpose(2, 1, 0).astype(bool)
    print(f"  {nx}x{ny}x{nz}, {len(faults)} faults: {int(z.minus.sum()):,} footwall + "
          f"{int(z.plus.sum()):,} hanging-wall cells ({int((z.cells & act).sum()):,} active) "
          f"within 100 m; {len(z.assumed)} near-vertical faults took + on the left of the chain "
          f"({len(caught)} warnings)")
    # a sample of E_01's cells re-measured by brute force against its pillar surface
    e = fault_zone(grid, {"E_01": faults["E_01"]}, 100.0)
    # its pillars are far from parallel: a face between two of them twists, and 4 flat
    # triangles a face would be ~20 m off here; the library measures the ruled surface exactly
    pts = surface_samples(grid, faults["E_01"], pad=110.0, n_depth=1200, n_along=80)
    rng = np.random.default_rng(1)
    pick = rng.choice(np.flatnonzero(np.isfinite(e.signed["E_01"].ravel())), 300, replace=False)
    lib = np.abs(e.signed["E_01"].ravel()[pick])
    brute = cKDTree(pts).query(centres(grid).reshape(-1, 3)[pick])[0]
    print(f"  E_01: 300 random zone cells re-measured by brute force ({len(pts):,} surface "
          f"points): the library is {np.max(brute - lib):.3f} m nearer at most, never farther "
          f"(by {max(np.max(lib - brute), 0.0):.0e} m)")
    assert np.all(lib <= brute + 1e-9) and np.max(brute - lib) < 0.05
    zone = fault_zone(grid, {"E_01": faults["E_01"]}, 400.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = dip_grid(grid, Dip("E_01", 60.0), faults=faults, zone=zone)
    mv = moved(grid, res.grid)
    a = audit(grid, res.grid)
    got = res.realized_dip("E_01")[:, 1]
    assert a.n_bad == 0 and not (mv & ~zone.cells).any() and np.allclose(got, 60.0, atol=1e-6)
    print(f"  E_01 -> 60 deg inside 400 m: exact on {len(got)} pillars; {int(mv.sum()):,} cells "
          f"moved ({a.n_moved:,} of them active), 0 of the {int((~zone.cells).sum()):,} outside;"
          f" {a.line()}")
    chain = res.chains["E_01"]
    mid = chain.ij[len(chain.ij) // 2]
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)
    cx, cy = coord[mid[1], mid[0], :2]
    lo, hi = chain.xy.min(0) - 700.0, chain.xy.max(0) + 700.0
    figure("K_norne_zones", "K: Norne. BEFORE: every fault's 100 m zone. AFTER: E_01 dipped "
           "60 deg inside its 400 m zone (zoomed)",
           [("BEFORE: every fault, D = 100 m", grid, z.side(), None),
            ("E_01's 400 m zone", grid, zone.side(), None),
            ("AFTER: E_01 at 60 deg inside it", res.grid, zone.side(), mv)],
           list(faults.values()), ("row", int(min(mid[1], ny - 1))),
           window=(cx - 900.0, cx + 900.0), after_window=((lo[0], hi[0]), (lo[1], hi[1])))


CASES = {}


def main(out, only=None):
    global OUT
    OUT = out
    for key, fn in sorted(CASES.items()):
        if only is None or key in only:
            fn()
    print("every selection re-measured by hand; every edit re-judged by lesson 32's checker")


if __name__ == "__main__":
    CASES.update({k[-1].upper(): v for k, v in globals().items() if k.startswith("case_")})
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, help="write the figures here")
    parser.add_argument("--case", help="comma-separated case letters to run (default: all)")
    args = parser.parse_args()
    only = set(args.case.upper().split(",")) if args.case else None
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        main(args.out, only)
    else:
        main(None, only)
