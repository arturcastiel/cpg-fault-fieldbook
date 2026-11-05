"""Stress the fault dip: 22 cases that mix it with curved traces, slides, growth, blind tips,
holes, junctions and field data, each judged by a cell checker that shares no code with the
library.

The checker (``audit`` below) rebuilds every cell from raw COORD and ZCORN and measures the
trilinear hexahedron itself: scaled Jacobians at 45 points per cell, volumes, and the plan
quadrilaterals of each column's four pillars at depths all four carry corners. A cell is
broken when it inverts, when its worst Jacobian falls below 0.05 *and* below what it was
(a 1 % wobble allowed), when it loses 90 % of its volume, or when its pillars fold in plan.
Only active cells are judged; a fold under a column with no active cell is reported apart.
Dips run with the default ``boundary="fixed"`` are also checked, from raw COORD, to hold the
grid's outline (no rim pillar moves across its side; the corners do not move). Every
accepted grid below has zero broken cells; every refusal is shown, and so is what the
checker says about the grid the gate refused.

- A angle sweep: one rim-to-rim fault set to 30, 35, 45, 60, 75 and 89 deg (boundary="free":
  the rigid hanging wall leans the east side of the grid too);
- B curved trace: an arc adapted into the pillars (adapt_grid), dipped 60 deg, slid 20 m;
- C S-shaped trace with an angle table 50 -> 80 deg and a throw table;
- D stair-stepped diagonal FAULTS (not adapted): dipped along their smoothed strike, slid;
- E fault with two tips: an angle per pillar (55 -> 80 -> 55) and a footwall-up throw;
- F graben and horst: four faults, two pairs dipping toward and away from each other; at 60
  deg the graben faults meet just below the cells and no slide fits (honest refusal);
- G domino: four parallel faults at 50, 55, 60 and 65 deg, each block stepped 10 m;
- H T-junction: a fault that ends on another; both dipped, the shared pillar compromises;
- I bedding dipping along strike: the rigid hanging wall keeps every cell volume
  (boundary="free");
- J pre-tilted pillars: the grid leans east, the fault is turned to dip west;
- K re-dip: dip 60, slide 15, re-dip 75, slide 10 more, re-dip 50; throws add up;
- L thin, pinched-out layers and ACTNUM holes;
- M skewed, uneven columns (sheared and stretched plan);
- N dome with a ring fault: dipped inward 70 deg and collapsed 25 m (a caldera);
- O shallow angles on bedding dipping toward the fault: with the rim fixed (default) 30 deg
  widens the fades and 15 deg is refused (the parallel east side cannot hold still); with
  boundary="free" 15 deg widens the fades, no widening raises, bisect stops short and says
  so, and 8 deg is refused outright;
- P random angles 40..85 deg: noise at every pillar folds the fault onto itself (refused);
  the same range drawn smoothly along strike is exact;
- Q the slide's cell gate: a throw on leaning pillars that inverts cells unless checked;
  dip_grid(..., slide=...) widens the fades until the full throw fits (cells change volume);
- R a conservative refusal: bedding steeper than the fault folds the plan quads but not the
  3-D cells; the gate refuses a grid the checker would pass;
- S growth fault: a dipped fault whose throw grows 1.0 -> 1.6 with depth;
- T blind fault: a dipped fault that cuts only the lower three layers;
- U seeded fuzz: random grids and faults, dip -> slide -> re-dip; the gate never accepts a
  grid the checker breaks, and every refusal is re-judged too;
- V a field deck (Norne, when it is on this machine; env CPG_NORNE_GRDECL as in lesson 30):
  E_01 and BC set to 60 deg; one accepted, one refused, both verdicts re-checked.

With matplotlib installed and --out DIR, every case writes one PNG figure. --case B,C runs
only those cases.
"""

import argparse
import functools
import importlib.util
import os
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import CornerPointGrid, read_grdecl
from fault_projections import (Dip, DipError, SlideError, Throw, adapt_grid, dip_grid,
                               fault_faces_for, slide_grid)
from fault_projections.dip import _deck_faults
from fault_projections.faults import FaultFace

from fieldbook.zslide import layer_cake, realized_throw

NX, NY, DX = 24, 12, 25.0                 # 600 m x 300 m in plan
DZ = [15.0, 25.0, 10.0, 30.0, 20.0]      # 100 m thick: cells about as tall as they are wide
TOP = 2000.0
PNG = importlib.util.find_spec("matplotlib") is not None
NORNE = "/mnt/c/projects/opm-machine/opm-tests/norne/INCLUDE/GRID/IRAP_1005.GRDECL"
MORPH = dict(boundary="slide", w0=0.05, polish="barrier", on_invalid="raise")
EAST, WEST = "east", "west"
OUT = None                       # set by main: where figures go (None: no figures)


def cake(top=None, dz=DZ, nx=NX, ny=NY, **kw):
    return layer_cake(nx, ny, len(dz), DX, dz,
                      top=top or (lambda x, y: TOP + 0.05 * x), **kw)


def x_fault(i, j0=1, j1=NY):
    return [FaultFace(i, j, "X") for j in range(j0, j1 + 1)]


def y_fault(j, i0=1, i1=NX):
    return [FaultFace(i, j, "Y") for i in range(i0, i1 + 1)]


def side(grid, faces, where):
    """A plan point on the east/west (or north/south) side of a straight X (or Y) fault."""
    f = faces[len(faces) // 2]
    nx, ny, _ = grid.cart_dims
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)
    if f.face == "X":
        p = coord[f.j, f.i, :2]
        return (p[0] + (DX if where == EAST else -DX), p[1])
    p = coord[f.j, f.i, :2]
    return (p[0], p[1] + (DX if where == "north" else -DX))


# ── the independent checker: raw COORD + ZCORN -> hexahedra; no library code ────────────────


def corners(grid):
    """(NX, NY, NZ, 2, 2, 2, 3) corner xyz, indexed [u, v, b] = [east, north, base]."""
    nx, ny, nz = grid.cart_dims
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6).transpose(1, 0, 2)
    z = np.asarray(grid.zcorn, float).reshape(nz, 2, ny, 2, nx, 2).transpose(4, 2, 0, 5, 3, 1)
    out = np.empty((nx, ny, nz, 2, 2, 2, 3))
    for u in (0, 1):
        for v in (0, 1):
            p = coord[u:u + nx, v:v + ny][:, :, None, None, :]
            zz = z[:, :, :, u, v, :]
            with np.errstate(divide="ignore", invalid="ignore"):   # zero-height pillars
                t = (zz - p[..., 2]) / (p[..., 5] - p[..., 2])
            out[..., u, v, :, 0] = p[..., 0] + t * (p[..., 3] - p[..., 0])
            out[..., u, v, :, 1] = p[..., 1] + t * (p[..., 4] - p[..., 1])
            out[..., u, v, :, 2] = zz
    return out


def scaled_jacobians(c):
    """(NX, NY, NZ, 45) scaled Jacobians at xi, eta in {0, .5, 1}, zeta in {0, .25, .., 1}."""
    di, dj, dk = c[..., 1, :, :, :] - c[..., 0, :, :, :], c[..., :, 1, :, :] - c[..., :, 0, :, :], \
        c[..., :, :, 1, :] - c[..., :, :, 0, :]
    out = []
    for xi in (0.0, 0.5, 1.0):
        for eta in (0.0, 0.5, 1.0):
            for zeta in (0.0, 0.25, 0.5, 0.75, 1.0):
                nu, nv, nb = [1 - xi, xi], [1 - eta, eta], [1 - zeta, zeta]
                e1 = np.einsum("...vbk,v,b->...k", di, nv, nb)
                e2 = np.einsum("...ubk,u,b->...k", dj, nu, nb)
                e3 = np.einsum("...uvk,u,v->...k", dk, nu, nv)
                det = np.einsum("...k,...k->...", np.cross(e1, e2), e3)
                n = np.prod([np.linalg.norm(e, axis=-1) for e in (e1, e2, e3)], axis=0)
                with np.errstate(divide="ignore", invalid="ignore"):
                    out.append(det / n)
    return np.stack(out, -1)


def hex_volumes(c):
    """Cell volumes: each face fanned into 4 triangles about its centroid."""
    faces = [((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)), ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)),
             ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)), ((0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0)),
             ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)), ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1))]
    centre = c.reshape(c.shape[:3] + (8, 3)).mean(-2)
    vol = np.zeros(c.shape[:3])
    for face in faces:
        p = [c[..., u, v, b, :] for u, v, b in face]
        fc = sum(p) / 4.0
        for m in range(4):
            vol += np.einsum("...k,...k->...", np.cross(p[m] - centre, p[(m + 1) % 4] - centre),
                             fc - centre) / 6.0
    return vol


def lattice_folds(before, after):
    """(NX, NY): columns whose pillars fold in plan (a corner turns over) at a depth all four
    carry corners; a corner folded before the edit is not counted again."""
    nx, ny, nz = after.cart_dims
    z = np.asarray(after.zcorn, float).reshape(nz, 2, ny, 2, nx, 2).transpose(4, 2, 0, 5, 3, 1)
    zlo, zhi = z.min(axis=(2, 5)).max(axis=(2, 3)), z.max(axis=(2, 5)).min(axis=(2, 3))
    live = zhi > zlo

    def quad(grid, depth):
        coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6).transpose(1, 0, 2)
        xy = []
        for u, v in ((0, 0), (1, 0), (1, 1), (0, 1)):
            p = coord[u:u + nx, v:v + ny]
            with np.errstate(divide="ignore", invalid="ignore"):   # zero-height pillars
                t = (depth - p[..., 2]) / (p[..., 5] - p[..., 2])
                xy.append(p[..., :2] + t[..., None] * (p[..., 3:5] - p[..., :2]))
        def cross(a, b):
            return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]

        return np.stack([cross(xy[(m + 1) % 4] - xy[m], xy[(m - 1) % 4] - xy[m])
                         for m in range(4)], -1)

    bad = np.zeros((nx, ny), bool)
    for s in np.linspace(0.0, 1.0, 9):
        depth = zlo + s * (np.where(live, zhi, zlo) - zlo)
        q0, q1 = quad(before, depth), quad(after, depth)
        sign = np.sign(q0.sum(-1, keepdims=True))
        bad |= ((q0 * sign > 0) & (q1 * sign <= 0)).any(-1) & live   # newly folded corners
    return bad


class Audit:
    """What the checker found in ``after`` relative to ``before``. Only cells with a corner
    that moved (by more than 1 micrometre: round-off in field coordinates is not a move) are
    measured: a cell whose eight corners stay where they were cannot break."""

    def __init__(self, before, after, min_jac=0.05):
        nx, ny, nz = before.cart_dims
        act = np.asarray(before.actnum).reshape(nz, ny, nx).transpose(2, 1, 0).astype(bool)
        c0, c1 = corners(before), corners(after)
        moved = ~np.isclose(c0, c1, rtol=0.0, atol=1e-6, equal_nan=True).all(axis=(3, 4, 5, 6))
        j0, j1 = scaled_jacobians(c0[moved]), scaled_jacobians(c1[moved])
        sign = np.sign(np.nansum(j0[act[moved]])) or 1.0
        j0, j1 = j0 * sign, j1 * sign
        live = act[moved][:, None] & np.isfinite(j0) & (j0 > 0)
        m0, m1 = np.where(live, j0, np.inf).min(-1), np.where(live, j1, np.inf).min(-1)
        v0, v1 = (hex_volumes(c[moved][:, None, None])[:, 0, 0] for c in (c0, c1))
        full = act[moved] & (v0 > 1e-9 * max(np.abs(v0).max(initial=0.0), 1e-300))
        shape = (nx, ny, nz)
        self.inverted, self.degenerate, self.crushed = (np.zeros(shape, bool) for _ in "idc")
        self.inverted[moved] = (live & ~(j1 > 0)).any(-1)
        self.degenerate[moved] = (m1 < min_jac) & (m1 < 0.99 * m0 - 1e-9) & ~(live & ~(j1 > 0)).any(-1)
        with np.errstate(divide="ignore", invalid="ignore"):
            self.crushed[moved] = full & ~(v1 / v0 > 0.1)
            self.dv = np.where(act, 0.0, np.nan)             # unmoved cells keep their volume
            self.dv[moved] = np.where(full, v1 / v0 - 1.0, np.nan)
        folds = lattice_folds(before, after)
        self.folds = folds & act.any(-1)                 # columns with a live cell
        self.dead_folds = folds & ~act.any(-1)           # reported, but no cell to break
        self.bad = (self.inverted | self.degenerate | self.crushed | self.folds[..., None]) & act
        self.n_moved = int((moved & act).sum())
        self.min_before = float(np.where(live, j0, np.inf).min(initial=np.inf))
        self.min_after = float(np.where(live, j1, np.inf).min(initial=np.inf))
        self._c1, self._act, self._sign = c1, act, sign

    @property
    def cell_min(self):
        """(NX, NY, NZ) worst 3-D scaled Jacobian of every active cell of ``after``."""
        j = scaled_jacobians(self._c1) * self._sign
        m = np.where(np.isfinite(j), j, np.inf).min(-1)
        return np.where(self._act & np.isfinite(m), m, np.nan)

    @property
    def n_bad(self):
        return int(self.bad.sum())

    def line(self):
        return (f"checker: {int(self.inverted.sum())} inverted, {int(self.degenerate.sum())} "
                f"degenerate, {int(self.crushed.sum())} crushed, {int(self.folds.sum())} plan folds"
                + (f" (+{int(self.dead_folds.sum())} in all-inactive columns)"
                   if self.dead_folds.any() else "")
                + f"; {self.n_moved} cells moved, their min 3-D Jacobian "
                  f"{self.min_before:.3g} -> {self.min_after:.3g}")


def audit(before, after):
    return Audit(before, after)


# ── checked steps ───────────────────────────────────────────────────────────────────────────


def coord_dip(grid, ij, normal):
    """Dip (deg) of pillar ij along ``normal``, straight from its six COORD numbers."""
    nx, ny, _ = grid.cart_dims
    x0, y0, z0, x1, y1, z1 = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)[ij[1], ij[0]]
    return float(np.degrees(np.arctan2(z1 - z0, (x1 - x0) * normal[0] + (y1 - y0) * normal[1])))


def angles(res, name):
    """Dips read from COORD at the chain pillars of ``name``, along the result's normals."""
    chain = res.chains[name]
    return np.array([coord_dip(res.grid, chain.ij[m], res.normals[name][m])
                     for m in range(len(chain.ij))])


def outline_shift(before, after):
    """Largest move of a rim pillar ACROSS its side of the grid (top and base points), and of
    a corner pillar in any direction, from raw COORD: what boundary="fixed" holds at 0."""
    nx, ny, _ = before.cart_dims
    c0 = np.asarray(before.coord, float).reshape(ny + 1, nx + 1, 6)
    c1 = np.asarray(after.coord, float).reshape(ny + 1, nx + 1, 6)
    worst = 0.0
    sides = ((np.s_[0, :], c0[0, :]), (np.s_[ny, :], c0[ny, :]),
             (np.s_[:, 0], c0[:, 0]), (np.s_[:, nx], c0[:, nx]))
    for idx, line in sides:
        d = c1[idx] - line
        for off in (0, 3):
            t = np.gradient(line[:, off:off + 2], axis=0)
            nrm = np.stack([-t[:, 1], t[:, 0]], -1) / np.linalg.norm(t, axis=1)[:, None]
            worst = max(worst, float(np.abs((d[:, off:off + 2] * nrm).sum(-1)).max()))
    for j, i in ((0, 0), (0, nx), (ny, 0), (ny, nx)):
        worst = max(worst, float(np.abs(c1[j, i, [0, 1, 3, 4]] - c0[j, i, [0, 1, 3, 4]]).max()))
    return worst


def dip(tag, grid, dips, *, exact=True, tol=1e-6, quiet=False, warn=False, **kw):
    """dip_grid, then assert the checker finds nothing and the angles hold (from COORD); with
    the default boundary="fixed" also that no rim pillar crossed its side of the grid.
    The "fades widened" warning is silenced unless ``warn``: the dip line prints the factor."""
    with warnings.catch_warnings():
        if not warn:
            warnings.filterwarnings("ignore", message="dip fades widened")
        res = dip_grid(grid, dips, **kw)
    a = audit(grid, res.grid)
    assert a.n_bad == 0, f"{tag}: the gate accepted a grid the checker breaks: {a.line()}"
    assert np.array_equal(res.grid.zcorn, grid.zcorn), "a dip never edits ZCORN"
    if kw.get("boundary", "fixed") == "fixed":
        rim = outline_shift(grid, res.grid)
        assert rim < 1e-6 * max(1.0, float(np.abs(np.asarray(grid.coord)).max())), (tag, rim)
    parts = []
    for name in res.chains:
        got, want = angles(res, name), res.requested[name]
        if exact and res.scale == 1.0:
            assert np.abs(got - want).max() < tol, (tag, name, got, want)
        parts.append(f"{name} {got.min():.2f}..{got.max():.2f} deg "
                     f"({res.modes[name]}, fades x{res.widen:g})")
    if not quiet:
        print(f"  dip    {'; '.join(parts)}")
        print(f"         {a.line()}")
    return res, a


def slide(tag, grid, throws, *, quiet=False, **kw):
    """slide_grid, then assert the checker finds nothing and COORD is untouched."""
    res = slide_grid(grid, throws, **kw)
    a = audit(grid, res.grid)
    assert a.n_bad == 0, f"{tag}: the slide gate accepted a grid the checker breaks: {a.line()}"
    assert np.array_equal(res.grid.coord, grid.coord), "a slide never edits COORD"
    if not quiet:
        print(f"  slide  {a.line()}")
    return res, a


def throw_of(before, after, faces, hanging):
    t = realized_throw(before, after, faces, hanging)[:, 1]
    return t


# ── figures ─────────────────────────────────────────────────────────────────────────────────


def face_pillars(f):
    return ((f.i, f.j - 1), (f.i, f.j)) if f.face == "X" else ((f.i - 1, f.j), (f.i, f.j))


def plan(ax, grid, a, faults, title):
    """Columns at the top coloured by their worst cell's 3-D scaled Jacobian; fault traces at
    the top (solid) and the base (dashed) of the pillars; broken cells outlined red."""
    from matplotlib.collections import PolyCollection

    nx, ny, nz = grid.cart_dims
    c = corners(grid)
    worst = np.nanmin(np.where(np.isnan(a.cell_min), np.inf, a.cell_min), axis=2)
    act = np.isfinite(worst)
    polys = [c[i, j, 0, [0, 1, 1, 0], [0, 0, 1, 1], 0, :2] for i in range(nx) for j in range(ny)]
    vals = np.where(act, worst, np.nan).ravel()
    pc = PolyCollection(polys, array=np.ma.masked_invalid(vals), cmap="viridis",
                        edgecolors="#555555", linewidths=0.25)
    pc.set_clim(0.0, 1.0)
    pc.cmap.set_bad("#dddddd")
    ax.add_collection(pc)
    bad = a.bad.any(-1)
    if bad.any():
        ax.add_collection(PolyCollection([polys[i * ny + j] for i, j in zip(*np.nonzero(bad))],
                                         facecolors="none", edgecolors="#d7301f", linewidths=1.6))
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)
    for faces in faults:
        for f in faces:
            p, q = face_pillars(f)
            for off, style in ((0, "-"), (3, "--")):
                ax.plot([coord[p[1], p[0], off], coord[q[1], q[0], off]],
                        [coord[p[1], p[0], off + 1], coord[q[1], q[0], off + 1]],
                        style, color="black", lw=1.8 if off == 0 else 1.0)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.tick_params(axis="x", labelrotation=30)
    return pc


def section(ax, grid, cut, title, a=None, faults=(), window=None):
    """Cell faces in a vertical section: ("row", j) shows the south faces of row j in x-z,
    ("col", i) the west faces of column i in y-z. Fault pillars in black; cells the checker
    finds broken (``a``) filled red."""
    from matplotlib.collections import PolyCollection

    nx, ny, nz = grid.cart_dims
    kind, n = cut
    c = corners(grid)
    act = np.asarray(grid.actnum).reshape(nz, ny, nx).transpose(2, 1, 0).astype(bool)
    polys, colours, live = [], [], []
    palette = ("#c6dbef", "#a1d99b", "#fdd0a2", "#dadaeb", "#fcbba1", "#d9f0a3", "#bcbddc")
    cells = [(i, n) for i in range(nx)] if kind == "row" else [(n, j) for j in range(ny)]
    axis = 0 if kind == "row" else 1
    for i, j in cells:
        for k in range(nz):
            q = c[i, j, k, :, 0, :, :] if kind == "row" else c[i, j, k, 0, :, :, :]
            pts = [q[0, 0], q[1, 0], q[1, 1], q[0, 1]]
            polys.append([(p[axis], p[2]) for p in pts])
            live.append(bool(act[i, j, k]))
            if a is not None and a.bad[i, j, k]:
                colours.append((0.89, 0.29, 0.2, 0.75))
            elif not act[i, j, k]:
                colours.append("#ffffff")
            else:
                colours.append(palette[k % len(palette)])
    ax.add_collection(PolyCollection(polys, facecolors=colours, edgecolors="#333333",
                                     linewidths=0.5))
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)
    seen = set()
    for faces in faults:
        for f in faces:
            for p in face_pillars(f):
                if (kind == "row" and p[1] == n or kind == "col" and p[0] == n) and p not in seen:
                    seen.add(p)
                    x0, y0, z0, x1, y1, z1 = coord[p[1], p[0]]
                    ax.plot([(x0, y0)[axis], (x1, y1)[axis]], [z0, z1], color="black", lw=2.0)
    pts = np.array([p for poly, on in zip(polys, live) for p in poly if on] or
                   [p for poly in polys for p in poly])
    if window:                                           # depths of the cells in view only
        inside = (pts[:, 0] >= window[0]) & (pts[:, 0] <= window[1])
        pts = pts[inside] if inside.any() else pts
    ax.set_xlim(window or (np.nanmin(pts[:, 0]) - 5, np.nanmax(pts[:, 0]) + 5))
    ax.set_ylim(np.nanmax(pts[:, 1]) + 3, np.nanmin(pts[:, 1]) - 3)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(("x" if kind == "row" else "y") + " (m)", fontsize=8, labelpad=1)
    ax.set_ylabel("depth (m)", fontsize=8)
    ax.tick_params(labelsize=7)


def figure(name, title, before, after, faults, cut, a=None, steps=None, window=None):
    """One PNG, BEFORE and AFTER side by side. Top: the plan of the original mesh and of the
    deformed mesh, each column coloured by its worst cell (the deformed one with broken
    cells outlined red). Below: true-scale sections of ``before``, of each grid in ``steps``
    (a list of (label, grid) or (label, grid, audit); with an audit its broken cells are
    filled red) and of ``after``."""
    if OUT is None or not PNG:
        return
    import matplotlib.pyplot as plt

    a = a or audit(before, after)
    shown = [("BEFORE (original mesh)", before, None)] \
        + [(*st[:2], st[2] if len(st) > 2 else None) for st in (steps or [])] \
        + [("AFTER (deformed mesh)", after, a)]
    fig = plt.figure(figsize=(13.0, 5.6 + 2.3 * len(shown)), layout="constrained")
    grid_spec = fig.add_gridspec(1 + len(shown), 2, height_ratios=[5.0] + [2.3] * len(shown))
    legend = "\nfault (black) at the top (solid) and base (dashed)"
    plan(fig.add_subplot(grid_spec[0, 0]), before, audit(before, before), faults,
         "BEFORE, plan: worst cell per column" + legend)
    pc = plan(fig.add_subplot(grid_spec[0, 1]), after, a, faults,
              "AFTER, plan: worst cell per column" + legend + "; red outline = broken cell")
    fig.colorbar(pc, ax=fig.axes[:2], fraction=0.03, pad=0.02).set_label(
        "min 3-D scaled Jacobian (1 = perfect box, <= 0 = inverted)", fontsize=7)
    for m, (label, g, ga) in enumerate(shown):
        section(fig.add_subplot(grid_spec[1 + m, :]), g, cut, f"{label}: {cut[0]} {cut[1]}",
                ga, faults, window)
    fig.suptitle(title, fontsize=11)
    fig.savefig(OUT / f"32_{name}.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote 32_{name}.png")


def plan_length(grid, faces):
    nx, ny, _ = grid.cart_dims
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)
    return float(sum(np.linalg.norm(coord[p[1], p[0], :2] - coord[q[1], q[0], :2])
                     for p, q in map(face_pillars, faces)))


def adapted(grid, trace, quiet=False):
    res = adapt_grid(grid, trace, **MORPH)
    a = audit(grid, res.grid)
    assert a.n_bad == 0, a.line()
    if not quiet:
        print(f"  adapt  {len(res.faults_after)} faces; mean distance to the trace "
              f"{res.distance_before.mean():.1f} -> {res.distance_after.mean():.1f} m; {a.line()}")
    return res.grid, list(res.faults_after)


# ── A angle sweep ────────────────────────────────────────────────────────────────────────────


def case_a():
    print("A angle sweep: rim-to-rim fault on x = 300 m, hanging wall east, rigid "
          "(boundary='free')")
    grid, faces = cake(), x_fault(12)
    hw = side(grid, faces, EAST)
    shown = {}
    for angle in (89, 75, 60, 45, 35, 30):
        res, a = dip("A", grid, Dip(faces, angle, hanging_side=hw, name="F"), quiet=True,
                     boundary="free")
        dv = np.nanmax(np.abs(a.dv[12:]))
        assert dv < 1e-9, "a rigid hanging wall keeps its cell volumes"
        got = angles(res, "F")
        print(f"  {angle:2d} deg: COORD says {got.min():.6f}..{got.max():.6f}; fades x{res.widen:g},"
              f" footwall fade {res.fade['F'][1]:5.0f} m; hanging-wall volumes {dv:.0e}; "
              f"min 3-D Jacobian {a.min_after:.3f}; {a.n_bad} broken")
        shown[angle] = res.grid
    figure("A_angle_sweep", "A: one fault set to 89, 60, 45 and 30 deg (boundary='free': rigid hanging wall)",
           grid, shown[30], [faces], ("row", 6),
           steps=[(f"{d} deg", shown[d]) for d in (89, 60, 45)], window=(-5, 605))


# ── B curved trace ───────────────────────────────────────────────────────────────────────────


def case_b():
    print("B curved trace: an arc adapted into the pillars, dipped 60 deg, slid 20 m")
    grid = cake()
    r, cx, cy = 200.0, 480.0, 150.0
    th = np.radians(np.linspace(-53.0, 53.0, 41))
    arc = np.c_[cx - r * np.cos(th), cy + r * np.sin(th)]
    grid, faces = adapted(grid, arc)
    hw = (cx, cy)
    res, _ = dip("B", grid, Dip(faces, 60.0, hanging_side=hw, name="ARC"))
    out, _ = slide("B", res.grid, Throw(faces, 20.0, hanging_side=hw, name="ARC"))
    t = throw_of(res.grid, out.grid, faces, hw)
    assert np.allclose(t, 20.0, atol=1e-6), t
    print(f"  throw  {t.min():.3f}..{t.max():.3f} m along {plan_length(grid, faces):.0f} m of arc")
    figure("B_curved_arc", "B: curved trace (adapted), dip 60 deg, throw 20 m", grid,
           out.grid, [faces], ("row", 6), steps=[("dipped", res.grid)])


# ── C S-shaped trace, angle table ────────────────────────────────────────────────────────────


def case_c():
    print("C S-shaped trace: angle table 50 -> 80 deg and throw table 8 -> 20 m along the arc")
    grid = cake()
    y = np.linspace(-10.0, 310.0, 65)
    trace = np.c_[300.0 + 50.0 * np.sin(2 * np.pi * y / 300.0), y]
    grid, faces = adapted(grid, trace)
    length = plan_length(grid, faces)
    hw = (560.0, 150.0)
    res, _ = dip("C", grid, Dip(faces, [(0.0, 50.0), (length, 80.0)], hanging_side=hw,
                                name="S"))
    got = res.realized_dip("S")
    start = "north" if res.chains["S"].xy[0, 1] > res.chains["S"].xy[-1, 1] else "south"
    print(f"  arc length runs from the chain end with the smaller (I, J), here the {start} end: "
          f"dip {got[0, 1]:.1f} deg there, {got[len(got) // 2, 1]:.1f} mid-way, "
          f"{got[-1, 1]:.1f} at the other end")
    # the slide reads arc length along the top of the grid it slides: the dipped one, whose
    # end pillars also lean along strike where the fault meets the fixed rim
    top = plan_length(res.grid, faces)
    print(f"  the top trace is {top:.2f} m long after the dip ({length:.2f} m before)")
    out, a = slide("C", res.grid, Throw(faces, [(0.0, 8.0), (top, 20.0)], hanging_side=hw,
                                        name="S"))
    t = throw_of(res.grid, out.grid, faces, hw)
    assert np.allclose(sorted((t[0], t[-1])), (8.0, 20.0), atol=1e-6), t
    shift = (corners(out.grid) - corners(res.grid))[:, :, 0, 0, 0, 0, 2]
    east = shift[13:, 10]
    print(f"  throw  {t[0]:.2f} m -> {t[-1]:.2f} m along the fault; along row 10 east of it the "
          f"hanging wall drops {east.min():.1f}..{east.max():.1f} m, largest step between "
          f"neighbours {np.abs(np.diff(east)).max():.2f} m (no tear inside the bend)")
    assert np.abs(np.diff(east)).max() < 1.0
    figure("C_s_trace_angle_table", "C: S-shaped trace, dip 50 -> 80 deg and throw 8 -> 20 m "
           "along the arc", grid, out.grid, [faces], ("row", 10), steps=[("dipped", res.grid)])


# ── D stair-stepped diagonal ─────────────────────────────────────────────────────────────────


def case_d():
    print("D stair-stepped diagonal FAULTS (not adapted): dip 50 deg, slide 15 m")
    grid = cake()
    faces = fault_faces_for(grid, np.array([[150.0, -1.0], [450.0, 301.0]]))
    hw = (520.0, 150.0)
    res, _ = dip("D", grid, Dip(faces, 50.0, hanging_side=hw, name="D"))
    out, _ = slide("D", res.grid, Throw(faces, 15.0, hanging_side=hw, name="D"))
    print(f"  {len(faces)} faces in steps; the dip follows the smoothed strike, so a step's "
          f"pillars lean the same way as the trace's mean direction")
    figure("D_staircase", "D: stair-stepped diagonal, dip 50 deg, throw 15 m", grid, out.grid,
           [faces], ("row", 6), steps=[("dipped", res.grid)])


# ── E two tips ───────────────────────────────────────────────────────────────────────────────


def case_e():
    print("E fault with two tips: an angle per pillar 55 -> 80 -> 55 deg, throw split 50/50")
    grid = cake()
    faces = x_fault(12, 3, 10)
    hw = side(grid, faces, EAST)
    m = len(faces) + 1
    table = 55.0 + 25.0 * np.sin(np.linspace(0.0, np.pi, m))
    res, _ = dip("E", grid, Dip(faces, table, hanging_side=hw, name="E"))
    out, _ = slide("E", res.grid, Throw(faces, 15.0, hanging_side=hw, hanging=0.5, name="E"))
    t = throw_of(res.grid, out.grid, faces, hw)
    print(f"  {res.modes['E']} mode (a fault with tips cannot be rigid), fade "
          f"{res.fade['E'][0]:.0f} m; throw 0 at the tips, {t.max():.2f} m mid-fault, the "
          f"footwall rises half of it")
    figure("E_tips_angle_per_pillar", "E: two tips, 55 -> 80 -> 55 deg, throw 15 m split "
           "hanging/footwall", grid, out.grid, [faces], ("row", 6),
           steps=[("dipped", res.grid)])


# ── F graben and horst ───────────────────────────────────────────────────────────────────────


def crossing_below(grid, row, i):
    """Depth where pillars i and i + 1 of pillar row ``row`` cross, and how far that is below
    the deepest corner of the cells between them (negative: never, they diverge)."""
    nx, ny, _ = grid.cart_dims
    coord = np.asarray(grid.coord, float).reshape(ny + 1, nx + 1, 6)[row]
    lean = (coord[:, 3] - coord[:, 0]) / (coord[:, 5] - coord[:, 2])
    gap, closing = coord[i + 1, 0] - coord[i, 0], lean[i] - lean[i + 1]
    if closing <= 0:
        return np.inf, np.inf
    depth = coord[i, 2] + gap / closing
    return depth, depth - corners(grid)[i, min(row, ny - 1), -1, :, :, 1, 2].max()


def case_f():
    print("F graben and horst: faults on x = 100, 225, 350, 475 m; a graben between the first "
          "two, a horst between the last two")
    grid = cake()
    faults = {"F1": x_fault(4), "F2": x_fault(9), "F3": x_fault(14), "F4": x_fault(19)}
    sides = {"F1": EAST, "F2": WEST, "F3": WEST, "F4": EAST}
    angle = {"F1": 60.0, "F2": 60.0, "F3": 65.0, "F4": 65.0}
    throws = {"F1": 15.0, "F2": 15.0, "F3": 10.0, "F4": 10.0}
    hws = {n: side(grid, faults[n], sides[n]) for n in faults}
    slides = [Throw(n, throws[n], hanging_side=hws[n]) for n in faults]

    res, _ = dip("F", grid, [Dip(n, angle[n], hanging_side=hws[n]) for n in faults],
                 faults=faults)
    depth, room = crossing_below(res.grid, 0, 6)
    print(f"  pivot 'middle': the facing faults meet just below the grid; graben pillars 6 and "
          f"7 cross at {depth:.0f} m, {room:.1f} m below the cells")
    try:
        slide("F", res.grid, slides, faults=faults)
        raise AssertionError("the slide gate should refuse this")
    except SlideError as exc:
        print(f"  slide  refused: ...{str(exc).split(' · ')[-1]}")
    blind = slide_grid(res.grid, slides, faults=faults, min_jacobian=None, on_invalid="ignore")
    a = audit(res.grid, blind.grid)
    print(f"         without the cell gate it would pass with {a.line()}")
    assert a.inverted.any(), "the refusal is real: the unguarded slide inverts cells"
    refused = blind.grid
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            dip_grid(grid, [Dip(n, angle[n], hanging_side=hws[n]) for n in faults],
                     faults=faults, slide=slides)
            raise AssertionError("no widening can make room below the meeting faults")
        except SlideError as exc:
            assert "converge" in str(exc)
    print("  dip_grid(..., slide=...) widens the fades x16 and still refuses: the faults "
          "themselves meet within the throw, which no change of the cells can fix")

    angle.update(F1=75.0, F2=75.0)
    res, _ = dip("F", grid, [Dip(n, angle[n], hanging_side=hws[n]) for n in faults],
                 faults=faults)
    depth, room = crossing_below(res.grid, 0, 6)
    print(f"  graben faults at 75 deg instead: pillars 6 and 7 now cross {room:.0f} m below the "
          f"cells")
    out, _ = slide("F", res.grid, slides, faults=faults)
    for n in faults:
        t = throw_of(res.grid, out.grid, faults[n], hws[n])
        assert np.allclose(t, throws[n], atol=1e-6), (n, t)
    print("  throws hold exactly: F1 and F2 15 m (graben), F3 and F4 10 m (horst)")
    figure("F_graben_horst", "F: graben and horst; with 60 deg graben faults the slide is refused"
           " (inverted cells in red), with 75 deg it passes", grid, out.grid,
           list(faults.values()), ("row", 6),
           steps=[("60 deg graben, slid WITHOUT the cell gate", refused, a),
                  ("75 deg graben, dipped", res.grid)])


# ── G domino ─────────────────────────────────────────────────────────────────────────────────


def case_g():
    print("G domino: four parallel faults dipping east at 50, 55, 60 and 65 deg, 10 m each")
    grid = cake()
    faults = {f"D{n + 1}": x_fault(i) for n, i in enumerate((4, 9, 14, 19))}
    hws = {n: side(grid, f, EAST) for n, f in faults.items()}
    res, _ = dip("G", grid, [Dip(n, 50.0 + 5 * m, hanging_side=hws[n])
                             for m, n in enumerate(faults)], faults=faults)
    out, _ = slide("G", res.grid, [Throw(n, 10.0, hanging_side=hws[n]) for n in faults],
                   faults=faults)
    base = corners(out.grid)[:, 6, 0, 0, 0, 0, 2] - corners(grid)[:, 6, 0, 0, 0, 0, 2]
    print(f"  the blocks step down 0, 10, 20, 30, 40 m: measured "
          f"{', '.join(f'{base[i]:.0f}' for i in (0, 5, 10, 15, 20))} m")
    assert np.allclose(base[[0, 5, 10, 15, 20]], [0, 10, 20, 30, 40], atol=1e-6)
    figure("G_domino", "G: domino, four faults at 50, 55, 60 and 65 deg, 10 m each", grid,
           out.grid, list(faults.values()), ("row", 6), steps=[("dipped", res.grid)])


# ── H T-junction ─────────────────────────────────────────────────────────────────────────────


def case_h():
    print("H T-junction: a north-dipping fault ends on an east-dipping one")
    grid = cake()
    main, branch = x_fault(12), y_fault(6, 1, 12)
    faults = {"MAIN": main, "BRANCH": branch}
    hw_main, hw_branch = side(grid, main, EAST), (150.0, 250.0)
    res, _ = dip("H", grid, [Dip("MAIN", 65.0, hanging_side=hw_main),
                             Dip("BRANCH", 70.0, hanging_side=hw_branch)], faults=faults)
    nx, ny, _ = grid.cart_dims
    x0, y0, z0, x1, y1, z1 = np.asarray(res.grid.coord, float).reshape(ny + 1, nx + 1, 6)[6, 12]
    east, north = np.degrees(np.arctan2(z1 - z0, x1 - x0)), np.degrees(np.arctan2(z1 - z0, y1 - y0))
    print(f"  the shared pillar (12, 6) leans along the line where the two planes meet: "
          f"{east:.4f} deg seen across MAIN, {north:.4f} deg across BRANCH")
    assert abs(east - 65.0) < 1e-6 and abs(north - 70.0) < 1e-6
    out, _ = slide("H", res.grid, [Throw("MAIN", 12.0, hanging_side=hw_main),
                                   Throw("BRANCH", 6.0, hanging_side=hw_branch)],
                   faults=faults)
    figure("H_t_junction", "H: T-junction, MAIN 65 deg east, BRANCH 70 deg north, throws 12 "
           "and 6 m", grid, out.grid, list(faults.values()), ("col", 6),
           steps=[("dipped", res.grid)])


# ── I bedding dipping along strike ───────────────────────────────────────────────────────────


def case_i():
    print("I bedding dipping 22 deg along strike: the rigid hanging wall keeps every volume")
    grid = cake(top=lambda x, y: TOP + 0.05 * x + 0.4 * y)
    faces = x_fault(12)
    hw = side(grid, faces, EAST)
    res, a = dip("I", grid, Dip(faces, 45.0, hanging_side=hw, decay=None, name="F"),
                 boundary="free")
    dv = np.nanmax(np.abs(a.dv[12:]))
    print(f"  hanging-wall cell volumes change by at most {dv:.1e} (the pivot depth varies 120 m "
          f"along the fault; the fit carries it off the fault unchanged along the normal)")
    assert dv < 1e-9
    figure("I_bedding_along_strike", "I: bedding dips 22 deg along strike; fault 45 deg, rigid (boundary='free')",
           grid, res.grid, [faces], ("row", 9), a=a)


# ── J pre-tilted pillars ─────────────────────────────────────────────────────────────────────


def case_j():
    print("J pre-tilted pillars: every pillar leans 35 m east; the fault is turned to dip west")
    grid = cake(lean=lambda x, y: (35.0, 0.0))
    faces = x_fault(12)
    before = coord_dip(grid, (12, 6), (1.0, 0.0))
    hw = side(grid, faces, WEST)
    res, _ = dip("J", grid, Dip(faces, 60.0, hanging_side=hw, name="F"))
    print(f"  the fault dipped {before:.1f} deg east before; now "
          f"{coord_dip(res.grid, (12, 6), (-1.0, 0.0)):.6f} deg west")
    out, _ = slide("J", res.grid, Throw(faces, 12.0, hanging_side=hw, name="F"))
    figure("J_pretilted_flip", "J: pillars leaning east, fault turned to dip 60 deg west, "
           "throw 12 m", grid, out.grid, [faces], ("row", 6), steps=[("dipped", res.grid)])


# ── K re-dip ─────────────────────────────────────────────────────────────────────────────────


def case_k():
    print("K re-dip: dip 60, slide 15, re-dip 75, slide 10, re-dip 50")
    grid0 = cake()
    faces = x_fault(12)
    hw = side(grid0, faces, EAST)
    steps, grid = [], grid0
    for what, value in (("dip", 60.0), ("slide", 15.0), ("dip", 75.0), ("slide", 10.0),
                        ("dip", 50.0)):
        if what == "dip":
            res, _ = dip("K", grid, Dip(faces, value, hanging_side=hw, name="F"))
        else:
            res, _ = slide("K", grid, Throw(faces, value, hanging_side=hw, name="F"))
        grid = res.grid
        steps.append((f"{what} {value:g}", grid))
    t = throw_of(grid0, grid, faces, hw)
    print(f"  total throw {t.min():.6f}..{t.max():.6f} m (a dip moves pillars, never depths)")
    assert np.allclose(t, 25.0, atol=1e-6)
    figure("K_redip_sequence", "K: dip 60, slide 15, re-dip 75, slide 10, re-dip 50", grid0,
           grid, [faces], ("row", 6), steps=steps[:-1])


# ── L thin layers, pinch-outs, holes ─────────────────────────────────────────────────────────


def stacked(horizons, actnum=None):
    """A grid on vertical pillars whose horizons are ``horizons[h](x, y)`` (depths)."""
    nz = len(horizons) - 1
    grid = cake(dz=[1.0] * nz, actnum=actnum)
    nx, ny, _ = grid.cart_dims
    px, py = np.meshgrid(DX * np.arange(nx + 1), DX * np.arange(ny + 1), indexing="ij")
    depth = np.array([np.broadcast_to(h(px, py), px.shape) for h in horizons], float)
    depth = np.maximum.accumulate(depth, axis=0)                 # horizons never cross
    slot = (np.arange(2 * nx) + 1) // 2, (np.arange(2 * ny) + 1) // 2
    z = depth[(np.arange(2 * nz) + 1) // 2][:, slot[0]][:, :, slot[1]]
    coord = np.asarray(grid.coord, float).reshape(-1, 6).copy()
    coord[:, 2], coord[:, 5] = depth.min(), depth.max()
    return CornerPointGrid(grid.cart_dims, coord.ravel(),
                           np.ascontiguousarray(z.transpose(0, 2, 1)).ravel(), grid.actnum)


def case_l():
    print("L thin and pinched-out layers, a zero-thickness layer and 15 % inactive cells")
    rng = np.random.default_rng(11)
    act = (rng.random((NX, NY, 6)) > 0.15).astype(int)
    top = lambda x, y: TOP + 0.05 * x                                   # noqa: E731
    h = [top,
         lambda x, y: top(x, y) + 0.5,                                  # 0.5 m
         lambda x, y: top(x, y) + 0.5 + 20.0 * np.clip(1 - x / 350.0, 0, 1),   # pinches out
         None, None, None, None]
    h[3] = lambda x, y: h[2](x, y)                                     # zero thickness
    h[4] = lambda x, y: h[3](x, y) + 30.0
    h[5] = lambda x, y: h[4](x, y) + 2.0
    h[6] = lambda x, y: h[5](x, y) + 25.0
    grid = stacked(h, act)
    faces = x_fault(12)
    hw = side(grid, faces, EAST)
    res, _ = dip("L", grid, Dip(faces, 55.0, hanging_side=hw, name="F"))
    out, _ = slide("L", res.grid, Throw(faces, 10.0, hanging_side=hw, name="F"))
    vol = hex_volumes(corners(grid))
    print(f"  {int((vol < 1e-9).sum())} cells of zero volume (the pinch-out and the empty layer) "
          f"and {int((act == 0).sum())} inactive cells: none judged, none broken")
    figure("L_pinchout_holes", "L: 0.5 m layer, pinch-out across the fault, empty layer, holes;"
           " 55 deg, 10 m", grid, out.grid, [faces], ("row", 6), steps=[("dipped", res.grid)])


# ── M skewed, uneven columns ─────────────────────────────────────────────────────────────────


def warped(grid, fn):
    coord = np.asarray(grid.coord, float).reshape(-1, 6).copy()
    for a in (0, 3):
        coord[:, a], coord[:, a + 1] = fn(coord[:, a], coord[:, a + 1])
    return CornerPointGrid(grid.cart_dims, coord.ravel(), grid.zcorn, grid.actnum)


def case_m():
    print("M skewed, uneven columns: x' = 600 (x/600)^1.6 + 0.4 y (5 to 39 m wide, sheared)")
    grid = warped(cake(), lambda x, y: (600.0 * (x / 600.0) ** 1.6 + 0.4 * y, y))
    widths = np.diff(600.0 * (DX * np.arange(NX + 1) / 600.0) ** 1.6)
    faces = x_fault(12)
    res, _ = dip("M", grid, Dip(faces, 55.0, hanging_side=(600.0, 150.0), name="F"))
    out, _ = slide("M", res.grid, Throw(faces, 12.0, hanging_side=(600.0, 150.0), name="F"))
    print(f"  columns {widths.min():.1f} to {widths.max():.1f} m wide; the fault strikes "
          f"{np.degrees(np.arctan(0.4)):.1f} deg off north and dips along its own normal")
    figure("M_skewed_columns", "M: sheared, uneven columns; 55 deg, 12 m", grid, out.grid,
           [faces], ("row", 6), steps=[("dipped", res.grid)])


# ── N dome and ring fault ────────────────────────────────────────────────────────────────────


def case_n():
    print("N dome with a ring fault: dipped 70 deg inward, the inside collapses 25 m (caldera)")
    grid = cake(top=lambda x, y: TOP + 0.03 * x
                - 45.0 * np.exp(-((x - 300.0) ** 2 + (y - 150.0) ** 2) / (2 * 110.0 ** 2)))
    th = np.linspace(0.0, 2.0 * np.pi, 33)
    ring = np.c_[300.0 + 170.0 * np.cos(th), 150.0 + 95.0 * np.sin(th)]
    grid, faces = adapted(grid, ring)
    centre = (300.0, 150.0)
    res, _ = dip("N", grid, Dip(faces, 70.0, hanging_side=centre, name="RING"))
    out, a = slide("N", res.grid, Throw(faces, 25.0, hanging_side=centre, name="RING"))
    t = throw_of(res.grid, out.grid, faces, centre)
    assert np.allclose(t, 25.0, atol=1e-6)
    print(f"  {res.modes['RING']} inside the ring; throw {t.min():.3f}..{t.max():.3f} m all round")
    figure("N_dome_ring_caldera", "N: dome, ring fault dipping 70 deg inward, 25 m collapse",
           grid, out.grid, [faces], ("row", 6), steps=[("dipped", res.grid)])


# ── O shallow angle ──────────────────────────────────────────────────────────────────────────


def case_o():
    print("O shallow angles on bedding dipping 14 deg toward the hanging wall")
    grid = cake(top=lambda x, y: TOP + 0.25 * x)
    faces = x_fault(12)
    hw = side(grid, faces, EAST)
    spec = Dip(faces, 15.0, hanging_side=hw, name="F")
    shift = 50.0 / np.tan(np.radians(15.0))
    try:
        dip_grid(grid, spec)
        raise AssertionError("15 deg against a fixed east side should be refused")
    except DipError:
        print(f"  15 deg, default boundary='fixed': refused even widened x16: the pillars must "
              f"shift +-{shift:.0f} m at the top and base, and the east side, parallel to the "
              f"fault and 300 m from it, holds still")
    held, _ = dip("O", grid, Dip(faces, 30.0, hanging_side=hw, name="F"), quiet=True)
    print(f"  30 deg, default boundary='fixed': fades widened x{held.widen:g}, every cell valid, "
          f"the outline held")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res, _ = dip("O", grid, spec, warn=True, boundary="free")
    print(f"  15 deg, boundary='free': fades widened x{res.widen:g} ({len(caught)} warning), "
          f"footwall fade {res.fade['F'][1]:.0f} m; every cell valid")
    assert res.widen > 1.0
    try:
        dip_grid(grid, spec, widen=0, boundary="free")
        raise AssertionError("widen=0 should be refused")
    except DipError as exc:
        print(f"  15 deg, free, widen=0: DipError ({str(exc).split(':')[0]})")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        part, _ = dip("O", grid, spec, widen=0, on_invalid="bisect", quiet=True, boundary="free")
    got = angles(part, "F")
    print(f"  15 deg, free, widen=0, on_invalid='bisect': {100 * part.scale:.1f} % of the change; "
          f"the fault dips {got.mean():.1f} deg, NOT 15 (scale {part.scale:.3f} says so)")
    assert part.scale < 1.0 and got.mean() > 16.0
    try:
        dip_grid(grid, Dip(faces, 8.0, hanging_side=hw, name="F"), boundary="free")
        raise AssertionError("8 deg should be refused")
    except DipError:
        print("  8 deg, free: refused even with the fades widened x16; the cells cannot take it")
    figure("O_shallow_widen_bisect", "O: bedding dipping 14 deg: 30 deg with the rim fixed; 15 deg "
           "with boundary='free' (after) vs no widening, bisected", grid, res.grid, [faces],
           ("row", 6), steps=[("30 deg, rim fixed (default)", held.grid),
                              (f"15 deg, free, widen=0, bisect: {got.mean():.1f} deg", part.grid)])


# ── P random angles ──────────────────────────────────────────────────────────────────────────


def base_reversals(res, name):
    """Fault segments whose base runs backwards along strike: the fault folds onto itself."""
    nx, ny, _ = res.grid.cart_dims
    coord = np.asarray(res.grid.coord, float).reshape(ny + 1, nx + 1, 6)
    line = np.array([coord[j, i] for i, j in res.chains[name].ij])
    top, base = np.diff(line[:, :2], axis=0), np.diff(line[:, 3:5], axis=0)
    return int(((top * base).sum(1) < 0).sum()), len(top)


def case_p():
    print("P random angles 40..85 deg along a curved trace: noise at every pillar vs a smooth draw")
    grid = cake()
    r, cx, cy = 260.0, 40.0, 150.0
    th = np.radians(np.linspace(-45.0, 45.0, 41))
    grid, faces = adapted(grid, np.c_[cx + r * np.cos(th), cy + r * np.sin(th)])
    hw = (560.0, 150.0)
    rng = np.random.default_rng(5)
    noisy = Dip(faces, rng.uniform(40.0, 85.0, len(faces) + 1), hanging_side=hw, name="P")
    try:
        dip_grid(grid, noisy)
        raise AssertionError("a fault twisted pillar by pillar should be refused")
    except DipError:
        pass
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw = dip_grid(grid, noisy, on_invalid="ignore")
        cut = dip_grid(grid, noisy, on_invalid="bisect")
    back, n = base_reversals(raw, "P")
    print(f"  noise at each pillar: REFUSED even widened x16. The fault itself is impossible: its "
          f"base runs backwards along strike on {back} of {n} segments (it folds onto itself). "
          f"bisect manages {100 * cut.scale:.0f} % of the change")
    assert back > 0
    length = plan_length(grid, faces)
    ctrl = np.c_[np.linspace(0.0, length, 5), rng.uniform(40.0, 85.0, 5)]
    res, _ = dip("P", grid, Dip(faces, ctrl, hanging_side=hw, name="P"))
    assert base_reversals(res, "P")[0] == 0
    print("  5 random control angles " + ", ".join(f"{v:.0f}" for v in ctrl[:, 1])
          + " deg, interpolated along strike: exact at every pillar, no widening")
    out, _ = slide("P", res.grid, Throw(faces, 12.0, hanging_side=hw, name="P"))
    figure("P_random_angles", "P: random 40..85 deg, smooth along strike (after) vs noise at "
           "every pillar (bisected)", grid, out.grid, [faces], ("row", 6),
           steps=[(f"noise per pillar, bisect {100 * cut.scale:.0f} %", cut.grid)])


# ── Q the slide's cell gate ──────────────────────────────────────────────────────────────────


def case_q():
    print("Q the slide's cell gate, and the cure: 18 m on a short 40 deg fault in 5 m layers")
    grid = layer_cake(12, 12, 4, 10.0, 5.0)
    tip = {"T": [FaultFace(6, j, "X") for j in range(4, 10)]}
    spec = Dip("T", 40.0, hanging_side="right")
    res = dip_grid(grid, spec, faults=tip)
    assert audit(grid, res.grid).n_bad == 0
    throw = Throw("T", 18.0, hanging_side="right")
    blind = slide_grid(res.grid, throw, faults=tip, min_jacobian=None)
    unchecked = audit(res.grid, blind.grid)
    print(f"  dip, then slide with min_jacobian=None (no cell gate): ok={blind.ok}, yet the "
          f"{unchecked.line()}")
    assert blind.ok and unchecked.inverted.any()
    try:
        slide_grid(res.grid, throw, faults=tip)
        raise AssertionError("the cell gate should refuse")
    except SlideError as exc:
        print(f"  dip, then slide (default): SlideError ...{str(exc).split(' · ')[-1]}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        part = slide_grid(res.grid, throw, faults=tip, on_invalid="bisect")
    print(f"  ... on_invalid='bisect': {100 * part.scale:.1f} % of the throw "
          f"({18.0 * part.scale:.1f} m) is the most these cells take")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        both = dip_grid(grid, spec, faults=tip, slide=throw)
    a = audit(grid, both.slide.grid)
    assert a.n_bad == 0 and both.slide.scale == 1.0
    assert np.allclose(angles(both, "T"), 40.0, atol=1e-6)
    dv = 100.0 * np.nanmax(np.abs(a.dv))
    print(f"  dip_grid(..., slide=throw): the fades widen x{both.widen:g} (to "
          f"{both.fade['T'][0]:.0f} m) so the pillars beside the fault lean with it; the full "
          f"18 m at exactly 40 deg, {a.line()}; cells change volume by up to {dv:.0f} %")
    figure("Q_slide_cell_gate", "Q: 18 m on a short 40 deg fault: slid without the cell gate "
           "(broken cells red) vs dip_grid(..., slide=...)", grid, both.slide.grid,
           tip.values(), ("row", 6), a=a,
           steps=[("dip, then 18 m WITHOUT the cell gate", blind.grid, unchecked),
                  (f"dip_grid(..., slide=...): fades x{both.widen:g}", both.grid)])


# ── R a conservative refusal ─────────────────────────────────────────────────────────────────


def case_r():
    print("R a conservative refusal: bedding rising 35 deg toward a fault dipping 30 deg")
    grid = cake(top=lambda x, y: TOP + 200.0 - 0.7 * x)
    faces = x_fault(12)
    spec = Dip(faces, 30.0, hanging_side=side(grid, faces, EAST), name="F")
    try:
        dip_grid(grid, spec)
        raise AssertionError("expected the gate to refuse")
    except DipError as exc:
        print(f"  dip_grid refuses: {str(exc)[str(exc).index('INVALID'):][:80]}...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        loose = dip_grid(grid, spec, on_invalid="ignore")
    a = audit(grid, loose.grid)
    print(f"  the refused grid: {a.line()}")
    print("  the plan view of a cell's top flips where bedding is steeper than the fault, but "
          "the 3-D cells stay sound: the gate is conservative here, not wrong")
    assert int(a.inverted.sum()) == 0
    figure("R_conservative_refusal", "R: refused, yet the independent 3-D check finds the cells "
           "sound (bedding steeper than the fault)", grid, loose.grid, [faces], ("row", 6), a=a)


# ── S growth fault ───────────────────────────────────────────────────────────────────────────


def case_s():
    print("S growth fault: dipped 55 deg, throw 10 m at the top growing to 16 m at the base")
    grid = cake()
    faces = x_fault(12)
    hw = side(grid, faces, EAST)
    res, _ = dip("S", grid, Dip(faces, 55.0, hanging_side=hw, name="G"))
    growth = [1.0, 1.15, 1.25, 1.35, 1.5, 1.6]
    out, _ = slide("S", res.grid, Throw(faces, 10.0, hanging_side=hw, growth=growth, name="G"))
    def thickness(g):                       # along the pillars, straight from ZCORN
        z = corners(g)[16:, ..., 2]
        return (z[..., 1] - z[..., 0]).mean(axis=(3, 4))

    got = 100.0 * (thickness(out.grid) / thickness(res.grid) - 1.0).mean(axis=(0, 1))
    want = 100.0 * 10.0 * np.diff(growth) / np.asarray(DZ)
    print("  hanging-wall layers thicken by " + ", ".join(f"{g:.1f}" for g in got)
          + " % (10 m x the growth step / layer thickness: "
          + ", ".join(f"{w:.1f}" for w in want) + ")")
    assert np.allclose(got, want, atol=0.05)
    figure("S_growth_fault", "S: growth fault, 55 deg, throw 10 m at the top to 16 m at the base",
           grid, out.grid, [faces], ("row", 6), steps=[("dipped", res.grid)])


# ── T blind fault ────────────────────────────────────────────────────────────────────────────


def case_t():
    print("T blind fault: dipped 55 deg, cuts only layers 3..5; layers 1..2 drape over it")
    grid = cake()
    faces = x_fault(12)
    hw = side(grid, faces, EAST)
    res, _ = dip("T", grid, Dip(faces, 55.0, hanging_side=hw, name="B"))
    out, _ = slide("T", res.grid, Throw(faces, 12.0, hanging_side=hw, layers=(3, 5), name="B"))
    c = corners(out.grid)[11:13, 6]          # the two columns either side of the fault
    gap = [float(c[1, k, 0, 0, 0, 2] - c[0, k, 1, 0, 0, 2]) for k in range(5)]
    gap.append(float(c[1, 4, 0, 0, 1, 2] - c[0, 4, 1, 0, 1, 2]))
    print("  offset at horizons 0..5: " + ", ".join(f"{g:.1f}" for g in gap) + " m: none above"
          " the blind tip (top of layer 3), rising on a cosine to the full 12 m at the base")
    assert np.allclose(gap[:3], 0.0, atol=1e-9) and abs(gap[5] - 12.0) < 1e-6
    assert np.all(np.diff(gap) >= -1e-9)
    figure("T_blind_fault", "T: blind fault, 55 deg, 12 m in layers 3..5 only", grid, out.grid,
           [faces], ("row", 6), steps=[("dipped", res.grid)])


# ── U seeded fuzz ────────────────────────────────────────────────────────────────────────────


def random_case(rng):
    """A random layer cake (tilted, domed, leaning, holes) and a random fault on it."""
    nx, ny, nz = (int(v) for v in (rng.integers(8, 18), rng.integers(6, 14), rng.integers(2, 7)))
    dx, dy = rng.uniform(20.0, 80.0, 2)
    dz = rng.uniform(0.5, 30.0, nz)
    sx, sy = np.tan(np.radians(rng.uniform(-25.0, 25.0, 2)))
    bump, cx, cy = rng.uniform(0.0, 40.0), nx * dx / 2, ny * dy / 2

    def top(x, y):
        return TOP + sx * x + sy * y - bump * np.exp(-((x - cx) ** 2 + (y - cy) ** 2)
                                                     / (0.3 * cx * cy + 1))
    lx, ly = rng.uniform(-0.3, 0.3, 2) * dz.sum()
    lean = (lambda x, y: (lx, ly)) if rng.random() < 0.4 else None
    act = (rng.random((nx, ny, nz)) > 0.15).astype(int) if rng.random() < 0.3 else None
    grid = layer_cake(nx, ny, nz, float(dx), dz, dy=float(dy), top=top, lean=lean, actnum=act)
    kind = str(rng.choice(["x", "y", "tip", "curve"]))
    if kind == "x":
        i = int(rng.integers(3, nx - 2))
        faces = [FaultFace(i, j, "X") for j in range(1, ny + 1)]
    elif kind == "y":
        j = int(rng.integers(3, ny - 2))
        faces = [FaultFace(i, j, "Y") for i in range(1, nx + 1)]
    elif kind == "tip":
        i, j0 = int(rng.integers(3, nx - 2)), int(rng.integers(1, ny // 2))
        faces = [FaultFace(i, j, "X") for j in range(j0, int(rng.integers(j0 + 2, ny)) + 1)]
    else:
        bend = rng.uniform(-0.25, 0.25) * nx * dx
        t = np.linspace(0.0, 1.0, 9)
        trace = np.c_[rng.uniform(0.35, 0.65) * nx * dx + bend * np.sin(np.pi * t),
                      -1.0 + t * (ny * dy + 2.0)]
        grid, faces = adapted(grid, trace, quiet=True)
    return grid, faces, kind


def judged(stats, before, call):
    """Run ``call`` (it raises on refusal); re-judge the verdict with the checker."""
    try:
        res = call(on_invalid="raise")
    except (DipError, SlideError):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            loose = call(on_invalid="ignore")
        a = audit(before, loose.grid)
        stats["refused" if a.n_bad else "conservative"] += 1
        return None
    a = audit(before, res.grid)
    stats["accepted"] += 1
    stats["disagree"] += a.n_bad > 0
    stats["points"].append((a.min_after, a.n_bad))
    return res


def case_u():
    print("U seeded fuzz: 2 seeds x 20 random grids; dip -> slide -> re-dip, every verdict re-judged")
    keys = ("accepted", "refused", "conservative", "disagree")
    total = {k: 0 for k in keys}
    dots = []
    for seed in (2024, 7):
        rng = np.random.default_rng(seed)
        stats = {k: 0 for k in keys} | {"points": []}
        for _ in range(20):
            grid, faces, kind = random_case(rng)
            angle, again = rng.uniform(25.0, 89.0, 2)
            pivot = str(rng.choice(["top", "middle", "bottom"]))
            side_ = str(rng.choice(["left", "right"]))
            nz = grid.cart_dims[2]
            faults = {"F": faces}
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                d1 = judged(stats, grid, functools.partial(
                    dip_grid, grid, Dip("F", angle, hanging_side=side_, pivot=pivot),
                    faults=faults))
                if d1 is None:
                    dots.append((angle, np.nan, kind))
                    continue
                dots.append((angle, stats["points"][-1][0], kind))
                throw = float(rng.uniform(0.05, 0.4)) * float(np.ptp(grid.zcorn)) / max(nz, 1)
                tip = kind == "tip"
                spec = Throw("F", throw, hanging_side=side_, decay="auto" if tip else None)
                s1 = judged(stats, d1.grid, functools.partial(
                    slide_grid, d1.grid, spec, faults=faults))
                if s1 is None:
                    continue
                judged(stats, s1.grid, functools.partial(
                    dip_grid, s1.grid, Dip("F", again, hanging_side=side_, pivot=pivot),
                    faults=faults))
        print(f"  seed {seed:4d}: " + ", ".join(f"{stats[k]} {k}" for k in keys))
        for k in keys:
            total[k] += stats[k]
    print(f"  all: {total['accepted']} steps accepted, every one passed by the checker; "
          f"{total['refused']} refused and the checker agrees they break cells; "
          f"{total['conservative']} refused where the checker would have passed (conservative)")
    assert total["disagree"] == 0, "the gate accepted a grid the checker breaks"
    if OUT is None or not PNG:
        return
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    for kind, mark in (("x", "o"), ("y", "s"), ("tip", "^"), ("curve", "D")):
        pts = np.array([(a, j) for a, j, k in dots if k == kind and np.isfinite(j)]).reshape(-1, 2)
        ax.scatter(pts[:, 0], pts[:, 1], marker=mark, label=f"{kind} fault, accepted")
        no = [a for a, j, k in dots if k == kind and not np.isfinite(j)]
        ax.scatter(no, np.zeros(len(no)), marker=mark, facecolors="none", edgecolors="red")
    ax.axhline(0.05, color="red", lw=0.8, ls="--")
    ax.text(26, 0.06, "floor 0.05", color="red", fontsize=8)
    ax.set_xlabel("requested dip (deg)")
    ax.set_ylabel("min 3-D scaled Jacobian after the dip")
    ax.set_title(f"U: random grids and faults, dip -> slide -> re-dip "
                 f"({total['accepted']} accepted, 0 broken; red rings at 0: first dip refused)",
                 fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "32_U_seeded_fuzz.png", dpi=110)
    plt.close(fig)
    print("  wrote 32_U_seeded_fuzz.png")


# ── V a field deck ───────────────────────────────────────────────────────────────────────────


def norne_deck(out):
    grid = Path(os.environ.get("CPG_NORNE_GRDECL", NORNE))
    files = [grid, grid.parent / "ACTNUM_0704.prop",
             grid.parent.parent / "FAULT" / "FAULT_JUN_05.INC"]
    if not all(f.is_file() for f in files):
        return None
    deck = out / "NORNE_GRID.GRDECL"
    deck.write_text("".join(f"INCLUDE\n  '{f.as_posix()}' /\n\n" for f in files))
    return deck


def case_v():
    print("V Norne: E_01 and BC set to 60 deg, every verdict re-checked")
    with TemporaryDirectory() as work:
        deck = norne_deck(Path(work))
        if deck is None:
            print("  skipped: Norne is not on this machine (set CPG_NORNE_GRDECL)")
            return
        grid = read_grdecl(deck)
        faults = _deck_faults(deck, grid.cart_dims[2])[0]
    nx, ny, nz = grid.cart_dims
    print(f"  {nx}x{ny}x{nz}, {int(np.asarray(grid.actnum).sum()):,} active cells, "
          f"{len(faults)} named faults")
    rows, shown = [], None
    for name, angle in (("E_01", 60.0), ("BC", 60.0)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                res = dip_grid(grid, Dip(name, angle), faults=faults, widen=2)
            except DipError as exc:
                loose = dip_grid(grid, Dip(name, angle), faults=faults, on_invalid="ignore")
                a = audit(grid, loose.grid)
                why = "the checker agrees" if a.n_bad else "the checker would pass it"
                rows.append(f"  {name:5s} {angle:4.0f} deg: refused "
                            f"({str(exc).split(': ')[-1].split(';')[0]}); {why}: {a.line()}")
                continue
        a = audit(grid, res.grid)
        assert a.n_bad == 0, f"{name} {angle}: {a.line()}"
        got = res.realized_dip(name)[:, 1]
        assert np.allclose(got, angle, atol=1e-6)
        rows.append(f"  {name:5s} {angle:4.0f} deg: exact on {len(got)} pillars, "
                    f"{res.modes[name]}, fades x{res.widen:g}; {a.line()}")
        if shown is None:
            shown = (name, angle, res, a)
    print("\n".join(rows))
    assert shown is not None, "no Norne fault was accepted"
    name, angle, res, a = shown
    cols = [ij for ij in res.chains[name].ij]
    mid = cols[len(cols) // 2]
    centre = float(np.asarray(grid.coord, float).reshape(ny + 1, -1, 6)[mid[1], mid[0], 0])
    figure("V_norne", f"V: Norne, fault {name} set to {angle:.0f} deg (sections zoomed on it)",
           grid, res.grid, [faults[name]], ("row", int(min(mid[1], ny - 1))), a=a,
           window=(centre - 700.0, centre + 700.0))


CASES = {}


def main(out, only=None):
    global OUT
    OUT = out
    for key, fn in sorted(CASES.items()):
        if only is None or key in only:
            fn()
    print("every case: the gate's verdict re-judged by an independent cell checker")


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
