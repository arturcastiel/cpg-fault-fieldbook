"""Fault-slide helpers shared by lessons 23-30.

The slide itself is :func:`fault_projections.slide_grid` (``slide_field`` and ``Throw`` here
only keep the lessons' names). What stays in this module is lesson scaffolding and the
independent checks the lessons use to verify the library:

- ``layer_cake`` and ``write_deck`` build and write synthetic cases;
- ``blocks``, ``unwelded_faces``, ``realized_throw`` and ``analytic_contacts`` recompute
  what a slide should give, without the library's code;
- ``process`` re-processes a slid grid with cpg2unstructured and reads its contacts.

Conventions: depths grow downward, so a positive throw moves the hanging wall down. Columns
and pillars are indexed ``[i, j]`` (0-based), as in ``pillar_lattice``. A slide only adds a
depth shift to ZCORN, so every corner moves along its own pillar and COORD never changes.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from cpg2unstructured import (CornerPointGrid, build_topology, from_grdecl, pillar_lattice,
                              read_grdecl, write_grid_include)
from cpg2unstructured.geometry import cell_geometry, face_geometry
from fault_projections.faults import FaultFace, faults_text
from fault_projections.from_faults import edge_face, face_edge, trace_faults
from fault_projections.slide import Throw as FpThrow
from fault_projections.slide import slide_grid

__all__ = [
    "layer_cake", "corner_cube", "with_cube", "faces_of", "face_columns",
    "blocks", "slide_blocks", "Throw", "slide_field", "unwelded_faces", "analytic_contacts",
    "Processed", "process", "crossed_faces", "contacts_by_layer", "write_deck", "realized_throw",
    "column_of_point", "cosine", "compare", "interior", "close_ring", "named_chains",
]


# ── synthetic grids ──────────────────────────────────────────────────────────────────────────


def layer_cake(nx, ny, nz, dx, dz, *, dy=None, top=None, origin=(0.0, 0.0), lean=None,
               actnum=None) -> CornerPointGrid:
    """An nx*ny*nz grid of stacked layers on straight pillars.

    ``dz`` is one thickness or one per layer. ``top(x, y)`` gives the depth of the top
    surface (default 0). ``lean(x, y)`` gives the ``(dx, dy)`` offset of a pillar's base
    relative to its top (default vertical). ``actnum`` is an optional ``[i, j, k]`` mask.
    """
    dy = dx if dy is None else dy
    thick = np.broadcast_to(np.asarray(dz, float), (nz,))
    x = origin[0] + dx * np.arange(nx + 1)
    y = origin[1] + dy * np.arange(ny + 1)
    px, py = np.meshgrid(x, y, indexing="ij")                       # (nx+1, ny+1)
    h0 = np.zeros_like(px) if top is None else np.vectorize(top)(px, py).astype(float)
    hz = h0[None] + np.concatenate([[0.0], np.cumsum(thick)])[:, None, None]
    z_top, z_bot = float(hz.min()), float(hz.max())
    off = np.zeros(px.shape + (2,))
    if lean is not None:
        off = np.array([[lean(a, b) for a, b in zip(ra, rb)] for ra, rb in zip(px, py)], float)
    # Pillars run from z_top to z_bot; the given xy is where the pillar crosses depth h0.
    frac0 = (h0 - z_top) / (z_bot - z_top)
    top_xy = np.stack([px, py], -1) - off * frac0[..., None]
    bot_xy = top_xy + off
    coord = np.concatenate([top_xy, np.full(px.shape + (1,), z_top), bot_xy,
                            np.full(px.shape + (1,), z_bot)], -1)
    coord = coord.transpose(1, 0, 2).reshape(-1, 6)                 # i fastest
    zc = hz[_slot_horizon(nz)][:, _slot_pillar(nx)][:, :, _slot_pillar(ny)]
    act = np.ones(nx * ny * nz, np.int8) if actnum is None else \
        np.asarray(actnum, np.int8).transpose(2, 1, 0).ravel()
    return CornerPointGrid((nx, ny, nz), coord, np.ascontiguousarray(
        zc.transpose(0, 2, 1)).ravel(), act)


def _slot_pillar(n):
    s = np.arange(2 * n)
    return s // 2 + s % 2


def _slot_column(n):
    return np.arange(2 * n) // 2


def _slot_horizon(nz):
    s = np.arange(2 * nz)
    return s // 2 + s % 2


def corner_cube(grid):
    """ZCORN as ``z[k-slot, j-slot, i-slot]``; slot 2c / 2c+1 = near / far corner of cell c."""
    nx, ny, nz = grid.cart_dims
    return grid.zcorn.reshape(2 * nz, 2 * ny, 2 * nx).copy()


def with_cube(grid, cube):
    return replace(grid, zcorn=np.ascontiguousarray(cube, dtype=float).ravel())


# ── faults and blocks ────────────────────────────────────────────────────────────────────────


def faces_of(faults):
    """A flat list of FaultFace from a face list, a dict of them, or a list of lists."""
    if isinstance(faults, dict):
        faults = list(faults.values())
    out = []
    for item in faults:
        out.extend(faces_of(item) if isinstance(item, (list, tuple, dict)) else [item])
    return out


def face_columns(f):
    """The 0-based (i, j) columns on either side of a 1-based FAULTS face."""
    if f.face == "X":
        return (f.i - 1, f.j - 1), (f.i, f.j - 1)
    return (f.i - 1, f.j - 1), (f.i - 1, f.j)


def interior(grid, faces):
    """Drop faces that sit on the grid rim (one side has no cell).

    Before fault-projections 0.2.0, `adapt_grid(...).faults_after` kept chain edges on the
    rim (invalid indices such as ``J = 0``); FAULTS needs both cells. Kept as a guard.
    """
    nx, ny, _ = grid.cart_dims
    keep = []
    for f in faces_of(faces):
        (i0, j0), (i1, j1) = face_columns(f)
        if 0 <= min(i0, i1) and max(i0, i1) < nx and 0 <= min(j0, j1) and max(j0, j1) < ny:
            keep.append(f)
    return keep


def close_ring(grid, faces):
    """Add the face that closes a ring trace.

    Before fault-projections 0.2.0, `adapt_grid(...).faults_after` left out the edge from
    the last chain pillar back to the first, so a closed trace came back one face short of
    a ring and cut off no block. Kept as a guard.
    """
    faces = faces_of(faces)
    for chain in trace_faults(grid, faces).chains:
        a, b = tuple(map(int, chain.ij[-1])), tuple(map(int, chain.ij[0]))
        if not chain.closed and abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1:
            faces.append(edge_face(a, b))
    return faces


def named_chains(grid, faces, prefix="F"):
    """Split faces into their connected chains: ``{"F1": faces, "F2": ...}`` west to east."""
    chains = trace_faults(grid, faces_of(faces)).chains
    chains = sorted(chains, key=lambda c: (float(np.mean(c.xy[:, 0])), float(np.mean(c.xy[:, 1]))))
    return {f"{prefix}{n + 1}": list(c.faces) for n, c in enumerate(chains)}


def _cut_pairs(faults):
    return {frozenset(face_columns(f)) for f in faces_of(faults)}


def blocks(grid, faults):
    """Label the columns ``[i, j]`` reachable from each other without crossing ``faults``."""
    nx, ny, _ = grid.cart_dims
    cut = _cut_pairs(faults)
    label = -np.ones((nx, ny), int)
    n = 0
    for start in np.ndindex(nx, ny):
        if label[start] >= 0:
            continue
        label[start] = n
        queue = deque([start])
        while queue:
            i, j = queue.popleft()
            for c in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
                if 0 <= c[0] < nx and 0 <= c[1] < ny and label[c] < 0 \
                        and frozenset({(i, j), c}) not in cut:
                    label[c] = n
                    queue.append(c)
        n += 1
    return label


def _column_shift(grid, per_column):
    """Add a per-column, per-horizon shift ``(nz+1, nx, ny)`` to every corner of the column."""
    nx, ny, nz = grid.cart_dims
    add = per_column[_slot_horizon(nz)][:, _slot_column(nx)][:, :, _slot_column(ny)]
    return with_cube(grid, corner_cube(grid) + add.transpose(0, 2, 1))


def slide_blocks(grid, faults, throws):
    """Move whole blocks along their pillars: ``throws = {(i, j) seed column: throw}``.

    A throw is one depth shift (positive = down) or one per horizon, ``nz + 1`` values from
    the top of layer 0 to the base of the last layer: unequal values thicken or thin the
    block's layers, as in a growth fault. Columns outside the seeded blocks stay put.
    """
    nx, ny, nz = grid.cart_dims
    label = blocks(grid, faults)
    shift = np.zeros((nz + 1, nx, ny))
    owner = {}
    for seed, throw in throws.items():
        lab = label[seed]
        if lab in owner:
            raise ValueError(f"columns {owner[lab]} and {seed} are in the same block: "
                             "the FAULTS do not separate them")
        owner[lab] = seed
        shift[:, label == lab] = np.broadcast_to(np.asarray(throw, float), (nz + 1,))[:, None]
    return _column_shift(grid, shift)


# ── throw fields for faults that end inside the grid ─────────────────────────────────────────


def cosine(u):
    """0 at u = 0, 1 at u >= 1, smooth at both ends."""
    return 0.5 * (1.0 - np.cos(np.pi * np.clip(u, 0.0, 1.0)))


@dataclass(frozen=True)
class Throw:
    """How one fault's throw is shaped, in the lessons' terms.

    ``faces`` must form one branch-free chain. ``throw`` is the peak throw (m, positive =
    hanging wall down) or a table ``((s, t), ...)`` along the chain's arc length, which must
    be 0 at every tip. ``hanging_side`` is a plan point on the hanging wall. The other
    fields are those of :class:`fault_projections.Throw`; ``snap`` defaults to ``True`` here
    so that the lessons' layer boundaries cross exactly at pillars.
    """

    faces: tuple
    throw: object = 10.0
    hanging_side: tuple = (0.0, 0.0)
    profile: str = "cosine"
    taper: float | None = None
    decay: object = "auto"
    hanging: float = 1.0
    snap: bool = True


def slide_field(grid, specs, *, faults=None, method="kernel"):
    """Slide faults with :func:`fault_projections.slide_grid`, editing ZCORN only.

    ``faults`` is the full FAULTS set (default: the faces of ``specs``). ``method`` is
    ``"kernel"`` (distance fade) or ``"laplacian"`` (harmonic). The lessons run their own
    checks on the result, so the library's gate and report are switched off here.
    """
    specs = [specs] if isinstance(specs, Throw) else list(specs)
    lib = [FpThrow(list(s.faces), s.throw, hanging_side=tuple(s.hanging_side),
                   hanging=s.hanging, profile=s.profile, taper=s.taper, decay=s.decay,
                   snap=s.snap) for s in specs]
    named = None
    if faults is not None:
        named = {k: faces_of(v) for k, v in faults.items()} if isinstance(faults, dict) \
            else {"FAULTS": faces_of(faults)}
    return slide_grid(grid, lib, faults=named, method=method, on_invalid="ignore",
                      report=False).grid


# ── independent checks (oracles for the library) ────────────────────────────────────────────


def _project(xy, pts):
    """Unsigned distance, arc length of the closest point, and side (+1 left, -1 right)."""
    seg = np.diff(xy, axis=0)
    length = np.linalg.norm(seg, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(length)])
    normal = np.stack([-seg[:, 1], seg[:, 0]], 1) / length[:, None]
    d_out, s_out, side = np.empty(len(pts)), np.empty(len(pts)), np.empty(len(pts))
    for n, p in enumerate(np.asarray(pts, float)):
        t = np.clip(((p - xy[:-1]) * seg).sum(1) / length ** 2, 0.0, 1.0)
        q = xy[:-1] + t[:, None] * seg
        dist = np.linalg.norm(p - q, axis=1)
        k = int(np.argmin(dist))
        nrm = normal[k]
        if t[k] <= 0.0 and k > 0:
            nrm = normal[k - 1] + normal[k]
        elif t[k] >= 1.0 and k < len(seg) - 1:
            nrm = normal[k] + normal[k + 1]
        d_out[n], s_out[n] = dist[k], cum[k] + t[k] * length[k]
        side[n] = np.sign(np.dot(p - q[k], nrm))
    return d_out, s_out, side


def _chain(grid, faces):
    trace = trace_faults(grid, list(faces))
    if not trace.is_simple or len(trace.chains) != 1:
        raise ValueError(f"a Throw needs one branch-free chain of faces: {trace.summary()}")
    return trace.chains[0]


def _ring(p):
    pi, pj = p
    return [(pi - 1, pj - 1), (pi, pj - 1), (pi, pj), (pi - 1, pj)]


def _chain_side(chain, m, cols):
    """+1 / -1: which side of the directed chain the group ``cols`` lies on at pillar m."""
    for k in (m - 1, m):
        if 0 <= k < len(chain.ij) - 1:
            a, b = chain.ij[k], chain.ij[k + 1]
            left = np.array([-(b[1] - a[1]), b[0] - a[0]])
            mid = (a + b) / 2.0
            for c in cols:
                off = np.array(c) + 0.5 - mid
                if abs(off @ (b - a)) < 1.0 - 1e-9 and abs(off @ left) < 1.0:
                    return float(np.sign(off @ left))
    return 0.0


# ── checks ───────────────────────────────────────────────────────────────────────────────────


def unwelded_faces(grid, tol=0.0):
    """The FaultFace of every lateral cell face whose two sides no longer share corners."""
    nx, ny, _ = grid.cart_dims
    z, out = corner_cube(grid), set()
    for i in range(1, nx):
        gap = np.abs(z[:, :, 2 * i] - z[:, :, 2 * i - 1]) > tol
        for j in np.flatnonzero(gap.any(0)) // 2:
            out.add(FaultFace(i, int(j) + 1, "X"))
    for j in range(1, ny):
        gap = np.abs(z[:, 2 * j, :] - z[:, 2 * j - 1, :]) > tol
        for i in np.flatnonzero(gap.any(0)) // 2:
            out.add(FaultFace(int(i) + 1, j, "Y"))
    return out


def _overlap(a0, a1, b0, b1, c0, c1, d0, d1):
    """Integral over t in [0, 1] of max(0, min(b, d) - max(a, c)), all four linear in t."""
    pts = {0.0, 1.0}
    for p0, p1, q0, q1 in ((a0, a1, c0, c1), (b0, b1, d0, d1)):
        den = (p1 - p0) - (q1 - q0)
        if den:
            t = (q0 - p0) / den
            if 0 < t < 1:
                pts.add(t)
    f = lambda t: min(b0 + (b1 - b0) * t, d0 + (d1 - d0) * t) - \
        max(a0 + (a1 - a0) * t, c0 + (c1 - c0) * t)  # noqa: E731
    total = 0.0
    ts = sorted(pts)
    for t0, t1 in zip(ts, ts[1:]):
        f0, f1 = f(t0), f(t1)
        if f0 >= 0 and f1 >= 0:
            total += 0.5 * (f0 + f1) * (t1 - t0)
        elif f0 > 0 or f1 > 0:
            pos = max(f0, f1)
            total += 0.5 * pos * (t1 - t0) * pos / (pos - min(f0, f1))
    return total


def analytic_contacts(grid, faults):
    """Exact contact areas across ``faults`` on vertical pillars, from the corner depths.

    Returns ``{(face, k on the first column, k on the second column): area}``, the first
    column being the one with the lower index across the face.
    """
    nx, ny, nz = grid.cart_dims
    coord = grid.coord.reshape(ny + 1, nx + 1, 6)
    if not np.allclose(coord[..., :2], coord[..., 3:5]):
        raise ValueError("analytic contacts need vertical pillars")
    z = corner_cube(grid)
    out = {}
    for f in faces_of(faults):
        (ia, ja), (ib, jb) = face_columns(f)
        if not (0 <= ia < nx and 0 <= ja < ny and 0 <= ib < nx and 0 <= jb < ny):
            continue
        p, q = face_edge(f)
        width = float(np.linalg.norm(coord[q[1], q[0], :2] - coord[p[1], p[0], :2]))
        if f.face == "X":    # columns (ia, ja) west and (ib, jb) east; pillars run along j
            sa = lambda k, e: z[k, 2 * ja + e, 2 * ia + 1]  # noqa: E731
            sb = lambda k, e: z[k, 2 * jb + e, 2 * ib]      # noqa: E731
        else:
            sa = lambda k, e: z[k, 2 * ja + 1, 2 * ia + e]  # noqa: E731
            sb = lambda k, e: z[k, 2 * jb, 2 * ib + e]      # noqa: E731
        for ka in range(nz):
            for kb in range(nz):
                a = _overlap(sa(2 * ka, 0), sa(2 * ka, 1), sa(2 * ka + 1, 0), sa(2 * ka + 1, 1),
                             sb(2 * kb, 0), sb(2 * kb, 1), sb(2 * kb + 1, 0), sb(2 * kb + 1, 1))
                if a * width > 1e-9:
                    out[(f, ka, kb)] = a * width
    return out


@dataclass
class Processed:
    """A re-processed grid: topology, face areas, cell volumes and (i, j, k) per cell."""

    topology: object
    areas: np.ndarray
    centroids: np.ndarray
    volumes: np.ndarray
    ijk: np.ndarray

    def volume_grid(self, cart_dims):
        """Cell volumes as ``[i, j, k]``; NaN for cells the topology dropped."""
        out = np.full(cart_dims, np.nan)
        out[tuple(self.ijk.T)] = self.volumes
        return out

    def contacts(self, faults):
        """``{(face, k first column, k second column): area}`` across ``faults``."""
        lookup = {}
        for f in faces_of(faults):
            a, b = face_columns(f)
            lookup[(a, b)] = f
        nb = self.topology.faces.neighbours
        out = {}
        for face in np.flatnonzero(self.topology.faces.direction <= 1):
            a, b = nb[face]
            if b < 0:
                continue
            ca, cb = self.ijk[a], self.ijk[b]
            key = None
            for first, second in ((ca, cb), (cb, ca)):
                f = lookup.get(((int(first[0]), int(first[1])), (int(second[0]), int(second[1]))))
                if f is not None:
                    key = (f, int(first[2]), int(second[2]))
            if key is not None:
                out[key] = out.get(key, 0.0) + float(self.areas[face])
        return out


def process(grid_or_path):
    """Build the topology (from a deck path or a grid) and its face areas and cell volumes."""
    grid = read_grdecl(Path(grid_or_path)) if isinstance(grid_or_path, (str, Path)) \
        else grid_or_path
    topology = build_topology(from_grdecl(grid))
    fg = face_geometry(topology.nodes, topology.faces, topology.cell_corners)
    volumes = cell_geometry(topology.nodes, topology.faces, topology.cell_faces, fg).volumes
    return Processed(topology, fg.areas, fg.centroids, volumes, topology.original_ijk)


def crossed_faces(topology) -> int:
    """Lateral faces with an intersection node: a layer boundary crosses inside them."""
    corners = topology.n_corner_nodes
    lateral = topology.faces.direction < 2
    return 0 if corners is None else int((topology.faces.nodes[lateral] >= corners).any(1).sum())


def contacts_by_layer(contacts, flip=()):
    """Sum contacts into ``{(k first side, k second side): area}``.

    Faces listed in ``flip`` have their two sides swapped first, so that every face can be
    read from the same geological side.
    """
    flip = set(flip)
    out = {}
    for (f, ka, kb), area in contacts.items():
        key = (kb, ka) if f in flip else (ka, kb)
        out[key] = out.get(key, 0.0) + area
    return out


def compare(found, expected, *, rtol=1e-9, atol=1e-6):
    """Assert two contact dictionaries agree, ignoring slivers below ``atol``."""
    f = {k: v for k, v in found.items() if v > atol}
    e = {k: v for k, v in expected.items() if v > atol}
    assert f.keys() == e.keys(), (sorted(set(f) ^ set(e), key=str)[:6])
    for key in e:
        assert np.isclose(f[key], e[key], rtol=rtol, atol=atol), (key, f[key], e[key])
    return sum(e.values())


# ── input and output ─────────────────────────────────────────────────────────────────────────


def write_deck(path, grid, faults, *, header="", extra=""):
    """Grid INCLUDE plus a FAULTS section (``{name: faces}``) and any extra keyword text."""
    path = Path(path)
    write_grid_include(grid, path, header=header)
    _, _, nz = grid.cart_dims
    with open(path, "a") as deck:
        if faults:
            deck.write("\n" + faults_text(list(faults.values()), nz, names=list(faults)))
        if extra:
            deck.write("\n" + extra)
    return path


def column_of_point(grid, point):
    """The (i, j) column whose top quad contains a plan point (by the pillar tops)."""
    top = pillar_lattice(grid).top
    nx, ny, _ = grid.cart_dims
    p = np.asarray(point, float)
    for i, j in np.ndindex(nx, ny):
        quad = np.array([top[i, j], top[i + 1, j], top[i + 1, j + 1], top[i, j + 1]])
        edge, rel = np.roll(quad, -1, 0) - quad, p - quad
        if np.all(edge[:, 0] * rel[:, 1] - edge[:, 1] * rel[:, 0] >= -1e-9):
            return i, j
    raise ValueError(f"point {point} is outside the grid")


def realized_throw(before, after, faces, hanging_side):
    """Arc length and relative throw (hanging wall down minus footwall) at each chain pillar."""
    chain = _chain(before, faces)
    nx, ny, _ = before.cart_dims
    shift = corner_cube(after)[0] - corner_cube(before)[0]           # top of layer 0
    _, _, hw = _project(chain.xy, [hanging_side])
    cum = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(chain.xy, axis=0), axis=1))])
    out = []
    for m, p in enumerate(map(tuple, chain.ij)):
        down = {1.0: [], -1.0: []}
        for c in _ring(p):
            if not (0 <= c[0] < nx and 0 <= c[1] < ny):
                continue
            sd = _chain_side(chain, m, [c])
            if sd:
                u, v = p[0] - c[0], p[1] - c[1]
                down[sd].append(shift[2 * c[1] + v, 2 * c[0] + u])
        if down[1.0] and down[-1.0]:
            hwv, fwv = np.mean(down[hw[0]]), np.mean(down[-hw[0]])
            out.append((cum[m], hwv - fwv))
        else:
            out.append((cum[m], 0.0))
    return np.array(out)
