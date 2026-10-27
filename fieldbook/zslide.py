"""Prototype of the planned ZCORN-only fault slide, shared by lessons 23-29.

This is teaching code, not a library call: the call itself is planned for fault-projections.
The defaults follow the decisions recorded for it:

- far field ``"auto"``: a fault that closes a block slides it rigidly; a fault that ends
  inside the grid tapers to zero at each tip and its motion decays away from it;
- decay by a distance kernel by default, or by a Laplacian (harmonic) solve on the pillars;
- the hanging wall takes the whole throw (``hanging=1.0``); a cosine profile toward each tip.

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
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

__all__ = [
    "layer_cake", "corner_cube", "with_cube", "faces_of", "face_columns",
    "blocks", "slide_blocks", "Throw", "slide_field", "unwelded_faces", "analytic_contacts",
    "Processed", "process", "contacts_by_layer", "write_deck", "realized_throw",
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

    `adapt_grid(...).faults_after` keeps a chain edge that runs along the rim, which
    comes out as an invalid index such as ``J = 0``; FAULTS needs both cells.
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

    `adapt_grid(...).faults_after` reads faces between consecutive chain pillars and
    leaves out the edge from the last pillar back to the first, so a closed trace comes
    back one face short of a ring and cuts off no block.
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
    """How one fault's throw is shaped; a prototype of the planned slide options.

    ``faces`` must form one branch-free chain. ``throw`` is the peak throw (m, positive =
    hanging wall down) or a table ``((s, t), ...)`` along the chain's arc length from its
    first pillar, which must be 0 at every tip. ``hanging_side`` is a plan point on the
    hanging wall. ``decay="auto"`` keeps the whole side moving when the fault has no tip,
    and otherwise fades the motion to zero at half the fault length. ``snap`` moves the
    throw at a fault pillar onto a layer boundary whenever the next face would otherwise
    carry that boundary across its inside (a contact polygon with five corners, which
    re-processing rejects); the crossing then happens exactly at a pillar.
    """

    faces: tuple
    throw: object = 10.0
    hanging_side: tuple = (0.0, 0.0)
    profile: str = "cosine"
    taper: float | None = None
    decay: object = "auto"
    hanging: float = 1.0
    snap: bool = True


class _Groups:
    """Around each pillar, the columns that stay welded together: one node per group."""

    def __init__(self, nx, ny, faults):
        cut = _cut_pairs(faults)
        self.nx, self.ny = nx, ny
        self.node_of = {}               # (pillar, column) -> node id
        self.pillar = []                # node id -> pillar (i, j)
        self.columns = []               # node id -> columns in the group
        for p in np.ndindex(nx + 1, ny + 1):
            pi, pj = p
            ring = [(pi - 1, pj - 1), (pi, pj - 1), (pi, pj), (pi - 1, pj)]
            present = [c for c in ring if 0 <= c[0] < nx and 0 <= c[1] < ny]
            parent = {c: c for c in present}

            def find(c):
                while parent[c] != c:
                    c = parent[c]
                return c
            for a, b in zip(ring, ring[1:] + ring[:1]):
                if a in parent and b in parent and frozenset({a, b}) not in cut:
                    parent[find(a)] = find(b)
            roots = {}
            for c in present:
                roots.setdefault(find(c), []).append(c)
            for cols in roots.values():
                node = len(self.pillar)
                self.pillar.append(p)
                self.columns.append(cols)
                for c in cols:
                    self.node_of[(p, c)] = node

    @staticmethod
    def _ring(p):
        pi, pj = p
        return [(pi - 1, pj - 1), (pi, pj - 1), (pi, pj), (pi - 1, pj)]


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


def _along(spec, s, s_tips, length):
    """Relative throw at arc lengths ``s`` along the chain."""
    s = np.asarray(s, float)
    if not np.isscalar(spec.throw) and np.ndim(spec.throw) == 2:
        table = np.asarray(spec.throw, float)
        out = np.interp(s, table[:, 0], table[:, 1])
        for tip in s_tips:
            if abs(np.interp(tip, table[:, 0], table[:, 1])) > 1e-9:
                raise ValueError("a throw table must be 0 at every tip")
        return out
    peak = float(spec.throw)
    if not s_tips:
        return np.full_like(s, peak)
    rise = spec.taper or (length / len(s_tips))
    u = np.min([np.abs(np.clip(s, 0.0, length) - t) for t in s_tips], axis=0) / rise
    shape = np.clip(u, 0.0, 1.0) if spec.profile == "linear" else cosine(u)
    return peak * shape


def _fade(spec, d, has_tips, length):
    decay = spec.decay
    if decay == "auto":
        decay = 0.5 * length if has_tips else None
    if decay is None:
        return np.ones_like(d), None
    return 1.0 - cosine(d / decay), float(decay)


def slide_field(grid, specs, *, faults=None, method="kernel"):
    """Apply one throw field per ``Throw``, editing ZCORN only.

    Every pillar carries one shift per group of columns that stay welded around it, so
    every face outside ``faults`` stays welded and the throw is zero wherever a fault ends.
    ``faults`` is the full FAULTS set (default: the faces of ``specs``); a field acts only on
    the blocks its fault touches. ``method`` is ``"kernel"`` (distance fade) or
    ``"laplacian"`` (harmonic between the fault and the decay distance, or the rim).
    """
    nx, ny, _ = grid.cart_dims
    specs = [specs] if isinstance(specs, Throw) else list(specs)
    allf = faces_of(faults) if faults is not None else faces_of([s.faces for s in specs])
    groups = _Groups(nx, ny, allf)
    label = blocks(grid, allf)
    xy = pillar_lattice(grid).top                                  # (nx+1, ny+1, 2)
    value = np.zeros(len(groups.pillar))
    for spec in specs:
        value += _field(grid, spec, groups, label, xy, method)
    add = np.zeros((nx, ny, 2, 2))
    for i, j in np.ndindex(nx, ny):
        for u, v in np.ndindex(2, 2):
            add[i, j, u, v] = value[groups.node_of[((i + u, j + v), (i, j))]]
    _, _, nz = grid.cart_dims
    cube = corner_cube(grid)
    per_slot = add.transpose(1, 3, 0, 2).reshape(2 * ny, 2 * nx)   # [j-slot, i-slot]
    return with_cube(grid, cube + per_slot[None])


def _field(grid, spec, groups, label, xy, method):
    chain = _chain(grid, spec.faces)
    cxy = chain.xy
    length = chain.length
    cum = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(cxy, axis=0), axis=1))])
    # A chain end is a tip when one welded group surrounds its pillar.
    ends = [(0, tuple(chain.ij[0])), (-1, tuple(chain.ij[-1]))]
    n_groups = lambda p: len({groups.node_of[(p, c)] for c in _Groups._ring(p)  # noqa: E731
                              if (p, c) in groups.node_of})
    s_tips = [0.0 if e == 0 else length for e, p in ends if n_groups(p) == 1]
    _, _, hw = _project(cxy, [spec.hanging_side])
    if hw[0] == 0:
        raise ValueError("hanging_side lies on the fault")
    touched = {label[c] for f in spec.faces for c in face_columns(f)
               if 0 <= c[0] < label.shape[0] and 0 <= c[1] < label.shape[1]}
    on_chain = {tuple(p): m for m, p in enumerate(chain.ij)}
    nodes = range(len(groups.pillar))
    d = np.zeros(len(groups.pillar))
    s = np.zeros(len(groups.pillar))
    side = np.zeros(len(groups.pillar))
    inside = np.zeros(len(groups.pillar), bool)
    for n in nodes:
        p, cols = groups.pillar[n], groups.columns[n]
        inside[n] = label[cols[0]] in touched
        if p in on_chain:
            m = on_chain[p]
            s[n] = cum[m]
            side[n] = _chain_side(chain, m, cols)
        else:
            d[n], s[n], side[n] = _project(cxy, [xy[p]])
    throw = _along(spec, s, s_tips, length)
    if spec.snap:
        order = [n for m in range(len(chain.ij)) for n in nodes
                 if groups.pillar[n] == tuple(chain.ij[m])]
        snapped = _snap(grid, chain, np.array([_along(spec, [cum[m]], s_tips, length)[0]
                                               for m in range(len(chain.ij))]), hw[0])
        for n in order:
            throw[n] = snapped[on_chain[groups.pillar[n]]]
    share = np.where(side == hw[0], spec.hanging, spec.hanging - 1.0)
    share[side == 0] = 0.0
    fade, decay = _fade(spec, d, bool(s_tips), length)
    if method == "kernel":
        return np.where(inside, share * throw * fade, 0.0)
    if method != "laplacian":
        raise ValueError(f"unknown method {method!r}")
    chain_node = np.array([groups.pillar[n] in on_chain for n in nodes])
    fixed = chain_node | ~inside
    if decay is not None:
        fixed |= d >= decay
    target = np.where(chain_node & inside, share * throw, 0.0)
    rows, cols_ = [], []
    for i, j in np.ndindex(groups.nx, groups.ny):
        corners = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
        ids = [groups.node_of[(p, (i, j))] for p in corners]
        for a, b in zip(ids, ids[1:] + ids[:1]):
            rows += [a, b]
            cols_ += [b, a]
    n = len(groups.pillar)
    w = coo_matrix((np.ones(len(rows)), (rows, cols_)), shape=(n, n)).tocsr()
    lap = (coo_matrix((np.asarray(w.sum(1)).ravel(), (range(n), range(n))), shape=(n, n))
           - w).tocsr()
    free = ~fixed
    out = target.copy()
    if free.any():
        rhs = -lap[free][:, fixed] @ target[fixed]
        out[free] = spsolve(lap[free][:, free].tocsc(), rhs)
    return out


def _snap(grid, chain, throw, hw):
    """Lower the throw at fault pillars until no face carries a layer boundary across.

    The horizons are read on each side of every fault pillar, so throw already on the
    fault (an earlier slide) counts: a hanging-wall horizon ``a`` meets a footwall horizon
    ``b`` when the added relative throw equals ``F_b - H_a``.
    """
    cube = corner_cube(grid)
    nx, ny, nz = grid.cart_dims

    def horizons(p, c):
        col = cube[:, 2 * c[1] + p[1] - c[1], 2 * c[0] + p[0] - c[0]]
        return np.append(col[0::2], col[-1])
    crit = []
    for m, p in enumerate(map(tuple, chain.ij)):
        sides = {}
        for c in _Groups._ring(p):
            if 0 <= c[0] < nx and 0 <= c[1] < ny:
                sides.setdefault(_chain_side(chain, m, [c]), c)
        hang, foot = sides.get(hw), sides.get(-hw)
        hang = hang if hang is not None else foot
        foot = foot if foot is not None else hang
        crit.append(np.subtract.outer(horizons(p, foot), horizons(p, hang)).ravel())
    out = throw.copy()
    tol = 1e-9 * max(1.0, float(np.abs(cube).max()))       # horizons meeting exactly stay put
    changed, passes = True, 0
    while changed:
        changed, passes = False, passes + 1
        if passes > 10 * len(out) * len(crit[0]):
            raise RuntimeError("snap did not settle; lower the throw")
        for m in range(1, len(out)):
            a, b = out[m - 1] - crit[m - 1], out[m] - crit[m]
            crossed = np.flatnonzero(((a < -tol) & (b > tol)) | ((a > tol) & (b < -tol)))
            if crossed.size:
                big, small = (m, m - 1) if abs(out[m]) >= abs(out[m - 1]) else (m - 1, m)
                options = crit[big][crossed]
                out[big] = options[np.argmin(np.abs(options - out[small]))]
                changed = True
    return out


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
        for c in _Groups._ring(p):
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
