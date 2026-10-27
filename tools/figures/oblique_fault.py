"""Example 08 figures: an oblique, dipping fault (staircase FAULTS), slid 25 m, regular vs adapted.

usage: python oblique_fault.py DIR
  DIR/08_oblique_fault/{A_regular,B_adapted}_{0_input,1_slid}.grdecl come from fieldbook
  lesson 21 (--out).  writes DIR/08_oblique_fault.png, DIR/08_oblique_fault_3d.html
"""

import sys
from collections import deque
from pathlib import Path


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from cpg2unstructured import build_topology, from_grdecl, read_grdecl
from cpg2unstructured.geometry import face_geometry
from fault_projections.from_faults import faults_from_deck

NX, NY, NZ, DX, DZ, THROW = 8, 6, 3, 100.0, 40.0, 25.0
TRACE = np.array([[300.0, 0.0], [600.0, 600.0]])
SEED = (0, NX - 1)
LEAN = 75.0  # hanging-wall pillars move 75 m east by the base
BASE = NZ * DZ
_n = np.cross(np.append(TRACE[1] - TRACE[0], 0.0), [LEAN, 0.0, BASE])
DIP = float(np.degrees(np.arccos(abs(_n[2]) / np.linalg.norm(_n))))
ROW = 0  # pillar row of the cross-section (south face)
VEX = 2.0  # vertical exaggeration of 3D views and section
OUT = Path(sys.argv[1])
DECKS = OUT / "08_oblique_fault"
KCOL = ["#9ecae1", "#a1d99b", "#fdd0a2"]
EXPOSED = "#d7191c"
PAIRCOL = {
    (0, 0): "#1f77b4",
    (1, 1): "#2ca02c",
    (2, 2): "#ff7f0e",
    (1, 0): "#9467bd",
    (2, 1): "#8c564b",
}
HEX_FACES = [[0, 2, 6, 4], [1, 5, 7, 3], [0, 4, 5, 1], [2, 3, 7, 6], [0, 1, 3, 2], [4, 6, 7, 5]]
HEX_EDGES = [
    (0, 1),
    (2, 3),
    (4, 5),
    (6, 7),
    (0, 2),
    (1, 3),
    (4, 6),
    (5, 7),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
]


def face_sides(f):
    return (f.j - 1, f.i - 1), (f.j - 1, f.i) if f.face == "X" else (f.j, f.i - 1)


def face_pillars(f):
    return ((f.i, f.j - 1), (f.i, f.j)) if f.face == "X" else ((f.i - 1, f.j), (f.i, f.j))


def moving_block(faces):
    cut = {frozenset(face_sides(f)) for f in faces}
    block = np.zeros((NY, NX), bool)
    block[SEED] = True
    q = deque([SEED])
    while q:
        j, i = q.popleft()
        for n in ((j, i + 1), (j, i - 1), (j + 1, i), (j - 1, i)):
            if (
                0 <= n[0] < NY
                and 0 <= n[1] < NX
                and not block[n]
                and frozenset({(j, i), n}) not in cut
            ):
                block[n] = True
                q.append(n)
    return block


def load(name):
    path = DECKS / name
    g = read_grdecl(path)
    t = build_topology(from_grdecl(g))
    fg = face_geometry(t.nodes, t.faces, t.cell_corners)
    faces = faults_from_deck(path)["OBLIQUE"]
    coord = g.coord.reshape(NY + 1, NX + 1, 6)
    return dict(
        g=g,
        t=t,
        areas=fg.areas,
        faces=faces,
        block=moving_block(faces),
        top=coord[:, :, :2],
        coord=coord,
        z=g.zcorn.reshape(2 * NZ, 2 * NY, 2 * NX),
    )


def pillar_xy(st, z):
    """(NY+1, NX+1, 2) pillar positions at depth z."""
    c = st["coord"]
    u = (z - c[..., 2:3]) / (c[..., 5:6] - c[..., 2:3])
    return c[..., :2] + u * (c[..., 3:5] - c[..., :2])


def side_of(st, cell, centroid, z):
    """Which column lies across the lateral face of `cell` whose centroid (at depth z) is given."""
    i, j, _ = st["t"].original_ijk[cell]
    top = pillar_xy(st, z)
    sides = {
        (i - 1, j): (top[j, i], top[j + 1, i]),
        (i + 1, j): (top[j, i + 1], top[j + 1, i + 1]),
        (i, j - 1): (top[j, i], top[j, i + 1]),
        (i, j + 1): (top[j + 1, i], top[j + 1, i + 1]),
    }

    def dist(seg):
        a, b = seg
        u = np.clip(np.dot(centroid - a, b - a) / np.dot(b - a, b - a), 0, 1)
        return np.linalg.norm(a + u * (b - a) - centroid)

    return min(sides, key=lambda n: dist(sides[n]))


def fault_pieces(st):
    """[(face, key)] on the fault; key (k fixed, k moving) or 'exposed'."""
    t, ijk, block, out = st["t"], st["t"].original_ijk, st["block"], []
    for f in np.flatnonzero(t.faces.direction <= 1):
        a, b = t.faces.neighbours[f]
        ca = ijk[a]
        if b >= 0:
            cb = ijk[b]
            if block[ca[1], ca[0]] == block[cb[1], cb[0]]:
                continue
            fixed, moving = (ca, cb) if block[cb[1], cb[0]] else (cb, ca)
            out.append((f, (int(fixed[2]), int(moving[2]))))
        else:
            c = t.nodes[t.faces.nodes[f]].mean(0)
            ni, nj = side_of(st, a, c[:2], c[2])
            if 0 <= ni < NX and 0 <= nj < NY and block[nj, ni] != block[ca[1], ca[0]]:
                out.append((f, "exposed"))
    return out


def colour(key):
    return EXPOSED if key == "exposed" else PAIRCOL[key]


def ordered(pts):
    pts = np.unique(np.round(pts, 9), axis=0)
    c = pts.mean(0)
    _, _, vt = np.linalg.svd(pts - c)
    uv = (pts - c) @ vt[:2].T
    return pts[np.argsort(np.arctan2(uv[:, 1], uv[:, 0]))]


def fault_length(st):
    """Fault area per metre of depth: sum of |edge x lean| over the FAULTS faces."""
    c = st["coord"]
    lean = (c[..., 3:] - c[..., :3]) / (c[..., 5:] - c[..., 2:3])
    return sum(
        float(np.linalg.norm(np.cross(c[b[1], b[0], :3] - c[a[1], a[0], :3], lean[a[1], a[0]])))
        for a, b in map(face_pillars, st["faces"])
    )


def plan_length(st):
    top = st["top"]
    return sum(
        float(np.linalg.norm(top[b[1], b[0]] - top[a[1], a[0]]))
        for a, b in map(face_pillars, st["faces"])
    )


def map_panel(ax, st, title):
    t, top = st["t"], st["top"]
    polys, vals = [], []
    for c in np.flatnonzero(t.original_ijk[:, 2] == 0):
        p = t.nodes[t.cell_corners[c]][[0, 1, 3, 2]]
        polys.append(p[:, :2])
        vals.append(p[:, 2].mean())
    pc = PolyCollection(
        polys, array=np.array(vals), cmap="viridis_r", edgecolors="0.35", linewidths=0.6
    )
    pc.set_clim(0, THROW)
    ax.add_collection(pc)
    segs = [[top[a[1], a[0]], top[b[1], b[0]]] for a, b in map(face_pillars, st["faces"])]
    ax.add_collection(LineCollection(segs, colors="k", linewidths=4.5, capstyle="round"))
    ax.plot(*TRACE.T, color="#d7191c", lw=2, ls="--")
    base = pillar_xy(st, NZ * DZ)
    segs_b = [[base[a[1], a[0]], base[b[1], b[0]]] for a, b in map(face_pillars, st["faces"])]
    ax.add_collection(LineCollection(segs_b, colors="w", linewidths=1.6, linestyles=":"))
    mid = TRACE.mean(0)
    nrm = np.array([2.0, -1.0]) / np.sqrt(5.0)
    ax.annotate(
        "", mid + 150 * nrm, mid + 20 * nrm, arrowprops=dict(arrowstyle="-|>", color="w", lw=2.2)
    )
    ax.text(
        *(mid + 160 * nrm),
        f"dips {DIP:.0f}°",
        color="w",
        fontsize=9.5,
        fontweight="bold",
        ha="left",
        va="center",
    )
    for n, f in enumerate(st["faces"]):
        (a, b) = face_pillars(f)
        m = (top[a[1], a[0]] + top[b[1], b[0]]) / 2
        off = np.array([14, -14]) if f.face == "X" else np.array([0, 16])
        ax.text(
            *(m + off),
            f"{f.face}",
            fontsize=7.5,
            ha="left" if f.face == "X" else "center",
            va="center",
            color="k",
            fontweight="bold",
        )
    ax.text(
        110,
        520,
        "fixed block\n(footwall)",
        fontsize=10,
        ha="center",
        bbox=dict(fc="white", ec="none", alpha=0.8),
    )
    ax.text(
        690,
        90,
        f"moving block (hanging wall)\n{st['block'].sum()} columns, -{THROW:g} m",
        fontsize=10,
        ha="center",
        color="k",
        bbox=dict(fc="white", ec="none", alpha=0.8),
    )
    ax.set_xlim(-10, NX * DX + LEAN + 10)
    ax.set_ylim(-10, NY * DX + 10)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(
        f"{title}\nFAULTS: {len(st['faces'])} faces, plan length {plan_length(st):.1f} m, "
        f"area {fault_length(st):.1f} m²/m",
        fontsize=10.5,
        loc="left",
    )
    return pc


def view3d_panel(ax, st, title):
    """Outer surface of both blocks, solid, plus the exposed fault pieces in red."""
    t = st["t"]
    exposed = {int(f) for f, key in fault_pieces(st) if key == "exposed"}
    shade = {0: 0.86, 1: 0.72, 2: 1.0}  # X, Y, top/bottom faces
    polys, fc = [], []
    for f in np.flatnonzero(t.faces.neighbours[:, 1] < 0):
        a = t.faces.neighbours[f, 0]
        k = t.original_ijk[a][2]
        rgb = np.array(to_rgba(EXPOSED if f in exposed else KCOL[k]))
        rgb[:3] *= shade[int(t.faces.direction[f])]
        polys.append(ordered(t.nodes[t.faces.nodes[f]]))
        fc.append(rgb)
    ax.add_collection3d(
        Poly3DCollection(polys, facecolors=fc, edgecolors=to_rgba("0.2", 0.7), linewidths=0.4)
    )
    ax.computed_zorder = False
    ax.plot(*TRACE.T, [0, 0], color="k", lw=1.5, ls="--", zorder=20)
    c = st["coord"]
    for i in {p[0] for f in st["faces"] for p in face_pillars(f) if p[1] == 0}:
        ax.plot([c[0, i, 0], c[0, i, 3]], [0, 0], [0, BASE], color="k", lw=2.5, zorder=21)
    ax.set_xlim(0, NX * DX + LEAN)
    ax.set_ylim(0, NY * DX)
    ax.set_zlim(BASE + THROW, 0)
    ax.set_box_aspect((NX * DX + LEAN, NY * DX, VEX * (BASE + THROW)), zoom=1.35)
    ax.view_init(elev=22, azim=-72)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("depth [m]")
    ax.set_title(title, fontsize=11)


def area_panel(ax, sts, labels):
    keys = list(PAIRCOL)
    w = 0.38
    for n, (st, lab, hatch) in enumerate(zip(sts, labels, ("", "//"))):
        tot = {}
        for f, key in fault_pieces(st):
            tot[key] = tot.get(key, 0.0) + st["areas"][f]
        L = fault_length(st)
        xs = np.arange(len(keys)) + (n - 0.5) * w
        vals = [tot.get(k, 0.0) for k in keys]
        ax.bar(
            xs,
            vals,
            w,
            color=[PAIRCOL[k] for k in keys],
            edgecolor="k",
            hatch=hatch,
            label=f"{lab}: A = {L:.1f} m²/m",
        )
        for x, k, v in zip(xs, keys, vals):
            ov = min((k[0] + 1) * DZ, (k[1] + 1) * DZ + THROW) - max(k[0] * DZ, k[1] * DZ + THROW)
            ax.plot([x - w / 2, x + w / 2], [L * ov] * 2, color="k", lw=2)
            ax.text(x, v + 400, f"{v:,.0f}", ha="center", fontsize=8)
        ex = tot.get("exposed", 0.0)
        ax.bar(len(keys) + (n - 0.5) * w, ex, w, color=EXPOSED, edgecolor="k", hatch=hatch)
        ax.text(len(keys) + (n - 0.5) * w, ex + 400, f"{ex:,.0f}", ha="center", fontsize=8)
    ax.set_xticks(
        range(len(keys) + 1),
        [
            f"k{a} | k{b}\n{min((a + 1) * DZ, (b + 1) * DZ + THROW) - max(a * DZ, b * DZ + THROW):g} m"
            for a, b in keys
        ]
        + ["exposed\n(both sides)"],
        fontsize=9,
    )
    ax.set_xlabel("fixed layer | moving layer, depth overlap")
    ax.set_ylabel("contact area [m²]")
    handles = [
        Patch(fc="white", ec="k", label=labels[0] + f" (A = {fault_length(sts[0]):.1f} m²/m)"),
        Patch(
            fc="white",
            ec="k",
            hatch="//",
            label=labels[1] + f" (A = {fault_length(sts[1]):.1f} m²/m)",
        ),
        Line2D([], [], color="k", lw=2, label="expected: A × depth overlap"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=9)
    ax.set_title(
        "Re-processed contact areas: same pairs on both decks;\neach is the fault area "
        "per metre of depth A × the depth overlap",
        fontsize=10.5,
        loc="left",
    )


def corner_xz(st, z, i, j, k, di, dk):
    """x, z of the south-side corner (di, dk) of cell (i, j, k) from COORD + ZCORN."""
    zz = z[2 * k + dk, 2 * j, 2 * i + di]
    c = st["coord"][j, i + di]
    u = (zz - c[2]) / (c[5] - c[2])
    return c[0] + u * (c[3] - c[0]), zz


def section_panel(ax, before, after):
    """South faces of cell row ROW (on pillar row ROW), before (dashed) and after the slide."""
    for st, solid in ((before, False), (after, True)):
        for i in range(NX):
            for k in range(NZ):
                q = [
                    corner_xz(st, st["z"], i, ROW, k, di, dk)
                    for di, dk in ((0, 0), (1, 0), (1, 1), (0, 1))
                ]
                if solid:
                    ax.add_patch(plt.Polygon(q, fc=KCOL[k], ec="0.25", lw=0.8))
                elif st["block"][ROW, i]:
                    ax.add_patch(plt.Polygon(q, fill=False, ec="0.25", lw=1.0, ls="--"))
    c = after["coord"]
    i0 = next(i for i in range(NX + 1) if any((i, ROW) in face_pillars(f) for f in after["faces"]))
    ax.plot(
        [c[ROW, i0, 0], c[ROW, i0, 3] + (c[ROW, i0, 3] - c[ROW, i0, 0]) * THROW / BASE],
        [0, BASE + THROW],
        color="k",
        lw=2.5,
    )
    for i, di, k, dk in ((i0, 0, 0, 0), (NX - 1, 1, 0, 0), (5, 0, NZ - 1, 1)):
        x0, z0 = corner_xz(before, before["z"], i, ROW, k, di, dk)
        x1, z1 = corner_xz(after, after["z"], i, ROW, k, di, dk)
        ax.annotate(
            "",
            (x1, z1),
            (x0, z0),
            arrowprops=dict(arrowstyle="-|>", color="#d7191c", lw=2, shrinkA=0, shrinkB=0),
        )
    ax.text(
        8,
        BASE + THROW + 1,
        f"red arrows: every moving corner slides {np.hypot(x1 - x0, z1 - z0):.1f} m"
        f"\nalong its pillar (Δx {x1 - x0:+.1f}, Δz {z1 - z0:+.1f}): the hanging\n"
        "wall moves as one piece, no cell changes volume",
        fontsize=9,
        color="#d7191c",
        va="bottom",
    )
    ax.set_xlim(-5, NX * DX + LEAN + 5)
    ax.set_ylim(BASE + THROW + 3, -3)
    ax.set_aspect(VEX)
    ax.set_xlabel(f"x [m]   (depth exaggerated {VEX:g}×)")
    ax.set_ylabel("depth [m]")
    box = dict(fc="white", ec="none", alpha=0.85, pad=1.5)
    ax.text(8, 6, "fixed block (footwall): vertical pillars", fontsize=9, va="top", bbox=box)
    ax.text(
        NX * DX - 25,
        2,
        "moving block (hanging wall): every pillar parallel to the "
        "fault\ndashed = before the slide",
        fontsize=9,
        ha="right",
        va="top",
        bbox=box,
    )
    ax.set_title(
        f"B, section along the south face (pillar row j={ROW}, the front of the 3D views):"
        f"\nthe fault pillar and every hanging-wall pillar lean {LEAN:g} m east over "
        f"{BASE:g} m of depth",
        fontsize=10.5,
        loc="left",
    )


def fig_png(a, b, b0):
    fig = plt.figure(figsize=(18, 17))
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.0, 0.62, 0.95],
        left=0.045,
        right=0.985,
        top=0.95,
        bottom=0.04,
        hspace=0.27,
    )
    top = gs[0].subgridspec(1, 2, wspace=0.0)
    view3d_panel(
        fig.add_subplot(top[0], projection="3d"),
        b0,
        f"B before the slide: the fault dips {DIP:.0f}° SE (depth exaggerated {VEX:g}×)",
    )
    view3d_panel(
        fig.add_subplot(top[1], projection="3d"),
        b,
        f"B after: the whole hanging wall slid {THROW:g} m down the fault;\n"
        "red = fault face left exposed above the moved block",
    )
    mid = gs[1].subgridspec(1, 2, width_ratios=[1.55, 1], wspace=0.1)
    section_panel(fig.add_subplot(mid[0]), b0, b)
    area_panel(fig.add_subplot(mid[1]), [a, b], ["A regular", "B adapted"])
    bot = gs[2].subgridspec(1, 2, wspace=0.08)
    ax_a, ax_b = fig.add_subplot(bot[0]), fig.add_subplot(bot[1])
    pc = map_panel(ax_a, a, "A  regular grid, staircase FAULTS")
    map_panel(ax_b, b, "B  same FAULTS, grid adapted to the trace")
    cb = fig.colorbar(
        pc, ax=[ax_a, ax_b], orientation="horizontal", fraction=0.035, pad=0.09, aspect=60
    )
    cb.set_label(
        f"top of layer k0 after the slide, mean depth [m]  ·  red dashes: trace at z = 0"
        f"  ·  white dots: fault at the base (z = {BASE:g})  ·  X / Y: FAULTS faces",
        fontsize=9.5,
    )
    fig.suptitle(
        f"A fault oblique to the grid (27° off J) that dips {DIP:.0f}° SE: staircase X/Y "
        f"FAULTS, slid {THROW:g} m by editing ZCORN only (COORD and FAULTS unchanged)",
        fontweight="bold",
        y=0.99,
    )
    path = OUT / "08_oblique_fault.png"
    path.unlink(missing_ok=True)
    fig.savefig(path, dpi=100)
    plt.close(fig)
    return path


def fig_html(states, titles):
    import plotly.graph_objects as go

    fig, groups, cells = go.Figure(), [], []
    for s, st in enumerate(states):
        t, start, last = st["t"], len(fig.data), s == len(states) - 1
        for k in range(NZ):
            xs, ys, zs, ii, jj, kk, txt = [], [], [], [], [], [], []
            for c in np.flatnonzero(t.original_ijk[:, 2] == k):
                i, j, _ = t.original_ijk[c]
                p, base = t.nodes[t.cell_corners[c]], len(xs)
                xs += list(p[:, 0])
                ys += list(p[:, 1])
                zs += list(p[:, 2])
                side = "moving" if st["block"][j, i] else "fixed"
                txt += [f"cell ({i},{j},{k}) {side}"] * 8
                for q in HEX_FACES:
                    ii += [base + q[0], base + q[0]]
                    jj += [base + q[1], base + q[2]]
                    kk += [base + q[2], base + q[3]]
            cells.append(len(fig.data))
            fig.add_trace(
                go.Mesh3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    i=ii,
                    j=jj,
                    k=kk,
                    color=KCOL[k],
                    opacity=1.0,
                    flatshading=True,
                    name=f"layer k{k}",
                    showlegend=s == 0,
                    hovertext=txt,
                    hoverinfo="text",
                    visible=last,
                )
            )
        ex, ey, ez = [], [], []
        for c in range(len(t.cell_corners)):
            p = t.nodes[t.cell_corners[c]]
            for a, b in HEX_EDGES:
                ex += [p[a, 0], p[b, 0], None]
                ey += [p[a, 1], p[b, 1], None]
                ez += [p[a, 2], p[b, 2], None]
        fig.add_trace(
            go.Scatter3d(
                x=ex,
                y=ey,
                z=ez,
                mode="lines",
                line=dict(color="#555", width=2),
                hoverinfo="skip",
                showlegend=False,
                visible=last,
            )
        )
        fig.add_trace(
            go.Scatter3d(
                x=TRACE[:, 0],
                y=TRACE[:, 1],
                z=[0, 0],
                mode="lines",
                line=dict(color=EXPOSED, width=6, dash="dash"),
                name="trace",
                showlegend=s == 0,
                visible=last,
            )
        )
        by_key = {}
        for f, key in fault_pieces(st):
            by_key.setdefault(str(key), (key, []))[1].append(f)
        for name, (key, items) in sorted(by_key.items()):
            xs, ys, zs, ii, jj, kk, txt = [], [], [], [], [], [], []
            for f in items:
                poly, base = ordered(t.nodes[t.faces.nodes[f]]), len(xs)
                xs += list(poly[:, 0])
                ys += list(poly[:, 1])
                zs += list(poly[:, 2])
                txt += [f"{name}<br>area {st['areas'][f]:.1f} m²"] * len(poly)
                for m in range(1, len(poly) - 1):
                    ii.append(base)
                    jj.append(base + m)
                    kk.append(base + m + 1)
            label = (
                "exposed fault piece" if key == "exposed" else f"fixed k{key[0]} ↔ moving k{key[1]}"
            )
            fig.add_trace(
                go.Mesh3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    i=ii,
                    j=jj,
                    k=kk,
                    color=colour(key),
                    flatshading=True,
                    name=label,
                    hovertext=txt,
                    hoverinfo="text",
                    visible=last,
                )
            )
        groups.append(range(start, len(fig.data)))
    n = len(fig.data)
    buttons = [
        dict(
            method="update",
            label=lab,
            args=[{"visible": [m in g for m in range(n)]}, {"title.text": full}],
        )
        for g, (lab, full) in zip(groups, titles)
    ]
    opacity = [
        dict(method="restyle", label=lab, args=[{"opacity": op}, cells])
        for lab, op in (("cells see-through", 0.25), ("cells solid", 1.0))
    ]
    fig.update_layout(
        title=dict(text=titles[-1][1] + f"  (depth exaggerated {VEX:g}×)", y=0.97),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                buttons=buttons,
                active=len(states) - 1,
                x=0.0,
                y=1.0,
                xanchor="left",
                yanchor="bottom",
            ),
            dict(
                type="buttons",
                direction="right",
                buttons=opacity,
                active=1,
                x=1.0,
                y=1.0,
                xanchor="right",
                yanchor="bottom",
            ),
        ],
        scene=dict(
            xaxis_title="x [m]",
            yaxis_title="y [m]",
            zaxis_title="depth [m]",
            zaxis=dict(autorange="reversed"),
            aspectmode="manual",
            aspectratio=dict(
                x=1.6,
                y=1.6 * NY * DX / (NX * DX + LEAN),
                z=1.6 * VEX * (BASE + THROW) / (NX * DX + LEAN),
            ),
            camera=dict(eye=dict(x=0.7, y=-1.9, z=0.6)),
        ),
        legend=dict(x=0.0, y=0.0, bgcolor="rgba(255,255,255,0.85)", font=dict(size=10)),
        height=850,
        margin=dict(l=0, r=0, t=110, b=0),
    )
    path = OUT / "08_oblique_fault_3d.html"
    fig.write_html(path, include_plotlyjs=True, full_html=True)
    return path


if __name__ == "__main__":
    st = {
        n: load(f"{n}.grdecl")
        for n in ("A_regular_0_input", "A_regular_1_slid", "B_adapted_0_input", "B_adapted_1_slid")
    }
    for n, s in st.items():
        tot = {}
        for f, key in fault_pieces(s):
            tot[key] = tot.get(key, 0.0) + s["areas"][f]
        print(
            n,
            f"L={fault_length(s):.1f}",
            {str(k): round(v, 1) for k, v in sorted(tot.items(), key=str)},
        )
    print(fig_png(st["A_regular_1_slid"], st["B_adapted_1_slid"], st["B_adapted_0_input"]))
    print(
        fig_html(
            list(st.values()),
            [
                (
                    "A input",
                    f"A input: regular grid, staircase FAULTS dipping {DIP:.0f}° SE, no throw",
                ),
                (
                    f"A slid {THROW:g} m",
                    f"A slid: the whole hanging wall {THROW:g} m down its pillars "
                    "(parallel to the fault); COORD and FAULTS unchanged",
                ),
                ("B input", "B input: grid adapted to the trace, same FAULTS, no throw"),
                (
                    f"B slid {THROW:g} m",
                    "B slid: same slide on the adapted grid; contact area follows the "
                    "shorter fault",
                ),
            ],
        )
    )
