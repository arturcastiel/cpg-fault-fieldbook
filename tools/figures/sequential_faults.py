"""Example 06: two sequential (stepped) normal faults, slid one after the other by ZCORN only.

usage: python sequential_faults.py OUTDIR
  OUTDIR/06_two_sequential_faults/{0_input,1_after_F1,2_after_F2}.grdecl + checks.txt
  OUTDIR/06_two_sequential_faults.png         2D cross-sections of the three stages
  OUTDIR/06_two_sequential_faults_3d.png      3D, final stage, true positions
  OUTDIR/06_two_sequential_faults_3d.html     interactive 3D, stage buttons

Deck: 9x4x3, DX=100 m, 10 m layers, both fault pillar lines tilted (dip east).
  F1 on pillar line i=3, F2 on pillar line i=6  ->  blocks W (i<3) | M (3..5) | E (i>=6)
  step 1: slide F1 by 7 m   (hanging wall = everything east of F1 = M + E)
  step 2: slide F2 by 12 m  (hanging wall = E)          ->  W 0, M 7, E 19 m
Helpers below are shaped like the planned cpg2python P1 API (zcorn_ops).
"""

import dataclasses
import sys
from collections import deque
from pathlib import Path


import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from cpg2unstructured import (
    CornerPointGrid,
    build_topology,
    from_grdecl,
    read_grdecl,
    write_grid_include,
)
from cpg2unstructured.cpg import _assemble_corner_xyz
from cpg2unstructured.geometry import cell_geometry, face_geometry
from fault_projections.faults import faults_text
from fault_projections.from_faults import edge_face, faults_from_deck

NX, NY, NZ = 9, 4, 3
DX, DZ = 100.0, 10.0
PILLAR_DEPTH = 80.0
TILT = 40.0  # fault pillar x runs F*DX -20 (top) -> +20 (bottom)
FAULT_I = {"F1": 3, "F2": 6}
THROW = {"F1": 7.0, "F2": 12.0}
OUT = Path(sys.argv[1])
DECKS = OUT / "06_two_sequential_faults"
KCOL = ["#9ecae1", "#a1d99b", "#fdd0a2"]
EXPOSED = "#d7191c"
PAIR_ORDER = [(0, 0), (1, 1), (2, 2), (1, 0), (2, 1), (2, 0), (0, 1), (1, 2), (0, 2)]
PAIR_COLS = [
    "#1f77b4",
    "#2ca02c",
    "#ff7f0e",
    "#9467bd",
    "#8c564b",
    "#17becf",
    "#e377c2",
    "#bcbd22",
    "#7f7f7f",
]
PAIRCOL = dict(zip(PAIR_ORDER, PAIR_COLS))


def BLOCK(i):
    return "W" if i < FAULT_I["F1"] else ("M" if i < FAULT_I["F2"] else "E")


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


# ── P1-shaped ZCORN helpers ───────────────────────────────────────────────────


def zcorn_3d(g):
    return np.asarray(g.zcorn, float).reshape(2 * NZ, 2 * NY, 2 * NX).transpose(2, 1, 0).copy()


def zcorn_flat(zc):
    return zc.transpose(2, 1, 0).ravel().copy()


def column_slot_mask(columns):
    return np.repeat(np.repeat(columns, 2, axis=0), 2, axis=1)


def apply_zcorn_displacement(zc, dz, columns):
    out = zc.copy()
    out[column_slot_mask(columns)] += dz
    return out


def with_zcorn(g, zc):
    return dataclasses.replace(
        g, zcorn=zcorn_flat(zc), props=dict(g.props), unrecognized=list(g.unrecognized)
    )


# ── deck + fault semantics (P2-shaped) ────────────────────────────────────────


def make_coord():
    rows = []
    for j in range(NY + 1):
        for i in range(NX + 1):
            x = i * DX
            dx = TILT / 2 if i in FAULT_I.values() else 0.0
            rows.append((x - dx, j * DX, 0.0, x + dx, j * DX, PILLAR_DEPTH))
    return np.asarray(rows, float)


def make_zc():
    zc = np.empty((2 * NX, 2 * NY, 2 * NZ))
    for k in range(NZ):
        zc[:, :, 2 * k], zc[:, :, 2 * k + 1] = k * DZ, (k + 1) * DZ
    return zc


def hanging_wall(faces):
    """(NX, NY) columns east of the fault: flood fill with this fault's faces cut."""
    cut = set()
    for f in faces:
        a, b = (
            ((f.i - 1, f.j - 1), (f.i, f.j - 1))
            if f.face == "X"
            else ((f.i - 1, f.j - 1), (f.i - 1, f.j))
        )
        cut.add(frozenset((a, b)))
    seed, other = (faces[0].i, faces[0].j - 1), (faces[0].i - 1, faces[0].j - 1)
    mask = np.zeros((NX, NY), bool)
    mask[seed] = True
    q = deque([seed])
    while q:
        i, j = q.popleft()
        for n in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
            if (
                0 <= n[0] < NX
                and 0 <= n[1] < NY
                and not mask[n]
                and frozenset(((i, j), n)) not in cut
            ):
                mask[n] = True
                q.append(n)
    assert not mask[other], "fault does not close a block"
    return mask


def write_deck(path, g, faults, header):
    write_grid_include(g, path, header=header)
    with open(path, "a") as fh:
        fh.write("\n" + faults_text(list(faults.values()), NZ, names=list(faults)))
    return path


# ── re-processing + measurements ──────────────────────────────────────────────


def process(path):
    g = read_grdecl(path)
    t = build_topology(from_grdecl(g))
    fg = face_geometry(t.nodes, t.faces, t.cell_corners)
    cg = cell_geometry(t.nodes, t.faces, t.cell_faces, fg)
    return dict(g=g, t=t, fg=fg, vol=cg.volumes, xyz=_assemble_corner_xyz(g), zc=zcorn_3d(g))


def fault_faces(st, F):
    """[(face, kw|None, ke|None)] on the fault plane between columns F-1 | F."""
    t = st["t"]
    ijk, nb = t.original_ijk, t.faces.neighbours
    out = []
    for f in np.where(t.faces.direction == 0)[0]:
        a, b = nb[f]
        ia = ijk[a][0]
        if b >= 0:
            cw, ce = sorted([ijk[a], ijk[b]], key=lambda v: v[0])
            if cw[0] == F - 1 and ce[0] == F:
                out.append((f, int(cw[2]), int(ce[2])))
        elif ia in (F - 1, F):
            out.append(
                (f, int(ijk[a][2]) if ia == F - 1 else None, int(ijk[a][2]) if ia == F else None)
            )
    return out


def expected_pairs(throw):
    pairs = set()
    for kw in range(NZ):
        for ke in range(NZ):
            if min((kw + 1) * DZ, (ke + 1) * DZ + throw) - max(kw * DZ, ke * DZ + throw) > 1e-9:
                pairs.add((kw, ke))
    return pairs


def throw_on(zc, F):
    d = zc[2 * F] - zc[2 * F - 1]
    return d.min(), d.max()


def pair_label(F, kw, ke):
    w, e = ("W", "M") if F == FAULT_I["F1"] else ("M", "E")
    if kw is not None and ke is not None:
        return f"{w} k{kw} ↔ {e} k{ke}"
    return f"{w} k{kw}: exposed" if kw is not None else f"{e} k{ke}: exposed"


# ── 2D figure ─────────────────────────────────────────────────────────────────


def front_segment(st, f):
    p = np.unique(st["t"].nodes[st["t"].faces.nodes[f]], axis=0)
    sel = p[np.isclose(p[:, 1], p[:, 1].min())]
    return sel[np.argsort(sel[:, 2])][:, [0, 2]]


def section(st, c):
    return st["t"].nodes[st["t"].cell_corners[c][[0, 1, 5, 4]]][:, [0, 2]]


def fig_2d(states, titles, moved):
    fig, axes = plt.subplots(3, 1, figsize=(15, 13.5))
    for n, (ax, st, title) in enumerate(zip(axes, states, titles)):
        t = st["t"]
        front = np.where(t.original_ijk[:, 1] == 0)[0]
        if n > 0:
            prev = states[n - 1]
            hw = moved[n - 1]
            ghost = [section(prev, c) for c in front if hw[t.original_ijk[c, 0], 0]]
            ax.add_collection(
                PolyCollection(
                    ghost,
                    facecolors="none",
                    edgecolors="0.4",
                    linestyles="--",
                    linewidths=0.8,
                    zorder=3,
                )
            )
        polys = [section(st, c) for c in front]
        ax.add_collection(
            PolyCollection(
                polys,
                facecolors=[KCOL[t.original_ijk[c, 2]] for c in front],
                edgecolors="0.35",
                linewidths=0.8,
            )
        )
        for c, p in zip(front, polys):
            i, _, k = t.original_ijk[c]
            ax.text(
                *p.mean(0), f"{BLOCK(i)} k{k}", ha="center", va="center", fontsize=7, color="0.3"
            )
        pp = st["g"].coord.reshape(NY + 1, NX + 1, 6)[0]
        for p in pp:
            s = np.array([-0.05, 0.75])
            ax.plot(
                p[0] + s * (p[3] - p[0]),
                p[2] + s * (p[5] - p[2]),
                color="0.55",
                lw=0.8,
                ls=":",
                zorder=0,
            )
        lines = []
        for name, F in FAULT_I.items():
            for f, kw, ke in fault_faces(st, F):
                seg = front_segment(st, f)
                if len(seg) < 2:
                    continue
                col = PAIRCOL[(kw, ke)] if kw is not None and ke is not None else EXPOSED
                ax.plot(*seg.T, color=col, lw=5, solid_capstyle="butt", zorder=4)
            lo, hi = throw_on(st["zc"], F)
            top = pp[F]
            ax.text(
                top[0],
                -3,
                f"{name}  throw {hi:.0f} m",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )
            contacts = {}
            for f, kw, ke in fault_faces(st, F):
                contacts[pair_label(F, kw, ke)] = (
                    contacts.get(pair_label(F, kw, ke), 0.0) + st["fg"].areas[f]
                )
            lines.append(
                f"{name}: " + ", ".join(f"{k} {v:,.0f}" for k, v in sorted(contacts.items()))
            )
        if n > 0:
            prev = states[n - 1]
            pillars = (3, 5, 6) if n == 1 else (6, 8, 9)
            for P in pillars:
                a = 2 * P if P < NX else 2 * P - 1
                p0, p1 = prev["xyz"][a, 0, 0][[0, 2]], st["xyz"][a, 0, 0][[0, 2]]
                on_fault = P in FAULT_I.values()
                ax.annotate(
                    "",
                    p1,
                    p0,
                    zorder=6,
                    arrowprops=dict(
                        arrowstyle="-|>,head_length=0.6,head_width=0.35",
                        lw=2.2,
                        shrinkA=0,
                        shrinkB=0,
                        color="#6a3d9a" if on_fault else "k",
                        path_effects=[pe.withStroke(linewidth=4.5, foreground="white")],
                    ),
                )
        ax.text(
            0.995,
            0.03,
            "fault contact / exposed area [m²], all rows j:\n" + "\n".join(lines),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7.5,
            bbox=dict(fc="white", ec="0.6", alpha=0.95),
            zorder=7,
        )
        disp = [st["zc"][2 * i, 0, 0] - states[0]["zc"][2 * i, 0, 0] for i in (0, 4, 7)]
        for x, (b, d) in zip((150, 450, 750), zip("WME", disp)):
            ax.text(x, -9, f"block {b}: {d:+.0f} m", ha="center", fontsize=9, color="0.25")
        ax.set_xlim(-10, NX * DX + 10)
        ax.set_ylim(65, -14)
        ax.set_aspect(4.0)
        ax.set_ylabel("depth z [m]  (vertical ×4)")
        ax.set_title(title, fontsize=11, loc="left")
    axes[-1].set_xlabel("x [m]  (cross-section j = 0)")
    seen = sorted(
        {
            (kw, ke)
            for st in states
            for F in FAULT_I.values()
            for _, kw, ke in fault_faces(st, F)
            if kw is not None and ke is not None
        },
        key=PAIR_ORDER.index,
    )
    handles = [Patch(fc=KCOL[k], ec="0.35", label=f"layer k{k}") for k in range(3)]
    handles += [
        Line2D([], [], color=PAIRCOL[p], lw=6, label=f"fault contact k{p[0]} ↔ k{p[1]}")
        for p in seen
    ]
    handles += [
        Line2D([], [], color=EXPOSED, lw=6, label="exposed fault piece → boundary"),
        Line2D([], [], color="0.4", ls="--", label="moved block, before this step"),
        Line2D([], [], color="#6a3d9a", lw=1.6, label="corner motion on a fault pillar (along it)"),
        Line2D([], [], color="k", lw=1.6, label="corner motion, interior (straight down)"),
        Line2D([], [], color="0.55", ls=":", label="pillar line (COORD, never changes)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8.5, frameon=False)
    fig.suptitle(
        "Two sequential faults: slide F1, then F2. Each slide is a ZCORN-only edit of one "
        "hanging wall; the deck is re-processed after each step",
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    path = OUT / "06_two_sequential_faults.png"
    fig.savefig(path, dpi=125)
    plt.close(fig)
    return path


# ── 3D figures ────────────────────────────────────────────────────────────────


def ordered(pts):
    pts = np.unique(np.round(pts, 9), axis=0)
    c = pts.mean(0)
    _, _, vt = np.linalg.svd(pts - c)
    uv = (pts - c) @ vt[:2].T
    return pts[np.argsort(np.arctan2(uv[:, 1], uv[:, 0]))]


def fault_polys(st):
    out = []
    for F in FAULT_I.values():
        for f, kw, ke in fault_faces(st, F):
            col = PAIRCOL[(kw, ke)] if kw is not None and ke is not None else EXPOSED
            out.append(
                (
                    ordered(st["t"].nodes[st["t"].faces.nodes[f]]),
                    col,
                    pair_label(F, kw, ke),
                    float(st["fg"].areas[f]),
                )
            )
    return out


def skip_face(i, n):
    return any((i == F - 1 and n == 1) or (i == F and n == 0) for F in FAULT_I.values())


def fig_3d_png(st):
    from matplotlib.colors import to_rgba

    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(projection="3d")
    t = st["t"]
    polys, fc, ec = [], [], []
    for c in range(len(t.cell_corners)):
        i, _, k = t.original_ijk[c]
        p = t.nodes[t.cell_corners[c]]
        solid = BLOCK(i) == "W"
        for n, q in enumerate(HEX_FACES):
            if skip_face(i, n):
                continue
            polys.append(p[q])
            fc.append(to_rgba(KCOL[k], 1.0 if solid else 0.10))
            ec.append(to_rgba("0.35", 1.0 if solid else 0.5))
    seen = {}
    for poly, col, lab, _ in fault_polys(st):
        polys.append(poly)
        fc.append(to_rgba(col))
        ec.append(to_rgba("k"))
        seen[
            (
                "exposed fault piece → boundary"
                if col == EXPOSED
                else "contact "
                + lab.split(": ")[0].replace("W ", "").replace("M ", "").replace("E ", "")
            )
        ] = col
    ax.computed_zorder = False
    ax.add_collection3d(Poly3DCollection(polys, facecolors=fc, edgecolors=ec, linewidths=0.4))
    pp = st["g"].coord.reshape(NY + 1, NX + 1, 6)
    segs = [
        [pp[j, F, :3], pp[j, F, :3] + 50.0 / PILLAR_DEPTH * (pp[j, F, 3:] - pp[j, F, :3])]
        for F in FAULT_I.values()
        for j in range(NY + 1)
    ]
    ax.add_collection3d(Line3DCollection(segs, colors="k", linestyles="--", linewidths=0.9))
    ax.set_xlim(0, NX * DX)
    ax.set_ylim(0, NY * DX)
    ax.set_zlim(50, 0)
    ax.set_box_aspect((2.2, 1.0, 0.55))
    ax.view_init(elev=22, azim=-40)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("depth [m]")
    handles = [Patch(fc=KCOL[k], ec="0.3", label=f"layer k{k}") for k in range(3)]
    handles += [Patch(fc=c, ec="k", label=lab) for lab, c in seen.items()]
    ax.legend(handles=handles, loc="upper left", fontsize=8, ncol=2)
    ax.set_title(
        "Two sequential faults after both slides (true positions). W block solid, M and E "
        "see-through;\nfault surfaces coloured by which layers touch (W|M on F1, M|E on F2)",
        fontsize=11,
    )
    fig.subplots_adjust(left=0, right=1, top=0.92, bottom=0)
    path = OUT / "06_two_sequential_faults_3d.png"
    fig.savefig(path, dpi=125)
    plt.close(fig)
    return path


def fig_3d_html(states, titles):
    import plotly.graph_objects as go

    fig = go.Figure()
    groups, cells = [], []
    for s, st in enumerate(states):
        t = st["t"]
        start = len(fig.data)
        for k in range(NZ):
            xs, ys, zs, ii, jj, kk, txt = [], [], [], [], [], [], []
            for c in np.where(t.original_ijk[:, 2] == k)[0]:
                i, j, _ = t.original_ijk[c]
                p = t.nodes[t.cell_corners[c]]
                base = len(xs)
                xs += list(p[:, 0])
                ys += list(p[:, 1])
                zs += list(p[:, 2])
                txt += [f"block {BLOCK(i)} cell ({i},{j},{k})"] * 8
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
                    flatshading=True,
                    opacity=0.25,
                    name=f"layer k{k}",
                    legendgroup=f"k{k}",
                    showlegend=s == 0,
                    hovertext=txt,
                    hoverinfo="text",
                    visible=s == len(states) - 1,
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
                visible=s == len(states) - 1,
            )
        )
        by_lab = {}
        for poly, col, lab, area in fault_polys(st):
            by_lab.setdefault((lab, col), []).append((poly, area))
        for (lab, col), items in sorted(by_lab.items()):
            xs, ys, zs, ii, jj, kk, txt = [], [], [], [], [], [], []
            for poly, area in items:
                base = len(xs)
                xs += list(poly[:, 0])
                ys += list(poly[:, 1])
                zs += list(poly[:, 2])
                txt += [f"{lab}<br>face area {area:.1f} m²"] * len(poly)
                for m in range(1, len(poly) - 1):
                    ii.append(base)
                    jj.append(base + m)
                    kk.append(base + m + 1)
            fig.add_trace(
                go.Mesh3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    i=ii,
                    j=jj,
                    k=kk,
                    color=col,
                    flatshading=True,
                    name=lab,
                    hovertext=txt,
                    hoverinfo="text",
                    visible=s == len(states) - 1,
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
        title=dict(text=titles[-1][1], y=0.97),
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
                active=0,
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
            aspectratio=dict(x=2.2, y=1.0, z=0.6),
            camera=dict(eye=dict(x=0.9, y=-2.0, z=0.9)),
        ),
        legend=dict(x=0.0, y=0.0, bgcolor="rgba(255,255,255,0.85)", font=dict(size=10)),
        height=850,
        margin=dict(l=0, r=0, t=110, b=0),
    )
    path = OUT / "06_two_sequential_faults_3d.html"
    fig.write_html(path, include_plotlyjs=True, full_html=True)
    return path


# ── main ──────────────────────────────────────────────────────────────────────


def main():
    DECKS.mkdir(parents=True, exist_ok=True)
    faults = {
        name: [edge_face((F, j), (F, j + 1)) for j in range(NY)] for name, F in FAULT_I.items()
    }
    coord, zc0 = make_coord(), make_zc()
    g0 = CornerPointGrid(
        cart_dims=(NX, NY, NZ),
        coord=coord,
        zcorn=zcorn_flat(zc0),
        actnum=np.ones(NX * NY * NZ, dtype=np.int8),
    )
    paths = [write_deck(DECKS / "0_input.grdecl", g0, faults, "example 06: input")]

    parsed = faults_from_deck(paths[0])
    checks = {"FAULTS read back from the deck (inline, both faults)": parsed == faults}
    g, zc, moved = read_grdecl(paths[0]), zcorn_3d(read_grdecl(paths[0])), []
    for step, name in enumerate(("F1", "F2"), start=1):
        hw = hanging_wall(parsed[name])
        moved.append(hw)
        zc = apply_zcorn_displacement(zc, THROW[name], hw)
        g = with_zcorn(g, zc)
        paths.append(
            write_deck(
                DECKS / f"{step}_after_{name}.grdecl",
                g,
                parsed,
                f"example 06: after sliding {name} by {THROW[name]:g} m",
            )
        )

    states = [process(p) for p in paths]
    zcs = [s["zc"] for s in states]
    checks["COORD bit-identical in all three decks"] = all(
        np.array_equal(s["g"].coord, states[0]["g"].coord) for s in states
    )
    checks["FAULTS unchanged in all three decks"] = all(
        faults_from_deck(p) == faults for p in paths
    )
    want = {0: (0, 0), 1: (7, 0), 2: (7, 12)}
    for s, zc in enumerate(zcs):
        got = tuple(throw_on(zc, F) for F in FAULT_I.values())
        checks[
            f"stage {s}: throw F1 {got[0][0]:g}..{got[0][1]:g}, F2 {got[1][0]:g}..{got[1][1]:g} m "
            f"(want {want[s][0]}, {want[s][1]})"
        ] = all(np.isclose(lo, w) and np.isclose(hi, w) for (lo, hi), w in zip(got, want[s]))
    blocks = np.array([[{"W": 0.0, "M": 7.0, "E": 19.0}[BLOCK(i)]] * NY for i in range(NX)])
    one_shot = zc0 + np.repeat(np.repeat(blocks, 2, 0), 2, 1)[..., None]
    rev = apply_zcorn_displacement(
        apply_zcorn_displacement(zc0, THROW["F2"], moved[1]), THROW["F1"], moved[0]
    )
    checks["order does not matter: F2 then F1 == F1 then F2 (bitwise)"] = np.array_equal(
        rev, zcs[2]
    )
    checks["sequential == one-shot block map W 0 / M 7 / E 19 m"] = np.allclose(one_shot, zcs[2])
    checks["sliding F1 carries F2 along (no throw created on F2)"] = np.allclose(
        zcs[1][2 * FAULT_I["F2"]], zcs[1][2 * FAULT_I["F2"] - 1]
    )

    def thick(zc):
        return zc[..., 1::2] - zc[..., 0::2]

    checks["every cell thickness preserved"] = all(
        np.array_equal(thick(z), thick(zc0)) for z in zcs
    )
    for s in (1, 2):
        for name, F in FAULT_I.items():
            th = throw_on(zcs[s], F)[1]
            got = sorted(
                {
                    (kw, ke)
                    for _, kw, ke in fault_faces(states[s], F)
                    if kw is not None and ke is not None
                }
            )
            checks[f"stage {s} {name}: re-processed contacts {got} == interval overlap"] = set(
                got
            ) == expected_pairs(th)

    ii = states[0]["t"].original_ijk[:, 0]
    v0 = states[0]["vol"]
    expect_cols = {1: {3, 5, 6}, 2: {3, 5, 6}}  # columns bounded by a moved tilted pillar
    for s in (1, 2):
        changed = {int(i) for i in np.unique(ii[~np.isclose(states[s]["vol"], v0, rtol=1e-12)])}
        checks[
            f"stage {s}: volume change confined to columns {sorted(changed)} "
            f"(bounded by a moved tilted pillar)"
        ] = changed == expect_cols[s]
    lines = ["Example 06: two sequential faults; checks", ""]
    lines += [f"[{'PASS' if v else 'FAIL'}] {k}" for k, v in checks.items()]
    lines += ["", "cell volume change vs input, per column i (min..max over j,k), %:"]
    v0 = states[0]["vol"]
    for s in (1, 2):
        dv = 100 * (states[s]["vol"] / v0 - 1)
        ii = states[0]["t"].original_ijk[:, 0]
        row = "  ".join(
            f"i{i}:{dv[ii == i].min():+.1f}..{dv[ii == i].max():+.1f}"
            if np.abs(dv[ii == i]).max() > 1e-9
            else f"i{i}:0"
            for i in range(NX)
        )
        lines.append(f"  stage {s}: {row}")
    for s in (1, 2):
        lines += ["", f"stage {s} fault contacts / exposed (area m², all rows):"]
        for name, F in FAULT_I.items():
            acc = {}
            for f, kw, ke in fault_faces(states[s], F):
                acc[pair_label(F, kw, ke)] = (
                    acc.get(pair_label(F, kw, ke), 0.0) + states[s]["fg"].areas[f]
                )
            lines.append(f"  {name}: " + ", ".join(f"{k} {v:,.0f}" for k, v in sorted(acc.items())))
    text = "\n".join(lines)
    (DECKS / "checks.txt").write_text(text + "\n")
    print(text)

    titles = [
        ("input", "Stage 0: input deck, two tilted faults, no throw"),
        (
            "after F1 (7 m)",
            "Stage 1: F1 slid 7 m, so M and E (its hanging wall) move down along their pillars",
        ),
        ("after F2 (12 m)", "Stage 2: F2 slid 12 m, so E moves a further 12 m (total 19 m)"),
    ]
    print(fig_2d(states, [t[1] for t in titles], moved))
    print(fig_3d_png(states[2]))
    print(fig_3d_html(states, titles))
    return all(checks.values())


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
