"""Gallery figures for fieldbook lessons 23-29: one PNG and one interactive HTML per lesson.

usage: python gallery.py DECK_DIR LESSON_FILE OUT_STEM "Figure title"
  DECK_DIR holds <case>_0_input.grdecl / <case>_1_slid.grdecl pairs (the lesson's --out).
  LESSON_FILE is the lesson source; its docstring bullets "- A name: text" give the subtitles.
  writes OUT_STEM.png (per case: plan map of the top shift + 3D view) and OUT_STEM_3d.html.
"""

import re
import sys
import textwrap
from pathlib import Path


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import TwoSlopeNorm, to_rgba
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from cpg2unstructured import build_topology, from_grdecl, read_grdecl
from fault_projections.from_faults import face_edge, faults_from_deck

KCOL = [
    "#9ecae1",
    "#a1d99b",
    "#fdd0a2",
    "#dadaeb",
    "#fcbba1",
    "#c7e9c0",
    "#fee391",
    "#bcbddc",
    "#9ecae1",
    "#a1d99b",
]
EXPOSED = "#d7191c"
SHADE = {0: 0.86, 1: 0.72, 2: 1.0}


def subtitles(lesson):
    doc = Path(lesson).read_text().split('"""')[1]
    out, key = {}, None
    for line in doc.splitlines():
        m = re.match(r"^- ([A-Z])[ _]?([^:]*):\s*(.*)", line.strip())
        if m:
            key = m.group(1)
            out[key] = (m.group(2).strip(), m.group(3).strip())
        elif key and line.startswith("  ") and line.strip():
            name, text = out[key]
            out[key] = (name, text + " " + line.strip())
        else:
            key = None
    return out


def load(deck_dir, name):
    g0 = read_grdecl(deck_dir / f"{name}_0_input.grdecl")
    g1 = read_grdecl(deck_dir / f"{name}_1_slid.grdecl")
    nx, ny, nz = g1.cart_dims
    faults = faults_from_deck(deck_dir / f"{name}_1_slid.grdecl") or {}
    try:
        topo = build_topology(from_grdecl(g1))
    except Exception as err:  # noqa: BLE001 - keep the plan panel
        print(f"  {name}: no 3D view ({type(err).__name__})")
        topo = None
    shape = (2 * nz, 2 * ny, 2 * nx)
    z0, z1 = g0.zcorn.reshape(shape), g1.zcorn.reshape(shape)
    coord = g1.coord.reshape(ny + 1, nx + 1, 6)
    act = (
        np.ones((nx, ny, nz), bool)
        if g1.actnum is None
        else np.asarray(g1.actnum).reshape(nz, ny, nx).transpose(2, 1, 0).astype(bool)
    )
    return dict(
        name=name,
        dims=(nx, ny, nz),
        shift=(z1 - z0),
        coord=coord,
        faults=faults,
        topo=topo,
        active=act,
    )


def column_shift(st, k):
    """Mean shift of the top of the first active cell in each column (k = 0 by default)."""
    nx, ny, nz = st["dims"]
    s = st["shift"]
    out = np.full((nx, ny), np.nan)
    for i, j in np.ndindex(nx, ny):
        ks = np.flatnonzero(st["active"][i, j])
        if ks.size:
            kk = ks[0] if k is None else k
            out[i, j] = s[2 * kk, 2 * j : 2 * j + 2, 2 * i : 2 * i + 2].mean()
    return out


def fault_segments(st, where="top"):
    c = st["coord"][..., :3] if where == "top" else st["coord"][..., 3:]
    segs = []
    for name, faces in st["faults"].items():
        for f in faces:
            (ia, ja), (ib, jb) = face_edge(f)
            segs.append((name, c[ja, ia], c[jb, ib]))
    return segs


def plan_panel(fig, ax, st, title):
    nx, ny, _ = st["dims"]
    top = st["coord"][..., :2]
    polys, vals = [], []
    sh = column_shift(st, None)
    for i, j in np.ndindex(nx, ny):
        if np.isnan(sh[i, j]):
            continue
        polys.append([top[j, i], top[j, i + 1], top[j + 1, i + 1], top[j + 1, i]])
        vals.append(sh[i, j])
    vals = np.array(vals)
    lim = max(float(np.abs(vals).max()), 1e-6)
    pc = PolyCollection(
        polys,
        array=vals,
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim),
        edgecolors="0.6",
        linewidths=0.25,
    )
    ax.add_collection(pc)
    for _, a, b in fault_segments(st):
        ax.plot([a[0], b[0]], [a[1], b[1]], color="k", lw=2.4, solid_capstyle="round")
    xy = top.reshape(-1, 2)
    ax.set_xlim(xy[:, 0].min(), xy[:, 0].max())
    ax.set_ylim(xy[:, 1].min(), xy[:, 1].max())
    ax.set_aspect("equal")
    ax.tick_params(labelsize=7)
    ax.set_title(title, fontsize=9, loc="left")
    cb = fig.colorbar(pc, ax=ax, shrink=0.8, pad=0.02)
    cb.ax.tick_params(labelsize=7)
    cb.set_label("top shift [m] (down +)", fontsize=7.5)


def ordered(p):
    c = p.mean(0)
    _, _, vt = np.linalg.svd(p - c)
    uv = (p - c) @ vt[:2].T
    return p[np.argsort(np.arctan2(uv[:, 1], uv[:, 0]))]


def boundary_polys(st):
    """Boundary faces of the slid topology: (polygon, rgba, is_fault_face, k)."""
    t = st["topo"]
    coord = st["coord"]
    pillars = []
    for name, faces in st["faults"].items():
        for f in faces:
            (ia, ja), (ib, jb) = face_edge(f)
            pillars.append((coord[ja, ia], coord[jb, ib]))
    span = np.ptp(coord[..., :2].reshape(-1, 2), axis=0).max()
    tol = 1e-4 * span
    out = []
    for f in np.flatnonzero(t.faces.neighbours[:, 1] < 0):
        a = t.faces.neighbours[f, 0]
        p = np.unique(np.round(t.nodes[t.faces.nodes[f]], 9), axis=0)
        if len(p) < 3:
            continue
        d = int(t.faces.direction[f])
        on_fault = False
        if d <= 1 and pillars:
            c = p.mean(0)
            for pa, pb in pillars:
                qa = pa[:3] + (pa[3:] - pa[:3]) * (c[2] - pa[2]) / (pa[5] - pa[2])
                qb = pb[:3] + (pb[3:] - pb[:3]) * (c[2] - pb[2]) / (pb[5] - pb[2])
                e = qb[:2] - qa[:2]
                u = np.dot(c[:2] - qa[:2], e) / np.dot(e, e)
                if -1e-6 <= u <= 1 + 1e-6 and np.linalg.norm(qa[:2] + u * e - c[:2]) < tol:
                    on_fault = True
                    break
        k = int(t.original_ijk[a][2])
        rgb = np.array(to_rgba(EXPOSED if on_fault else KCOL[k % len(KCOL)]))
        rgb[:3] *= SHADE[d]
        q = ordered(p)
        normal = np.zeros(3)
        if on_fault:
            area = 0.5 * np.cross(q[1:-1] - q[0], q[2:] - q[0]).sum(0)
            out_dir = q.mean(0) - t.nodes[t.cell_corners[a]].mean(0)
            normal = area if np.dot(area, out_dir) > 0 else -area
            zs = t.nodes[:, 2]
            normal = normal * (0.2 + (zs.max() - q[:, 2].mean()) / np.ptp(zs))  # favour high faces
        out.append((q, rgb, on_fault, k, normal))
    return out


def best_azimuth(polys):
    """The view (from the south-east or south-west) that faces most exposed fault area."""
    best, score = -58.0, -1.0
    for az in (-58.0, -122.0, 58.0, 122.0):
        cam = np.array([np.cos(np.radians(az)), np.sin(np.radians(az)), 0.0])
        sc = sum(max(0.0, float(np.dot(n, cam))) for *_, n in polys)
        if sc > score + 1e-9:
            best, score = az, sc
    return best


def exaggeration(st):
    nodes = st["topo"].nodes
    xr, yr = np.ptp(nodes[:, 0]), np.ptp(nodes[:, 1])
    zr = max(np.ptp(nodes[:, 2]), 1e-6)
    raw = 0.3 * max(xr, yr) / zr
    nice = [1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]
    return min(nice, key=lambda v: abs(np.log(v / raw))), xr, yr, zr


def view3d(ax, st, title):
    polys = boundary_polys(st)
    ax.add_collection3d(
        Poly3DCollection(
            [p for p, *_ in polys],
            facecolors=[c for _, c, *_ in polys],
            edgecolors=to_rgba("0.25", 0.45),
            linewidths=0.2,
        )
    )
    nodes = st["topo"].nodes
    vex, xr, yr, zr = exaggeration(st)
    ax.set_xlim(nodes[:, 0].min(), nodes[:, 0].max())
    ax.set_ylim(nodes[:, 1].min(), nodes[:, 1].max())
    ax.set_zlim(nodes[:, 2].max(), nodes[:, 2].min())
    ax.set_box_aspect((xr, yr, vex * zr), zoom=1.4)
    ax.view_init(elev=26, azim=best_azimuth(polys))
    ax.tick_params(labelsize=6)
    ax.set_xlabel("x", fontsize=7, labelpad=-4)
    ax.set_ylabel("y", fontsize=7, labelpad=-4)
    ax.set_zlabel("depth", fontsize=7, labelpad=-2)
    n_fault = sum(1 for _, _, f, _, _ in polys if f)
    ax.set_title(
        f"{title.splitlines()[0]}  ·  3D, depth ×{vex:g}; red = exposed fault face ({n_fault})",
        fontsize=8.5,
    )


def figure_png(states, subs, out_stem, title):
    n = len(states)
    cols = min(3, n)
    rows = int(np.ceil(n / cols))
    fig = plt.figure(figsize=(6.4 * cols, 8.6 * rows))
    gs = fig.add_gridspec(
        2 * rows,
        cols,
        height_ratios=[1, 1.15] * rows,
        left=0.04,
        right=0.98,
        top=1 - 0.135 / rows,
        bottom=0.01,
        hspace=0.18,
        wspace=0.18,
    )
    for m, st in enumerate(states):
        r, c = divmod(m, cols)
        letter = st["name"][0]
        name, text = subs.get(letter, (st["name"][2:].replace("_", " "), ""))
        if sum(s["name"][0] == letter for s in states) > 1:  # two decks for one case
            name = st["name"][2:].replace("_", " ")
        head = f"{letter}  {name or st['name'][2:].replace('_', ' ')}"
        lines = textwrap.wrap(text.rstrip(";. "), 66)
        body = "\n".join(lines[:4]) + (" ..." if len(lines) > 4 else "")
        plan_panel(fig, fig.add_subplot(gs[2 * r, c]), st, head + ("\n" + body if body else ""))
        if st["topo"] is not None:
            view3d(fig.add_subplot(gs[2 * r + 1, c], projection="3d"), st, head)
    fig.suptitle(
        title + "\nplan: shift of the top of each column, FAULTS in black (COORD never "
        "changes)  ·  3D: slid grid, red = fault face left exposed by the slide",
        fontsize=13,
        fontweight="bold",
        y=0.997,
    )
    path = Path(f"{out_stem}.png")
    path.unlink(missing_ok=True)
    fig.savefig(path, dpi=90)
    plt.close(fig)
    return path


def figure_html(states, subs, out_stem, title):
    import plotly.graph_objects as go

    fig, groups, labels = go.Figure(), [], []
    for s, st in enumerate(states):
        if st["topo"] is None:
            continue
        start = len(fig.data)
        polys = boundary_polys(st)
        vex, xr, yr, zr = exaggeration(st)
        keys = sorted({("fault" if f else k) for _, _, f, k, _ in polys}, key=str)
        for key in keys:
            xs, ys, zs, ii, jj, kk = [], [], [], [], [], []
            for p, _, f, k, _ in polys:
                if ("fault" if f else k) != key:
                    continue
                base = len(xs)
                xs += list(p[:, 0])
                ys += list(p[:, 1])
                zs += list(p[:, 2])
                for q in range(1, len(p) - 1):
                    ii.append(base)
                    jj.append(base + q)
                    kk.append(base + q + 1)
            color = EXPOSED if key == "fault" else KCOL[key % len(KCOL)]
            fig.add_trace(
                go.Mesh3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    i=ii,
                    j=jj,
                    k=kk,
                    color=color,
                    flatshading=True,
                    visible=s == 0,
                    showlegend=False,
                    name="exposed fault" if key == "fault" else f"layer k{key}",
                )
            )
        groups.append((range(start, len(fig.data)), vex, xr, yr, zr, best_azimuth(polys)))
        letter = st["name"][0]
        labels.append((letter, f"{letter}  {subs.get(letter, (st['name'][2:], ''))[0]}"))
    n = len(fig.data)

    def scene(vex, xr, yr, zr, az):
        big = max(xr, yr)
        eye = dict(x=1.9 * np.cos(np.radians(az)), y=1.9 * np.sin(np.radians(az)), z=0.9)
        return dict(
            xaxis_title="x [m]",
            yaxis_title="y [m]",
            zaxis_title="depth [m]",
            zaxis=dict(autorange="reversed"),
            aspectmode="manual",
            aspectratio=dict(x=1.6 * xr / big, y=1.6 * yr / big, z=1.6 * vex * zr / big),
            camera=dict(eye=eye),
        )

    buttons = [
        dict(
            method="update",
            label=letter,
            args=[
                {"visible": [m in g for m in range(n)]},
                {
                    "title.text": f"{text}  (depth ×{vex:g}; red = exposed fault face)",
                    "scene": scene(vex, xr, yr, zr, az),
                },
            ],
        )
        for (g, vex, xr, yr, zr, az), (letter, text) in zip(groups, labels)
    ]
    g0 = groups[0]
    fig.update_layout(
        title=dict(
            text=f"{title}<br>{labels[0][1]}  (depth ×{g0[1]:g}; red = exposed fault face)", y=0.97
        ),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                buttons=buttons,
                active=0,
                x=0.0,
                y=1.0,
                xanchor="left",
                yanchor="bottom",
            )
        ],
        scene=scene(*g0[1:]),
        height=850,
        margin=dict(l=0, r=0, t=110, b=0),
    )
    path = Path(f"{out_stem}_3d.html")
    path.unlink(missing_ok=True)
    fig.write_html(path, include_plotlyjs=True, full_html=True)
    return path


if __name__ == "__main__":
    deck_dir, lesson, out_stem, title = Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
    subs = subtitles(lesson)
    names = sorted(p.name[: -len("_1_slid.grdecl")] for p in deck_dir.glob("*_1_slid.grdecl"))
    states = [load(deck_dir, name) for name in names]
    print(figure_png(states, subs, out_stem, title))
    print(figure_html(states, subs, out_stem, title))
