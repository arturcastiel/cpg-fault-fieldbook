"""Example 07 figures: a fault that ends inside the grid, slid with a tapered throw.

usage: python partial_fault.py DIR
  DIR/07_partial_fault/{input,tapered_8m}.grdecl come from fieldbook lesson 20 (--out).
  writes DIR/07_partial_fault.png, DIR/07_partial_fault_3d.html
"""

import sys
from pathlib import Path


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from cpg2unstructured import build_topology, from_grdecl, read_grdecl
from cpg2unstructured.geometry import face_geometry

NX, NY, NZ, DX = 6, 6, 3, 100.0
LINE, TIP = 3, 4
OUT = Path(sys.argv[1])
DECKS = OUT / "07_partial_fault"
KCOL = ["#9ecae1", "#a1d99b", "#fdd0a2"]
EXPOSED, WELDED = "#d7191c", "#d1d1d1"
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


def load(name):
    g = read_grdecl(DECKS / name)
    t = build_topology(from_grdecl(g))
    fg = face_geometry(t.nodes, t.faces, t.cell_corners)
    return dict(g=g, t=t, areas=fg.areas, z=g.zcorn.reshape(2 * NZ, 2 * NY, 2 * NX))


def line_faces(st):
    """[(face, row j, key)] on pillar line LINE; key (kw, ke) | 'exposed' | 'welded'."""
    t, ijk, out = st["t"], st["t"].original_ijk, []
    for f in np.flatnonzero(t.faces.direction == 0):
        a, b = t.faces.neighbours[f]
        cells = sorted([ijk[a]] + ([ijk[b]] if b >= 0 else []), key=lambda c: c[0])
        j = int(cells[0][1])
        if b >= 0 and (cells[0][0], cells[1][0]) == (LINE - 1, LINE):
            key = "welded" if j >= TIP else (int(cells[0][2]), int(cells[1][2]))
        elif b < 0 and cells[0][0] in (LINE - 1, LINE):
            key = "exposed"
        else:
            continue
        out.append((f, j, key))
    return out


def ordered(pts, axes=None):
    pts = np.unique(np.round(pts, 9), axis=0)
    c = pts.mean(0)
    if axes is None:
        _, _, vt = np.linalg.svd(pts - c)
        uv = (pts - c) @ vt[:2].T
    else:
        uv = (pts - c)[:, axes]
    return pts[np.argsort(np.arctan2(uv[:, 1], uv[:, 0]))]


def colour(key):
    return EXPOSED if key == "exposed" else (WELDED if key == "welded" else PAIRCOL[key])


def map_panel(ax, st):
    t, z = st["t"], st["z"]
    polys, vals = [], []
    for c in np.flatnonzero(t.original_ijk[:, 2] == 0):
        i, j, _ = t.original_ijk[c]
        p = t.nodes[t.cell_corners[c]][[0, 1, 3, 2]]
        polys.append(p[:, :2])
        vals.append(p[:, 2].mean())
    pc = PolyCollection(
        polys, array=np.array(vals), cmap="viridis_r", edgecolors="0.3", linewidths=0.6
    )
    pc.set_clim(0, 8)
    ax.add_collection(pc)
    cb = plt.colorbar(pc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("top of layer k0, mean depth [m]")
    x_top = LINE * DX - 20
    ax.plot([x_top, x_top], [0, TIP * DX], color="k", lw=4, solid_capstyle="butt")
    ax.plot([x_top], [TIP * DX], marker="*", ms=18, color="gold", mec="k", zorder=5)
    ax.annotate(
        "fault tip\n(throw 0)",
        (x_top, TIP * DX),
        (x_top - 150, TIP * DX + 90),
        arrowprops=dict(arrowstyle="->"),
        fontsize=9,
    )
    ax.text(
        x_top - 12,
        TIP * DX / 2,
        "FAULTS P",
        rotation=90,
        ha="right",
        va="center",
        fontsize=10,
        fontweight="bold",
    )
    for row in range(NY + 1):
        s = min(2 * row, 2 * NY - 1)
        th = z[0, s, 2 * LINE] - z[0, s, 2 * LINE - 1]
        ax.text(
            x_top + 14,
            row * DX,
            f"{th:g} m",
            fontsize=8.5,
            va="center",
            ha="left",
            color="w" if row < 3 else "k",
            bbox=dict(fc="#6a3d9a" if th > 0 else "0.6", ec="none", alpha=0.9, pad=1.5),
        )
    ax.text(150, 300, "W side\nnot moved", ha="center", va="center", fontsize=9, color="k")
    ax.text(
        470,
        120,
        "E side: corners move\ndown their pillars by\nthe row's throw",
        ha="center",
        va="center",
        fontsize=9,
        color="k",
    )
    ax.text(
        470,
        520,
        "past the tip: throw 0,\nfaces stay welded",
        ha="center",
        va="center",
        fontsize=9,
        color="k",
    )
    ax.set_xlim(-10, NX * DX + 10)
    ax.set_ylim(-10, NY * DX + 10)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]  (strike)")
    ax.set_title(
        "Map view after the tapered slide (purple tags: throw on each fault pillar)",
        fontsize=10.5,
        loc="left",
    )


def allan_panel(ax, st):
    t = st["t"]
    polys, cols = [], []
    for f, j, key in line_faces(st):
        p = ordered(t.nodes[t.faces.nodes[f]], axes=[1, 2])[:, [1, 2]]
        polys.append(p)
        cols.append(colour(key))
    ax.add_collection(PolyCollection(polys, facecolors=cols, edgecolors="k", linewidths=0.5))
    z = st["z"]
    for k in range(NZ):
        for side, slot, ls in (("W", 2 * LINE - 1, "-"), ("E", 2 * LINE, "--")):
            ys = np.arange(NY + 1) * DX
            sl = [min(2 * r, 2 * NY - 1) for r in range(NY + 1)]
            ax.plot(ys, z[2 * k, sl, slot], color="k" if side == "W" else "#6a3d9a", ls=ls, lw=1)
    ax.axvline(TIP * DX, color="k", lw=1.2, ls=":")
    ax.plot([TIP * DX], [0], marker="*", ms=16, color="gold", mec="k", clip_on=False, zorder=5)
    ax.text(
        (TIP + 1) * DX,
        15,
        "past the tip: no fault,\nthe two columns share corners\n(throw 0, faces welded)",
        fontsize=9,
        ha="center",
        va="center",
        bbox=dict(fc="white", ec="0.5", alpha=0.9),
    )
    for j in range(TIP):
        tot = {}
        for f, jj, key in line_faces(st):
            if jj == j:
                tot[key] = tot.get(key, 0.0) + st["areas"][f]
        ax.text(
            j * DX + 50,
            42,
            f"row {j}\nW k1|E k0\n{tot[(1, 0)]:,.0f} m²",
            ha="center",
            va="top",
            fontsize=8,
        )
    ax.set_xlim(0, NY * DX)
    ax.set_ylim(52, -3)
    ax.set_xlabel("y [m] along strike (fault plane seen from the east, projected)")
    ax.set_ylabel("depth [m]")
    ax.set_title(
        "Fault plane after re-processing: who touches whom (from build_topology faces)",
        fontsize=10.5,
        loc="left",
    )
    handles = [Patch(fc=PAIRCOL[p], ec="k", label=f"W k{p[0]} ↔ E k{p[1]}") for p in PAIRCOL]
    handles += [
        Patch(fc=EXPOSED, ec="k", label="exposed fault piece"),
        Patch(fc=WELDED, ec="k", label="no fault: welded face"),
        Line2D([], [], color="k", lw=1, label="W layer tops"),
        Line2D([], [], color="#6a3d9a", ls="--", lw=1, label="E layer tops"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8, ncol=3, framealpha=0.95)


def view3d_panel(ax, st):
    t = st["t"]
    polys, fc, ec = [], [], []
    for c in range(len(t.cell_corners)):
        i, j, k = t.original_ijk[c]
        p = t.nodes[t.cell_corners[c]]
        west = i < LINE
        for n, q in enumerate(HEX_FACES):
            if (i == LINE - 1 and n == 1) or (i == LINE and n == 0):
                continue
            polys.append(p[q])
            fc.append(to_rgba(KCOL[k], 1.0 if west else 0.12))
            ec.append(to_rgba("0.35", 1.0 if west else 0.45))
    for f, j, key in line_faces(st):
        polys.append(ordered(t.nodes[t.faces.nodes[f]]))
        fc.append(to_rgba(colour(key)))
        ec.append(to_rgba("k"))
    ax.computed_zorder = False
    ax.add_collection3d(Poly3DCollection(polys, facecolors=fc, edgecolors=ec, linewidths=0.4))
    pp = st["g"].coord.reshape(NY + 1, NX + 1, 6)
    ax.add_collection3d(
        Line3DCollection(
            [
                [pp[j, LINE, :3], pp[j, LINE, :3] + 45 / 80 * (pp[j, LINE, 3:] - pp[j, LINE, :3])]
                for j in range(NY + 1)
            ],
            colors="k",
            linestyles="--",
            linewidths=0.8,
        )
    )
    tip = pp[TIP, LINE, :3]
    ax.scatter(*tip, marker="*", s=220, color="gold", edgecolors="k", depthshade=False, zorder=10)
    ax.set_xlim(0, NX * DX)
    ax.set_ylim(0, NY * DX)
    ax.set_zlim(45, 0)
    ax.set_box_aspect((1.3, 1.3, 0.45))
    ax.view_init(elev=24, azim=-58)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("depth [m]")
    ax.set_title(
        "3D, true positions: W solid, E see-through; fault coloured by contact", fontsize=10.5
    )


def fig_png(st):
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1, 1.15])
    map_panel(fig.add_subplot(gs[0, 0]), st)
    view3d_panel(fig.add_subplot(gs[0, 1], projection="3d"), st)
    allan_panel(fig.add_subplot(gs[1, :]), st)
    fig.suptitle(
        "A fault that ends inside the grid: the throw tapers from 8 m at the edge to 0 at "
        "the tip (ZCORN-only; COORD and FAULTS unchanged)",
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    path = OUT / "07_partial_fault.png"
    path.unlink(missing_ok=True)
    fig.savefig(path, dpi=110)
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
                txt += [f"cell ({i},{j},{k}) {'W' if i < LINE else 'E'} side"] * 8
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
                    opacity=0.25,
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
        by_key = {}
        for f, j, key in line_faces(st):
            by_key.setdefault(str(key), (key, []))[1].append((f, j))
        for name, (key, items) in sorted(by_key.items()):
            xs, ys, zs, ii, jj, kk, txt = [], [], [], [], [], [], []
            for f, j in items:
                poly, base = ordered(t.nodes[t.faces.nodes[f]]), len(xs)
                xs += list(poly[:, 0])
                ys += list(poly[:, 1])
                zs += list(poly[:, 2])
                txt += [f"row {j}: {name}<br>area {st['areas'][f]:.1f} m²"] * len(poly)
                for m in range(1, len(poly) - 1):
                    ii.append(base)
                    jj.append(base + m)
                    kk.append(base + m + 1)
            label = name if isinstance(key, str) else f"W k{key[0]} ↔ E k{key[1]}"
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
            aspectratio=dict(x=1.4, y=1.4, z=0.5),
            camera=dict(eye=dict(x=1.1, y=-1.7, z=0.9)),
        ),
        legend=dict(x=0.0, y=0.0, bgcolor="rgba(255,255,255,0.85)", font=dict(size=10)),
        height=850,
        margin=dict(l=0, r=0, t=110, b=0),
    )
    path = OUT / "07_partial_fault_3d.html"
    fig.write_html(path, include_plotlyjs=True, full_html=True)
    return path


if __name__ == "__main__":
    before, after = load("input.grdecl"), load("tapered_8m.grdecl")
    print(fig_png(after))
    print(
        fig_html(
            [before, after],
            [
                ("input", "Input: fault P ends at pillar (3, 4); no throw"),
                (
                    "tapered 8 m",
                    "Tapered slide: throw 8, 6, 4, 2, 0 m on the fault pillars; welded past the tip",
                ),
            ],
        )
    )
