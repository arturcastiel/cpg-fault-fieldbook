"""Example 09 figures: control the throw of a fault that ends inside the grid (lesson 22).

usage: python throw_control.py DIR
  DIR/09_throw_control/<case>_{0_input,1_slid}.grdecl come from fieldbook lesson 22 (--out).
  writes DIR/09_throw_control.png (maps + profiles), DIR/09_throw_control_3d.png,
  DIR/09_throw_control_3d.html
"""

import sys
from pathlib import Path


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import TwoSlopeNorm, to_rgba
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from cpg2unstructured import build_topology, from_grdecl, read_grdecl
from fault_projections.from_faults import faults_from_deck

NX, NY, NZ, DX, DZ, LINE, TMAX = 16, 12, 3, 50.0, 20.0, 8, 16.0
VEX = 4.0
OUT = Path(sys.argv[1])
DECKS = OUT / "09_throw_control"
KCOL = ["#9ecae1", "#a1d99b", "#fdd0a2"]
EXPOSED = "#d7191c"
CASES = [
    (
        "A_linear_whole_side",
        "A  linear taper, whole side moves",
        "= lesson 20: the east side tilts like a hinged plank",
    ),
    (
        "B_short_taper",
        "B  cosine taper over the last 150 m only",
        "almost a rigid slide, with a smooth hinge at the tip",
    ),
    (
        "C_cosine_whole_side",
        "C  cosine taper over the full length",
        "smooth along strike; the east rim still moves",
    ),
    (
        "D_decay",
        "D  as C, plus decay over 300 m",
        "motion fades away from the fault; rims stay put",
    ),
    ("E_decay_split", "E  as D, split 50 / 50", "hanging wall down 8 m, footwall up 8 m"),
    (
        "F_two_tips",
        "F  two tips, decay 300 m, split 70 / 30",
        "an isolated fault: bell-shaped throw",
    ),
]
COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6761d"]


def load(name):
    g0, g1 = (
        read_grdecl(DECKS / f"{name}_0_input.grdecl"),
        read_grdecl(DECKS / f"{name}_1_slid.grdecl"),
    )
    z0 = g0.zcorn.reshape(2 * NZ, 2 * NY, 2 * NX)
    z1 = g1.zcorn.reshape(2 * NZ, 2 * NY, 2 * NX)
    rows = sorted(f.j - 1 for f in faults_from_deck(DECKS / f"{name}_1_slid.grdecl")["P"])
    return dict(
        name=name,
        shift=(z1 - z0)[0],
        z1=z1,
        rows=(rows[0], rows[-1] + 1),
        t=build_topology(from_grdecl(g1)),
    )


def pillar_shift(st, j, side):
    """Top shift at pillar row j along x; side 'east' or 'west' picks the cells' corner."""
    s, js = st["shift"], min(2 * j, 2 * NY - 1)
    xs, vals = [], []
    for i in range(NX + 1):
        if side == "west" and i > LINE or side == "east" and i < LINE:
            continue
        if side == "west":
            islot = 2 * i - 1 if i > 0 else 0
        else:
            islot = 2 * i if i < NX else 2 * NX - 1
        xs.append(i * DX)
        vals.append(s[js, islot])
    return np.array(xs), np.array(vals)


def throw_along(st):
    s = st["shift"]
    ys, ts = [], []
    for j in range(NY + 1):
        js = min(2 * j, 2 * NY - 1)
        ys.append(j * DX)
        ts.append(s[js, 2 * LINE] - s[js, 2 * LINE - 1])
    return np.array(ys), np.array(ts)


def tips(st):
    return [r for r in st["rows"] if 0 < r < NY]


def map_panel(ax, st, title, sub, norm):
    s = st["shift"]
    polys, vals = [], []
    for j in range(NY):
        for i in range(NX):
            polys.append(
                [
                    (i * DX, j * DX),
                    ((i + 1) * DX, j * DX),
                    ((i + 1) * DX, (j + 1) * DX),
                    (i * DX, (j + 1) * DX),
                ]
            )
            vals.append(s[2 * j : 2 * j + 2, 2 * i : 2 * i + 2].mean())
    pc = PolyCollection(
        polys, array=np.array(vals), cmap="RdBu_r", norm=norm, edgecolors="0.55", linewidths=0.3
    )
    ax.add_collection(pc)
    for side, lo, hi in (("west", 0, LINE + 1), ("east", LINE, NX + 1)):
        grid = np.array([pillar_shift(st, j, side)[1] for j in range(NY + 1)])
        X, Y = np.meshgrid(np.arange(lo, hi) * DX, np.arange(NY + 1) * DX)
        levels = [v for v in (-12, -8, -4, -1, 1, 4, 8, 12) if grid.min() < v < grid.max()]
        if levels:
            cs = ax.contour(X, Y, grid, levels=levels, colors="k", linewidths=0.6, alpha=0.6)
            ax.clabel(cs, fontsize=6.5, fmt="%g")
    r0, r1 = st["rows"]
    ax.plot([LINE * DX] * 2, [r0 * DX, r1 * DX], color="k", lw=3.5, solid_capstyle="butt")
    for r in tips(st):
        ax.plot(LINE * DX, r * DX, marker="*", ms=15, mfc="gold", mec="k", zorder=5)
    ax.set_xlim(0, NX * DX)
    ax.set_ylim(0, NY * DX)
    ax.set_aspect("equal")
    ax.set_title(f"{title}\n{sub}", fontsize=10, loc="left")
    ax.tick_params(labelsize=8)
    return pc


def profiles(ax_t, ax_x, states):
    for st, (name, title, _), c in zip(states, CASES, COLORS):
        ys, ts = throw_along(st)
        ax_t.plot(
            ys,
            ts,
            color=c,
            lw=2,
            marker="o",
            ms=3.5,
            label=title.split("  ")[0] + "  " + title.split("  ")[1],
        )
        peak = int(np.argmax(ts))
        for side in ("west", "east"):
            xs, v = pillar_shift(st, peak, side)
            ax_x.plot(
                xs,
                v,
                color=c,
                lw=2,
                label=f"{name[0]}  (row y = {peak * DX:g} m)" if side == "west" else None,
            )
    for st in states:
        for r in tips(st):
            ax_t.axvline(r * DX, color="0.6", ls=":", lw=1)
    ax_t.set_xlabel("y along the fault [m]")
    ax_t.set_ylabel("relative throw [m]")
    ax_t.set_title(
        "Throw along the fault (east minus west, on the fault pillars);\n"
        "C, D, E share one curve (decay and split act away from the fault); dotted: tip rows",
        fontsize=10,
        loc="left",
    )
    ax_t.legend(fontsize=8, loc="upper right")
    ax_t.set_ylim(-0.5, TMAX + 4)
    ax_x.axvline(LINE * DX, color="k", lw=2)
    ax_x.axhline(0, color="0.5", lw=0.8)
    ax_x.set_ylim(TMAX + 10, -13)
    ax_x.set_xlabel("x across the fault [m]")
    ax_x.set_ylabel("vertical shift of the top [m]  (down +)")
    ax_x.set_title(
        "Top of layer k0 across the fault, on the row of largest throw:\n"
        "whole side (A–C) vs decay (D–F); split moves the footwall up",
        fontsize=10,
        loc="left",
    )
    ax_x.text(LINE * DX - 15, -11, "footwall (west)", ha="right", fontsize=9)
    ax_x.text(LINE * DX + 15, -11, "hanging wall (east)", ha="left", fontsize=9)
    ax_x.legend(fontsize=8, loc="lower left", ncol=2)


def fig_png(states):
    fig = plt.figure(figsize=(17, 15.5))
    gs = fig.add_gridspec(
        3,
        3,
        height_ratios=[1, 1, 0.95],
        left=0.04,
        right=0.93,
        top=0.93,
        bottom=0.05,
        hspace=0.32,
        wspace=0.12,
    )
    norm = TwoSlopeNorm(vmin=-TMAX, vcenter=0.0, vmax=TMAX)
    for n, (st, (_, title, sub)) in enumerate(zip(states, CASES)):
        pc = map_panel(fig.add_subplot(gs[n // 3, n % 3]), st, title, sub, norm)
    cax = fig.add_axes([0.945, 0.40, 0.012, 0.5])
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(
        "vertical shift of the top of layer k0 [m]\nred: down   ·   blue: up", fontsize=9.5
    )
    profiles(fig.add_subplot(gs[2, :2]), fig.add_subplot(gs[2, 2]), states)
    fig.suptitle(
        "A fault that ends inside the grid (black line, gold star = tip): six ways to "
        f"shape a {TMAX:g} m throw by editing ZCORN only\nevery case: COORD and FAULTS "
        "unchanged, only the FAULTS faces open, no cell changes volume, contacts match "
        "the analytic overlaps",
        fontsize=12.5,
        fontweight="bold",
        y=0.985,
    )
    path = OUT / "09_throw_control.png"
    path.unlink(missing_ok=True)
    fig.savefig(path, dpi=100)
    plt.close(fig)
    return path


def ordered(p):
    c = p.mean(0)
    _, _, vt = np.linalg.svd(p - c)
    uv = (p - c) @ vt[:2].T
    return p[np.argsort(np.arctan2(uv[:, 1], uv[:, 0]))]


def view3d(ax, st, title, sub):
    t = st["t"]
    shade = {0: 0.86, 1: 0.72, 2: 1.0}
    polys, fc = [], []
    for f in np.flatnonzero(t.faces.neighbours[:, 1] < 0):
        a = t.faces.neighbours[f, 0]
        p = t.nodes[t.faces.nodes[f]]
        rim = np.allclose(p[:, 0], 0) or np.allclose(p[:, 0], NX * DX)
        fault = t.faces.direction[f] == 0 and not rim
        rgb = np.array(to_rgba(EXPOSED if fault else KCOL[t.original_ijk[a][2]]))
        rgb[:3] *= shade[int(t.faces.direction[f])]
        polys.append(ordered(np.unique(np.round(p, 9), axis=0)))
        fc.append(rgb)
    ax.add_collection3d(
        Poly3DCollection(polys, facecolors=fc, edgecolors=to_rgba("0.25", 0.5), linewidths=0.25)
    )
    ax.computed_zorder = False
    r0, r1 = st["rows"]
    for r in tips(st):
        ax.scatter(
            [LINE * DX],
            [r * DX],
            [st["z1"][0, min(2 * r, 2 * NY - 1), 2 * LINE]],
            marker="*",
            s=160,
            c="gold",
            edgecolors="k",
            zorder=30,
            depthshade=False,
        )
    ax.set_xlim(0, NX * DX)
    ax.set_ylim(0, NY * DX)
    ax.set_zlim(NZ * DZ + TMAX, -TMAX)
    ax.set_box_aspect((NX * DX, NY * DX, VEX * (NZ * DZ + 2 * TMAX)), zoom=1.25)
    ax.view_init(elev=24, azim=-58)
    ax.set_xlabel("x [m]", fontsize=8)
    ax.set_ylabel("y [m]", fontsize=8)
    ax.set_zlabel("depth [m]", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(f"{title}\n{sub}", fontsize=10)


def fig_3d(states):
    fig = plt.figure(figsize=(18, 12.5))
    for n, (st, (_, title, sub)) in enumerate(zip(states, CASES)):
        view3d(fig.add_subplot(2, 3, n + 1, projection="3d"), st, title, sub)
    fig.subplots_adjust(left=0.0, right=1.0, top=0.9, bottom=0.02, wspace=0.0, hspace=0.08)
    fig.suptitle(
        f"The same six throw fields in 3D (depth exaggerated {VEX:g}×): red = fault face "
        "left exposed; gold star = tip, where the throw is 0",
        fontsize=13,
        fontweight="bold",
        y=0.985,
    )
    path = OUT / "09_throw_control_3d.png"
    path.unlink(missing_ok=True)
    fig.savefig(path, dpi=100)
    plt.close(fig)
    return path


def fig_html(states):
    import plotly.graph_objects as go

    fig, groups = go.Figure(), []
    for s, st in enumerate(states):
        t, start = st["t"], len(fig.data)
        for k in range(NZ):
            xs, ys, zs, ii, jj, kk = [], [], [], [], [], []
            for c in np.flatnonzero(t.original_ijk[:, 2] == k):
                p, base = t.nodes[t.cell_corners[c]], len(xs)
                xs += list(p[:, 0])
                ys += list(p[:, 1])
                zs += list(p[:, 2])
                for q in (
                    [0, 2, 6, 4],
                    [1, 5, 7, 3],
                    [0, 4, 5, 1],
                    [2, 3, 7, 6],
                    [0, 1, 3, 2],
                    [4, 6, 7, 5],
                ):
                    ii += [base + q[0], base + q[0]]
                    jj += [base + q[1], base + q[2]]
                    kk += [base + q[2], base + q[3]]
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
                    name=f"layer k{k}",
                    showlegend=s == 0,
                    visible=s == 0,
                )
            )
        r0, r1 = st["rows"]
        fig.add_trace(
            go.Scatter3d(
                x=[LINE * DX] * 2,
                y=[r0 * DX, r1 * DX],
                z=[-TMAX] * 2,
                mode="lines",
                line=dict(color="black", width=8),
                name="fault trace (above the grid)",
                showlegend=s == 0,
                visible=s == 0,
            )
        )
        groups.append(range(start, len(fig.data)))
    n = len(fig.data)
    buttons = [
        dict(
            method="update",
            label=title.split("  ")[0],
            args=[
                {"visible": [m in g for m in range(n)]},
                {"title.text": f"{title}: {sub}  (depth exaggerated {VEX:g}×)"},
            ],
        )
        for g, (_, title, sub) in zip(groups, CASES)
    ]
    fig.update_layout(
        title=dict(text=f"{CASES[0][1]}: {CASES[0][2]}  (depth exaggerated {VEX:g}×)", y=0.97),
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
        scene=dict(
            xaxis_title="x [m]",
            yaxis_title="y [m]",
            zaxis_title="depth [m]",
            zaxis=dict(autorange="reversed"),
            aspectmode="manual",
            aspectratio=dict(x=1.6, y=1.2, z=1.6 * VEX * (NZ * DZ + 2 * TMAX) / (NX * DX)),
            camera=dict(eye=dict(x=1.2, y=-1.6, z=0.8)),
        ),
        height=850,
        margin=dict(l=0, r=0, t=110, b=0),
    )
    path = OUT / "09_throw_control_3d.html"
    fig.write_html(path, include_plotlyjs=True, full_html=True)
    return path


if __name__ == "__main__":
    states = [load(name) for name, _, _ in CASES]
    print(fig_png(states))
    print(fig_3d(states))
    print(fig_html(states))
