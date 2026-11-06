"""Draw the flowcharts and illustrations of the wiki's Workflows section.

usage: python tools/figures/workflows.py [ASSETS_DIR]

Writes PNGs under ASSETS_DIR/workflows (by default docs/assets/workflows):

- chooser.png and flow-N-*.png: flowcharts. Each box names the call to make and the lesson
  that runs it. They are drawn, not computed.
- concepts.png and why-pillars.png: SCHEMATICS, drawn by hand in code: one idea per panel,
  not a grid any call returned.
- ill-N-*.png and at-a-glance.png: BEFORE / AFTER illustrations COMPUTED with the libraries
  on one synthetic layer cake (16 x 8 columns of 50 m, six 15 m layers, top at 2000 m). Every
  panel is a grid an API call returned, drawn from its own COORD and ZCORN, and the numbers
  in the titles are measured on it.

Needs matplotlib; about a minute.
"""

import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DPI = 150

# ── flowchart engine ─────────────────────────────────────────────────────────────────────────

KINDS = {
    "start": dict(fc="#deebf7", ec="#3182bd", label="what you start from"),
    "data": dict(fc="#fdf1dd", ec="#d9822b", label="data or file"),
    "step": dict(fc="#ffffff", ec="#404040", label="library call"),
    "check": dict(fc="#efedf5", ec="#6a51a3", label="check or gate"),
    "decision": dict(fc="#fff7bc", ec="#b8860b", label="decision"),
    "end": dict(fc="#e5f5e0", ec="#31a354", label="result"),
    "stop": dict(fc="#fde0dd", ec="#b2182b", label="refused: no grid written"),
    "note": dict(fc="#f5f5f5", ec="#9e9e9e", label="where to read more"),
}
TITLE_KW = dict(fontsize=9.5, family="DejaVu Sans", weight="bold", color="#1a1a1a")
CODE_KW = dict(fontsize=8.0, family="DejaVu Sans Mono", color="#08519c")
REF_KW = dict(fontsize=7.5, family="DejaVu Sans", style="italic", color="#555555")
PAD_X, PAD_Y, GAP = 9.0, 7.0, 3.0


@dataclass
class Box:
    """One flowchart node at grid position (col, row); text is title, code and a reference."""

    key: str
    kind: str
    col: float
    row: float
    title: str
    code: str = ""
    ref: str = ""
    width: float = 0.0          # minimum width in points


@dataclass
class Edge:
    """An arrow from box ``a`` to box ``b``.

    ``route``: "v" leaves a's bottom (or top) and enters b's top (or bottom), with one
    horizontal jog at ``via`` (a row number; default halfway); "h" is the same sideways;
    "vh" goes down then into b's side; "hv" goes sideways then into b's top or bottom.
    The label sits on the segment that enters b unless ``seg`` says otherwise.
    """

    a: str
    b: str
    label: str = ""
    route: str = "v"
    via: float | None = None
    seg: int | None = None      # which segment carries the label


def _measure(fig, renderer, text, kw):
    if not text:
        return 0.0, 0.0
    t = fig.text(0, 0, text, linespacing=1.25, **kw)
    bb = t.get_window_extent(renderer)
    t.remove()
    return bb.width * 72 / fig.dpi, bb.height * 72 / fig.dpi


def _sizes(boxes):
    import matplotlib.pyplot as plt

    fig = plt.figure(dpi=DPI)
    renderer = fig.canvas.get_renderer()
    out = {}
    for b in boxes:
        parts = [_measure(fig, renderer, txt, kw)
                 for txt, kw in ((b.title, TITLE_KW), (b.code, CODE_KW), (b.ref, REF_KW))]
        w = max(p[0] for p in parts) + 2 * PAD_X
        h = sum(p[1] for p in parts) + GAP * (sum(1 for p in parts if p[1]) - 1) + 2 * PAD_Y
        w = max(w, b.width)
        if b.kind == "decision":
            w, h = w * 1.42, h * 1.62
        elif b.kind == "data":
            w += 16
        out[b.key] = (w, h, parts)
    plt.close(fig)
    return out


def flowchart(path, title, subtitle, boxes, edges, xs=230.0, ys=92.0):
    """Lay the boxes out on a (col, row) grid of ``xs`` by ``ys`` points and draw the chart."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Polygon

    size = _sizes(boxes)
    pos = {b.key: (b.col * xs, b.row * ys) for b in boxes}
    by_key = {b.key: b for b in boxes}
    x0 = min(pos[k][0] - size[k][0] / 2 for k in pos) - 24
    x1 = max(pos[k][0] + size[k][0] / 2 for k in pos) + 24
    y0 = min(pos[k][1] - size[k][1] / 2 for k in pos) - 58
    y1 = max(pos[k][1] + size[k][1] / 2 for k in pos) + 48
    fig = plt.figure(figsize=((x1 - x0) / 72, (y1 - y0) / 72), dpi=DPI)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    def anchor(k, side):
        (x, y), (w, h, _) = pos[k], size[k]
        return {"top": (x, y - h / 2), "bottom": (x, y + h / 2),
                "left": (x - w / 2, y), "right": (x + w / 2, y)}[side]

    for e in edges:
        (ax_, ay), (bx, by) = pos[e.a], pos[e.b]
        down, right = by > ay, bx > ax_
        if e.route == "v":
            s, t = anchor(e.a, "bottom" if down else "top"), anchor(e.b, "top" if down else "bottom")
            ym = e.via * ys if e.via is not None else (s[1] + t[1]) / 2
            pts = [s, t] if abs(s[0] - t[0]) < 0.5 else [s, (s[0], ym), (t[0], ym), t]
        elif e.route == "h":
            s, t = anchor(e.a, "right" if right else "left"), anchor(e.b, "left" if right else "right")
            xm = e.via * xs if e.via is not None else (s[0] + t[0]) / 2
            pts = [s, t] if abs(s[1] - t[1]) < 0.5 else [s, (xm, s[1]), (xm, t[1]), t]
        elif e.route == "vh":
            s, t = anchor(e.a, "bottom" if down else "top"), anchor(e.b, "left" if right else "right")
            pts = [s, (s[0], t[1]), t]
        else:  # "hv"
            s, t = anchor(e.a, "right" if right else "left"), anchor(e.b, "top" if down else "bottom")
            pts = [s, (t[0], s[1]), t]
        pts = np.array(pts, float)
        ax.plot(pts[:-1, 0], pts[:-1, 1], color="#555555", lw=1.1, solid_capstyle="round",
                zorder=1)
        ax.annotate("", xy=pts[-1], xytext=pts[-2], zorder=1,
                    arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.1, shrinkA=0,
                                    shrinkB=0, mutation_scale=11))
        if e.label:
            # default: on the segment that enters b (a branch's own drop), so labels on a
            # shared bus are never ambiguous; "hv" labels its first leg
            seg = e.seg if e.seg is not None else (0 if e.route == "hv" else len(pts) - 2)
            p, q = pts[seg], pts[seg + 1]
            ax.text((p[0] + q[0]) / 2, (p[1] + q[1]) / 2, e.label, ha="center", va="center",
                    fontsize=7.5, style="italic", color="#333333", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none"))

    for k, (x, y) in pos.items():
        b, (w, h, parts) = by_key[k], size[k]
        st = KINDS[b.kind]
        common = dict(fc=st["fc"], ec=st["ec"], lw=1.3, zorder=2)
        if b.kind in ("start", "end"):
            patch = FancyBboxPatch((x - w / 2, y - h / 2), w, h, **common,
                                   boxstyle="round,pad=0,rounding_size=11")
        elif b.kind == "decision":
            patch = Polygon([(x, y - h / 2), (x + w / 2, y), (x, y + h / 2), (x - w / 2, y)],
                            **common)
        elif b.kind == "data":
            s = 8.0
            patch = Polygon([(x - w / 2 + s, y - h / 2), (x + w / 2, y - h / 2),
                             (x + w / 2 - s, y + h / 2), (x - w / 2, y + h / 2)], **common)
        elif b.kind == "check":
            c = 7.0
            patch = Polygon([(x - w / 2 + c, y - h / 2), (x + w / 2 - c, y - h / 2),
                             (x + w / 2, y - h / 2 + c), (x + w / 2, y + h / 2 - c),
                             (x + w / 2 - c, y + h / 2), (x - w / 2 + c, y + h / 2),
                             (x - w / 2, y + h / 2 - c), (x - w / 2, y - h / 2 + c)], **common)
        else:
            patch = FancyBboxPatch((x - w / 2, y - h / 2), w, h, **common,
                                   boxstyle="square,pad=0",
                                   linestyle="--" if b.kind == "note" else "-")
        ax.add_patch(patch)
        text_h = sum(p[1] for p in parts) + GAP * (sum(1 for p in parts if p[1]) - 1)
        top = y - text_h / 2
        for txt, kw, (_, th) in zip((b.title, b.code, b.ref), (TITLE_KW, CODE_KW, REF_KW), parts):
            if not txt:
                continue
            ax.text(x, top, txt, ha="center", va="top", linespacing=1.25, zorder=4,
                    multialignment="left" if kw is CODE_KW else "center", **kw)
            top += th + GAP

    ax.text(x0 + 14, y0 + 12, title, ha="left", va="top", fontsize=13, weight="bold",
            color="#1a1a1a")
    ax.text(x0 + 14, y0 + 32, subtitle, ha="left", va="top", fontsize=9, color="#444444")
    used = [k for k in KINDS if any(b.kind == k for b in boxes)]
    lx = x0 + 14
    for k in used:
        st = KINDS[k]
        ax.add_patch(FancyBboxPatch((lx, y1 - 28), 18, 11, boxstyle="square,pad=0", fc=st["fc"],
                                    ec=st["ec"], lw=1.1,
                                    linestyle="--" if k == "note" else "-"))
        t = ax.text(lx + 24, y1 - 22.5, st["label"], ha="left", va="center", fontsize=7.5,
                    color="#444444")
        renderer = fig.canvas.get_renderer()
        lx += 24 + t.get_window_extent(renderer).width * 72 / DPI + 18
    fig.savefig(path, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(path)


# ── the flowcharts ───────────────────────────────────────────────────────────────────────────


def chooser(out):
    B, E = Box, Edge
    boxes = [
        B("s", "start", 2.5, 0, "A corner-point grid and a fault",
          ref="grid: GRDECL or EGRID (lessons 01, 02, 14)\nfault: a map trace, a 3D polygon,"
              " facets, a quad mesh,\nor the FAULTS already in the deck"),
        B("d1", "decision", 2.5, 1.05, "Should the grid's\ngeometry change?"),
        B("d2", "decision", 0.75, 2.05, "What do you\nhave or want?"),
        B("d3", "decision", 4.0, 2.05, "What should move?"),
        B("w1", "end", 0.0, 3.5, "1 · Stair-step projection",
          ref="map trace → FAULTS through all layers\ngrid unchanged"),
        B("w2", "end", 1.0, 3.5, "2 · Full 3D projection",
          ref="3D surface → X, Y and Z faces\nwith real K ranges; grid unchanged;\n"
              "to slide it: 9"),
        B("w7", "end", 2.05, 3.5, "7 · Fault zones",
          ref="the cells within ±D of a fault;\noptionally keep 4, 5, 6 inside them"),
        B("w3", "end", 3.1, 3.5, "3 · 2D / 2.5D deformation",
          ref="pillars move in plan onto the fault\n(stage 1) or tilt into it (stage 2)"),
        B("w4", "end", 4.1, 3.5, "4 · Fault slide",
          ref="corners move down their pillars\n(ZCORN only): the throw"),
        B("w5", "end", 5.1, 3.5, "5 · Slide + pillar shift",
          ref="fit or tilt the pillars, then\nslide along them, in one call"),
        B("w6", "end", 6.1, 3.5, "6 · Set a dip angle",
          ref="tilt a fault's pillars to an angle\n(COORD only)"),
        B("w9", "end", 7.1, 3.5, "9 · Slide a 3D surface",
          ref="vertical: slide its 3D FAULTS;\ndipping: move the pillars onto it\n"
              "first; horizontal: not possible"),
        B("w8", "note", 3.1, 4.55, "8 · On a field deck: wrap any of these in the deck workflow",
          ref="intake (read only) → edit → the four gates → INCLUDE files; the master deck"
              " is never rewritten"),
    ]
    edges = [
        E("s", "d1"),
        E("d1", "d2", "no", route="hv"),
        E("d1", "d3", "yes", route="hv"),
        E("d2", "w1", "a map trace", via=2.75),
        E("d2", "w2", "a 3D surface", via=2.75),
        E("d2", "w7", "only the cells near it", via=2.75),
        E("d3", "w3", "pillars, in plan", via=2.75),
        E("d3", "w4", "layers, in depth", via=2.75),
        E("d3", "w5", "both", via=2.75),
        E("d3", "w6", "the dip angle", via=2.75),
        E("d3", "w9", "slip on a\n3D surface", via=2.75),
    ]
    flowchart(out / "chooser.png", "Which workflow?",
              "Answer two questions; each result is a page of this section with its own "
              "flowchart, the calls, the lessons and pictures.", boxes, edges, xs=205, ys=90)


def flow_1(out):
    B, E = Box, Edge
    boxes = [
        B("g", "data", 0, 0, "The grid", "grid = read_grdecl(\"FIELD.GRDECL\")",
          "or read_egrid · lessons 01, 02, 17 (MAPAXES)"),
        B("t", "data", 1.2, 0, "A fault trace in map view", "line = np.array([[x0, y0], [x1, y1]])",
          "(N, 2) XY, in the grid's frame"),
        B("m", "step", 0, 1, "Build the topology", "mesh = build_topology(from_grdecl(grid))",
          "lesson 03 · Methods: Topology and quality"),
        B("p", "step", 0, 2, "Project it into map view",
          "top = projected_topology(mesh, \"xy\")\ngraph = build_projected_graph(top)",
          "lesson 04"),
        B("s", "step", 0, 3.05, "Snap the trace to grid edges",
          "paths = stitch_paths(\n    project_polylines([line], graph), graph)", "lesson 04"),
        B("f", "step", 0, 4.1, "Turn the paths into fault faces",
          "faces = [path_to_faults(p, top) for p in paths]", "lesson 04"),
        B("o", "end", 0, 5.1, "FAULTS through all layers",
          "faults_text(faces, nz=nz)  ·  write_faults(path, faces, nz)",
          "the grid is unchanged: the fault is a staircase"),
        B("short", "note", 1.2, 4.1, "Shortcut: all of the above in one call",
          "faces = fault_faces_for(grid, [line])",
          "used by lessons 21, 25, 30, 31; returns FaultFace(i, j, 'X'|'Y')"),
        B("c", "check", 1.2, 5.1, "Is the staircase good enough?",
          ref="it is longer than the trace: about 1.35x at 30 deg,\nat any resolution (lesson"
              " 28 A, B); to put the fault\non the trace, deform the grid: workflow 3"),
    ]
    edges = [E("g", "m"), E("m", "p"), E("p", "s"), E("s", "f"), E("f", "o"),
             E("t", "s", route="vh"), E("short", "f", "same faces", route="h"),
             E("o", "c", route="h")]
    flowchart(out / "flow-1-stair-step.png", "Workflow 1 · Stair-step projection",
              "Put a map-view fault trace on the grid as FAULTS without moving a single "
              "corner. Learn: Fault projection (lesson 04).", boxes, edges, xs=300, ys=80)


def flow_2(out):
    B, E = Box, Edge
    boxes = [
        B("g", "data", 1.5, 0, "The grid, as cell faces",
          "mesh = build_topology(from_grdecl(read_grdecl(path)))", "lessons 01, 03"),
        B("d", "decision", 1.5, 1.25, "What is the fault\nsurface?"),
        B("a", "step", 0, 2.75, "A planar polygon", "surface = project_fault_surface(\n"
          "    polygon, mesh)", "(N, 3) coplanar XYZ · lesson 05"),
        B("b", "step", 1.0, 2.75, "Planar facets", "surface = project_faceted_surface(\n"
          "    quads, mesh)", "list of (4, 3) quads · lesson 25 C"),
        B("c", "step", 2.0, 2.75, "A quad mesh (shared seams)",
          "surface = project_quad_mesh_surface(\n    QuadMesh(verts, quads), mesh)",
          "lessons 13, 25 D"),
        B("e", "data", 3.3, 1.25, "Sticks or a point cloud?", "p0, n, _, _ = fit_plane(pts)\n"
          "quads = facet_surface(cloud, ...)", "lesson 25 B, C"),
        B("r", "end", 1.0, 4.05, "A FaultSurface", "surface.cell_ijk, .direction, .area",
          "X, Y and Z cell faces (Z: a layer top or base) with real K ranges"),
        B("o", "end", 1.0, 5.1, "3D FAULTS",
          "write_faults_3d(path, {\"F1\": surface}, mesh)", "faults_text_3d returns the text"),
        B("n", "check", 2.55, 4.6, "The grid is unchanged",
          ref="a dipping fault is a staircase in depth too.\nTo make cells follow it: "
              "workflow 3 (stage 2)\nor 6 (dip angle). The staircase cannot be slid\n"
              "(SlideError); to slide the fault: workflow 9"),
    ]
    edges = [E("g", "d"), E("d", "a", "planar", via=2.05), E("d", "b", "curved", via=2.05),
             E("d", "c", "a mesh", via=2.05),
             E("e", "d", "first make a\npolygon or facets", route="h"),
             E("a", "r", via=3.45), E("b", "r", via=3.45), E("c", "r", via=3.45), E("r", "o"),
             E("r", "n", route="hv")]
    flowchart(out / "flow-2-full-3d.png", "Workflow 2 · Full 3D projection",
              "Put a 3D fault surface on the grid's cell faces in X, Y and Z, without moving the "
              "grid. Learn: Fault projection; Methods: Surface projection.", boxes, edges,
              xs=245, ys=86)


def flow_3(out):
    B, E = Box, Edge
    boxes = [
        B("g", "data", 1.05, 0, "The grid", "grid = read_grdecl(path)", "lessons 01, 14"),
        B("d", "decision", 1.05, 1.25, "What do you know\nabout the fault?"),
        B("s1", "step", 0, 3.1, "Stage 1 (2D): a map trace",
          "res = adapt_grid(grid, trace,\n    boundary=\"pinned\", w0=0.05,\n"
          "    polish=\"barrier\", on_invalid=\"raise\")",
          "pillars move rigidly in plan onto the trace\nseveral traces: pass them together"
          " (07)\nlessons 06, 07"),
        B("s2", "step", 1.05, 3.1, "Stage 2 (2.5D): a dipping plane",
          "res = adapt_grid(grid, polygon,\n    regenerate_faults=False)",
          "(N, 3) polygon: top and base move apart,\nso the pillars tilt into the plane\n"
          "lesson 08; facets: lesson 16"),
        B("s12", "step", 2.15, 3.1, "Both stages in one call",
          "res = adapt_grid(grid, polygon,\n    top_trace=\"auto\")",
          "stage 1 snaps a chain to the fault's trace\non the top (its FAULTS); stage 2 tilts\n"
          "lesson 30 C"),
        B("k", "check", 1.05, 4.6, "Judge it: the four gates",
          "evaluate(grid, res.grid, trace, targets, min_sj=0.3)",
          "validity, quality floor, corridor, chain, rim (lesson 11);\nevaluate needs the "
          "targets of the low-level path:\ntargets_from_polyline + morph_pillars (lesson 06)"),
        B("q", "decision", 1.05, 5.85, "Pass?"),
        B("fix", "note", 2.75, 5.85, "Change one setting at a time",
          ref="boundary pinned or slide, w0, polish (lesson 09);\nsheared outline (lesson 10);"
              " Decisions page"),
        B("r", "end", 1.05, 7.05, "AdaptResult",
          "res.grid (new COORD) · res.faults_after (FAULTS)",
          "ZCORN depths unchanged; export: adapt_deck (lesson 12) or workflow 8"),
    ]
    edges = [E("g", "d"), E("d", "s1", "a map trace", via=1.95), E("d", "s2", "a dipping plane",
             via=1.95), E("d", "s12", "a plane, and\nI want FAULTS", via=1.95),
             E("s1", "k", via=4.0), E("s2", "k", via=4.0), E("s12", "k", via=4.0),
             E("k", "q"), E("q", "fix", "no", route="h"), E("q", "r", "yes"),
             E("fix", "d", "retry", route="vh")]
    flowchart(out / "flow-3-deformation.png", "Workflow 3 · 2D / 2.5D deformation",
              "Move the pillars so the grid follows the fault: in plan (stage 1) or tilted "
              "into a dipping plane (stage 2). Learn: Deformation; Validation.", boxes, edges,
              xs=265, ys=84)


def flow_4(out):
    B, E = Box, Edge
    boxes = [
        B("g", "data", 1, 0, "A grid and its FAULTS",
          "faults = faults_from_deck(deck)  # by name\nfaces = fault_faces_for(grid, [line])\n"
          "faces = adapt_grid(...).faults_after",
          "from the deck, from a trace (workflow 1), or after\nadapting (workflow 3)"),
        B("t", "step", 1, 1.35, "Define each throw",
          "Throw(\"F1\", 16.0, hanging_side=EAST,\n      decay=\"auto\")",
          "throw: metres, one per horizon, or a distance table\n(lessons 22, 26); layers=(k0,"
          " k1); growth; snap=True\nfor dipping layers (25 F); profile and taper at tips (20)"),
        B("d", "decision", 1, 2.75, "Does the fault cut\nthe grid in two?"),
        B("rg", "note", 0, 2.75, "Rigid block",
          ref="rim to rim or a closed ring:\nthe hanging wall drops as one\n(lesson 19)"),
        B("dc", "note", 2, 2.75, "Fades to its tips",
          ref="ends inside the grid: the throw\ntapers at each tip and decays\naway (lessons"
              " 20, 22)"),
        B("s", "step", 1, 4.1, "Slide",
          "res = slide_grid(grid, [Throw(...), ...],\n    faults=faults)",
          "decay=\"auto\" picks rigid or fade per fault;\nseveral faults: one list; in"
          " sequence: call again (19)"),
        B("k", "check", 1, 5.25, "Gates and re-processing",
          "min_jacobian=0.05 · on_invalid=\"raise\"",
          "re-read the deck: read_grdecl → from_grdecl → build_topology;\nthe opened faces "
          "must be FAULTS faces (lessons 19 to 22)"),
        B("r", "end", 1, 6.4, "SlideResult",
          "res.grid · res.opened · res.modes · res.volume_change",
          "only ZCORN changes: COORD, ACTNUM and FAULTS stay"),
        B("deck", "note", 2.25, 6.4, "On a deck file",
          "slide_deck(deck, Throw(\"E_01\", 5.0), out_dir)", "writes a ZCORN INCLUDE "
          "(lesson 30 B)"),
    ]
    edges = [E("g", "t"), E("t", "d"), E("d", "rg", "yes", route="h"),
             E("d", "dc", "no", route="h"), E("rg", "s", route="vh"),
             E("dc", "s", route="vh"), E("s", "k"), E("k", "r"), E("r", "deck", route="h")]
    flowchart(out / "flow-4-slide.png", "Workflow 4 · Fault slide",
              "Offset the layers across a fault by moving corners down their own pillars "
              "(ZCORN only). Learn: Slide a fault in Z; Gallery.", boxes, edges, xs=255,
              ys=84)


def flow_5(out):
    B, E = Box, Edge
    boxes = [
        B("g", "data", 1, 0, "A grid (and FAULTS, or a fault plane)",
          ref="the deck on disk or a grid in memory"),
        B("d", "decision", 1, 1.2, "How should the\npillars move?"),
        B("a", "step", 0, 3.0, "Fit them to a fault, then slide",
          "adapt_deck(deck, polygon, out_dir,\n    top_trace=\"auto\", fault_names=\"F\",\n"
          "    slide=Throw(\"F\", 15.0,\n        hanging_side=(x, y)),\n"
          "    boundary=\"slide\")",
          "rim to rim, 65 deg: lesson 30 C; two traces: 30 A;\nin memory: adapt_grid, then"
          " slide_grid\non res.grid with res.faults_after"),
        B("b", "step", 2, 3.0, "Tilt them to an angle, then slide",
          "dip_grid(grid, Dip(\"F\", 60.0),\n    faults=faults,\n"
          "    slide=Throw(\"F\", 15.0))", "on a deck: dip_deck(..., slide=...)\nlesson 31 I"),
        B("o", "check", 1, 4.5, "Order: pillars first, then ZCORN",
          ref="the slide runs on the moved pillars, so the hanging wall\nmoves along the "
              "tilted fault (dip-slip); lesson 21 B shows\nthe same slide on a regular and on"
              " an adapted grid"),
        B("v", "check", 1, 5.55, "Volumes change, by the predicted wedge",
          ref="a slide along tilted pillars changes cell volumes:\nthrow x thickness x "
              "cot(dip) per metre of fault\n(lessons 25 B, 30 C)"),
        B("r", "end", 1, 6.6, "One result, both edits",
          "adapt_deck: result.grid  ·  dip_grid: res.slide.grid",
          "COORD and ZCORN changed; res.slide is the SlideResult (dip_grid's\nres.grid is "
          "the dip alone); adapt_deck regenerates FAULTS, dip_grid keeps them"),
        B("z", "note", 2.3, 6.05, "Keep it local", "dip_grid(..., slide=..., zone=z)",
          "workflow 7, lesson 33 I"),
    ]
    edges = [E("g", "d"), E("d", "a", "onto a trace or plane", via=1.85),
             E("d", "b", "to a dip angle", via=1.85), E("a", "o", via=3.95),
             E("b", "o", via=3.95), E("o", "v"), E("v", "r"), E("z", "r", route="vh")]
    flowchart(out / "flow-5-slide-pillars.png", "Workflow 5 · Slide + pillar shift",
              "Move the pillars (COORD) and slide the layers along them (ZCORN) in one call. "
              "Learn: lessons 21, 30, 31 I.", boxes, edges, xs=280, ys=86)


def flow_6(out):
    B, E = Box, Edge
    boxes = [
        B("g", "data", 1, 0, "A grid with named FAULTS",
          "faults = faults_from_deck(deck)  # or {\"F\": faces}", "lesson 31 J reads Norne's"),
        B("p", "step", 1, 1.25, "Say what angle", "Dip(\"F\", 60.0, hanging_side=EAST,\n"
          "    pivot=\"middle\", decay=\"auto\")",
          "angle: one value, a table along strike, one per pillar (31 E);\npivot: top, middle"
          " or bottom (31 C); decay: rigid or a fade (31 A, D)"),
        B("s", "step", 1, 2.45, "Tilt the pillars",
          "res = dip_grid(grid, Dip(...), faults=faults)",
          "default boundary=\"fixed\": the outline holds still; \"free\" lets it move (31 A)"),
        B("k", "check", 1, 3.6, "The gate",
          ref="every cell must stay valid: the footwall fade widens\n(with a warning), or the"
              " call refuses or bisects (31 H)"),
        B("q", "decision", 1, 4.75, "Passed?"),
        B("fix", "note", 2.4, 4.75, "Lower the angle, pivot=\"top\",\nor a wider fade",
          ref="Norne BC at 70 deg fails, pivot=\"top\" passes (31 J);\nshallow angles and the "
              "fixed outline (32 O)"),
        B("r", "end", 1, 5.9, "DipResult", "res.grid (COORD only) · res.report",
          "ZCORN, ACTNUM and FAULTS are byte-identical"),
        B("v", "check", 0, 5.9, "Check it yourself",
          ref="the angle from COORD's six numbers per pillar;\n22 stress cases judged "
              "independently (lesson 32)"),
        B("deck", "note", 2.4, 5.9, "On a deck file",
          "dip_deck(deck, Dip(\"E_01\", 60.0), out_dir)", "add slide= for workflow 5 (31 I)"),
    ]
    edges = [E("g", "p"), E("p", "s"), E("s", "k"), E("k", "q"), E("q", "fix", "no", route="h"),
             E("fix", "p", "retry", route="vh"), E("q", "r", "yes"), E("r", "v", route="h"),
             E("r", "deck", route="h")]
    flowchart(out / "flow-6-dip.png", "Workflow 6 · Set a fault's dip angle",
              "Tilt the pillars of a fault to an angle, and move their neighbours so no cell "
              "folds (COORD only). Learn: lessons 31, 32.", boxes, edges, xs=290, ys=84)


def flow_7(out):
    B, E = Box, Edge
    boxes = [
        B("g", "data", 1, 0, "A grid with named faults", "faults = {\"F\": faces}",
          "or faults_from_deck(deck)"),
        B("z", "step", 1, 1.2, "Find the cells within ±D",
          "z = fault_zone(grid, faults, 100.0,\n    hanging_side=EAST)",
          "D: one distance or (footwall, hanging wall) (33 C);\nmeasured along the fault's "
          "normal, in 3D for dipping faults (33 B)"),
        B("d", "decision", 1, 2.55, "What do you need?"),
        B("sel", "step", 0, 4.0, "The cells only",
          "z.cells · z.minus · z.plus · z.side(\"F\")\nwrite_zone(z, \"zone.grdecl\")",
          "FLTZONE: 0 outside, 1 footwall, 2 hanging wall;\nrename it to MULTNUM or FIPNUM "
          "(33 F)"),
        B("con", "step", 2, 4.0, "Keep an edit inside them",
          "dip_grid(grid, Dip(...), zone=z)\nslide_grid(grid, Throw(...), zone=z)\n"
          "dip_grid(grid, Dip(...), slide=Throw(...), zone=z)", "lessons 33 G, H, I"),
        B("k", "check", 2, 5.3, "Outside the zone nothing moves",
          ref="every corner outside is bit-identical; the edit\nfades to zero at the zone's "
              "edge; a zone too\nsmall for the edit is refused (33 J)"),
        B("dft", "note", 0.95, 5.3, "Without zone= (the default)",
          ref="the edit may move the whole domain,\nwith the grid's outline fixed"),
        B("r1", "end", 0, 5.3, "A region keyword", ref="a cell mask, or an INCLUDE file"),
        B("r2", "end", 2, 6.4, "A local edit", ref="the same angle or throw at the fault"),
    ]
    edges = [E("g", "z"), E("z", "d"), E("d", "sel", "select", via=3.15),
             E("d", "con", "select and confine", via=3.15), E("sel", "r1"), E("con", "k"),
             E("k", "r2")]
    flowchart(out / "flow-7-zones.png", "Workflow 7 · Fault zones",
              "Find the cells within -D and +D of a fault, and optionally keep a dip or a "
              "slide inside them. Learn: lesson 33.", boxes, edges, xs=290, ys=84)


def flow_8(out):
    B, E = Box, Edge
    boxes = [
        B("i", "start", 1.33, 0, "A field deck", "run 14_field_deck_intake --deck FIELD.GRDECL",
          "read only: dimensions, ACTNUM, MAPAXES, cell quality (lessons 14, 17)"),
        B("f", "data", 1.33, 1.15, "Its faults", "faults = faults_from_deck(deck)",
          "FAULTS through INCLUDE, by name (lesson 30 B)"),
        B("d", "decision", 1.33, 2.4, "Which edit?"),
        B("a", "step", 0, 3.85, "Adapt", "adapt_deck(deck, traces, out_dir)",
          "workflow 3 · lesson 12"),
        B("s", "step", 0.9, 3.85, "Slide", "slide_deck(deck, throws, out_dir)",
          "workflow 4 · lesson 30 B"),
        B("p", "step", 1.8, 3.85, "Dip", "dip_deck(deck, dips, out_dir)",
          "workflow 6 · lesson 31 I, J"),
        B("zn", "step", 2.65, 3.85, "Zone", "write_zone(z, path)", "workflow 7 · lesson 33 K"),
        B("k", "check", 1.33, 5.15, "Check before you keep it",
          ref="four gates (11) · same-scale SVGs (15) · a JSON dossier (18)\n"
              "read the written files back (29 D) · bit-exact bundle (29 A)"),
        B("o", "end", 1.33, 6.25, "Separate INCLUDE files",
          "<prefix>_grid.grdecl · <prefix>_faults.grdecl · report",
          "the master deck is never rewritten"),
        B("u", "note", 2.75, 6.25, "Your part",
          ref="INCLUDE order, ACTNUM, MAPAXES,\nsimulator parsing and flow"),
    ]
    edges = [E("i", "f"), E("f", "d"), E("d", "a", via=3.15), E("d", "s", via=3.15),
             E("d", "p", via=3.15), E("d", "zn", via=3.15), E("a", "k", via=4.55),
             E("s", "k", via=4.55), E("p", "k", via=4.55), E("zn", "k", via=4.55), E("k", "o"),
             E("o", "u", route="h")]
    flowchart(out / "flow-8-field-deck.png", "Workflow 8 · A field deck, end to end",
              "Read a deck without changing it, edit it with one of workflows 3 to 7, check, "
              "and write separate INCLUDE files. Learn: Validation and export.", boxes, edges,
              xs=250, ys=84)


def flow_9(out):
    B, E = Box, Edge
    boxes = [
        B("s", "start", 1.55, 0, "A 3D fault surface",
          ref="a polygon, facets or a quad mesh (workflow 2);\nsticks or a point cloud: "
              "fit_plane / facet_surface (lesson 25)"),
        B("d", "decision", 1.55, 1.3, "How does it lie\nagainst the pillars?"),
        B("v1", "step", 0, 2.95, "Project it in 3D",
          "surface = project_fault_surface(\n    wall, mesh)",
          "X and Y faces only, with a K range\nworkflow 2 · lesson 05"),
        B("v2", "data", 0, 4.05, "Write it into the deck",
          "write_grid_include(grid, deck)\nfaults_text_3d({\"F9\": surface}, mesh)",
          "or write_faults_3d + INCLUDE"),
        B("v3", "step", 0, 5.2, "Slide the deck",
          "slide_deck(deck, Throw(\"F9\", 10.0,\n    hanging_side=EAST), \"out\")",
          "a K range becomes Throw.layers:\nblind tips stay welded · workflow 4"),
        B("k", "check", 1.62, 2.95, "Move the pillars onto it first",
          ref="its 3D staircase cannot be slid: the X faces\nfall on several pillar lines and "
              "the Z steps\nhave no pillar to slide along (SlideError)"),
        B("a", "step", 1.1, 4.4, "A · Fit the pillars to it",
          "a = adapt_grid(grid, plane,\n    top_trace=\"auto\",\n    boundary=\"slide\")\n"
          "slide_grid(a.grid, Throw(\n    a.faults_after, 15.0, ...))",
          "a polygon or facets · lesson 30 C\nworkflow 5 A"),
        B("b", "step", 2.2, 4.4, "B · Tilt a pillar line to its angle",
          "angle = dip of fit_plane(plane)\nfaces = fault_faces_for(grid,\n"
          "  [surface_top_trace(grid, plane)])\ndip_grid(grid, Dip(faces, angle,\n"
          "    pivot=\"top\"), slide=Throw(...))",
          "lesson 31 I · workflows 5 B, 6"),
        B("h", "stop", 3.35, 2.95, "Not possible",
          ref="slide and dip both refuse a fault\nwhose faces are all Z: pillars run\n"
              "unbroken from top to base.\nKeep its 3D FAULTS as a barrier"),
        B("o", "end", 1.1, 6.3, "A slid grid",
          "res.grid · res.realized_throw(\"F\")",
          "check res.ok and res.unexpected (empty);\nfor a deck: workflow 8"),
    ]
    edges = [
        E("s", "d"),
        E("d", "v1", "along them\n(vertical)", via=2.2),
        E("d", "k", "across them\n(dipping)", via=2.2),
        E("d", "h", "flat\n(horizontal)", via=2.2),
        E("v1", "v2"), E("v2", "v3"),
        E("k", "a", "fit them", via=3.7), E("k", "b", "tilt them", via=3.7),
        E("v3", "o", via=5.85), E("a", "o", via=5.85), E("b", "o", via=5.85),
    ]
    flowchart(out / "flow-9-slide-3d.png", "Workflow 9 · Slide a 3D fault surface",
              "A slide moves corners along their pillars. A surface that lies along them "
              "slides as projected; one that crosses them needs the pillars moved first.",
              boxes, edges, xs=235, ys=86)


# ── computed illustrations ───────────────────────────────────────────────────────────────────

NX, NY, NZ, DX, DZ, TOP = 16, 8, 6, 50.0, 15.0, 2000.0
EAST = (700.0, 200.0)
ROW = 3                                      # the section: south faces of cell row 3 (y = 150)
PALETTE = ("#c6dbef", "#a1d99b", "#fdd0a2", "#dadaeb", "#fcbba1", "#d9f0a3")
RED, BLUE = "#b2182b", "#2166ac"
ZLIM = (TOP + NZ * DZ + 25.0, TOP - 5.0)    # one depth range for every section, so all compare


def cake():
    from fieldbook.zslide import layer_cake

    return layer_cake(NX, NY, NZ, DX, [DZ] * NZ, top=lambda x, y: TOP)


def x_fault(i, j0=1, j1=NY):
    from fault_projections import FaultFace

    return [FaultFace(i, j, "X") for j in range(j0, j1 + 1)]


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
            t = (zz - p[..., 2]) / (p[..., 5] - p[..., 2])
            out[..., u, v, :, 0] = p[..., 0] + t * (p[..., 3] - p[..., 0])
            out[..., u, v, :, 1] = p[..., 1] + t * (p[..., 4] - p[..., 1])
            out[..., u, v, :, 2] = zz
    return out


def moved(before, after):
    return ~(corners(before) == corners(after)).all(axis=(3, 4, 5, 6))


def face_segment(c, f):
    """Top-surface end points of FAULTS face ``f`` (1-based, Eclipse style)."""
    i, j = f.i - 1, f.j - 1
    if f.face == "X":
        return c[i, j, 0, 1, 0, 0, :2], c[i, j, 0, 1, 1, 0, :2]
    return c[i, j, 0, 0, 1, 0, :2], c[i, j, 0, 1, 1, 0, :2]


def plan(ax, grid, title, faults=(), trace=None, values=None, cmap=None, clim=None):
    """Top surface in plan: one quad per column, FAULTS in red, an authored trace dashed."""
    from matplotlib.collections import LineCollection, PolyCollection

    c = corners(grid)
    polys = [c[i, j, 0, [0, 1, 1, 0], [0, 0, 1, 1], 0, :2] for i in range(NX) for j in range(NY)]
    if values is None:
        pc = PolyCollection(polys, facecolors="#f0f0f0", edgecolors="#555555", linewidths=0.5)
    else:
        pc = PolyCollection(polys, array=np.asarray(values, float).ravel(), cmap=cmap,
                            edgecolors="#555555", linewidths=0.5)
        pc.set_clim(*clim)
    ax.add_collection(pc)
    if len(faults):
        ax.add_collection(LineCollection([face_segment(c, f) for f in faults], colors=RED,
                                         linewidths=2.6, capstyle="round"))
    if trace is not None:
        ax.plot(trace[:, 0], trace[:, 1], "--", color=BLUE, lw=1.6)
    ax.set_xlim(-15, NX * DX + 15)
    ax.set_ylim(-15, NY * DX + 15)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x (m)", fontsize=8, labelpad=1)
    ax.set_ylabel("y (m)", fontsize=8)
    ax.tick_params(labelsize=7)
    return pc


def section(ax, grid, title, *, pillars=(), plane=None, faces3d=None, zone=None, hit=None,
            window=(150.0, 650.0)):
    """South faces of cell row ROW in x-z, at true scale. ``pillars``: pillar lines i drawn red
    (a fault); ``plane``: an authored fault ((x, z), (x, z)) dashed blue; ``faces3d``: a
    FaultSurface whose X and Z faces in this row are drawn red; ``zone``: cells drawn pale
    outside it; ``hit``: cells hatched (they moved)."""
    from matplotlib.collections import LineCollection, PolyCollection

    c = corners(grid)
    polys, colours, hatched = [], [], []
    for i in range(NX):
        for k in range(NZ):
            q = c[i, ROW, k, :, 0, :, :]
            poly = [(p[0], p[2]) for p in (q[0, 0], q[1, 0], q[1, 1], q[0, 1])]
            polys.append(poly)
            colours.append(PALETTE[k] if zone is None or zone[i, ROW, k] else "#f4f4f4")
            if hit is not None and hit[i, ROW, k]:
                hatched.append(poly)
    ax.add_collection(PolyCollection(polys, facecolors=colours, edgecolors="#333333",
                                     linewidths=0.5))
    if hatched:
        ax.add_collection(PolyCollection(hatched, facecolors="none", edgecolors="#737373",
                                         linewidths=0, hatch="////"))
    segs = [[(c[i - 1, ROW, k, 1, 0, 0, 0], c[i - 1, ROW, k, 1, 0, 0, 2]),
             (c[i - 1, ROW, k, 1, 0, 1, 0], c[i - 1, ROW, k, 1, 0, 1, 2])]
            for i in pillars for k in range(NZ)]
    if faces3d is not None:
        for (a, _), d in zip(faces3d.cell_ijk, faces3d.direction):
            i, j, k = a
            if j != ROW:
                continue
            if d == 0:
                segs.append([(c[i, j, k, 1, 0, 0, 0], c[i, j, k, 1, 0, 0, 2]),
                             (c[i, j, k, 1, 0, 1, 0], c[i, j, k, 1, 0, 1, 2])])
            elif d == 2:
                segs.append([(c[i, j, k, 0, 0, 1, 0], c[i, j, k, 0, 0, 1, 2]),
                             (c[i, j, k, 1, 0, 1, 0], c[i, j, k, 1, 0, 1, 2])])
    if segs:
        ax.add_collection(LineCollection(segs, colors=RED, linewidths=2.6, capstyle="round"))
    if plane is not None:
        (xa, za), (xb, zb) = plane
        ax.plot([xa, xb], [za, zb], "--", color=BLUE, lw=1.6)
    assert np.nanmax(c[:, ROW, ..., 2]) < ZLIM[0]
    ax.set_xlim(*window)
    ax.set_ylim(*ZLIM)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x (m)", fontsize=8, labelpad=1)
    ax.set_ylabel("depth (m)", fontsize=8)
    ax.tick_params(labelsize=7)


def save(fig, path, title):
    import matplotlib.pyplot as plt

    fig.suptitle(title, fontsize=11)
    fig.savefig(path, dpi=110, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(path)


def compute():
    """Run every case once; the per-workflow figures and the overview draw from this."""
    from cpg2unstructured import (build_topology, check_validity, from_grdecl, pillar_lattice,
                                  write_grid_include)
    from fault_projections import (Dip, SlideError, Throw, adapt_grid, dip_grid,
                                   fault_faces_for, fault_zone, faults_text_3d,
                                   project_fault_surface, slide_deck, slide_grid)
    from fault_projections.from_faults import faults_from_deck

    grid = cake()
    r = {"grid": grid}
    # 1 + 3: one oblique interior trace
    trace = np.array([[170.0, 70.0], [610.0, 330.0]])
    stair = fault_faces_for(grid, [trace])
    s1 = adapt_grid(grid, trace, boundary="pinned", w0=0.05, polish="barrier",
                    on_invalid="raise")
    assert s1.ok and check_validity(pillar_lattice(s1.grid), pillar_lattice(grid)).ok
    r.update(trace=trace, stair=stair, stage1=s1)
    r["trace_len"] = float(np.linalg.norm(trace[1] - trace[0]))
    r["stair_len"] = DX * len(stair)
    c0, c1 = corners(grid), corners(s1.grid)
    r["stage1_plan_shift"] = float(np.abs(c1[..., :2] - c0[..., :2]).max())
    r["stage1_dz"] = float(np.abs(c1[..., 2] - c0[..., 2]).max())
    # 2: a 45 degree plane, projected as its top trace (2D) and as a surface (3D)
    xa, xb = 310.0, 410.0
    poly = np.array([[xa, -10.0, TOP - 5], [xa, NY * DX + 10, TOP - 5],
                     [xb, NY * DX + 10, TOP + 95], [xb, -10.0, TOP + 95]])
    mesh = build_topology(from_grdecl(grid))
    surf = project_fault_surface(poly, mesh)
    x_top = xa + 5.0                                     # where the plane meets the top
    flat = fault_faces_for(grid, [np.array([[x_top, -10.0], [x_top, NY * DX + 10]])])
    r.update(plane3d=((xa, TOP - 5), (xb, TOP + 95)), surf=surf, flat=flat,
             surf_counts=np.bincount(surf.direction, minlength=3))
    # 2, the limit: a 3D staircase cannot be slid; dip the pillars, then slide (workflow 6)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d45 = dip_grid(grid, Dip(x_fault(8), 45.0, hanging_side=EAST, name="F"),
                       slide=Throw(x_fault(8), 15.0, hanging_side=EAST, name="F"))
    assert d45.ok and d45.slide is not None and d45.slide.ok
    x0, _, z0, x1, _, z1 = np.asarray(d45.grid.coord, float).reshape(NY + 1, NX + 1, 6)[ROW, 8]
    ext = 5.0 / (z1 - z0)                                # reach 5 m past the top and base
    xa45, za45 = x0 - ext * (x1 - x0), z0 - 5.0
    xb45, zb45 = x1 + ext * (x1 - x0), z1 + 5.0
    poly45 = np.array([[xa45, -10.0, za45], [xa45, NY * DX + 10, za45],
                       [xb45, NY * DX + 10, zb45], [xb45, -10.0, zb45]])
    surf45 = project_fault_surface(poly45, mesh)
    with TemporaryDirectory() as tmp:
        deck = Path(tmp) / "deck.grdecl"
        write_grid_include(grid, deck, header="workflows: full 3D, then slide")
        with open(deck, "a") as fh:
            fh.write("\n" + faults_text_3d({"F": surf45}, mesh))
        seen = sorted({f.i for f in faults_from_deck(deck)["F"]})
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                slide_deck(deck, Throw("F", 15.0, hanging_side=EAST), Path(tmp) / "out")
            refused = None
        except SlideError as exc:
            refused = str(exc)
    assert refused is not None and len(seen) > 1
    poly45a = np.array([[xa45, -2000.0, za45], [xa45, 2000.0, za45], [xb45, 2000.0, zb45],
                        [xb45, -2000.0, zb45]])
    a45 = adapt_grid(grid, poly45a, top_trace="auto", boundary="slide", w0=0.05,
                     polish="barrier", on_invalid="raise")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a45s = slide_grid(a45.grid, Throw(a45.faults_after, 15.0, hanging_side=EAST, name="F"))
    assert a45.ok and a45s.ok and not a45s.unexpected
    r.update(plane45=((xa45, za45), (xb45, zb45)), surf45=surf45, dip45=d45,
             surf45_counts=np.bincount(surf45.direction, minlength=3), seen45=seen,
             refused45=refused, throw45=float(np.nanmedian(d45.slide.realized_throw("F")[:, 1])),
             adapt45=a45, adapt45_slid=a45s,
             adapt45_lines=sorted({f.i for f in a45.faults_after if f.face == "X"}),
             adapt45_throw=float(np.nanmedian(a45s.realized_throw("F")[:, 1])))
    r.update(compute_9(grid, mesh))
    # 3 stage 2: an interior dipping plane, top x = 300, base x = 380
    poly2 = np.array([[300.0, 40.0, TOP - 5], [300.0, 360.0, TOP - 5],
                      [380.0, 360.0, TOP + 95], [380.0, 40.0, TOP + 95]])
    s2 = adapt_grid(grid, poly2, boundary="pinned", w0=0.05, polish="barrier",
                    regenerate_faults=False, on_invalid="raise")
    assert s2.ok
    lat0, lat2 = pillar_lattice(grid), pillar_lattice(s2.grid)
    r.update(stage2=s2, plane2=((300.0, TOP - 5), (380.0, TOP + 95)),
             stage2_tilt=float(np.abs((lat2.top - lat0.top) - (lat2.bot - lat0.bot)).max()))
    # 4: a vertical fault on pillar line 8 (x = 400), rim to rim; and one that ends inside
    faces = x_fault(8)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sl = slide_grid(grid, Throw(faces, 20.0, hanging_side=EAST, name="F"))
        part = slide_grid(grid, Throw(x_fault(8, 1, 5), 20.0, hanging_side=EAST, name="P"))
    key = {(f.i, f.j, f.face) for f in faces}
    assert sl.ok and part.ok and {(f.i, f.j, f.face) for f in sl.opened} <= key
    r.update(faces=faces, slide=sl, partial=part, partial_faces=x_fault(8, 1, 5))
    r["partial_max"] = float((corners(part.grid)[..., 2] - corners(grid)[..., 2]).max())
    r["slide_coord_same"] = bool(np.array_equal(np.asarray(sl.grid.coord),
                                                np.asarray(grid.coord)))
    # 5a: adapt to a rim-to-rim 60 degree plane (top trace auto) and slide 15 m
    tan60 = np.tan(np.radians(60.0))
    za, zb, x0 = TOP - 10, TOP + 100, 330.0
    xa5, xb5 = x0 + (za - TOP) / tan60, x0 + (zb - TOP) / tan60
    poly5 = np.array([[xa5, -2000.0, za], [xa5, 2000.0, za], [xb5, 2000.0, zb],
                      [xb5, -2000.0, zb]])
    a5 = adapt_grid(grid, poly5, top_trace="auto", boundary="slide", w0=0.05, polish="barrier",
                    on_invalid="raise")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a5s = slide_grid(a5.grid, Throw(a5.faults_after, 15.0, hanging_side=EAST, name="F"))
        # 5b and 6: dip the vertical fault to 60 degrees, with and without a slide
        d6 = dip_grid(grid, Dip(faces, 60.0, hanging_side=EAST, name="F"))
        d6f = dip_grid(grid, Dip(faces, 60.0, hanging_side=EAST, name="F", decay=150.0))
        d5 = dip_grid(grid, Dip(faces, 60.0, hanging_side=EAST, name="F"),
                      slide=Throw(faces, 15.0, hanging_side=EAST, name="F"))
    assert a5.ok and a5s.ok and d6.ok and d6f.ok and d5.ok and d5.slide is not None
    r.update(plane5=((xa5, za), (xb5, zb)), adapt5=a5, adapt5_slid=a5s, dip=d6, dip_fade=d6f,
             dip_slide=d5)
    # 7: the cells within 100 m, and a 60 degree dip kept inside them
    z = fault_zone(grid, {"F": faces}, 100.0, hanging_side=EAST)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dz = dip_grid(grid, Dip(faces, 60.0, hanging_side=EAST, name="F"), zone=z)
    assert dz.ok
    r.update(zone=z, dip_zone=dz, moved_all=moved(grid, d6.grid), moved_zone=moved(grid, dz.grid))
    assert not (r["moved_zone"] & ~z.cells).any()
    return r


def compute_9(grid, mesh):
    """Workflow 9: slide a 3D fault surface, by how it lies against the (vertical) pillars."""
    from cpg2unstructured import fit_plane, write_grid_include
    from fault_projections import (Dip, DipError, SlideError, Throw, adapt_grid, dip_deck,
                                   dip_grid, fault_faces_for, faults_text_3d,
                                   project_fault_surface, slide_deck, slide_grid,
                                   surface_top_trace)
    from fault_projections.from_faults import faults_from_deck

    r = {}
    span = (-10.0, NY * DX + 10)

    def wall(x0, z0, x1, z1, y=span):
        return np.array([[x0, y[0], z0], [x0, y[1], z0], [x1, y[1], z1], [x1, y[0], z1]])

    def deck_with(tmp, name, surface):
        path = Path(tmp) / f"{name}.grdecl"
        write_grid_include(grid, path, header=f"workflow 9: {name}")
        with open(path, "a") as fh:
            fh.write("\n" + faults_text_3d({name: surface}, mesh))
        return path

    def refusal(call, *args, exc):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                call(*args)
        except exc as e:
            return str(e)
        return None

    with TemporaryDirectory() as tmp:
        # along the pillars: a vertical wall through every layer, and one through layers 1-3
        for name, zb in (("V", TOP + 95), ("B", TOP + 45)):
            s = project_fault_surface(wall(410.0, TOP - 5, 410.0, zb), mesh)
            n = np.bincount(s.direction, minlength=3)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = slide_deck(deck_with(tmp, name, s), Throw(name, 10.0, hanging_side=EAST),
                                 Path(tmp) / name)
            assert n[2] == 0 and res.ok and not res.unexpected
            prof = [float(np.nanmedian(res.realized_throw(name, h)[:, 1])) for h in range(NZ + 1)]
            r[f"v9_{name}"] = dict(plane=((410.0, TOP - 5), (410.0, zb)), surf=s, counts=n,
                                   res=res, prof=prof,
                                   k=sorted({int(c[0][2]) + 1 for c in s.cell_ijk}))
        assert r["v9_V"]["prof"][0] == r["v9_V"]["prof"][-1] == 10.0
        assert r["v9_B"]["k"] == [1, 2, 3] and r["v9_B"]["prof"][0] == 10.0
        assert r["v9_B"]["prof"][3] == r["v9_B"]["prof"][-1] == 0.0
        # across the pillars: a 45 degree plane, top x = 350 (on pillar line 7), base x = 440
        za, zb, xa = TOP - 5, TOP + 95, 345.0
        xb = xa + (zb - za)
        s45 = project_fault_surface(wall(xa, za, xb, zb), mesh)
        d = deck_with(tmp, "D", s45)
        seen = sorted({f.i for f in faults_from_deck(d)["D"]})
        msg = refusal(slide_deck, d, Throw("D", 15.0, hanging_side=EAST), Path(tmp) / "D",
                      exc=SlideError)
        assert msg and len(seen) > 1
        r.update(d9_plane=((xa, za), (xb, zb)), d9_surf=s45, d9_seen=seen, d9_refused=msg,
                 d9_counts=np.bincount(s45.direction, minlength=3))
        # flat: a horizontal surface at the base of layer 3, x 200-600: Z faces only
        sh = project_fault_surface(np.array([[200.0, span[0], TOP + 45], [600.0, span[0], TOP + 45],
                                             [600.0, span[1], TOP + 45],
                                             [200.0, span[1], TOP + 45]]), mesh)
        d = deck_with(tmp, "H", sh)
        h_slide = refusal(slide_deck, d, Throw("H", 10.0, hanging_side=EAST), Path(tmp) / "H",
                          exc=SlideError)
        h_dip = refusal(dip_deck, d, Dip("H", 45.0, hanging_side=EAST), Path(tmp) / "H",
                        exc=DipError)
        nh = np.bincount(sh.direction, minlength=3)
        assert h_slide and h_dip and nh[0] == nh[1] == 0 and nh[2] > 0
        r.update(h9_surf=sh, h9_plane=((200.0, TOP + 45), (600.0, TOP + 45)), h9_counts=nh,
                 h9_slide=h_slide, h9_dip=h_dip)
    # route A: fit the pillars to the plane (top trace auto), then slide 15 m
    plane = wall(xa, za, xb, zb, y=(-2000.0, 2000.0))
    a = adapt_grid(grid, plane, top_trace="auto", boundary="slide", w0=0.05, polish="barrier",
                   on_invalid="raise")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a_s = slide_grid(a.grid, Throw(a.faults_after, 15.0, hanging_side=EAST, name="F"))
    lines = sorted({f.i for f in a.faults_after if f.face == "X"})
    assert a.ok and a_s.ok and not a_s.unexpected
    assert abs(angle_of(a_s.grid, lines[0]) - 45.0) < 0.5   # needs cpg2python PR #50
    # route B: the plane's angle (fit_plane) on the pillar line under its top trace
    _, normal, _, _ = fit_plane(plane)
    angle = float(np.degrees(np.arccos(abs(normal[2]))))
    faces = fault_faces_for(grid, [surface_top_trace(grid, plane)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = dip_grid(grid, Dip(faces, angle, hanging_side=EAST, pivot="top", name="F"),
                     slide=Throw(faces, 15.0, hanging_side=EAST, name="F"))
    bi = faces[0].i
    x0, _, z0, x1, _, z1 = np.asarray(b.grid.coord, float).reshape(NY + 1, NX + 1, 6)[ROW, bi]
    assert b.ok and b.slide is not None and b.slide.ok and {f.i for f in faces} == {bi}
    assert abs(x0 - (xa + z0 - za)) < 1e-6 and abs(x1 - (xa + z1 - za)) < 1e-6
    # route B snaps: a plane meeting the top at x = 365 tilts pillar line 7 (x = 350)
    off = wall(360.0, za, 460.0, zb, y=(-2000.0, 2000.0))
    off_faces = fault_faces_for(grid, [surface_top_trace(grid, off)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bo = dip_grid(grid, Dip(off_faces, 45.0, hanging_side=EAST, pivot="top", name="F"))
    xo = np.asarray(bo.grid.coord, float).reshape(NY + 1, NX + 1, 6)[ROW, off_faces[0].i]
    # route A also takes facets: 60 degrees down to z = 2045 m, then 45 degrees
    t60 = np.tan(np.radians(60.0))
    xm = 330.0 + (TOP + 45 - (TOP - 10)) / t60
    fy = (-2000.0, 2000.0)
    facets = [wall(330.0, TOP - 10, xm, TOP + 45, y=fy),
              wall(xm, TOP + 45, xm + 55.0, TOP + 100, y=fy)]
    af = adapt_grid(grid, facets, top_trace="auto", boundary="slide", w0=0.05,
                    polish="barrier", on_invalid="raise")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        afs = slide_grid(af.grid, Throw(af.faults_after, 15.0, hanging_side=EAST, name="F"))
    assert af.ok and afs.ok and not afs.unexpected
    r.update(a9=a, a9_slid=a_s, a9_lines=lines, a9_angle=angle_of(a_s.grid, lines[0]),
             a9_throw=float(np.nanmedian(a_s.realized_throw("F")[:, 1])),
             b9=b, b9_line=bi, b9_angle=angle, b9_top=(x0, z0), b9_base=(x1, z1),
             b9_throw=float(np.nanmedian(b.slide.realized_throw("F")[:, 1])),
             b9_offset=float(365.0 - xo[0]), b9_offset_line=off_faces[0].i,
             f9_throw=float(np.nanmedian(afs.realized_throw("F")[:, 1])))
    return r


def angle_of(grid, i=8):
    """Dip of fault pillar (i, ROW) from its COORD numbers."""
    x0, _, z0, x1, _, z1 = np.asarray(grid.coord, float).reshape(NY + 1, NX + 1, 6)[ROW, i]
    return float(np.degrees(np.arctan2(z1 - z0, x1 - x0)))


def illustrations(out, r):
    import matplotlib.pyplot as plt

    g = r["grid"]
    # 1 stair-step
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 3.9), layout="constrained")
    plan(axs[0], g, "BEFORE: the grid and the authored trace (dashed)", trace=r["trace"])
    plan(axs[1], g, f"AFTER: {len(r['stair'])} FAULTS faces (red) on the SAME grid\n"
         f"staircase {r['stair_len']:.0f} m for a {r['trace_len']:.0f} m trace "
         f"({r['stair_len'] / r['trace_len']:.2f}x)", faults=r["stair"], trace=r["trace"])
    save(fig, out / "ill-1-stair-step.png", "1 · Stair-step projection: "
         "fault_faces_for(grid, [trace]); no corner moves")
    # 2 full 3D
    fig, axs = plt.subplots(3, 1, figsize=(9.0, 7.6), layout="constrained")
    section(axs[0], g, "BEFORE: a 45 deg fault plane (dashed), section along y = 150 m",
            plane=r["plane3d"])
    section(axs[1], g, f"2D: its top trace projected: {len(r['flat'])} FAULTS faces, each "
            "through all layers (vertical)", pillars=sorted({f.i for f in r["flat"]}),
            plane=r["plane3d"])
    n = r["surf_counts"]
    section(axs[2], g, f"3D: project_fault_surface: {n[0]} X + {n[1]} Y + {n[2]} Z faces; "
            "the staircase follows the dip", faces3d=r["surf"], plane=r["plane3d"])
    save(fig, out / "ill-2-full-3d.png", "2 · Full 3D projection: the grid is unchanged, "
         "only the chosen faces differ")
    # 2, the limit: staircase vs dip + slide
    fig, axs = plt.subplots(4, 1, figsize=(9.0, 10.0), layout="constrained")
    n = r["surf45_counts"]
    section(axs[0], g, f"3D projection of a 45 deg plane (dashed): {n[0]} X + {n[2]} Z faces; "
            "the grid is unchanged", faces3d=r["surf45"], plane=r["plane45"])
    seen = r["seen45"]
    section(axs[1], g, f"slide_deck on that FAULTS: Z faces dropped, the X faces become "
            f"{len(seen)} full-height lines\n(i = {', '.join(map(str, seen))}): SlideError, "
            "no grid written", pillars=seen, plane=r["plane45"])
    axs[1].text(0.985, 0.08, "REFUSED", transform=axs[1].transAxes, ha="right",
                fontsize=13, weight="bold", color=RED,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=RED, lw=1.2))
    lines, a45s = r["adapt45_lines"], r["adapt45_slid"]
    section(axs[2], a45s.grid, f"workflow 5 A: adapt_grid(polygon, top_trace=\"auto\") + "
            f"slide_grid(Throw 15):\nthe pillars fitted to the plane "
            f"({angle_of(a45s.grid, lines[0]):.1f} deg), throw {r['adapt45_throw']:.1f} m",
            pillars=lines, plane=r["plane45"])
    d45 = r["dip45"]
    section(axs[3], d45.slide.grid, f"workflow 5 B / 6: dip_grid(Dip 45, slide=Throw 15): "
            f"{angle_of(d45.slide.grid):.1f} deg, throw {r['throw45']:.1f} m,\nthe hanging "
            "wall slides along the tilted pillars", pillars=[8], plane=r["plane45"])
    save(fig, out / "ill-2-slide-limit.png", "2 · Limit: a 3D staircase cannot be slid; "
         "move the pillars onto the fault, then slide")
    # 3 deformation
    fig = plt.figure(figsize=(12.5, 5.5), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[3.2, 1.45])
    plan(fig.add_subplot(gs[0, 0]), g, "BEFORE, plan: the grid and a map trace (dashed)",
         trace=r["trace"])
    s1 = r["stage1"]
    plan(fig.add_subplot(gs[0, 1]), s1.grid, f"AFTER stage 1, plan: pillars moved up to "
         f"{r['stage1_plan_shift']:.0f} m onto the trace;\n{len(s1.faults_after)} FAULTS "
         f"faces regenerated on it; depths change by {r['stage1_dz']:.0f} m",
         faults=s1.faults_after, trace=r["trace"])
    section(fig.add_subplot(gs[1, 0]), g, "BEFORE, section: a dipping plane (dashed)",
            plane=r["plane2"])
    section(fig.add_subplot(gs[1, 1]), r["stage2"].grid, "AFTER stage 2, section: pillar tops "
            f"and bases moved apart\n(up to {r['stage2_tilt']:.0f} m); the pillars tilt into "
            "the plane",
            plane=r["plane2"])
    save(fig, out / "ill-3-deformation.png", "3 · 2D / 2.5D deformation: adapt_grid moves "
         "pillars (COORD); ZCORN depths stay")
    # 4 slide
    fig = plt.figure(figsize=(12.5, 5.6), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 3.2])
    section(fig.add_subplot(gs[0, 0]), g, "BEFORE: a vertical fault at x = 400 m (red), "
            "rim to rim", pillars=[8])
    section(fig.add_subplot(gs[0, 1]), r["slide"].grid, "AFTER slide_grid, 20 m: the east "
            "block drops as one (rigid)", pillars=[8])
    top0 = corners(g)[:, :, 0, :, :, 0, 2].mean(axis=(2, 3))
    top1 = corners(r["partial"].grid)[:, :, 0, :, :, 0, 2].mean(axis=(2, 3))
    ax = fig.add_subplot(gs[1, 0])
    pc = plan(ax, r["partial"].grid, "AFTER, a fault that ends inside the grid (red): how far "
              "each column's top dropped", faults=r["partial_faces"], values=top1 - top0,
              cmap="Reds", clim=(0, 20))
    cb = fig.colorbar(pc, cax=ax.inset_axes((1.02, 0.0, 0.03, 1.0)))
    cb.set_label("drop (m)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.text(0.02, 0.95, "What changed (measured)\n\n"
            f"rigid: COORD identical: {r['slide_coord_same']}\n"
            f"rigid: faces opened: {len(r['slide'].opened)}, all FAULTS faces\n"
            f"partial: the largest corner drop is {r['partial_max']:.0f} m, at the\n"
            f"fault; it fades to 0 at the tip and away from it\n(decay=\"auto\"; mode "
            f"{r['partial'].modes['P']})", va="top",
            fontsize=9, family="DejaVu Sans Mono")
    save(fig, out / "ill-4-slide.png", "4 · Fault slide: corners move down their pillars "
         "(ZCORN only)")
    # 5 slide + pillar shift
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 4.6), layout="constrained")
    section(axs[0, 0], g, "BEFORE A: no FAULTS; a 60 deg plane (dashed)", plane=r["plane5"])
    a5 = r["adapt5"]
    section(axs[0, 1], r["adapt5_slid"].grid, "AFTER A: adapt_grid(top_trace=\"auto\") + "
            "15 m slide", pillars=sorted({f.i for f in a5.faults_after if f.face == "X"}),
            plane=r["plane5"])
    section(axs[1, 0], g, "BEFORE B: a vertical fault (red)", pillars=[8])
    d5 = r["dip_slide"]
    section(axs[1, 1], d5.slide.grid, f"AFTER B: dip_grid(Dip 60, slide=Throw 15): "
            f"{angle_of(d5.slide.grid):.1f} deg", pillars=[8])
    save(fig, out / "ill-5-slide-pillars.png", "5 · Slide + pillar shift: the pillars move "
         "(COORD), then the hanging wall slides along them (ZCORN)")
    # 6 dip
    fig, axs = plt.subplots(3, 1, figsize=(9.0, 6.0), layout="constrained")
    w = (50.0, 750.0)
    section(axs[0], g, "BEFORE: a vertical fault (red)", pillars=[8], window=w)
    section(axs[1], r["dip"].grid, f"AFTER dip_grid(Dip(60)): {angle_of(r['dip'].grid):.1f}"
            " deg; rigid, the outline fixed (default)", pillars=[8], window=w)
    section(axs[2], r["dip_fade"].grid, f"AFTER dip_grid(Dip(60, decay=150)): "
            f"{angle_of(r['dip_fade'].grid):.1f} deg; the tilt fades out", pillars=[8],
            window=w)
    save(fig, out / "ill-6-dip.png", "6 · Set a dip angle: COORD only; ZCORN, ACTNUM and "
         "FAULTS unchanged")
    # 7 zones
    z = r["zone"]
    fig = plt.figure(figsize=(12.5, 9.0), layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[3.2, 2.2, 2.2])
    from matplotlib.colors import ListedColormap

    cm = ListedColormap(["#92c5de", "#f0f0f0", "#f4a582"])
    plan(fig.add_subplot(gs[0, 0]), g, f"SELECT: fault_zone(D = 100 m): "
         f"{int(z.minus.sum())} footwall + {int(z.plus.sum())} hanging-wall cells",
         faults=r["faces"], values=z.side()[:, :, 0], cmap=cm, clim=(-1.5, 1.5))
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    ma, mz = r["moved_all"], r["moved_zone"]
    ax.text(0.02, 0.95, "What moved (measured, every corner)\n\n"
            f"whole domain (default): {int(ma.sum())} of {ma.size} cells\n"
            f"inside the zone:        {int(mz.sum())} cells, "
            f"{int((mz & ~z.cells).sum())} outside it\n\n"
            "blue: footwall band (1 in FLTZONE)\norange: hanging-wall band (2)\n"
            "hatched cells in the sections moved", va="top", fontsize=9,
            family="DejaVu Sans Mono")
    w = (50.0, 750.0)
    section(fig.add_subplot(gs[1, :]), r["dip"].grid, "60 deg dip, whole domain (default): "
            "cells far from the fault move (hatched)", pillars=[8], hit=ma, window=w)
    section(fig.add_subplot(gs[2, :]), r["dip_zone"].grid, f"60 deg dip, zone=z: "
            f"{angle_of(r['dip_zone'].grid):.1f} deg; outside the zone (grey) nothing moves",
            pillars=[8], zone=z.cells, hit=mz, window=w)
    save(fig, out / "ill-7-zones.png", "7 · Fault zones: select the cells within ±D, and "
         "optionally keep an edit inside them")


def stamp(ax, text="REFUSED"):
    ax.text(0.985, 0.08, text, transform=ax.transAxes, ha="right", fontsize=13, weight="bold",
            color=RED, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=RED, lw=1.2))


def illustration_9(out, r):
    import matplotlib.pyplot as plt

    g = r["grid"]
    fig, axs = plt.subplots(5, 2, figsize=(12.5, 11.0), layout="constrained")
    v, b = r["v9_V"], r["v9_B"]
    section(axs[0, 0], g, f"BEFORE: a vertical wall at x = 410 m through every layer,\n"
            f"projected in 3D: {v['counts'][0]} X faces (red) at x = 400 m, no Z faces",
            faces3d=v["surf"], plane=v["plane"])
    section(axs[0, 1], v["res"].grid, f"AFTER slide_deck(Throw 10): the hanging wall drops "
            f"{v['prof'][0]:.1f} m\nat the top and {v['prof'][-1]:.1f} m at the base",
            faces3d=v["surf"])
    section(axs[1, 0], g, f"BEFORE: the same wall through layers 1–3 only (blind):\n"
            f"{b['counts'][0]} X faces, K = {b['k'][0]}..{b['k'][-1]} in its FAULTS",
            faces3d=b["surf"], plane=b["plane"])
    section(axs[1, 1], b["res"].grid, f"AFTER: the K range becomes Throw.layers = (1, 3): "
            f"{b['prof'][0]:.1f} m at the top,\n{b['prof'][3]:.1f} m from the base of layer 3 "
            "down (the blind tip stays welded)", faces3d=b["surf"])
    n, seen = r["d9_counts"], r["d9_seen"]
    section(axs[2, 0], g, f"BEFORE: a 45 deg plane (dashed), projected in 3D: {n[0]} X + "
            f"{n[2]} Z faces,\na staircase across the pillars", faces3d=r["d9_surf"],
            plane=r["d9_plane"])
    section(axs[2, 1], g, f"slide_deck on those FAULTS: the X faces fall on {len(seen)} "
            f"pillar lines\n(i = {', '.join(map(str, seen))}): SlideError, no grid written",
            pillars=seen, plane=r["d9_plane"])
    stamp(axs[2, 1])
    a = r["a9_slid"]
    section(axs[3, 0], a.grid, f"AFTER A: adapt_grid(plane, top_trace=\"auto\") + slide_grid"
            f"(Throw 15):\npillars fitted to the plane, {r['a9_angle']:.1f} deg; throw "
            f"{r['a9_throw']:.1f} m", pillars=r["a9_lines"], plane=r["d9_plane"])
    bg = r["b9"].slide.grid
    section(axs[3, 1], bg, f"AFTER B: dip_grid(Dip({r['b9_angle']:.1f}, pivot=\"top\"), slide="
            f"Throw 15):\n{angle_of(bg, r['b9_line']):.1f} deg; throw {r['b9_throw']:.1f} m",
            pillars=[r["b9_line"]], plane=r["d9_plane"])
    nh = r["h9_counts"]
    section(axs[4, 0], g, f"BEFORE: a horizontal surface at 2045 m (x 200–600 m), "
            f"projected in 3D:\n{nh[2]} Z faces (red), no X or Y faces", faces3d=r["h9_surf"],
            plane=r["h9_plane"])
    section(axs[4, 1], g, "slide_deck: SlideError; dip_deck: DipError (all its faces are Z);\n"
            "the grid is unchanged", faces3d=r["h9_surf"], plane=r["h9_plane"])
    stamp(axs[4, 1])
    save(fig, out / "ill-9-slide-3d.png", "9 · Slide a 3D fault surface: vertical slides as "
         "projected; dipping, move the pillars first; flat cannot slide")


# ── drawn schematics (not computed) ──────────────────────────────────────────────────────────

FW, HW = "#c6dbef", "#fdd0a2"
OK_C, NO_C = "#238b45", RED


def _poly(ax, pts, **kw):
    from matplotlib.patches import Polygon

    kw.setdefault("ec", "#333333")
    kw.setdefault("lw", 0.8)
    ax.add_patch(Polygon(pts, closed=True, **kw))


def _px(p, y, H):
    """x on pillar ``p`` = (x at the top y = 0, x at the base y = -H) at depth ``y``."""
    return p[0] + (p[1] - p[0]) * (-y) / H


def _cell(ax, pl, pr, k, h, H, s=0.0, **kw):
    """Cell of layer k between pillars pl and pr, moved down its pillars by s."""
    y0, y1 = -k * h - s, -(k + 1) * h - s
    _poly(ax, [(_px(pl, y0, H), y0), (_px(pr, y0, H), y0), (_px(pr, y1, H), y1),
               (_px(pl, y1, H), y1)], **kw)


def _pillar(ax, p, H, ext=0.25, **kw):
    kw.setdefault("color", "#7f7f7f")
    kw.setdefault("lw", 1.2)
    ax.plot([_px(p, ext, H), _px(p, -H - ext, H)], [ext, -H - ext], zorder=0, **kw)


def _arrow(ax, a, b, color=BLUE, lw=1.8):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                                    mutation_scale=13))


def _frame(ax, title, caption, colour="#1a1a1a"):
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=10, weight="bold", color=colour)
    ax.text(0.5, -0.02, caption, transform=ax.transAxes, ha="center", va="top", fontsize=8.3,
            color="#333333", linespacing=1.3)


def _block(ax, pillars, nz, h, hw_from, s, layer_colours=True):
    """Columns between consecutive pillars; columns >= hw_from moved down by s."""
    H = nz * h
    for c in range(len(pillars) - 1):
        for k in range(nz):
            moved = c >= hw_from
            fc = PALETTE[k] if layer_colours else (HW if moved else FW)
            _cell(ax, pillars[c], pillars[c + 1], k, h, H, s if moved else 0.0, fc=fc)
    for p in pillars:
        _pillar(ax, p, H)
    return H


def _staircase(ax, xline, nz, h, lw=2.6, colour=RED):
    """X faces at the pillar nearest the plane in each layer, joined by Z steps."""
    xs = [round(xline(-(k + 0.5) * h)) for k in range(nz)]
    for k, x in enumerate(xs):
        ax.plot([x, x], [-k * h, -(k + 1) * h], color=colour, lw=lw, solid_capstyle="round")
        if k + 1 < nz and xs[k + 1] != x:
            ax.plot([x, xs[k + 1]], [-(k + 1) * h] * 2, color=colour, lw=lw,
                    solid_capstyle="round")
    return xs


def why_pillars(out):
    """Why a corner-point slide needs the fault to lie along the pillars (drawn)."""
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 4, figsize=(15.0, 4.3), layout="constrained")
    h, s = 0.6, 0.3
    # a: vertical fault, vertical pillars
    ax = axs[0]
    P = [(0, 0), (1, 1), (2, 2)]
    H = _block(ax, P, 3, h, 1, s)
    ax.plot([1, 1], [0, -H - s], color=RED, lw=2.6)
    _arrow(ax, (1.18, -0.35), (1.18, -0.35 - s - 0.25))
    _frame(ax, "✓  A vertical fault", "corners move down their own pillar;\non a vertical "
           "pillar that is slip\nalong the fault (workflows 4, 9)", OK_C)
    ax.set_xlim(-1.3, 3.3)
    ax.set_ylim(-H - 0.75, 0.4)
    # b: the pillars tilted into the fault
    ax = axs[1]
    t = 0.9
    P = [(0, t), (1, 1 + t), (2, 2 + t)]
    H = _block(ax, P, 3, h, 1, s)
    ax.plot([_px(P[1], 0, H), _px(P[1], -H - s, H)], [0, -H - s], color=RED, lw=2.6)
    y0, y1 = -0.3, -0.3 - s - 0.3
    _arrow(ax, (_px(P[1], y0, H) + 0.2, y0), (_px(P[1], y1, H) + 0.2, y1))
    _frame(ax, "✓  Pillars tilted into the fault", "the pillar lies in the fault, so down "
           "the\npillar is down-dip slip (workflows 5, 6,\n9 A and B)", OK_C)
    ax.set_xlim(-0.85, 3.75)
    ax.set_ylim(-H - 0.75, 0.4)
    # c: a dipping staircase on vertical pillars
    ax = axs[2]
    nz, nx = 3, 4
    H = nz * h
    fault = [1, 2, 3]                            # the X face's pillar in each layer
    for c in range(nx):
        for k in range(nz):
            _cell(ax, (c, c), (c + 1, c + 1), k, h, H, fc=HW if c >= fault[k] else FW)
    for c in range(nx + 1):
        _pillar(ax, (c, c), H)

    def plane(y):
        return 0.75 + (-y) / h

    ax.plot([plane(0.1), plane(-H - 0.1)], [0.1, -H - 0.1], "--", color=BLUE, lw=1.4)
    _staircase(ax, plane, nz, h)
    _cell(ax, (1, 1), (2, 2), 0, h, H, s, fc="none", ec=NO_C, lw=1.3, ls="--")
    _poly(ax, [(1, -h), (2, -h), (2, -h - s), (1, -h - s)], fc="#fcbba1", ec=NO_C,
          hatch="xxx", lw=0, alpha=0.9)
    ax.text(2.06, -h - s / 2, "overlap", color=NO_C, fontsize=8.5, va="center", weight="bold")
    _arrow(ax, (1.5, -0.12), (1.5, -0.12 - s), color=NO_C)
    _frame(ax, "✗  A dipping staircase", "on a Z step, down the pillar runs into the\n"
           "cell below; the X faces sit on 3 pillar\nlines: SlideError (move the pillars)",
           NO_C)
    ax.set_xlim(-0.3, 4.3)
    ax.set_ylim(-H - 0.75, 0.4)
    # d: a horizontal fault
    ax = axs[3]
    nz, H = 2, 2 * h
    for c in range(3):
        for k in range(nz):
            _cell(ax, (c, c), (c + 1, c + 1), k, h, H, fc=HW if k == 0 else FW)
    for c in range(4):
        _pillar(ax, (c, c), H, lw=2.4, color="#525252")
    ax.plot([0, 3], [-h, -h], color=RED, lw=2.6)
    _arrow(ax, (0.6, -0.3), (2.4, -0.3), color=NO_C)
    ax.text(1.5, -0.18, "the slip", color=NO_C, fontsize=8.5, ha="center", va="bottom")
    for c in (1, 2):
        ax.text(c, -h, "✗", color=NO_C, fontsize=13, ha="center", va="center", weight="bold",
                bbox=dict(boxstyle="circle,pad=0.05", fc="white", ec="none"))
    _frame(ax, "✗  A horizontal fault", "the slip is sideways, but corners only move\nalong "
           "their pillars and a pillar is one\nline from top to base: no grid can do it", NO_C)
    ax.set_xlim(-0.8, 3.8)
    ax.set_ylim(-H - 1.35, 0.4)
    fig.suptitle("Why a slide needs the fault along the pillars (schematic; sections, depth "
                 "down)", fontsize=11)
    fig.savefig(out / "why-pillars.png", dpi=110, facecolor="white")
    plt.close(fig)
    print(out / "why-pillars.png")


def concepts(out):
    """One small schematic per workflow (drawn, not computed)."""
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(3, 3, figsize=(13.5, 10.6), layout="constrained")
    h = 0.6
    # 1 stair-step, plan
    ax = axs[0, 0]
    for i in range(5):
        for j in range(4):
            _poly(ax, [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)], fc="#f0f0f0", lw=0.6)
    stair = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (3, 2), (4, 2), (4, 3), (5, 3), (5, 4)]
    ax.plot(*zip(*stair), color=RED, lw=2.6)
    ax.plot([0, 5], [0.1, 3.9], "--", color=BLUE, lw=1.5)
    _frame(ax, "1 · Stair-step projection", "a map trace (dashed) becomes FAULTS on cell "
           "edges;\nthe grid does not move")
    # 2 full 3D, section
    ax = axs[0, 1]
    nz, H = 4, 4 * h
    for c in range(5):
        for k in range(nz):
            _cell(ax, (c, c), (c + 1, c + 1), k, h, H, fc=PALETTE[k])

    def plane2(y):
        return 0.9 + 0.75 * (-y) / h

    ax.plot([plane2(0.1), plane2(-H - 0.1)], [0.1, -H - 0.1], "--", color=BLUE, lw=1.5)
    _staircase(ax, plane2, nz, h)
    _frame(ax, "2 · Full 3D projection", "a 3D surface becomes X, Y and Z faces with\n"
           "real K ranges; the grid does not move")
    # 3 deformation, plan: the staircase nodes moved onto the trace
    ax = axs[0, 2]
    a, b = np.array([0.0, 0.1]), np.array([5.0, 3.9])
    nodes = np.array([[(i, j) for j in range(5)] for i in range(6)], float)
    for i, j in stair:
        p = nodes[i, j]
        t = np.clip((p - a) @ (b - a) / ((b - a) @ (b - a)), 0, 1)
        nodes[i, j] = a + t * (b - a)
    for i in range(5):
        for j in range(4):
            _poly(ax, [nodes[i, j], nodes[i + 1, j], nodes[i + 1, j + 1], nodes[i, j + 1]],
                  fc="#f0f0f0", lw=0.6)
    ax.plot(*zip(*[nodes[i, j] for i, j in stair]), color=RED, lw=2.6)
    _frame(ax, "3 · 2D / 2.5D deformation", "pillars move onto the trace (stage 1)\nor tilt "
           "into a dipping plane (stage 2)")
    # 4 slide
    ax = axs[1, 0]
    P = [(c, c) for c in range(5)]
    H = _block(ax, P, 3, h, 2, 0.35)
    ax.plot([2, 2], [0, -H - 0.35], color=RED, lw=2.6)
    _arrow(ax, (2.25, -0.3), (2.25, -0.85))
    _frame(ax, "4 · Fault slide", "corners move down their pillars (ZCORN):\nthe hanging "
           "wall drops by the throw")
    # 5 slide + pillar shift
    ax = axs[1, 1]
    t = 0.8
    P = [(c, c + t) for c in range(5)]
    H = _block(ax, P, 3, h, 2, 0.35)
    ax.plot([_px(P[2], 0, H), _px(P[2], -H - 0.35, H)], [0, -H - 0.35], color=RED, lw=2.6)
    _arrow(ax, (_px(P[2], -0.3, H) + 0.25, -0.3), (_px(P[2], -0.85, H) + 0.25, -0.85))
    _frame(ax, "5 · Slide + pillar shift", "fit or tilt the pillars (COORD), then slide\nthe "
           "hanging wall along them (ZCORN)")
    # 6 dip angle
    ax = axs[1, 2]
    lean = [0.0, 0.25, 0.8, 0.25, 0.0]
    P = [(c, c + lean[c]) for c in range(5)]
    H = _block(ax, P, 3, h, 9, 0.0)
    ax.plot([2, 2], [0, -H], ":", color="#555555", lw=1.4)
    ax.plot([2, 2 + lean[2]], [0, -H], color=RED, lw=2.6)
    from matplotlib.patches import Arc

    theta = np.degrees(np.arctan2(H, lean[2]))
    ax.add_patch(Arc((2, 0), 1.0, 1.0, theta1=-90, theta2=-theta, color=BLUE, lw=1.4))
    ax.text(2.3, -0.75, "θ", color=BLUE, fontsize=12, ha="center", va="center")
    _frame(ax, "6 · Set a dip angle", "tilt a fault's pillars to an angle (COORD only);\n"
           "ZCORN, ACTNUM and FAULTS stay")
    # 7 zones, plan
    ax = axs[2, 0]
    for i in range(8):
        for j in range(4):
            side = "#92c5de" if 2 <= i < 4 else ("#f4a582" if 4 <= i < 6 else "#f0f0f0")
            _poly(ax, [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)], fc=side, lw=0.6)
    ax.plot([4, 4], [0, 4], color=RED, lw=2.6)
    for x0, x1, lab in ((4, 2, "−D"), (4, 6, "+D")):
        ax.annotate("", xy=(x1, 4.3), xytext=(x0, 4.3),
                    arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.2))
        ax.text((x0 + x1) / 2, 4.55, lab, ha="center", va="bottom", fontsize=10)
    ax.set_ylim(-0.3, 5.2)
    _frame(ax, "7 · Fault zones", "select the cells within ±D of a fault (FLTZONE);\n"
           "optionally keep a dip or slide inside them")
    # 8 field deck
    ax = axs[2, 1]
    box = dict(boxstyle="round,pad=0.35", lw=1.2)
    items = [(0.5, 3.6, "FIELD.DATA + GRDECL\n(read only)", "#deebf7", "#3182bd"),
             (0.5, 2.55, "edit: workflows 3–7, 9", "#ffffff", "#404040"),
             (0.5, 1.5, "four gates · SVGs · dossier", "#efedf5", "#6a51a3"),
             (0.5, 0.45, "<prefix>_grid.grdecl\n<prefix>_faults.grdecl", "#e5f5e0", "#31a354")]
    for x, y, txt, fc, ec in items:
        ax.text(x, y, txt, ha="center", va="center", fontsize=9, family="DejaVu Sans",
                bbox=dict(fc=fc, ec=ec, **box))
    for y0, y1 in ((3.25, 2.85), (2.25, 1.8), (1.2, 0.85)):
        _arrow(ax, (0.5, y0), (0.5, y1), color="#555555", lw=1.2)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-0.1, 4.2)
    _frame(ax, "8 · A field deck, end to end", "the master deck is never rewritten:\nthe "
           "edit goes into separate INCLUDE files")
    # 9 slide a 3D surface: along, across, flat
    ax = axs[2, 2]
    H = 2.0
    for dx in (-0.45, 0.45):                         # vertical pillars around the wall
        ax.plot([0.4 + dx, 0.4 + dx], [0.2, -H - 0.2], color="#7f7f7f", lw=1.0)
    ax.plot([0.4, 0.4], [0, -H], color=RED, lw=2.6)
    _arrow(ax, (0.6, -0.4), (0.6, -1.2))
    ax.plot([1.75, 2.45], [0, -H], color=RED, lw=2.6)
    for dx in (-0.45, 0.45):
        ax.plot([1.75 + dx, 2.45 + dx], [0.2, -H - 0.2], color="#7f7f7f", lw=1.0, ls="--")
    _arrow(ax, (2.25, -0.45), (2.25 + 0.7 * 0.4, -0.45 - 0.8))
    for c in (3.3, 3.9, 4.5):
        ax.plot([c, c], [0.2, -H - 0.2], color="#7f7f7f", lw=1.0)
    ax.plot([3.1, 4.7], [-H / 2, -H / 2], color=RED, lw=2.6)
    for x, mark, col, lab in ((0.4, "✓", OK_C, "vertical"), (2.1, "✓", OK_C, "dipping:\npillars "
                              "moved"), (3.9, "✗", NO_C, "horizontal")):
        ax.text(x, 0.45, mark, ha="center", va="bottom", fontsize=15, color=col, weight="bold")
        ax.text(x, -H - 0.35, lab, ha="center", va="top", fontsize=8.5, color="#333333")
    ax.set_xlim(-0.3, 5.0)
    ax.set_ylim(-H - 1.0, 1.0)
    _frame(ax, "9 · Slide a 3D fault surface", "a slide moves corners along the pillars: it "
           "needs\nthe fault along them (see 'Why' on page 9)")
    for ax in axs.flat:
        ax.margins(0.06)
    fig.suptitle("The workflows in one picture each (schematic; plan views for 1, 3 and 7, "
                 "sections for the rest)", fontsize=11.5)
    fig.savefig(out / "concepts.png", dpi=110, facecolor="white")
    plt.close(fig)
    print(out / "concepts.png")


def overview(out, r):
    import matplotlib.pyplot as plt

    g = r["grid"]
    fig, axs = plt.subplots(5, 2, figsize=(12.5, 12.4), layout="constrained")
    plan(axs[0, 0], g, "1 · Stair-step: FAULTS on the unchanged grid", faults=r["stair"],
         trace=r["trace"])
    plan(axs[0, 1], r["stage1"].grid, "3 · Stage 1: pillars moved onto the trace",
         faults=r["stage1"].faults_after, trace=r["trace"])
    section(axs[1, 0], g, "2 · Full 3D: X and Z faces follow a 45 deg plane",
            faces3d=r["surf"], plane=r["plane3d"])
    section(axs[1, 1], r["stage2"].grid, "3 · Stage 2: pillars tilted into the plane",
            plane=r["plane2"])
    section(axs[2, 0], r["slide"].grid, "4 · Slide: 20 m, ZCORN only", pillars=[8])
    section(axs[2, 1], r["dip_slide"].slide.grid, "5 · Slide + pillar shift: 60 deg, 15 m",
            pillars=[8])
    section(axs[3, 0], r["dip"].grid, "6 · Dip angle: 60 deg, COORD only", pillars=[8])
    section(axs[3, 1], r["dip_zone"].grid, "7 · Zone: the same dip kept within 100 m",
            pillars=[8], zone=r["zone"].cells, hit=r["moved_zone"])
    b = r["v9_B"]
    section(axs[4, 0], b["res"].grid, "9 · A blind vertical 3D fault (layers 1–3) slid 10 m",
            faces3d=b["surf"])
    section(axs[4, 1], r["a9_slid"].grid, "9 · A dipping 3D plane: pillars fitted, slid 15 m",
            pillars=r["a9_lines"], plane=r["d9_plane"])
    save(fig, out / "at-a-glance.png", "The workflows at a glance: one synthetic grid "
         "(16 x 8 x 6 cells), the grids the calls return")


def main(assets: Path) -> int:
    import matplotlib

    matplotlib.use("Agg")
    out = assets / "workflows"
    out.mkdir(parents=True, exist_ok=True)
    for chart in (chooser, flow_1, flow_2, flow_3, flow_4, flow_5, flow_6, flow_7, flow_8,
                  flow_9):
        chart(out)
    why_pillars(out)
    concepts(out)
    r = compute()
    illustrations(out, r)
    illustration_9(out, r)
    overview(out, r)
    return 0


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "docs" / "assets"
    raise SystemExit(main(target))
