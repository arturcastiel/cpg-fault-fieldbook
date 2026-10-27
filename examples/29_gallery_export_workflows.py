"""Gallery 7: end-to-end deck workflows: adapt on disk, slide, hand to a simulator or a viewer.

Each case starts from a layer cake (24x16x4, 50 m columns, layers 15/25/20/30 m), adapts it to
fault traces, slides it by editing ZCORN only, and then follows the files the way a project
would: written, hashed, read back, viewed and written again:

- A simulator bundle: a graben grid delivered as a directory of INCLUDE files (GRID, FAULTS,
  MULTFLT per named fault, PORO/PERMX/NTG per layer), a VTK view of the processed topology and
  a JSON manifest with sha256 hashes, dimensions, throws and the pinned library revisions;
- B adapt_deck end to end: a deck on disk is adapted by adapt_deck into INCLUDE files, read
  back through a master deck, slid and written; the COORD diff stays inside the morph band
  around the fault, and the slide leaves COORD bit-identical;
- C VTK growth series: a fault that ends inside the grid grows from 0 to 30 m in six steps;
  every step is written as VTK (grid and fault faces) with a ParaView .pvd collection and
  .vtk.series files whose time is the throw, for animating fault growth;
- D round trip: a growth-fault deck on an undulating top is written, read, and written again
  byte for byte; its FAULTS and the re-processed contact areas survive unchanged
  (reproducibility).

Pass --out DIR to keep the decks (<case>_0_input.grdecl, <case>_1_slid.grdecl).
"""

import argparse
import hashlib
import importlib.metadata
import json
import platform
import re
import xml.etree.ElementTree as ElementTree
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
from cpg2unstructured import format_keyword, pillar_lattice, read_grdecl, write_grid_include
from cpg2unstructured.morph import support_mask
from fault_projections import adapt_deck, adapt_grid
from fault_projections.adapt import fault_faces_for
from fault_projections.faults import faults_text
from fault_projections.from_faults import faults_from_deck
from fault_projections.targets import targets_from_polylines
from fault_projections.vtk import write_grid_vtk, write_surface_vtk

from fieldbook.cli import upstream_revisions
from fieldbook.zslide import (Throw, analytic_contacts, blocks, column_of_point, compare,
                              face_columns, interior, layer_cake, named_chains, process,
                              realized_throw, slide_blocks, slide_field, unwelded_faces,
                              write_deck)

NX, NY, DX = 24, 16, 50.0
DZ = [15.0, 25.0, 20.0, 30.0]
NZ = len(DZ)
MORPH = dict(boundary="pinned", w0=0.05, polish="barrier", on_invalid="raise")
LAYERS = {"PORO": [0.18, 0.25, 0.12, 0.22], "PERMX": [50.0, 400.0, 5.0, 150.0],
          "NTG": [0.8, 0.95, 0.6, 0.9]}


def base_grid(top=lambda x, y: 1000.0 + 0.02 * x):
    return layer_cake(NX, NY, NZ, DX, DZ, top=top)


def adapted(grid, traces):
    res = adapt_grid(grid, traces, **MORPH)
    return res.grid, named_chains(res.grid, interior(res.grid, res.faults_after))


def props_text():
    """PORO, PERMX and NTG, one value per layer, in Eclipse order (i fastest, k slowest)."""
    return "".join(format_keyword(kw, np.repeat(v, NX * NY)) for kw, v in LAYERS.items())


def write_pair(name, grid, slid, faults, out, extra=""):
    for stage, g in (("0_input", grid), ("1_slid", slid)):
        write_deck(out / f"{name}_{stage}.grdecl", g, faults, header=f"{name}: {stage[2:]}",
                   extra=extra)
    return out / f"{name}_1_slid.grdecl"


def check(grid, slid, faces, volume_scale=1.0):
    """Assert the slide invariants; return both re-processed grids and the fault contact."""
    assert np.array_equal(slid.coord, grid.coord), "COORD must not change"
    assert unwelded_faces(slid) == set(faces), "exactly the FAULTS faces must open"
    before, after = process(grid), process(slid)
    dims = grid.cart_dims
    assert np.allclose(after.volume_grid(dims), before.volume_grid(dims) * volume_scale)
    return before, after, compare(after.contacts(faces), analytic_contacts(slid, faces))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def trace_distance(points, line):
    """Plan distance from each point to a polyline."""
    a, ab = line[:-1], np.diff(line, axis=0)
    t = np.clip(((points[:, None] - a) * ab).sum(-1) / (ab ** 2).sum(-1), 0.0, 1.0)
    return np.linalg.norm(points[:, None] - (a + t[..., None] * ab), axis=-1).min(1)


def vtk_counts(path):
    """Points and cell count of a legacy ASCII VTK file."""
    lines = Path(path).read_text().splitlines()
    head = next(n for n, line in enumerate(lines) if line.startswith("POINTS"))
    points = np.loadtxt(lines[head + 1:head + 1 + int(lines[head].split()[1])], ndmin=2)
    keyword = next(line.split() for line in lines if line.startswith(("CELLS", "POLYGONS")))
    return points, int(keyword[1])


def fault_face_ids(processed, faces):
    """Topology faces whose two cells sit across one of the FAULTS faces."""
    pairs = {frozenset(face_columns(f)) for f in faces}
    nb, ij = processed.topology.faces.neighbours, list(map(tuple, processed.ijk[:, :2].tolist()))
    return np.array([n for n in np.flatnonzero(processed.topology.faces.direction <= 1)
                     if nb[n, 1] >= 0 and frozenset({ij[nb[n, 0]], ij[nb[n, 1]]}) in pairs])


def cell_depth(processed):
    """Mean depth of each cell's eight corners."""
    return processed.topology.nodes[processed.topology.cell_corners, 2].mean(1)


def write_vtu(path, topology, cell_data):
    """XML VTK UnstructuredGrid (ASCII, full precision) of the hexahedra, with cell data."""
    hexes = topology.cell_corners[:, [0, 1, 3, 2, 4, 5, 7, 6]]     # (i, j, k) bits -> VTK order

    def array(name, values, kind=None, components=1):
        values = np.asarray(values)
        kind = kind or ("Int64" if values.dtype.kind in "iu" else "Float64")
        text = " ".join(map(repr, values.ravel().tolist()))
        return (f'<DataArray type="{kind}" Name="{name}" NumberOfComponents="{components}" '
                f'format="ascii">{text}</DataArray>\n')
    n = len(hexes)
    Path(path).write_text(
        '<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1">\n'
        f'<UnstructuredGrid><Piece NumberOfPoints="{len(topology.nodes)}" NumberOfCells="{n}">\n'
        f"<Points>{array('Points', topology.nodes, components=3)}</Points>\n<Cells>\n"
        + array("connectivity", hexes) + array("offsets", 8 * np.arange(1, n + 1))
        + array("types", np.full(n, 12), "UInt8") + "</Cells>\n<CellData>\n"
        + "".join(array(k, v) for k, v in cell_data.items())
        + "</CellData>\n</Piece></UnstructuredGrid>\n</VTKFile>\n")


def case_a(out):
    name = "A_simulator_bundle"
    west = np.array([[380.0, 0.0], [420.0, 400.0], [400.0, 800.0]])
    east = np.array([[780.0, 0.0], [760.0, 400.0], [800.0, 800.0]])
    grid, faults = adapted(base_grid(), [west, east])
    faces = [f for chain in faults.values() for f in chain]
    seed, throw = column_of_point(grid, (600.0, 400.0)), 25.0
    slid = slide_blocks(grid, faces, {seed: throw})
    _, after, contact = check(grid, slid, faces)
    multflt = {"F1": 0.1, "F2": 0.01}
    mult_text = "MULTFLT\n" + "".join(f"  '{k}' {m} /\n" for k, m in multflt.items()) + "/\n"
    write_pair(name, grid, slid, faults, out, extra=props_text() + "\n" + mult_text)

    bundle = out / name
    bundle.mkdir(exist_ok=True)
    write_grid_include(slid, bundle / "grid.grdecl", header=f"{name}: adapted and slid grid")
    (bundle / "faults.inc").write_text(faults_text(list(faults.values()), NZ, names=list(faults)))
    (bundle / "multflt.inc").write_text(mult_text)
    (bundle / "props.inc").write_text(props_text())
    parts = ("grid.grdecl", "faults.inc", "multflt.inc", "props.inc")
    (bundle / "case.grdecl").write_text("".join(f"INCLUDE\n  '{p}' /\n" for p in parts))
    write_grid_vtk(str(bundle / "topology.vtk"), after.topology)
    files = sorted(p.name for p in bundle.iterdir() if p.name != "manifest.json")
    label = blocks(grid, faces)
    manifest = {
        "case": name, "cart_dims": list(grid.cart_dims), "active_cells": int(slid.actnum.sum()),
        "throws": [{"seed_column": list(seed), "block": int(label[seed]), "throw_m": throw}],
        "faults": {k: {"faces": len(v), "multflt": multflt[k]} for k, v in faults.items()},
        "layers": {"dz_m": DZ, **LAYERS},
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "cpg2unstructured": importlib.metadata.version("cpg2unstructured"),
                     "fault-projections": importlib.metadata.version("fault-projections"),
                     "revisions": upstream_revisions()},
        "files": {f: {"bytes": (bundle / f).stat().st_size, "sha256": sha256(bundle / f)}
                  for f in files},
    }
    (bundle / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    # A receiver checks the bundle from disk only.
    manifest = json.loads((bundle / "manifest.json").read_text())
    assert set(manifest["files"]) == set(files) == set(parts) | {"case.grdecl", "topology.vtk"}
    assert all(sha256(bundle / f) == m["sha256"] for f, m in manifest["files"].items())
    back = read_grdecl(bundle / "grid.grdecl")
    assert np.array_equal(back.coord, slid.coord) and np.array_equal(back.zcorn, slid.zcorn)
    deck = read_grdecl(bundle / "case.grdecl")                      # follows the INCLUDEs
    assert np.array_equal(deck.zcorn, slid.zcorn)
    assert {"FAULTS", "MULTFLT"} <= set(deck.unrecognized)
    assert all(np.array_equal(deck.props[k], np.repeat(v, NX * NY)) for k, v in LAYERS.items())
    # faults_from_deck reads one file and does not follow INCLUDE: read the FAULTS include.
    named = faults_from_deck(bundle / "faults.inc")
    assert named == faults
    assert re.findall(r"^\s*'(\w+)'", (bundle / "multflt.inc").read_text(), re.M) == list(named)
    points, cells = vtk_counts(bundle / "topology.vtk")
    nodes = after.topology.nodes
    assert cells == len(after.topology.cell_corners) == manifest["active_cells"]
    error = np.abs(points - nodes)
    assert np.all(error <= 5e-6 * np.abs(nodes) + 1e-12)              # 6 significant digits
    size = sum(m["bytes"] for m in manifest["files"].values())
    counts = ", ".join(f"{k} {len(v)}" for k, v in named.items())
    mults = ", ".join(f"{k} {m:g}" for k, m in multflt.items())
    print(f"{name:20s} graben {throw:+g} m, contact {contact:,.0f} m2; {len(files)} files + "
          f"manifest ({size / 1e6:.2f} MB), sha256 verified; FAULTS {counts}; MULTFLT {mults}")
    print(f"{'':20s} read back: COORD/ZCORN bit-exact, PORO/PERMX/NTG per layer; VTK {cells} "
          f"hexahedra, points rounded by <= {error.max() * 100:.2f} cm (view only)")


def case_b(out, rings=3):
    name = "B_adapt_deck_e2e"
    work = out / name
    work.mkdir(exist_ok=True)
    raw = base_grid()
    trace = np.array([[430.0, 0.0], [520.0, 300.0], [560.0, 560.0], [540.0, 800.0]])
    stair = interior(raw, fault_faces_for(raw, trace))
    deck = write_deck(work / "raw.grdecl", raw, {"F1": stair}, header=f"{name}: staircase FAULTS")
    res = adapt_deck(deck, trace, work, prefix="adapted", fault_names=["F1"],
                     support_rings=rings, **MORPH)
    assert res.faces_invariant
    includes = (res.grid_path.name, res.faults_path.name)
    (work / "case.grdecl").write_text("".join(f"INCLUDE\n  '{p}' /\n" for p in includes))
    grid = read_grdecl(work / "case.grdecl")
    faces = interior(grid, faults_from_deck(res.faults_path)["F1"])
    assert faces == interior(res.grid, res.faults_after)
    assert np.array_equal(grid.coord, res.grid.coord) and np.array_equal(grid.zcorn, raw.zcorn)

    # COORD diff of the adaptation: only pillars in the morph band around the trace moved.
    moved = np.any(grid.coord != raw.coord, axis=1).reshape(NY + 1, NX + 1).T     # [i, j]
    band = support_mask(moved.shape, targets_from_polylines(raw, trace).ij, rings)
    assert moved.any() and not (moved & ~band).any()
    reach = trace_distance(pillar_lattice(raw).top[moved], trace).max()

    seed, throw = column_of_point(grid, (1000.0, 400.0)), 20.0
    slid = slide_blocks(grid, faces, {seed: throw})
    *_, contact = check(grid, slid, faces)
    back = read_grdecl(write_pair(name, grid, slid, {"F1": faces}, out))
    assert np.array_equal(back.coord, grid.coord) and np.array_equal(back.zcorn, slid.zcorn)
    changed = int(np.count_nonzero(back.zcorn != grid.zcorn))
    print(f"{name:20s} adapt_deck moved {moved.sum()} of {moved.size} pillars (COORD only), "
          f"all in the {rings}-ring band, <= {reach:.0f} m from the trace; FAULTS {len(stair)} "
          f"staircase -> {len(faces)} on the trace")
    print(f"{'':20s} slide east block {throw:+g} m: COORD bit-identical on disk, ZCORN "
          f"{changed:,} of {grid.zcorn.size:,} values changed; contact {contact:,.0f} m2")


def case_c(out):
    name = "C_vtk_growth_series"
    work = out / name
    work.mkdir(exist_ok=True)
    trace = np.array([[620.0, 0.0], [660.0, 250.0], [640.0, 520.0]])
    grid, faults = adapted(base_grid(), trace)
    faces, hanging = faults["F1"], (1000.0, 100.0)
    start = process(grid)
    v0, depth0 = start.volume_grid(grid.cart_dims), cell_depth(start)
    steps, contacts = [], []
    for n, throw in enumerate(np.linspace(0.0, 30.0, 6)):
        slid = slide_field(grid, Throw(faces=tuple(faces), throw=throw, hanging_side=hanging),
                           faults=faces)
        proc = process(slid)
        assert np.array_equal(slid.coord, grid.coord)
        assert unwelded_faces(slid) == (set(faces) if throw else set())
        assert np.allclose(proc.volume_grid(grid.cart_dims), v0), "volumes changed"
        assert np.isclose(realized_throw(grid, slid, faces, hanging)[:, 1].max(), throw)
        contacts.append(compare(proc.contacts(faces), analytic_contacts(slid, faces)))
        ids = fault_face_ids(proc, faces)
        write_grid_vtk(str(work / f"grid_{n}.vtk"), proc.topology)
        write_surface_vtk(str(work / f"fault_{n}.vtk"), proc.topology,
                          SimpleNamespace(face_ids=ids))      # only face_ids is read
        assert np.array_equal(proc.ijk, start.ijk)
        write_vtu(work / f"grid_{n}.vtu", proc.topology, {
            "layer": proc.ijk[:, 2] + 1, "volume": proc.volumes,
            "subsidence": cell_depth(proc) - depth0})
        steps.append((float(throw), n, len(proc.topology.cell_corners), len(ids)))
    write_pair(name, grid, slid, faults, out)

    # A .pvd may only list XML files (.vtu, ...); legacy .vtk series go in a .vtk.series file.
    rows = "".join(f'    <DataSet timestep="{t:g}" part="0" file="grid_{n}.vtu"/>\n'
                   for t, n, _, _ in steps)
    (work / "growth.pvd").write_text('<?xml version="1.0"?>\n<VTKFile type="Collection" '
                                     f'version="0.1">\n  <Collection>\n{rows}  </Collection>\n'
                                     "</VTKFile>\n")
    for kind in ("grid", "fault"):
        series = [{"name": f"{kind}_{n}.vtk", "time": t} for t, n, _, _ in steps]
        (work / f"{kind}.vtk.series").write_text(
            json.dumps({"file-series-version": "1.0", "files": series}, indent=1) + "\n")

    collection = ElementTree.parse(work / "growth.pvd").getroot()
    assert collection.get("type") == "Collection"
    listed = collection.find("Collection").findall("DataSet")
    assert [float(d.get("timestep")) for d in listed] == [s[0] for s in steps]
    for (t, n, cells, n_faces), item in zip(steps, listed):
        piece = ElementTree.parse(work / item.get("file")).getroot().find(".//Piece")
        assert int(piece.get("NumberOfCells")) == cells == vtk_counts(work / f"grid_{n}.vtk")[1]
        assert vtk_counts(work / f"fault_{n}.vtk")[1] == n_faces
    for kind in ("grid", "fault"):
        series = json.loads((work / f"{kind}.vtk.series").read_text())["files"]
        assert [s["time"] for s in series] == [s[0] for s in steps]
        assert all((work / s["name"]).is_file() for s in series)
    print(f"{name:20s} throw {' -> '.join(f'{s[0]:g}' for s in steps)} m: fault contact "
          f"{contacts[0]:,.0f} -> {contacts[-1]:,.0f} m2, volumes unchanged at every step")
    print(f"{'':20s} {len(steps)} steps x (grid .vtk, fault .vtk with {steps[-1][3]} faces, grid "
          f".vtu) + growth.pvd + 2 .vtk.series, time = throw")


def case_d(out):
    name = "D_round_trip"
    def top(x, y):
        return 1000.0 + 0.02 * x + 4.0 * np.sin(x / 170.0) * np.cos(y / 130.0)
    grid, faults = adapted(base_grid(top), np.array([[260.0, 0.0], [720.0, 800.0]]))
    faces = faults["F1"]
    seed = column_of_point(grid, (1100.0, 400.0))
    growth = np.array([6.0, 9.5, 13.0, 16.5, 21.0])     # per horizon: hanging layers thicken
    slid = slide_blocks(grid, faces, {seed: growth})
    label = blocks(grid, faces)
    scale = np.where((label == label[seed])[..., None], 1.0 + np.diff(growth) / DZ, 1.0)
    before, after, contact = check(grid, slid, faces, volume_scale=scale)
    first = write_pair(name, grid, slid, faults, out, extra=props_text())

    back = read_grdecl(first)
    for key in ("coord", "zcorn", "actnum"):
        assert np.array_equal(getattr(back, key), getattr(slid, key)), key
    assert all(np.array_equal(back.props[k], np.repeat(v, NX * NY)) for k, v in LAYERS.items())
    named = faults_from_deck(first)
    assert named == faults
    again = write_deck(out / f"{name}_2_rewritten.grdecl", back, named, header=f"{name}: slid",
                       extra="".join(format_keyword(k, back.props[k]) for k in LAYERS))
    assert again.read_bytes() == first.read_bytes(), "the second generation must be identical"
    assert faults_from_deck(again) == faults
    reread = process(again)
    assert reread.contacts(faces) == after.contacts(faces), "contacts must survive exactly"
    assert np.array_equal(reread.volumes, after.volumes)
    dv = after.volumes.sum() - before.volumes.sum()
    print(f"{name:20s} write -> read -> write: byte-identical ({sha256(first)[:12]}...), "
          f"COORD/ZCORN/ACTNUM/props bit-exact, FAULTS F1 {len(named['F1'])} faces")
    print(f"{'':20s} growth fault {growth[0]:g} -> {growth[-1]:g} m: hanging-wall volume "
          f"{dv:+,.0f} m3 as predicted; {len(after.contacts(faces))} contact pairs "
          f"({contact:,.0f} m2) identical after the round trip")


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    case_a(out)
    case_b(out)
    case_c(out)
    case_d(out)
    print("every case: files read back from disk reproduce the grid, FAULTS and contacts; "
          "COORD never changes in a slide")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
