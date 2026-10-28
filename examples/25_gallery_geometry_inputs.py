"""Gallery 3: from fault data to a slid grid - traces, sticks, clouds and meshes adapt the grid.

A fault reaches a model as a map trace, as sticks picked on seismic sections, as a point cloud,
or as a triangulated / quad surface. Each case turns one of these inputs into FAULTS on a
24x16x4 layer cake (50 m columns, layers 15/25/20/30 m), with adapt_grid moving the pillars
onto the fault, and then slides a block by editing ZCORN only:

- A arc trace: the same curved trace as a staircase on the regular grid and as an adapted
  pillar chain; the staircase overstates the fault length (and its contact area) by 16 %;
- B fault sticks: noisy XYZ sticks on a 65 degree plane -> fit_plane -> a planar polygon that
  tilts the pillars into the fault; the slide along leaning pillars changes cell volumes, by
  the wedge the fault geometry predicts (within 1 %);
- C point cloud: a curved, east-dipping fault sampled as 3000 noisy points, faceted into
  planar quads that tilt the pillars; the projected 3D FAULTS cover the slid columns;
- D quad mesh: a kinked two-panel fault mesh; on the regular grid it projects to a 3D staircase
  (X, Y and Z faces) that no column slide can open, on the adapted grid to whole columns;
- E dome with a crestal graben: two traces over the crest of a dome; the graben drops 20 m and
  the crest moves to the footwall shoulders (an exploration trap that the fault splits);
- F dipping layers: a partial normal fault in layers dipping 8.5 degrees (7.5 m per column);
  a throw field that dies out at both tips. Unsnapped, horizons cross inside fault faces and
  the faces are cut at intersection nodes; the snap keeps layer boundaries on the pillars.

Each case asserts that COORD stays put and exactly the FAULTS faces open, plus what the input
kind lets us check: lengths and contact areas, plane fit and distance, volumes, coverage of the
projected surface, or the realized throw. Pass --out DIR to keep the decks
(<case>_0_input.grdecl, <case>_1_slid.grdecl).
"""

import argparse
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import (build_topology, cell_quality, fit_plane, from_grdecl,
                              pillar_lattice, surface_boundary_loops)
from fault_projections import (QuadMesh, adapt_grid, facet_surface, faults_text_3d,
                               project_faceted_surface, project_quad_mesh_surface,
                               surface_to_fault_faces, validate_quad_mesh)
from fault_projections.adapt import fault_faces_for
from fault_projections.faceting import sample_cloud_grid
from fault_projections.from_faults import face_edge, trace_faults

from fieldbook.zslide import (Throw, analytic_contacts, blocks, column_of_point, compare,
                              corner_cube, crossed_faces, face_columns, faces_of, interior,
                              layer_cake, named_chains, process, realized_throw, slide_blocks,
                              slide_field, unwelded_faces, with_cube, write_deck)

NX, NY, DX = 24, 16, 50.0
DZ = [15.0, 25.0, 20.0, 30.0]
H = sum(DZ)
# "slide" lets rim pillars move along the rim, so the chain reaches the fault at y = 0 and 800 m.
MORPH = dict(boundary="slide", w0=0.05, polish="barrier", on_invalid="raise")


def cake(top):
    return layer_cake(NX, NY, len(DZ), DX, DZ, top=top)


def flat(x, y):
    return 1000.0


def run(name, grid, faults, throws, out, extra=""):
    """Slide (seeded blocks or a Throw field), write both decks, assert exactly the FAULTS open."""
    faces = faces_of(faults)
    if isinstance(throws, Throw):
        slid = slide_field(grid, throws, faults=faces)
    else:
        slid = slide_blocks(grid, faces, throws)
    assert np.array_equal(slid.coord, grid.coord), "COORD must not move"
    assert unwelded_faces(slid) == set(faces), "exactly the FAULTS faces must open"
    for n, (g, tag) in enumerate(((grid, "input"), (slid, "slid"))):
        write_deck(out / f"{name}_{n}_{tag}.grdecl", g, faults, header=f"{name}: {tag}",
                   extra=extra)
    return slid


def volumes(grid):
    return process(grid).volume_grid(grid.cart_dims)


def vertical_checks(grid, slid, faces):
    """Vertical pillars: no cell changes volume, contacts equal the analytic overlaps."""
    assert np.allclose(volumes(grid), volumes(slid)), "volumes changed"
    return compare(process(slid).contacts(faces), analytic_contacts(slid, faces))


def plan_length(grid, faces):
    top = pillar_lattice(grid).top
    return sum(float(np.linalg.norm(top[face_edge(f)[0]] - top[face_edge(f)[1]])) for f in faces)


def pillar_ends(grid, ij):
    """(top xyz, base xyz) of the pillars ij."""
    lat = pillar_lattice(grid)
    i, j = np.asarray(ij).T
    return [np.c_[xy[i, j], z[i, j]] for xy, z in ((lat.top, lat.z_top), (lat.bot, lat.z_bot))]


def two_stage(grid, trace, surface):
    """Stage 1 snaps a pillar chain to the top trace (and gives the FAULTS), stage 2 tilts it."""
    r1 = adapt_grid(grid, trace, **MORPH)
    faces = interior(r1.grid, r1.faults_after)
    r2 = adapt_grid(r1.grid, surface, taper_cells=0.0, **MORPH)
    return r1.grid, r2.grid, faces, trace_faults(r1.grid, faces).chains[0].ij


def sj_min(grid):
    lat = pillar_lattice(grid)
    return min(cell_quality(lat).sj_min, cell_quality(lat, plane="bot").sj_min)


def hanging_wall_change(grid, slid, faces, seed):
    """Relative volume change per cell; asserts the footwall is untouched."""
    v0, v1 = volumes(grid), volumes(slid)
    label = blocks(grid, faces)
    hw = (label == label[seed])[:, :, None] & np.ones(len(DZ), bool)
    assert np.allclose(v0[~hw], v1[~hw], rtol=1e-12, atol=0.0), "footwall volumes changed"
    return (v1 - v0) / v0, float(v1.sum() - v0.sum())


def drape(grid, top):
    """Re-sample the horizons at the moved pillars: adapt_grid moves COORD, not ZCORN."""
    nx, ny, _ = grid.cart_dims
    xy = pillar_lattice(grid).top
    cube = corner_cube(grid)
    pi, pj = (np.arange(2 * nx) + 1) // 2, (np.arange(2 * ny) + 1) // 2
    now = np.empty((nx + 1, ny + 1))
    now[pi[:, None], pj[None, :]] = cube[0].T
    shift = np.vectorize(top)(xy[..., 0], xy[..., 1]) - now
    return with_cube(grid, cube + shift[pi][:, pj].T[None]), float(np.abs(shift).max())


def case_a(out):
    print("A: a curved trace - staircase FAULTS vs an adapted pillar chain")
    grid = cake(lambda x, y: 1000.0 + 0.02 * x)
    r, cx, cy = 650.0, 1100.0, 400.0
    th = np.radians(np.linspace(-40.0, 40.0, 41))
    arc = np.c_[cx - r * np.cos(th), cy + r * np.sin(th)]
    true = 2.0 * r * np.arcsin(400.0 / r)                 # arc length between y = 0 and 800 m
    np.savetxt(out / "A_arc_trace.csv", arc, delimiter=",", header="x,y", comments="")
    stair = interior(grid, fault_faces_for(grid, arc))
    label = blocks(grid, stair)
    cut = [f for f in stair if label[face_columns(f)[0]] != label[face_columns(f)[1]]]
    print(f"  staircase: {len(stair)} faces, {len(stair) - len(cut)} dead-end spur dropped "
          "(both its sides are in one block, it cannot open)")
    res = adapt_grid(grid, arc, **MORPH)
    t, lengths = 15.0, {}
    for name, g, faces in (("A_arc_staircase", grid, cut),
                           ("A_arc_trace", res.grid, interior(res.grid, res.faults_after))):
        slid = run(name, g, {"ARC": faces}, {column_of_point(g, (cx, cy)): t}, out)
        total = vertical_checks(g, slid, faces)
        lengths[name] = length = plan_length(g, faces)
        assert np.isclose(total, length * (H - t)), "contact = fault length x overlap"
        print(f"  {name:16s} {len(faces)} faces, {length:6.1f} m long (true arc {true:.1f} m, "
              f"{100 * (length / true - 1):+5.1f} %); east side -{t:g} m; "
              f"contact {total:,.0f} m2")
    assert lengths["A_arc_staircase"] > 1.1 * true, "the staircase overstates the length"
    assert abs(lengths["A_arc_trace"] - true) < 0.02 * true


def case_b(out):
    print("B: fault sticks -> fit_plane -> planar polygon -> pillars tilted into the plane")
    grid = cake(flat)
    dip, strike = np.radians(65.0), np.radians(20.0)
    along = np.array([np.sin(strike), np.cos(strike), 0.0])
    down = np.array([np.cos(strike) / np.tan(dip), -np.sin(strike) / np.tan(dip), 1.0])
    rng = np.random.default_rng(7)
    z = np.linspace(-10.0, 100.0, 6)
    sticks = [np.array([560.0, 400.0, 1000.0]) + s * along + np.outer(z, down)
              + rng.normal(0.0, 0.5, (6, 3)) for s in (-380.0, -190.0, 0.0, 190.0, 380.0)]
    pts = np.vstack(sticks)
    np.savetxt(out / "B_fault_sticks.csv", np.c_[np.repeat(np.arange(5), 6), pts],
               delimiter=",", header="stick,x,y,z", comments="", fmt="%g")
    p0, n, _, _ = fit_plane(pts)
    rms = float(np.sqrt(np.mean(((pts - p0) @ n) ** 2)))
    fit_dip = float(np.degrees(np.arccos(abs(n[2]))))
    assert rms < 1.0 and abs(fit_dip - 65.0) < 1.0, (rms, fit_dip)
    st = np.cross(n, [0.0, 0.0, 1.0])
    st /= np.linalg.norm(st)
    dd = np.cross(st, n)
    dd /= dd[2]                                              # down-dip, per metre of depth

    def on_plane(depth, s):
        return p0 + (depth - p0[2]) * dd + s * st

    polygon = np.array([on_plane(990.0, -2000.0), on_plane(990.0, 2000.0),
                        on_plane(1100.0, 2000.0), on_plane(1100.0, -2000.0)])
    trace = np.array([on_plane(1000.0, -600.0)[:2], on_plane(1000.0, 600.0)[:2]])
    g1, g, faces, chain = two_stage(grid, trace, polygon)
    gap1, gap = (np.abs((np.vstack(pillar_ends(x, chain)) - p0) @ n) for x in (g1, g))
    assert gap1.max() > 30.0 and gap.max() < 0.15 * DX, (gap1.max(), gap.max())
    print(f"  {len(sticks)} sticks, {len(pts)} points: plane rms {rms:.2f} m, dip {fit_dip:.1f} "
          f"(true 65); chain pillar ends to the plane {gap1.max():.0f} m (vertical) -> "
          f"max {gap.max():.1f}, mean {gap.mean():.1f} m (tilted), sj_min {sj_min(g):.2f}")
    seed, t = column_of_point(g, (1100.0, 400.0)), 20.0
    slid = run("B_fault_sticks", g, {"STICKS": faces}, {seed: t}, out)
    rel, dv = hanging_wall_change(g, slid, faces, seed)
    # The block slides down a wedge: it loses t * H * cot(dip) per metre of strike.
    expected = -t * H / np.tan(np.radians(fit_dip)) * (NY * DX / abs(st[1]))
    assert 0.0 < np.abs(rel).max() and abs(dv - expected) < 0.03 * abs(expected), (dv, expected)
    changed = int(np.sum(np.abs(rel) > 1e-9))
    print(f"  hanging wall -{t:g} m along leaning pillars: {changed} cells change "
          f"volume (max {100 * np.abs(rel).max():.1f} %), total {dv:,.0f} m3 vs the wedge "
          f"{expected:,.0f} m3")


def case_c(out):
    print("C: point cloud -> planar facets -> pillars tilted onto a curved fault")
    grid = cake(flat)
    tan70 = np.tan(np.radians(70.0))

    def strike_x(y):
        return 480.0 + 150.0 * ((y - 400.0) / 400.0) ** 2

    rng = np.random.default_rng(3)
    xy = np.c_[rng.uniform(440.0, 700.0, 3000), rng.uniform(-60.0, 860.0, 3000)]
    cloud = np.c_[xy, 1000.0 + (xy[:, 0] - strike_x(xy[:, 1])) * tan70
                  + rng.normal(0.0, 0.2, len(xy))]
    np.savetxt(out / "C_point_cloud.csv", cloud, delimiter=",", header="x,y,z", comments="",
               fmt="%.2f")
    box = dict(nu=4, nv=17, xlim=(440.0, 700.0), ylim=(-25.0, 825.0))
    quads = facet_surface(cloud, **box)
    flatness = max(float(np.abs((q - fit_plane(q)[0]) @ fit_plane(q)[1]).max()) for q in quads)
    assert len(quads) == 4 * 17 and flatness < 1e-6
    # The top trace: where each west-east row of the resampled cloud crosses z = 1000 m.
    rows = sample_cloud_grid(cloud, **box).transpose(1, 0, 2)
    trace = np.array([[np.interp(1000.0, r[:, 2], r[:, 0]), r[0, 1]] for r in rows])
    _, g, faces, chain = two_stage(grid, trace, quads)
    ends = np.vstack(pillar_ends(g, chain))
    grad = np.c_[np.ones(len(ends)), -300.0 * (ends[:, 1] - 400.0) / 400.0 ** 2,
                 np.full(len(ends), -1.0 / tan70)]
    gap = np.abs(ends[:, 0] - strike_x(ends[:, 1]) - (ends[:, 2] - 1000.0) / tan70)
    gap /= np.linalg.norm(grad, axis=1)
    rim = np.isin(np.tile(chain[:, 1], 2), (0, NY))      # rim pillars only slide along the rim
    assert gap[~rim].max() < 0.15 * DX, gap[~rim].max()
    ug = build_topology(from_grdecl(g))
    surf = project_faceted_surface(quads, ug)
    found = {(f.i, f.j, f.k, f.face) for f in surface_to_fault_faces(ug, surf)}
    columns = {(f.i, f.j, k, f.face) for f in faces for k in range(1, len(DZ) + 1)}
    assert columns <= found, "the projected surface must cover every chain column"
    extra = sorted(found - columns)
    assert all(f[1] in (1, NY) for f in extra), "partial faces only at the rim rows"
    (out / "C_point_cloud_faults3d.inc").write_text(faults_text_3d({"CLOUD": surf}, ug))
    print(f"  {len(cloud)} points -> {len(quads)} facets (flat to {flatness:.0e} m); chain pillar "
          f"ends to the true surface: max {gap[~rim].max():.1f}, mean {gap[~rim].mean():.1f} m "
          f"inside, max {gap[rim].max():.1f} m at the rim")
    print(f"  projected {len(found)} cell faces = {len(faces)} columns x {len(DZ)} layers + "
          f"{len(extra)} partial (i, j, k, face) {extra} where a rim pillar misses the surface")
    seed, t = column_of_point(g, (1100.0, 400.0)), 20.0
    slid = run("C_point_cloud", g, {"CLOUD": faces}, {seed: t}, out)
    rel, dv = hanging_wall_change(g, slid, faces, seed)
    print(f"  hanging wall -{t:g} m on the column faces (a column slide cannot open the partial "
          f"rim faces): max volume change {100 * np.abs(rel).max():.1f} %, total {dv:,.0f} m3")


def case_d(out):
    print("D: a structured quad mesh fault -> 3D FAULTS -> slide")
    grid = cake(flat)
    cot = 1.0 / np.tan(np.radians(70.0))
    top = np.array([[420.0, 0.0], [560.0, 400.0], [440.0, 800.0]])
    normals = []
    for a, b in zip(top[:-1], top[1:]):
        s = (b - a) / np.linalg.norm(b - a)
        normals.append(np.array([s[1], -s[0]]))             # both panels dip to the east
    shift = (1090.0 - 1000.0) * cot
    lines = np.array(normals)
    kink = np.linalg.solve(lines, [normals[0] @ top[0] + shift, normals[1] @ top[1] + shift])
    ends = top[[0, 2]] + shift * lines
    # Slide the base ends along their panel's strike back onto the rims y = 0 and y = 800 m.
    south = ends[0] + (0.0 - ends[0][1]) / (kink[1] - ends[0][1]) * (kink - ends[0])
    north = kink + (800.0 - kink[1]) / (ends[1][1] - kink[1]) * (ends[1] - kink)
    base = np.array([south, kink, north])
    verts = np.r_[np.c_[top, np.full(3, 1000.0)], np.c_[base, np.full(3, 1090.0)]]
    mesh = QuadMesh(verts, np.array([[0, 1, 4, 3], [1, 2, 5, 4]]))
    report = validate_quad_mesh(mesh)
    assert (report.n_quads, report.n_interior_edges, report.n_boundary_edges) == (2, 1, 6)
    np.savetxt(out / "D_quad_mesh_vertices.csv", verts, delimiter=",", header="x,y,z",
               comments="", fmt="%.3f")

    def project(g):
        ug = build_topology(from_grdecl(g))
        surf = project_quad_mesh_surface(mesh, ug)
        assert len(surface_boundary_loops(ug, surf.face_ids)) == 1, "one rim: no holes"
        return ug, surf, surface_to_fault_faces(ug, surf)

    _, _, regular = project(grid)
    kinds = Counter(f.face for f in regular)
    stacks = Counter((f.i, f.j, f.face) for f in regular)
    assert kinds["Z"] > 0, "on vertical columns the dipping mesh needs Z faces"
    print(f"  mesh: {report.n_quads} quads, planar to {report.max_planar_deviation:.0e} m; "
          f"regular grid: {len(regular)} faces (X {kinds['X']}, Y {kinds['Y']}, Z {kinds['Z']}), "
          f"{sum(v == len(DZ) for v in stacks.values())} of {len(stacks)} columns full height")
    _, g, faces, _ = two_stage(grid, top, [verts[q] for q in mesh.quads])
    ug, surf, adapted = project(g)
    columns = {(f.i, f.j, k, f.face) for f in faces for k in range(1, len(DZ) + 1)}
    assert {(f.i, f.j, f.k, f.face) for f in adapted} == columns, "whole columns, nothing else"
    (out / "D_quad_mesh_faults3d.inc").write_text(faults_text_3d({"MESH": surf}, ug))
    seed, t = column_of_point(g, (1100.0, 400.0)), 20.0
    slid = run("D_quad_mesh", g, {"MESH": faces}, {seed: t}, out)
    rel, dv = hanging_wall_change(g, slid, faces, seed)
    print(f"  adapted grid: {len(adapted)} faces = {len(faces)} columns x {len(DZ)} layers, "
          f"sj_min {sj_min(g):.2f}; hanging wall -{t:g} m: max volume change "
          f"{100 * np.abs(rel).max():.1f} %, total {dv:,.0f} m3")


def case_e(out):
    print("E: a dome with a crestal graben")

    def dome(x, y):
        return 1000.0 + 1e-4 * ((x - 600.0) ** 2 + (y - 400.0) ** 2)

    west = np.array([[470.0, -10.0], [500.0, 400.0], [470.0, 810.0]])
    east = np.array([[700.0, -10.0], [710.0, 400.0], [750.0, 810.0]])
    res = adapt_grid(cake(dome), [west, east], **MORPH)
    g, fix = drape(res.grid, dome)
    faults = named_chains(g, interior(g, res.faults_after))
    faces = faces_of(faults)
    seed, t = column_of_point(g, (600.0, 400.0)), 20.0
    slid = run("E_dome_graben", g, faults, {seed: t}, out)
    total = vertical_checks(g, slid, faces)
    label = blocks(g, faces)
    graben = (label == label[seed])[np.arange(2 * NX) // 2][:, np.arange(2 * NY) // 2].T
    top0, top1 = corner_cube(g)[0], corner_cube(slid)[0]
    assert np.isclose(top1.min(), min(top0[~graben].min(), top0[graben].min() + t))
    assert top1.min() > top0.min() and not graben.flat[np.argmin(top1)], "crest on a shoulder"
    print(f"  {len(faults)} faults, {len(faces)} faces, horizons re-draped (up to {fix:.2f} m); "
          f"graben -{t:g} m: crest {top0.min():.1f} -> {top1.min():.1f} m (now on a shoulder), "
          f"graben top {top0[graben].min() + t:.1f} m; contact {total:,.0f} m2")


def case_f(out):
    print("F: a partial normal fault in dipping layers (7.5 m per column)")

    def tilt(x, y):
        return 1000.0 + 0.15 * x

    trace = np.array([[330.0, 180.0], [600.0, 420.0], [880.0, 610.0]])
    res = adapt_grid(cake(tilt), trace, **MORPH)
    g, fix = drape(res.grid, tilt)
    faces = interior(g, res.faults_after)
    hanging, peak = (800.0, 300.0), 25.0
    raw = slide_field(g, Throw(faces=tuple(faces), throw=peak, hanging_side=hanging, snap=False),
                      faults=faces)
    raw_total = vertical_checks(g, raw, faces)   # horizons cross inside faces: still exact
    crossed = crossed_faces(process(raw).topology)
    assert crossed > 0, "an unsnapped field should leave a horizon inside a fault face"
    spec = Throw(faces=tuple(faces), throw=peak, hanging_side=hanging)
    slid = run("F_dipping_layers", g, {"F1": faces}, spec, out)
    total = vertical_checks(g, slid, faces)
    throw = realized_throw(g, slid, faces, hanging)
    unsnapped = realized_throw(g, raw, faces, hanging)[:, 1]
    assert np.all(throw[:, 1] <= unsnapped + 1e-9), "the snap only lowers the throw"
    snapped = int(np.sum(~np.isclose(throw[:, 1], unsnapped)))
    assert throw[0, 1] == 0.0 and throw[-1, 1] == 0.0 and 0.9 * peak < throw[:, 1].max() <= peak
    np.savetxt(out / "F_realized_throw.csv", throw, delimiter=",", header="s,throw",
               comments="", fmt="%.3f")
    print(f"  {len(faces)} faces, {plan_length(g, faces):.0f} m, horizons re-draped (up to "
          f"{fix:.2f} m); unsnapped: {crossed} faces cut at a horizon crossing, contact "
          f"{raw_total:,.0f} m2; snapped: {snapped} pillars lowered")
    print(f"  realized throw (m): {' '.join(f'{v:.0f}' for v in throw[:, 1])}; "
          f"contact {total:,.0f} m2")


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    for case in (case_a, case_b, case_c, case_d, case_e, case_f):
        case(out)
    print("every case: COORD unchanged, exactly the FAULTS faces open; on vertical pillars (A, E, "
          "F) volumes and analytic contacts hold, on tilted ones (B-D) only the hanging wall "
          "changes volume")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
