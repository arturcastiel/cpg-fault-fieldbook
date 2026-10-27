"""Gallery 4: ways to define a throw — tables, growth, large throws, restoration and inversion.

A throw can be one number, a measured profile along the fault, one value per horizon, or a sum
of slides of either sign. Each case starts from a layer cake (24x16x4, 50 m columns, layers
25/30/25/35 m, top dipping 1/64 to the east), adapts it with adapt_grid to an oblique fault
that ends inside the grid or to a curved fault that crosses it, and slides along the pillars
by editing ZCORN only:

- A throw table: a measured throw-distance profile (0, 18, 10, 22, 0 m) on the partial fault;
  the realized throw at every fault pillar is the table (well data, seismic picks);
- B growth fault: the hanging wall of the crossing fault slides 5 m at the top and 27 m at the
  base, so its layers thicken, as in a syn-sedimentary fault; volume grows by exactly the
  extra thickness times the plan area;
- C large throw: a 60 m peak throw on 25-35 m layers; without snapping, re-processing rejects
  the grid, and snapping puts every horizon crossing on a pillar (throw beyond the layering);
- D restoration: the opposite throw restores ZCORN to round-off (structural restoration, undo);
  this holds while the snap is inactive, as the snap is not linear in the throw;
- E inversion: +20 m then -8 m equals one +12 m slide, for blocks and for fields (tectonic
  inversion of a normal fault);
- F kernel vs laplacian: the same partial fault faded by a distance kernel or a harmonic solve;
  only the far field differs, the throw on the fault is identical.

Each case asserts that COORD stays put, exactly the FAULTS faces open, and the re-processed
contacts match the analytic overlaps. D, E and F also write a third deck (<case>_2_*.grdecl:
restored, field, laplacian), A and C a CSV of the throw at the fault pillars. Pass --out DIR to
keep the decks (<case>_0_input.grdecl, <case>_1_slid.grdecl).
"""

import argparse
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from cpg2unstructured import NonMatchingZError, pillar_lattice
from fault_projections import adapt_grid
from fault_projections.from_faults import trace_faults

from fieldbook.zslide import (Throw, analytic_contacts, blocks, column_of_point, compare,
                              corner_cube, interior, layer_cake, process, realized_throw,
                              slide_blocks, slide_field, unwelded_faces, write_deck)

NX, NY, DX = 24, 16, 50.0
DZ = [25.0, 30.0, 25.0, 35.0]
MORPH = dict(boundary="pinned", w0=0.05, polish="barrier", on_invalid="raise")
PARTIAL = np.array([[330.0, 170.0], [560.0, 400.0], [880.0, 610.0]])     # both tips inside
CROSSING = np.array([[540.0, 0.0], [640.0, 400.0], [580.0, 800.0]])      # rim to rim
HW_PARTIAL, HW_CROSSING = (750.0, 250.0), (1000.0, 400.0)              # on the hanging walls


def base_grid():
    # A dyadic dip (1/64) keeps every depth exact in binary, so sums of throws are exact too.
    return layer_cake(NX, NY, len(DZ), DX, DZ, top=lambda x, y: 1000.0 + x / 64.0)


def adapted(grid, trace, label):
    res = adapt_grid(grid, trace, **MORPH)
    faces = interior(res.grid, res.faults_after)
    chain = trace_faults(res.grid, faces).chains[0]
    shift = np.linalg.norm(pillar_lattice(res.grid).top - pillar_lattice(grid).top, axis=-1)
    print(f"{label}: {len(faces)} FAULTS faces, {chain.length:.0f} m, pillar "
          f"{tuple(map(int, chain.ij[0]))} to {tuple(map(int, chain.ij[-1]))}; pillars moved "
          f"<= {shift.max():.1f} m")
    return res.grid, faces, chain


def column_areas(grid):
    """Plan area of every column [i, j] (half the cross product of its diagonals)."""
    top = pillar_lattice(grid).top
    d1, d2 = top[1:, 1:] - top[:-1, :-1], top[:-1, 1:] - top[1:, :-1]
    return 0.5 * np.abs(d1[..., 0] * d2[..., 1] - d1[..., 1] * d2[..., 0])


def checked(name, grid, slid, faces, out, *, tag="1_slid", volumes=True):
    """Write the decks; assert COORD, the open faces, volumes and contacts; return volumes."""
    if tag == "1_slid":
        write_deck(out / f"{name}_0_input.grdecl", grid, {"F1": faces}, header=f"{name}: input")
    write_deck(out / f"{name}_{tag}.grdecl", slid, {"F1": faces}, header=f"{name}: {tag}")
    assert np.array_equal(slid.coord, grid.coord), "COORD changed"
    assert unwelded_faces(slid) == set(faces), "exactly the FAULTS faces must open"
    v0 = process(grid).volume_grid(grid.cart_dims)
    after = process(slid)
    v1 = after.volume_grid(grid.cart_dims)
    if volumes:
        assert np.allclose(v1, v0, rtol=1e-12, atol=0.0), "a cell changed volume"
    compare(after.contacts(faces), analytic_contacts(slid, faces))
    return v0, v1


def top_shift(grid, slid):
    """Depth shift of the top of layer 0 at every corner, as [j-slot, i-slot]."""
    return corner_cube(slid)[0] - corner_cube(grid)[0]


def case_a(grid, faces, chain, out):
    s = chain.length * np.array([0.0, 0.2, 0.45, 0.75, 1.0])
    table = tuple(zip(s, (0.0, 18.0, 10.0, 22.0, 0.0)))
    spec = Throw(tuple(faces), throw=table, hanging_side=HW_PARTIAL)
    slid = slide_field(grid, spec)
    # Every value is below the thinnest layer (25 m), so the snap has nothing to do.
    assert np.array_equal(slid.zcorn, slide_field(grid, replace(spec, snap=False)).zcorn)
    checked("A_throw_table", grid, slid, faces, out)
    got = realized_throw(grid, slid, faces, HW_PARTIAL)
    want = np.interp(got[:, 0], s, [t for _, t in table])
    assert np.allclose(got[:, 1], want, rtol=0.0, atol=1e-9), np.abs(got[:, 1] - want).max()
    rows = "\n".join(f"{a:.2f},{b:.6f},{c:.6f}" for a, b, c in zip(got[:, 0], want, got[:, 1]))
    (out / "A_throw_table.csv").write_text("s_m,table_m,realized_m\n" + rows + "\n")
    print("A_throw_table      measured s -> throw (m): "
          + ", ".join(f"{a:.0f} -> {t:g}" for a, t in table))
    print(f"  realized at {len(got)} fault pillars = the table within "
          f"{np.abs(got[:, 1] - want).max():.1e} m; snap inactive (all throws < 25 m)")


def case_b(grid, faces, out):
    growth = np.array([5.0, 9.0, 14.0, 20.0, 27.0])     # per horizon, top of layer 0 to base
    seed = column_of_point(grid, HW_CROSSING)
    hw = blocks(grid, faces) == blocks(grid, faces)[seed]
    slid = slide_blocks(grid, faces, {seed: growth})
    v0, v1 = checked("B_growth_fault", grid, slid, faces, out, volumes=False)
    area = column_areas(grid)[hw]
    thicker = np.diff(growth)
    for k, dt in enumerate(thicker):
        dv = v1[..., k][hw] - v0[..., k][hw]
        assert np.allclose(dv, dt * area, rtol=1e-12, atol=0.0), f"layer {k}"
    assert np.allclose(v1[~hw], v0[~hw], rtol=1e-12, atol=0.0), "the footwall changed"
    gain = [float((v1[..., k] - v0[..., k])[hw].sum()) for k in range(len(DZ))]
    print(f"B_growth_fault     hanging wall ({hw.sum()} columns, {area.sum() / 1e6:.3f} km2) "
          f"down {'/'.join(f'{t:g}' for t in growth)} m, top to base")
    print(f"  layers {'/'.join(f'+{t:g}' for t in thicker)} m thicker "
          f"({'/'.join(f'{t / dz:.0%}' for t, dz in zip(thicker, DZ))}), volume "
          f"{'/'.join(f'+{g / 1e6:.3f}' for g in gain)} Mm3")
    print("  = extra thickness x plan area in every hanging-wall cell; footwall volumes unchanged")


def case_c(grid, faces, out):
    spec = Throw(tuple(faces), throw=60.0, hanging_side=HW_PARTIAL, snap=False)
    raw = slide_field(grid, spec)
    try:
        process(raw)
        rejected = False
    except NonMatchingZError:
        rejected = True
    assert rejected, "a 60 m throw on 25-35 m layers must cross a horizon inside a face"
    snapped = slide_field(grid, replace(spec, snap=True))
    checked("C_large_throw", grid, snapped, faces, out)
    t_raw = realized_throw(grid, raw, faces, HW_PARTIAL)[:, 1]
    t_snap = realized_throw(grid, snapped, faces, HW_PARTIAL)[:, 1]
    offsets = np.unique(np.round(np.subtract.outer(np.cumsum([0.0] + DZ),
                                                   np.cumsum([0.0] + DZ)), 9))
    moved = ~np.isclose(t_snap, t_raw)
    assert moved.any() and np.all(np.abs(t_snap) <= np.abs(t_raw) + 1e-9)
    assert all(np.isclose(offsets, t).any() for t in t_snap[moved]), "snap onto a boundary"
    rows = "\n".join(f"{m},{a:.6f},{b:.6f}" for m, (a, b) in enumerate(zip(t_raw, t_snap)))
    (out / "C_large_throw.csv").write_text("pillar,unsnapped_m,snapped_m\n" + rows + "\n")
    print("C_large_throw      peak 60 m on 25-35 m layers: snap=False rejected, snap=True "
          "re-processes")
    print("  unsnapped: " + " ".join(f"{t:.0f}" for t in t_raw) + " m")
    print("  snapped:   " + " ".join(f"{t:.0f}" for t in t_snap)
          + f" m; {moved.sum()} pillars moved, max {t_snap.max():.1f} m")


def case_d(grid, faces, out):
    spec = Throw(tuple(faces), throw=20.0, hanging_side=HW_PARTIAL, decay=300.0)
    slid = slide_field(grid, spec)
    checked("D_restoration", grid, slid, faces, out)
    restored = slide_field(slid, replace(spec, throw=-spec.throw))
    write_deck(out / "D_restoration_2_restored.grdecl", restored, {"F1": faces},
               header="D_restoration: slid back")
    err = float(np.abs(restored.zcorn - grid.zcorn).max())
    assert err <= 1e-9 and not unwelded_faces(restored, tol=1e-9), err
    print(f"D_restoration      +20 m then -20 m: ZCORN back within {err:.1e} m on "
          f"{np.count_nonzero(slid.zcorn != grid.zcorn)} moved corners; all faces welded")


def case_e(grid, faces, out):
    seed = column_of_point(grid, HW_CROSSING)
    twice = slide_blocks(slide_blocks(grid, faces, {seed: 20.0}), faces, {seed: -8.0})
    once = slide_blocks(grid, faces, {seed: 12.0})
    assert np.array_equal(twice.zcorn, once.zcorn), "block slides must add exactly"
    checked("E_inversion", grid, twice, faces, out)
    spec = Throw(tuple(faces), throw=20.0, hanging_side=HW_CROSSING, decay=400.0)
    field2 = slide_field(slide_field(grid, spec), replace(spec, throw=-8.0))
    field1 = slide_field(grid, replace(spec, throw=12.0))
    checked("E_inversion", grid, field2, faces, out, tag="2_field")
    err = float(np.abs(field2.zcorn - field1.zcorn).max())
    assert err <= 1e-9, err
    print(f"E_inversion        +20 m then -8 m == +12 m: blocks bit for bit, field (decay "
          f"400 m) to {err:.1e} m")


def case_f(grid, faces, out):
    spec = Throw(tuple(faces), throw=18.0, hanging_side=HW_PARTIAL, decay=350.0)
    kernel = slide_field(grid, spec, method="kernel")
    harmonic = slide_field(grid, spec, method="laplacian")
    checked("F_kernel_vs_laplacian", grid, kernel, faces, out)
    checked("F_kernel_vs_laplacian", grid, harmonic, faces, out, tag="2_laplacian")
    area = column_areas(grid)
    lines = []
    for name, slid in (("kernel", kernel), ("laplacian", harmonic)):
        shift = np.abs(top_shift(grid, slid))
        per_col = shift.reshape(NY, 2, NX, 2).mean((1, 3)).T          # [i, j]
        lines.append(f"{name}: max {shift.max():.1f} m, |shift| > 1 m on "
                     f"{area[per_col > 1.0].sum() / 1e6:.3f} km2")
    t_k = realized_throw(grid, kernel, faces, HW_PARTIAL)
    t_l = realized_throw(grid, harmonic, faces, HW_PARTIAL)
    assert np.allclose(t_k, t_l, rtol=0.0, atol=1e-12), "the throw on the fault must agree"
    diff = np.abs(kernel.zcorn - harmonic.zcorn).max()
    print(f"F_kernel_vs_laplacian  {lines[0]}")
    print(f"  {lines[1]}; shifts differ by up to {diff:.1f} m off the fault")
    print(f"  the throw on the fault is identical for both (peak {t_k[:, 1].max():.1f} m)")


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    grid = base_grid()
    g_part, f_part, chain = adapted(grid, PARTIAL, "partial oblique fault")
    g_cross, f_cross, _ = adapted(grid, CROSSING, "crossing curved fault")
    case_a(g_part, f_part, chain, out)
    case_b(g_cross, f_cross, out)
    case_c(g_part, f_part, out)
    case_d(g_part, f_part, out)
    case_e(g_cross, f_cross, out)
    case_f(g_part, f_part, out)
    print("every case: COORD unchanged, exactly the FAULTS faces open, contacts match the "
          "analytic overlaps")


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--out", type=Path, help="keep the decks in this directory")
args = parser.parse_args()
if args.out is not None:
    main(args.out)
else:
    with TemporaryDirectory() as temporary:
        main(Path(temporary))
