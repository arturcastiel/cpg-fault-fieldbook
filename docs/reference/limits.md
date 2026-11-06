# Scope and known failures

- **Sheared rim:** `boundary="slide"` uses index-space tangents,
  so physical outline motion can be nonzero despite `result.ok`.
  Reproduce with `10_sheared_rim`; prefer `pinned`.
- **Field topology:** since cpg2unstructured 0.4.0, `build_topology`
  processes the documented Norne deck and every cell closes. Beside
  skewed fault pillars, cell volumes deviate from the trilinear ones
  (99th percentile 0.09 %, worst 4.4 %). `regenerate_faults=False`
  adapts the grid only and **does not output regenerated FAULTS**.
- **3D input:** dipping polygons/facets produce independent top/bottom
  targets, but `adapt_grid` does not regenerate FAULTS for them.
- **Sliding a 3D projection:** the staircase that `project_fault_surface`
  gives a dipping fault cannot be slid. A slide moves corners along their
  pillars and takes column faces over one K range; `slide_deck` drops the
  Z steps, sees parallel lines and raises `SlideError`. Move the pillars
  onto the fault first (`adapt_grid(top_trace="auto")` or `dip_grid`, then
  slide: workflows 5 and 9). A vertical fault projected in 3D slides (its K range
  becomes `Throw.layers`). A horizontal fault (Z faces only) can be neither
  slid nor dipped: pillars cannot be offset across it.
- **Dip by `dip_grid` from a 3D surface (workflow 9, route B):** the faces
  come from the plane's top trace, snapped to the nearest pillar line, so
  the tilted line can lie up to half a cell from the plane (15 m measured
  for a plane meeting the top 15 m off a line). `adapt_grid` moves the
  pillars onto the plane itself.
- **Near-flat cells:** validity alone allows a low scaled Jacobian;
  add a project-specific floor and inspect both planes.
- **Deck semantics:** generating an INCLUDE does not rewrite or
  validate a simulator deck, transmissibility or flow.
- **Z slide (lessons 19–22):** these lessons edit ZCORN directly to show
  the mechanics; the library call is `fault_projections.slide_grid`
  (lessons 23–30).
  - A rigid slide needs FAULTS that close a block. A fault ending
    inside the grid needs a throw that tapers to zero at its tip.
    Lesson 22 shapes that taper with a profile, a decay distance and
    a footwall share; the decay is a distance kernel, not a solved
    elastic field.
  - If the throw passes a multiple of the layer thickness inside one
    cell face, re-processing gives that face extra nodes where the
    layer edges cross (lesson 20). Tools that need 4-corner faces
    should set `snap=True`.
  - On tilted fault pillars, cells bounded by a moved pillar change
    volume, while their thickness is preserved, unless every pillar of
    the moving block is parallel to the fault (lesson 21): then the
    block translates rigidly.
  - A dipping fault that ends on the grid rim cannot be adapted with a
    stage-2 polygon under `boundary="pinned"`, because the rim pillars
    stay vertical. Lesson 21 leans the hanging-wall pillars first,
    then adapts with a 2D trace (stage 1 keeps the lean). With
    `boundary="slide"` the rim pillars can tilt: lesson 30 C adapts a
    rim-to-rim 65° polygon with `top_trace="auto"`. To set the angle of
    a fault that is already in `FAULTS`, use `dip_grid` (lesson 31): it
    leans the fault's pillars without adapting, those on the rim too
    (along the rim, with the default fixed outline).
  - Slides across K faces are not covered: k-direction faults still
    raise `NonMatchingZError`.
- **Slide gallery (lessons 23–29):** the helper `fieldbook/zslide.py`
  calls `fault_projections.slide_grid`, with `snap=True` by default.
  - The snap lowers a throw to the nearest offset at which layer
    boundaries meet on a pillar, so realized throws differ slightly
    from a requested table (lesson 27 B: 20.9–49.1 m for 20–50 m).
    Restoration and inversion are exact only while the snap is
    inactive.
  - The snap only covers the fault whose `Throw` is being applied. A
    splay field also changes the main fault's throw near the junction,
    so layer boundaries may cross inside faces there (lesson 24).
  - Slides along tilted pillars (stage-2 polygons, facets, meshes)
    change cell volumes by the fault wedge (lesson 25 B–D).
- **Fixed in cpg2unstructured 0.4.0:** inactive cells and pinch-outs
  beside a fault get boundary faces instead of raising (lesson 28 E),
  and six `__all__` names are now importable from the package root.
- **Fixed upstream in fault-projections 0.2.0** (found by the gallery):
  rim and spur faces in `adapt_grid(...).faults_after`, rings left
  open, chains merged across traces (now `faults_by_trace`), traces
  with different vertex counts, facet tapers at inner facet tips,
  `BoundaryDoesNotBound` on stage-1 grids, fault rim pillars left behind
  under `boundary="pinned"` (now `slide_rim_ends=True`), rim faces from
  steep traces, `%.6g` VTK points and the missing `.vtu`/`.pvd`
  writers, INCLUDE in `faults_from_deck`, and `adapt_deck` merging
  every fault into one group. The lessons no longer call the old
  guards `interior` and `close_ring`; `named_chains` and `taper_cells=0`
  remain where a lesson needs them.
- **Fixed upstream in fault-projections 0.4.0:** a stage-2 polygon or
  facets alone give no connected chain. `adapt_grid(..., top_trace=...)`
  now runs stage 1 on the top trace and then stage 2 (`"auto"` derives the
  trace), and lesson 25 uses it. `faces_invariant` no longer reports a change
  for a trace through inactive columns (an L-shaped active region): only
  faces with active cells on both sides are compared. `slide_deck(png=True)`
  writes a figure.
- **Fault dip (lessons 31 and 32):** `fault_projections.dip_grid` changes
  COORD only.
  - The default `boundary="fixed"` holds the grid outline. A hanging wall
    that reaches the rim is then not rigid next to it (case A: cells there
    change up to 21.7 %), and a shallow fault near a parallel side can be
    refused (lesson 32 O: 15° refused, 30° passes). `boundary="free"` keeps
    the hanging wall rigid but leans the rim and any blocks beyond other
    faults; a `decay` distance fades the tilt out instead.
  - Cell volumes are not conserved one by one. A rigid hanging wall away from
    the rim keeps them to rounding (6.2e-15 in case I); footwall and
    fade-zone cells change
    (case A: up to 21.65 %; the tipped fault of case D: up to 29.4 %). The
    middle pivot keeps column totals only.
  - Bisection (`on_invalid="bisect"`) keeps the grid valid but misses
    the angle (case H: 45° for 35° requested).
  - The gate is conservative where bedding is steeper than the fault: the
    plan view of a cell's top flips while the 3-D cell is sound, and the
    dip is refused (lesson 32 R). In lesson 32's fuzz, 4 of 105 verdicts
    were such refusals; none was an accepted grid with broken cells.
  - **Field result (lesson 31 J, Norne):** E_01 → 60° passes. BC → 70°
    fails the gate with the default middle pivot, even with the fades
    widened ×16: its present dip jumps from 59° to 84° within four
    pillars, and top-layer cells fold. `pivot="top"` passes at 65° and
    70°; 65° with the middle pivot now fails on three already-poor cells
    (Jacobian 0.021 → 0.0207, beyond the gate's 1 % allowance). A smarter
    spread for faults with sharply varying dip is not implemented.
- **Fault zones (lesson 33):** `fault_zone` selects cells by their centre's
  distance from the fault surface; it never changes the grid. `zone=`
  confines a dip or slide, but the whole change is then taken up inside the
  zone, so its cells deform more, and a zone too narrow for the change is
  refused (`ZoneError`, or the gate). The zone is measured on the input
  grid.
- **Fixed upstream after fault-projections 0.4.0** (pinned revision):
  - `adapt_grid` and `adapt_deck` default to `polish="auto"`: when the
    morph leaves an active column with a scaled Jacobian below 0.05 and
    below half its old value, the barrier polish runs. On upstream's
    65° stage-2 test polygon (24 × 16 × 4 cells, `top_trace="auto"`,
    `boundary="slide"`), the minimum scaled Jacobian is 3.7e-4 with
    `polish="none"` and 0.277 with the default. `polish="none"` restores
    the old results.
  - `surface_top_trace` follows the longest path of the points' minimum
    spanning tree, so U-shaped and hooked faults no longer fold.
  - `Throw.layers=(k1, k2)` slides a blind fault; `slide_deck` reads it
    from the `FAULTS` K range (earlier versions slid every layer).
- **Fixed upstream in cpg2unstructured after 0.4.0** (09affc9, found by
  workflow 9): `morph_pillars` dropped a pillar's whole move when its top
  target was within `write_tol` of its current top, so the base move was
  lost too. `adapt_grid(polygon, top_trace="auto")` on a dipping plane
  whose top trace lies exactly on a pillar line left those pillars vertical
  and still returned `ok` (90° instead of 45°). Each end is now zeroed on
  its own; the plane reaches 45.1°, and off-line cases are unchanged.
