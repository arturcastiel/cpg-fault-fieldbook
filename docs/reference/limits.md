# Scope and known failures

- **Sheared rim:** `boundary="slide"` uses index-space tangents,
  so physical outline motion can be nonzero despite `result.ok`.
  Reproduce with `10_sheared_rim`; prefer `pinned`.
- **Field topology:** full `build_topology` currently fails on the
  documented Norne deck with `NonMatchingZError` (a fault polygon at
  i=4–5, j=2 needs more than four corners);
  `regenerate_faults=True` is unsafe to assume there. `False`
  adapts the grid only and **does not output regenerated FAULTS**.
- **3D input:** dipping polygons/facets produce independent top/bottom
  targets, but `adapt_grid` does not regenerate FAULTS for them.
- **Near-flat cells:** validity alone allows a low scaled Jacobian;
  add a project-specific floor and inspect both planes.
- **Deck semantics:** generating an INCLUDE does not rewrite or
  validate a simulator deck, transmissibility or flow.
- **Z slide (lessons 19–22):** this is not yet a library call; the
  lessons edit ZCORN directly.
  - A rigid slide needs FAULTS that close a block. A fault ending
    inside the grid needs a throw that tapers to zero at its tip.
    Lesson 22 shapes that taper with a profile, a decay distance and
    a footwall share; the decay is a distance kernel, not a solved
    elastic field.
  - If the throw passes a multiple of the layer thickness inside one
    cell face, re-processing raises `NonMatchingZError`.
  - On tilted fault pillars, cells bounded by a moved pillar change
    volume, while their thickness is preserved, unless every pillar of
    the moving block is parallel to the fault (lesson 21): then the
    block translates rigidly.
  - A dipping fault that ends on the grid rim cannot be adapted with a
    stage-2 polygon under `boundary="pinned"`, because the rim pillars
    stay vertical. Lesson 21 leans the hanging-wall pillars first,
    then adapts with a 2D trace (stage 1 keeps the lean).
  - Inactive cells and slides across K faces are not covered.
- **Slide gallery (lessons 23–29):** the helper `fieldbook/zslide.py`
  is a prototype, not a library call.
  - The snap lowers a throw to the nearest offset at which layer
    boundaries meet on a pillar, so realized throws differ slightly
    from a requested table (lesson 27 B: 20.9–49.1 m for 20–50 m).
    Restoration and inversion are exact only while the snap is
    inactive.
  - The snap only covers the fault whose `Throw` is being applied. A
    splay field also changes the main fault's throw near the junction;
    keep such throws below the layer spacing (lesson 24).
  - A cell beside a slid fault that is inactive, or a pinch-out there,
    makes re-processing raise `NonMatchingZError` (lesson 28 E).
  - Slides along tilted pillars (stage-2 polygons, facets, meshes)
    change cell volumes by the fault wedge (lesson 25 B–D).
- **Upstream gaps found by the gallery:**
  - `adapt_grid(...).faults_after` can include rim faces (`I=0` or
    `J=0`), omits the closing face of a closed trace, and merges the
    chains of all traces into one list. The helper adds `interior`,
    `close_ring` and `named_chains`.
  - A path edge on the south or west rim can emit a spur face on the
    far side; arc ends can leave a dead-end spur. Drop faces whose two
    columns lie in the same block.
  - `adapt_grid` with several traces needs the same vertex count per
    trace (they are stacked into one array).
  - `targets_from_facets` and `targets_from_polygon` taper at every
    facet tip; use `taper_cells=0` for multi-facet faults. Stage 2
    alone gives no connected chain: run stage 1 with the top trace
    first. `project_quad_mesh_surface` raises `BoundaryDoesNotBound`
    on a grid adapted by stage 1 only.
  - `boundary="pinned"` leaves fault rim pillars where they were
    (25–50 m off the trace). Without the barrier polish the
    minimum scaled Jacobian can fall to about 2e-4.
  - `fault_projections.vtk` writes points with `%.6g`, which rounds
    UTM coordinates by up to centimetres; it has no `.vtu` or `.pvd`
    writer. `faults_from_deck` does not follow INCLUDE, and
    `adapt_deck` merges all FAULTS into one group.
- **Steep traces with `adapt_grid`:** on lesson 21's 8×6 grid, a
  trace from the south rim to the north rim at about 30° or more off
  J makes the regenerated FAULTS include a rim face (`j=0`), so
  `faces_invariant` is False. Use `fault_faces_for` on the unadapted
  grid, or a gentler trace.
- **Packaging caveat:** six `cpg2unstructured.__all__` names are
  unbound at package root in the pinned revision: `depth_crossings`,
  `corner_folds`, `proximity_violations`, `elastic_stiffness`,
  `stiffness_from_quality`, `support_mask`. Use their source
  submodules rather than `from cpg2unstructured import ...`.
