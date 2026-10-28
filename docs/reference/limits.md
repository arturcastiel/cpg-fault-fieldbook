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
    then adapts with a 2D trace (stage 1 keeps the lean).
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
- **Still open upstream:**
  - Stage 2 alone gives no connected chain: run stage 1 with the top
    trace first.
  - Without the barrier polish, the minimum scaled Jacobian can fall
    to about 2e-4.
  - A trace that runs into inactive columns of an irregular active
    region can still give `faces_invariant` False.
