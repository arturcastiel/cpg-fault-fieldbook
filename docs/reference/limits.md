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
