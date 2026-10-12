# Scope and known failures

- **Sheared rim:** `boundary="slide"` uses index-space tangents,
  so physical outline motion can be nonzero despite `result.ok`.
  Reproduce with `09_sheared_rim`; prefer `pinned`.
- **Field topology:** full `build_topology` currently fails on the
  documented Norne deck (face with more than two neighbours);
  `regenerate_faults=True` is unsafe to assume there. `False`
  adapts the grid only and **does not output regenerated FAULTS**.
- **3D input:** dipping polygons/facets produce independent top/bottom
  targets, but `adapt_grid` does not regenerate FAULTS for them.
- **Near-flat cells:** validity alone allows a low scaled Jacobian;
  add a project-specific floor and inspect both planes.
- **Deck semantics:** generating an INCLUDE does not rewrite or
  validate a simulator deck, transmissibility or flow.
- **Packaging caveat:** six `cpg2unstructured.__all__` names are
  unbound at package root in the pinned revision: `depth_crossings`,
  `corner_folds`, `proximity_violations`, `elastic_stiffness`,
  `stiffness_from_quality`, `support_mask`. Use their source
  submodules rather than `from cpg2unstructured import ...`.
