# API coverage: measured, not assumed

The [function reference](../functions/index.md) targets the main
learning surface first. The legacy [full API directory](../legacy-field-guide.html#api)
names the complete declared export list. **A mention in that directory
is not a MATLAB-style reference page or an executed example.**
The [generated public-symbol inventory](api-inventory.md) records each
name in both packages' `__all__`, whether the root attribute exists and
whether it has a detailed function page. Its test detects upstream
export drift.

| Area | Function-page coverage | Runnable examples | Status |
|---|---|---|---|
| Text grid, topology, graph | Readers, pivot, topology | 01, 02 | Bundled fixture |
| 2D projection | Core graph/path-to-FAULTS | 03 | Bundled fixture |
| 3D projection | Planar polygon; coherent quad reference | 04, 12 | Quad example needs sibling checkout |
| Targets and deformation | Single, crossing, dipping, faceted adaptation | 05–09, 15 | Synthetic only |
| Validation | Quality, validity, corridor | 05, 08–10 | Project floors explicit |
| Export / read-only deck intake | `adapt_deck`, readers | 11, 13 | Grid INCLUDE only; no simulator check |
| EGRID/MAPAXES | Reader reference | 02_egrid, 17_mapaxes_guard | Bundled upstream reader fixture; no simulator round trip |
| ZCORN-only Z slide | `cpg2unstructured.zcorn_ops`, `fault_projections.slide_grid`; lessons 19–22 edit the corner cube directly | 19_z_slide_sequential, 20_z_slide_partial_fault, 21_z_slide_oblique_fault, 22_z_slide_throw_control | Synthetic only |
| Slide gallery | `fieldbook/zslide.py` over `fault_projections.slide_grid` | 23–29 gallery lessons, 38 cases | Synthetic only; exercises `adapt_grid`, `adapt_deck`, stage-2 polygons, faceted and quad-mesh projection, `cell_corridor`, VTK writers |
| Shapefile and VTK plotting | Legacy description | None yet | Not verified end-to-end here |
| Other root exports/types/errors | Legacy directory | Upstream tests only | Reference expansion pending |

**Definition of done for a full reference:** each public symbol has a
syntax/shape/output/error page and at least one checked minimal
example or a clear reason it is not executable. This initial
release does not meet that condition; do not call it complete.
