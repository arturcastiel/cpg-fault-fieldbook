# Function reference

Source-backed reference pages for the grid and fault-projection APIs used by
this fieldbook. The package import names are `cpg2unstructured` (from the
`cpg2python` sibling checkout) and `fault_projections` (from
`fault-projections`). Examples use explicit submodule imports to make ownership
clear; they have not been run as part of writing these pages.

| Page | Functions | Main use |
|---|---|---|
| [Grid input](grid-input.md) | `read_grdecl`, `read_egrid` | Load Eclipse text or binary corner-point data |
| [Topology and quality](topology-quality.md) | `from_grdecl`, `build_topology`, `check_validity`, `cell_quality` | Build a corner cube/mesh and assess geometry |
| [Pillar morphing](morph-pillars.md) | `morph_pillars` | Move pillars toward prescribed targets |
| [Targets and corridors](targets-corridors.md) | `targets_from_polyline`, `targets_from_polylines`, `cell_corridor` | Turn traces into pillar targets and post-morph checks |
| [Grid adaptation](adapt-grid.md) | `adapt_grid`, `adapt_deck` | Run an end-to-end fault-driven adaptation |
| [3-D surface projection](surface-projection.md) | `project_fault_surface`, `project_quad_mesh_surface` | Map planar or quad-mesh faults onto cell faces |
| [Fault slide](fault-slide.md) | `slide_grid`, `slide_deck`, `zcorn_ops` | Give existing FAULTS a throw by editing ZCORN only |

## Reading these references

- Shapes are NumPy shapes; `(NX, NY, NZ)` follows the logical grid dimensions.
  Pillar indices are zero-based unless explicitly described as Eclipse output.
- “Source-backed example” means the call pattern follows implementation,
  examples, or tests in the linked checkout. It does **not** mean this page
  executed the snippet.
- Maturity describes the documented implementation and evidence, not a
  production guarantee. Read each function’s limits before exporting a deck.
- `Grdecl` is an alias of `CornerPointGrid`; EGRID input does not create a
  different legacy type.

