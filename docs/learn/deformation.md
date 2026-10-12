# 3. Deform a grid to an authored fault

**First exercise:** inspect claimed pillar targets and solve directly.
This exposes the exact targets needed later for a chain check. Run
`python -m fieldbook.cli run 05_one_trace`.

{{ example:05_one_trace.py }}

Two crossing faults should be passed **together**: a pillar claimed
by both retains two target residuals rather than silently giving
ownership to the last solver call. Run
`python -m fieldbook.cli run 06_crossing_traces`.

{{ example:06_crossing_traces.py }}

A dipping planar 3D polygon supplies distinct top and bottom XY
targets; it does **not** regenerate a 3D fault-surface FAULTS file.
Run `python -m fieldbook.cli run 07_dipping_fault`.

{{ example:07_dipping_fault.py }}

For a fault represented by adjacent planar facets, each facet
contributes its own target rows. This is **not** the seam-preserving
quad-surface FAULTS projection path: here only grid pillars move.
Run `python -m fieldbook.cli run 15_faceted_fault`:

{{ example:15_faceted_fault.py }}

Compare one variable at a time against a common grid and fault with
`08_settings`. Unpolished solutions can be valid but contain cells
below a project-defined scaled-Jacobian floor:

{{ example:08_settings.py }}

For a sheared outline, `slide` may report `ok=True` yet move the
physical boundary; `09_sheared_rim` measures the defect:

{{ example:09_sheared_rim.py }}

These outcomes are measured on a small synthetic fixture only.
See the [conditional recommendation](../decisions/settings.md).
