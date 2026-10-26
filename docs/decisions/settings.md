# Which deformation settings should I start with?

**Conditional starting point:** elastic smoothing, `w0=0.05`, barrier
polish, active-column boundary/validity scopes and `on_invalid="raise"`.
Use `pinned` for interior faults or non-axis-aligned outlines; consider
`slide` for a rim-touching fault on an axis-aligned outline. Do not
choose `free` merely to make the trace fit: its apparent gain may
come from moving the reservoir footprint.

| Priority | Start with | Trade-off |
|---|---|---|
| Preserve an interior outline | `pinned`, `elastic`, `barrier` | The corrected study measured 8.39 m interior residual vs 8.36 m for `slide`; avoids rim motion. |
| Rim-touching, axis-aligned | `slide`, `elastic`, `barrier` | 16.97 m touching residual, no newly sub-0.3 cells in the 51 matching study runs; only five non-degenerate touching faults informed the mean. |
| Sheared/rotated outline | `pinned` | `slide` moved some sheared edges by up to 47.4 m despite `ok=True`. |
| More fidelity, accept quality risk | `polish_fidelity=20.0` | 4.51 m interior vs 8.36 m standard barrier, but 10/51 runs gained cells below SJ 0.3. |
| TPFA orthogonality | `polish_ortho=2.0` | 11.72 m interior; non-orthogonality p95 delta 10.5° vs 15.8°. |

Those numbers come from the corrected local boundary study
(`/home/arturcastiel/projects/results-faults/boundary-study/RECOMMENDATION.md`):
408 configurations on a 12 × 12 two-layer Cartesian grid, plus
supporting sheared/graded probes. The earlier benchmark's
`lap_free_polish` **composite-score winner is not a fixed-extent
recommendation**: `free` moves the outline. The main sweep did not
settle `support_rings`, `poisson`, field grids, 3D fault surfaces or
irregular ACTNUM. The [local source-linked snapshot](../legacy-field-guide.html)
has additional limitations and citations.

Use `09_settings` and `10_sheared_rim` to reproduce the *mechanism*
on a different, tiny fixture; their unitless numbers are not the
metre-valued boundary-study numbers. Compare authored-fault distance,
top and bottom quality, fault corridor, chain and outline before
retaining a candidate. An initial weight is a starting point, not an
universal optimum.
