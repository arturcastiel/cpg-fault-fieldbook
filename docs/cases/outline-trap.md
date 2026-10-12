# Case: a valid `slide` morph escapes a sheared physical outline

**Fixture:** the bundled 5 × 5 × 3 grid, with both endpoints
sheared by `x_world = x_local + 0.5*y`; an authored rim-to-rim
trace transformed into the same frame. Run
`python -m fieldbook.cli run 09_sheared_rim`.

| Boundary | `result.ok` | Max physical outline deviation |
|---|---|---:|
| `pinned` | True | 0 |
| `slide` | True | ~0.1576 grid units |
| `free` | True | ~0.4966 grid units |

**Verdict:** `slide` fails this outline-preservation requirement
even though its cells pass the library's validity gate. The measured
deviation is *perpendicular to the original physical edges*, not
merely the permissible tangential rim displacement.
`free` intentionally changes the footprint. On a sheared/rotated
outline use `pinned` pending a fix and explicitly compare the
original and deformed rim in the deck's coordinate frame.
These numbers illustrate the defect on a tiny synthetic grid;
they are not field measurements.
