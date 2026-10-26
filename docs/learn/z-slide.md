# 5. Slide a fault in Z by editing ZCORN

A Z slide moves cell corners **down their own pillars**, so only ZCORN
changes. COORD, ACTNUM and the FAULTS keyword stay the same, which
leaves the pillar geometry and the fault surface untouched. Everything
derived from depth has to be **re-processed** from the new deck with
`read_grdecl` → `from_grdecl` → `build_topology`:

- which layers touch across the fault;
- the exposed fault pieces that become boundary;
- face areas and cell volumes.

The lessons below write each deck, read it back, and assert what the
rebuilt topology reports.

The lessons index corners with the corner cube
`z = zcorn.reshape(2*NZ, 2*NY, 2*NX)`. Slots `2c` and `2c+1` are the
near and far corners of cell `c` along each axis. On fault pillar
line `i`, slot `2i-1` belongs to the west column and slot `2i` to the
east column, so the throw is `z[:, :, 2*i] - z[:, :, 2*i-1]`.

> **Status.** No library call wraps the slide yet; the lessons write
> out the corner-cube edit in numpy. A planned ZCORN module in
> `cpg2unstructured` and a slide function in `fault_projections` will
> replace these helpers. The lessons need `cpg2unstructured` 0.3.0,
> whose `build_topology` pairs cells across a fault by depth overlap.

## Two faults, one after the other

Run `python -m fieldbook.cli run 19_z_slide_sequential`. Add
`--out DIR` to keep the three decks.

A rigid slide needs a **hanging wall**: the columns you can reach
from the east side of the fault's FAULTS faces without crossing them.
That flood fill cuts only the fault being slid. Sliding F1 therefore
moves blocks M and E together, carrying F2 with them, and F2 has no
throw until it is slid itself.

![Two sequential faults, three stages](../assets/z-slide/sequential-faults.png)

**What to expect:**

- **Throws:** F1 is 7 m, then F2 is 12 m, both exact on every pillar.
- **Unchanged keywords:** COORD and FAULTS are unchanged in all three
  decks.
- **Slides commute:** F2-then-F1 is bit-identical to F1-then-F2,
  which equals one block map (W 0, M 7, E 19 m).
- **Contacts:** after re-processing, the contact pairs match the
  depth-interval overlaps. For example, on F2 after 12 m the contact
  W k1 | E k0 is 3,578 m².
- **Volume side effect on tilted faults:** a corner that moves down a
  tilted pillar also moves in x.
  - Cells bounded by a moved tilted pillar change volume: columns 3,
    5 and 6 here, from −8.8 % to +4.2 %.
  - Cell thickness never changes.
  - Column 5 widens: it is the far side of block M, and F2's pillar
    moved with M.

{{ example:19_z_slide_sequential.py }}

## A fault that ends inside the grid

Run `python -m fieldbook.cli run 20_z_slide_partial_fault`. Add
`--out DIR` to keep the decks.

FAULTS `P` covers only rows 0–3 of pillar line 3. Around the tip
pillar, the east side is still joined to the west side, so the flood
fill finds no closed block and a rigid slide is impossible.

Instead, the throw is a **field on pillars**: 8, 6, 4, 2 and 0 m along
the fault, applied to every corner east of the line.

- **Zero at the tip:** the four cells around the tip pillar share its
  corners, so the throw there must be zero.
- **Welded past the tip:** the throw stays zero beyond the tip, so
  those faces stay welded.

![Partial fault with a tapered throw](../assets/z-slide/partial-fault.png)

**What to expect:**

- **Faces that open:** only the four FAULTS faces open. Every other
  face, including those past the tip, keeps shared corners.
- **Contact shapes:** contacts become trapezoids, and triangles in the
  last row before the tip.
- **Contact areas:** the rebuilt areas match the analytic overlaps.
  For example, W k1 | E k0 goes 783 → 559 → 335 → 112 m² towards the
  tip.
- **Current limit:** a 15 m taper passes the 10 m layer thickness
  inside one cell face. The rebuilt fault polygon would need five
  corners, so `build_topology` raises `NonMatchingZError`. Keep each
  face's throw range clear of layer-thickness multiples until that
  limit is lifted.

{{ example:20_z_slide_partial_fault.py }}

## An oblique, dipping fault

Run `python -m fieldbook.cli run 21_z_slide_oblique_fault`. Add
`--out DIR` to keep the four decks.

At the top (z = 0), the fault trace runs from (300, 0) to (600, 600),
27° off the J axis, so no single pillar line carries it.
`fault_faces_for` turns it into a **staircase** of 6 X faces and 3 Y
faces. The flood fill and the throw check therefore treat both
directions:

- an X face at pillar line `i` compares slots `2i-1 | 2i` along I;
- a Y face at pillar row `j` compares slots `2j-1 | 2j` along J.

The staircase pillars keep their tops on the trace and lean 40 m east
at 80 m depth. Every fault pillar leans by the same vector, so the
fault dips 65.9° to the south-east. Each FAULTS face is a planar
parallelogram, with area `|edge × lean|` per metre of depth.

The lesson runs the same slide on two decks with identical FAULTS:

- **A, regular grid:** the staircase follows the grid lines, so the
  fault is 900 m long in plan and 970.8 m² per metre of depth.
- **B, adapted grid:** `adapt_grid(..., boundary="pinned")` takes a
  2D trace (stage 1). It moves pillars rigidly, so the fault pillars
  land on the trace and keep their lean. The FAULTS faces are
  unchanged (`faces_invariant`). The fault is 711.3 m long in plan and
  775.1 m² per metre of depth.

A stage-2 dipping polygon is not used here because `boundary="pinned"`
keeps the rim pillars vertical, and this fault ends on the rim.

![Oblique dipping fault, regular and adapted grids](../assets/z-slide/oblique-fault.png)

**What to expect:**

- **The slide:** in both decks the south-east hanging wall of 21
  columns slides 7 m. COORD and FAULTS are unchanged.
- **Slip along the fault:** a corner on a fault pillar slides 7.83 m
  along it (Δx +3.5 m, Δz +7 m), so it stays on the fault surface.
  The slip follows the pillar lean, which is not exactly the dip
  direction.
- **Faces that open:** exactly the 9 FAULTS faces lose shared
  corners, X and Y alike.
- **Contacts:** the same five layer pairs touch in both decks. Each
  area is the fault area per metre of depth times the depth overlap,
  so A carries 25 % more contact area than B. For example, fixed k1 |
  moving k0 is 6,796 m² in A and 5,426 m² in B.
- **Volumes:** only moving cells that touch a fault pillar change
  volume (27 cells, −4.0 % to −1.6 %). The fixed block is unchanged.
- **Steep traces:** keep `adapt_grid` traces gentle on small grids.
  At about 30° or more off the axis here, the regenerated FAULTS
  include a face on the grid rim and `faces_invariant` is False. The
  staircase from `fault_faces_for` on the regular grid is still
  usable.

{{ example:21_z_slide_oblique_fault.py }}

The figures are drawn from the decks that `--out` writes. The plotted
numbers are the ones the lessons print. All cases are synthetic layer
cakes: nothing here checks simulator acceptance, transmissibility or
flow. See also the [known limits](../reference/limits.md).
