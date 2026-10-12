# Grid input: `read_grdecl` and `read_egrid`

Both readers produce the same format-agnostic `CornerPointGrid`. The text
reader also retains selected cell properties and records names of unrecognized
keywords; the EGRID reader reads only the root corner-point grid geometry and
`ACTNUM`.

## `read_grdecl`

### Syntax

```python
read_grdecl(path: str | pathlib.Path) -> Grdecl
```

`Grdecl` is an alias for `CornerPointGrid`, not a separate text-only class.
The argument is the path to a GRDECL text file.

### Input and output

The parser requires `DIMENS` or `SPECGRID`, plus `COORD` and `ZCORN`.
`DIMENS` and `SPECGRID`, if both present, must agree. `ACTNUM` is optional and
defaults to all cells active. `INCLUDE 'file.inc' /` is followed recursively;
the included path is resolved relative to the file containing that directive.

The returned object carries:

| Field | Shape and meaning |
|---|---|
| `cart_dims` | `(NX, NY, NZ)` |
| `coord` | `( (NX+1)*(NY+1), 6 )`, columns `(xt, yt, zt, xb, yb, zb)`; text input is `float64` |
| `zcorn` | flat `(8*NX*NY*NZ,)`, `float64` |
| `actnum` | flat `(NX*NY*NZ,)`, `int8` |
| `props` | dict of supported cell-property keywords to flat arrays |
| `unrecognized` | names of skipped, unsupported keywords |
| `source` | `"grdecl"` |

Supported property names include `PORO`, `PERMX/PERMY/PERMZ`, `NTG`,
`MULTX/MULTY/MULTZ`, `PRESSURE`, saturation keywords, and the integer region
keywords `SATNUM`, `PVTNUM`, `EQLNUM`, `FIPNUM`. Property arrays are checked
against `NX*NY*NZ`; region identifiers use `int32`.

### Source-based variations

```python
from cpg2unstructured.grdecl_reader import read_grdecl

grid = read_grdecl("FIELD.GRDECL")
print(grid.cart_dims, grid.coord.shape, grid.zcorn.size)
poro = grid.props.get("PORO")           # present only if supplied
unknown_names = grid.unrecognized
```

For an include-based model, point `read_grdecl` at the master text deck. The
reader follows geometry/property includes; it is not a complete Eclipse deck
interpreter and does not preserve arbitrary simulator sections for writing.

### Failure modes and limits

- Missing dimensions/geometry, conflicting dimensions, or wrong array lengths
  raise `GrdeclConsistencyError`; malformed blocks, unsupported token syntax,
  and missing `/` terminators raise `GrdeclSyntaxError`.
- Include cycles raise `GrdeclIncludeError`. Whitespace within a quoted
  include filename is unsupported by the whitespace tokenizer.
- Duplicate supported property keywords warn and use the last value.
- Unknown uppercase keywords are skipped to the next `/` and listed by name.
  Their data are not retained; this is not a safe read-modify-write path for a
  complete simulator deck.

### Maturity

Implemented and covered by parser tests for dimensions, shapes, optional
`ACTNUM`, properties, includes, run-length values, and malformed inputs.
Coverage does not imply support for every Eclipse keyword or deck section.

## `read_egrid`

### Syntax

```python
read_egrid(
    path: str | pathlib.Path,
    *,
    apply_mapaxes: bool = True,
    dtype: numpy.dtype | type = numpy.float32,
) -> CornerPointGrid
```

Accepts binary `.EGRID` and formatted `.FEGRID`; `resfo` detects the encoding.

### Input and output

The reader streams the file and retains `GRIDHEAD`, `COORD`, `ZCORN`, optional
`ACTNUM`, and optional `MAPAXES`. The returned `CornerPointGrid` has the same
field shapes as the text reader. `COORD` and `ZCORN` default to `float32`;
`dtype=numpy.float64` widens them. `ACTNUM` is `int8`; stored raw `mapaxes` is
`float64`. `props` and `unrecognized` are empty because this reader handles
grid geometry rather than cell properties.

When `MAPAXES` exists, the default applies its planar rotation to the X/Y values
of both endpoints in `COORD`. `mapaxes` retains the original six values and
`mapaxes_applied` records whether the transform was performed.

### Source-based variations

```python
import numpy as np
from cpg2unstructured.egrid_reader import read_egrid

field = read_egrid("FIELD.EGRID")  # MAPAXES is applied when present
wide = read_egrid("FIELD.EGRID", dtype=np.float64)
raw_local = read_egrid("FIELD.EGRID", apply_mapaxes=False)
# raw_local warns if MAPAXES exists; coords then remain in the local frame.
```

The same pivot can be handed to `from_grdecl`:

```python
from cpg2unstructured.cpg import from_grdecl

cube = from_grdecl(field)
```

This is the hand-off shown by the sibling `read_egrid_demo.py`; it does not
mean that the EGRID reader itself builds mesh topology.

### Failure modes and limits

- Missing required EGRID keywords raise `GrdeclSyntaxError`; malformed
  dimensions or lengths raise `GrdeclConsistencyError`.
- Only corner-point grids (`GRIDHEAD[0] == 1`) are supported. LGR grids are not
  exposed; the reader stops at `ENDGRID` and first-occurrence handling selects
  root-grid records.
- If `MAPAXES` is present and `apply_mapaxes=False`, a `RuntimeWarning` makes
  the unrotated coordinate choice explicit. Keep fault inputs in the same
  coordinate frame as the grid.
- The text reader and EGRID reader are not interchangeable for cell-property
  ingestion: EGRID `props` is empty.

### Maturity

Implemented; the reader has a property-based round-trip test and a source
example for MAPAXES and dtype choices. EGRID fixture tests and the demo are
linked below. Downstream mesh construction has its own topology limitations.

## Source and tests

- [Text reader source](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/cpg2unstructured/grdecl_reader.py)
  · [grid data type and `Grdecl` alias](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/cpg2unstructured/corner_point.py)
- [GRDECL parser tests](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/tests/test_grdecl_reader.py)
- [EGRID reader source](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/cpg2unstructured/egrid_reader.py)
  · [EGRID demo](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/examples/read_egrid_demo.py)
  · [EGRID property tests](https://github.com/arturcastiel/cpg2python/blob/e4637ccedf82c82f5c081bb909b354a7b919c6be/tests/test_egrid_property.py)
- [Broader fieldbook context](../legacy-field-guide.html)

