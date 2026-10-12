# Offline fixture

`fieldbook/data/quickstart.grdecl` comes from
`cpg2python/examples/quickstart.grdecl` at
commit `e4637ccedf82c82f5c081bb909b354a7b919c6be`; upstream is
GPL-3.0-or-later. It is 5 x 5 x 3 with 70 active cells, five inactive cells,
and a PORO property. XY spans 0–5 in arbitrary consistent units; **do not
interpret these numbers as metres**.

`fieldbook/data/cart_2x2x2_no_actnum.EGRID` and
`fieldbook/data/cart_2x2x2_rotated.EGRID` are also verbatim upstream
fixtures from `cpg2python/tests/data/` at the same revision, under
GPL-3.0-or-later. They demonstrate binary reading, ACTNUM defaults and
MAPAXES frames. They are reader fixtures, not a validated reservoir model.

| Bundled file under `fieldbook/data/` | SHA-256 |
|---|---|
| `quickstart.grdecl` | `66a2af4640f66ab89f228d9258218c7fc9bba8cf8b9b91b19995aa5bc4813550` |
| `cart_2x2x2_no_actnum.EGRID` | `51f8a92578b12de9cb76fbaf7a9885061ca9587ee6f3e7d219c73e462350965b` |
| `cart_2x2x2_rotated.EGRID` | `313bd0ee6220a92393b35328963c68a241c61812ea55ece59ed0a29f86cab140` |

The fixture is copied unchanged so all essential lessons work offline and
without relying on a sibling checkout at runtime. Never commit private
reservoir models to this repository.
