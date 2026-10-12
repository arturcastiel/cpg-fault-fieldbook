# Data, downloads and coordinate frames

The bundled toy GRDECL is licensed and works offline. No external
download is required for any lesson except the optional use of an
existing sibling repository's quad example. No field deck is
redistributed.

Optional downloads must be individually licensed and attributed by
the operator. `fieldbook.downloads.fetch(url, destination, sha256,
max_bytes=...)` requires HTTPS, a known SHA-256 and a byte limit,
downloads to a temporary file, verifies the digest, then atomically
replaces the destination. It is **never** invoked by setup, wiki build
or CI. The opt-in CLI additionally requires license and attribution
and writes a `.source.json` record beside the file:

```bash
python -m fieldbook.cli fetch 'https://example.org/path/to/approved.deck' \
  --sha256 '<known-64-character-hex-digest>' --max-bytes 104857600 \
  --license '<verified-license>' --attribution '<source and author>' \
  --out downloads/approved.deck
```

The URL and hash above are **placeholders**, not a real dataset. Obtain
them from the data publisher and inspect the permission to redistribute
before fetching. There is no approved public field-data catalog in this
initial release; inventing a URL/hash or redistributing a private
deck would make the guide less trustworthy.

GRDECL `COORD` pillars, EGRID MAPAXES, fault trace XY and exported
INCLUDE files must share one coordinate frame and one unit system.
The toy XY range 0–5 has no asserted physical unit. Convert fault
coordinates **before** projection/adaptation, and inspect the grid
writer's MAPAXES guard before output. `read_grdecl` and `read_egrid`
return the same corner-grid data type; `from_grdecl` pivots either
into a cube for topology.
