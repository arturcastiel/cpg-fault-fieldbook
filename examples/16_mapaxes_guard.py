"""Show why applied MAPAXES coordinates cannot silently round-trip."""

from tempfile import TemporaryDirectory
from pathlib import Path

from cpg2unstructured import MapaxesNotInvertedError, read_egrid, write_grid_include
from fieldbook.fixtures import egrid

grid = read_egrid(egrid(rotated=True))
assert grid.mapaxes_applied
with TemporaryDirectory() as temporary:
    path = Path(temporary) / "unsafe.grdecl"
    try:
        write_grid_include(grid, path)
    except MapaxesNotInvertedError as error:
        print("Expected MAPAXES frame guard:", error)
    else:
        raise AssertionError("writer accepted transformed coordinates without opt-in")
    assert not path.exists()
print("To override, explicitly set allow_mapaxes_applied=True after reviewing deck frame.")
