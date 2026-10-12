import hashlib
from pathlib import Path

from fieldbook.fixtures import egrid, quickstart

HASHES = {
    "quickstart.grdecl": "66a2af4640f66ab89f228d9258218c7fc9bba8cf8b9b91b19995aa5bc4813550",
    "cart_2x2x2_no_actnum.EGRID": "51f8a92578b12de9cb76fbaf7a9885061ca9587ee6f3e7d219c73e462350965b",
    "cart_2x2x2_rotated.EGRID": "313bd0ee6220a92393b35328963c68a241c61812ea55ece59ed0a29f86cab140",
}


def test_bundled_fixture_hashes():
    for path in (quickstart(), egrid(), egrid(rotated=True)):
        assert hashlib.sha256(path.read_bytes()).hexdigest() == HASHES[path.name]


def test_upstream_fixture_matches_when_checkout_is_present():
    upstream = Path(__file__).resolve().parents[2] / "cpg2python"
    if upstream.exists():
        for local, source in (
            (quickstart(), upstream / "examples" / "quickstart.grdecl"),
            (egrid(), upstream / "tests" / "data" / egrid().name),
            (egrid(rotated=True), upstream / "tests" / "data" / egrid(rotated=True).name),
        ):
            assert local.read_bytes() == source.read_bytes()
