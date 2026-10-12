from pathlib import Path

from fieldbook.inventory import render


def test_declared_upstream_symbols_are_tracked():
    root = Path(__file__).resolve().parents[1]
    assert (root / "docs/reference/api-inventory.md").read_text(encoding="utf-8") == render()
