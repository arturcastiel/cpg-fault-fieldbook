from pathlib import Path

from fieldbook.cli import lessons
from fieldbook.docs import check
from fieldbook.docs_hook import expand

ROOT = Path(__file__).resolve().parents[1]


def test_all_examples_are_referenced_in_wiki():
    pages = "\n".join(page.read_text(encoding="utf-8")
                      for page in (ROOT / "docs").rglob("*.md"))
    for script in lessons():
        assert f"{{{{ example:{script.name} }}}}" in pages


def test_snippets_and_links():
    check()
    generated = expand("{{ example:01_read_grid.py }}")
    assert "```python" in generated and "read_grdecl" in generated
