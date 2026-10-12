"""MkDocs hook: pull snippets verbatim from the scripts pytest executes."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = re.compile(r"\{\{ example:([0-9][0-9]_[a-z_]+\.py) \}\}")


def expand(markdown: str) -> str:
    def replace(match: re.Match[str]) -> str:
        path = ROOT / "examples" / match.group(1)
        return "```python\n" + path.read_text(encoding="utf-8").rstrip() + "\n```"

    return EXAMPLE.sub(replace, markdown)


def on_page_markdown(markdown, page, config, files):
    return expand(markdown)
