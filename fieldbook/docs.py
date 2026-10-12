"""Check documented snippets and build a static MkDocs wiki on demand."""

import argparse
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlsplit, unquote

from fieldbook.docs_hook import EXAMPLE, ROOT, expand


LINK = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)")
DIRECTIVE = re.compile(r"\{\{\s*example:[^}]+\}\}")


def check() -> None:
    docs = ROOT / "docs"
    pages = list(docs.rglob("*.md"))
    for page in pages:
        text = page.read_text(encoding="utf-8")
        for directive in DIRECTIVE.findall(text):
            if not EXAMPLE.fullmatch(directive):
                raise ValueError(f"{page}: invalid example directive: {directive}")
        expand(text)  # raises for missing scripts
        for target in LINK.findall(text):
            parsed = urlsplit(target)
            if parsed.scheme or not parsed.path:
                continue
            path = (page.parent / unquote(parsed.path)).resolve()
            if not path.exists():
                raise FileNotFoundError(f"{page}: broken link {target}")
    print(f"Checked {len(pages)} wiki pages; source examples and local links resolve")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="check only; do not build HTML")
    args = parser.parse_args(argv)
    check()
    if not args.check:
        subprocess.run(
            [sys.executable, "-m", "mkdocs", "build", "--strict",
             "--config-file", str(ROOT / "mkdocs.yml")],
            check=True,
            cwd=ROOT,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
