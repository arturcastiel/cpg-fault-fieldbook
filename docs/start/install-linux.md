# Install on Linux or WSL

Requirements: Python 3.10+, Git, and the two upstream checkouts. Clone
these projects **beside** the fieldbook. The tested commits are recorded
in the repository-root `upstream-lock.toml`; do not assume a
same-version branch is a same-revision implementation.

```bash
cd /path/to/library-development/cpg-fault-fieldbook
python3 -m venv .venv
.venv/bin/python -m pip install -e "../cpg2python[fast]"
.venv/bin/python -m pip install -e "../fault-projections"
.venv/bin/python -m pip install -e ".[test,docs]"
.venv/bin/python -m fieldbook.cli doctor
.venv/bin/python -m fieldbook.cli run 01_read_grid
.venv/bin/python -m pytest -q
.venv/bin/python -m fieldbook.docs
```

If using preinstalled dependencies, omit installs and invoke your active
interpreter. The beginner lessons run offline with the bundled licensed
fixture; the quad demo (`13_quad_surface`) additionally requires the
sibling upstream source checkout. No download occurs without an explicit
call to [`fetch`](../reference/data.md).

**Revisions:** `git -C ../cpg2python rev-parse HEAD` and
`git -C ../fault-projections rev-parse HEAD` should match the lock file.
Building static HTML creates `site/index.html`; serve it locally or
open it directly. Nothing here installs a reservoir simulator.
