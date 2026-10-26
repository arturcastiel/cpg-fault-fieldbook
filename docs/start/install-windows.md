# Native Windows installation

Open PowerShell in the `cpg-fault-fieldbook` directory. Install Python
3.10–3.12, Git and the two upstream repositories as sibling directories.
Use the commits in `upstream-lock.toml`; a version string alone does not
prove the checkout matches.

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e "../cpg2python[fast]"
.\.venv\Scripts\python.exe -m pip install -e "../fault-projections"
.\.venv\Scripts\python.exe -m pip install -e ".[test,docs]"
.\.venv\Scripts\python.exe -m fieldbook.cli doctor
.\.venv\Scripts\python.exe -m fieldbook.cli run 01_read_grid
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m fieldbook.docs
```

Open `site\index.html` for the generated wiki. Use `--out output\demo`
with export lessons: `run 12_export_includes --out output\demo`.
The toy deck uses **arbitrary consistent units**; do not silently assume
metres or project a trace defined in a different MAPAXES frame. WSL
users follow the [Linux directions](install-linux.md) and use
`/mnt/c/...` paths instead. Windows CI validates installation and
the lesson suite; a reservoir simulator is not part of CI.
