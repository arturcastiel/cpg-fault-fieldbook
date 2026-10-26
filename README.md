# CPG + Fault Fieldbook

A runnable, **experimental** companion to [cpg2python](../cpg2python/README.md)
(`cpg2unstructured`) and [fault-projections](../fault-projections/README.md)
(`fault_projections`). The libraries do the work; this repository teaches their
calls and independently checks what each example actually establishes. Grid
adaptation is not fracture propagation, transmissibility modeling, or a
simulator-ready deck.

## Start here

Python 3.10+; from a checkout beside both upstream repositories:

```bash
python -m venv .venv
# Linux/WSL: .venv/bin/python ; Windows PowerShell: .venv\Scripts\python.exe
.venv/bin/python -m pip install -e "../cpg2python[fast]"
.venv/bin/python -m pip install -e "../fault-projections"
.venv/bin/python -m pip install -e ".[test]"
.venv/bin/python -m fieldbook.cli doctor
.venv/bin/python -m fieldbook.cli list
.venv/bin/python -m fieldbook.cli run 01_read_grid
.venv/bin/python -m pytest -q
```

On Windows replace `.venv/bin/python` with `.venv\Scripts\python.exe`;
see [Windows instructions](docs/start/install-windows.md). The copies in
`../cpg2python` and `../fault-projections` should match
[`upstream-lock.toml`](upstream-lock.toml). Verify revisions with
`git -C ../cpg2python rev-parse HEAD` and
`git -C ../fault-projections rev-parse HEAD`. The editable installs above are
for local learning; changing either checkout changes results. Do not silently
assume package-version equality implies a matching commit.

The numbered scripts in [`examples/`](examples/) run offline against a
redistributable [five-by-five grid](fixtures/README.md); no private deck is
bundled or fetched. Scripts that emit files require an explicit `--out`
directory. Start with [the wiki](docs/index.md) or [the original detailed
standalone guide](docs/legacy-field-guide.html). Use
`python -m fieldbook.docs` to build static HTML into `site/` after installing
`.[docs]` and run `python -m fieldbook.docs --check` to verify that example
inclusions and local links resolve.

**CI.** `fault-projections` is a private repository, so the workflow needs a
repository secret `FAULT_PROJECTIONS_TOKEN`: a fine-grained personal access
token with *Contents: read-only* on `arturcastiel/fault-projections` only.
Set it with `gh secret set FAULT_PROJECTIONS_TOKEN -R
arturcastiel/cpg-fault-fieldbook`. Pull requests from forks do not receive
secrets, so their CI stops at that step.

## What is covered now

| Track | Runnable source | What you learn |
|---|---|---|
| Grid | `01_read_grid`, `02_egrid`, `03_mesh_graph` | Text/binary input, MAPAXES, pivot, topology and properties |
| Faults | `04_project_trace`, `05_project_surface` | 2.5D FAULTS and planar 3D faces |
| Adapt | `06_one_trace`, `07_crossing_traces`, `08_dipping_fault` | Targets, high-level and low-level solves |
| Judge | `09_settings`, `10_sheared_rim`, `11_four_gates` | Fidelity, quality, rim displacement and corridor |
| Export | `12_export_includes` | Separate grid INCLUDE, no master-deck rewrite |
| Real-deck intake | `14_field_deck_intake --deck PATH` | Read-only dimensions, ACTNUM, frame and top/bottom quality |
| Figures | `15_visualize --out PATH` | Same-scale before/after SVGs plus numerical gates |
| Curved fault | `16_faceted_fault` | Two joined planar facets drive independent top/bottom targets |
| Coordinate frame | `17_mapaxes_guard` | Expected writer error for transformed EGRID geometry |
| Reproducibility | `18_case_dossier --out PATH` | JSON case report with input hash, revisions, settings and gates |

The [function reference](docs/functions/index.md) documents high-impact calls;
the [API coverage inventory](docs/reference/coverage.md) distinguishes
documented workflows from the rest of the exported surface. The existing
full-length HTML field guide is retained as a dated snapshot, not an
automatically updated reference. An example passing on the bundled fixture
does **not** qualify an arbitrary field deck. Optional download support is
checksum-checked, requires an explicitly supplied source/license/attribution,
and never runs during tests or docs builds. See the
[opt-in data guide](docs/reference/data.md); there is no preapproved
public field-deck catalog in this release.

## Licensing and data

This companion and its upstream-derived quickstart fixture are distributed
under GPL-3.0-or-later. See [LICENSE](LICENSE) and
[fixture provenance](fixtures/README.md). Neither upstream implementation
is copied into this package. Do not commit private decks or generated reports.
