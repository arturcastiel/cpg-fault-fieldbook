"""The numbered scripts shown in the wiki must actually finish successfully."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from fieldbook.cli import lessons, main

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("script", lessons(), ids=lambda script: script.stem)
def test_lesson(script, tmp_path):
    if script.stem == "13_quad_surface" and not (
        ROOT.parent / "fault-projections" / "examples" / "quad_mesh_fault.py"
    ).is_file():
        pytest.skip("quad lesson requires the pinned sibling source checkout")
    command = [sys.executable, str(script)]
    if script.stem in ("12_export_includes", "13_quad_surface", "15_visualize",
                       "18_case_dossier"):
        command.extend(["--out", str(tmp_path / "outputs")])
    if script.stem == "14_field_deck_intake":
        command.extend(["--deck", str(ROOT / "fieldbook" / "data" / "quickstart.grdecl")])
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True,
                            timeout=120, env=os.environ.copy())
    assert result.returncode == 0, result.stdout + result.stderr
    if script.stem == "12_export_includes":
        assert (tmp_path / "outputs" / "fault_adapted_grid.grdecl").is_file()
        assert not (tmp_path / "outputs" / "fault_adapted_faults.grdecl").exists()
    if script.stem == "15_visualize":
        assert (tmp_path / "outputs" / "before.svg").is_file()
        assert (tmp_path / "outputs" / "after.svg").is_file()
    if script.stem == "18_case_dossier":
        import json

        case = json.loads((tmp_path / "outputs" / "one-trace.json").read_text())
        assert case["verdict"] == "synthetic_geometry_gates_pass"
        assert case["simulator_validated"] is False
        assert len(case["input"]["sha256"]) == 64
        assert case["gates"]["quality_floor_ok"]


def test_cli_unknown_lesson_is_error(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["run", "not_a_lesson"])
    assert exc.value.code == 2
    assert "unknown lesson" in capsys.readouterr().err
