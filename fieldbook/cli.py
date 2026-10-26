"""List and run independent lessons; preserve script exit codes."""

import argparse
import importlib.metadata
import json
import re
import subprocess
import sys
from pathlib import Path

from fieldbook.downloads import fetch


EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def lessons() -> list[Path]:
    return sorted(EXAMPLES.glob("[0-9][0-9]_*.py"))


def upstream_revisions() -> dict[str, str]:
    text = (EXAMPLES.parent / "upstream-lock.toml").read_text(encoding="utf-8")
    revisions = {}
    for package in ("cpg2unstructured", "fault_projections"):
        match = re.search(
            rf"\[{package}\]\s*.*?^revision = \"([0-9a-f]{{40}})\"",
            text, flags=re.MULTILINE | re.DOTALL,
        )
        if match is None:
            raise ValueError(f"missing pinned 40-character revision for {package}")
        revisions[package] = match.group(1)
    return revisions


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the offline CPG/fault lessons")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list")
    commands.add_parser("doctor")
    download = commands.add_parser("fetch", help="explicit checksum-verified HTTPS download")
    download.add_argument("url")
    download.add_argument("--sha256", required=True)
    download.add_argument("--max-bytes", required=True, type=int)
    download.add_argument("--license", required=True, dest="license_name")
    download.add_argument("--attribution", required=True)
    download.add_argument("--out", required=True, type=Path)
    run = commands.add_parser("run")
    run.add_argument("lesson", help="lesson stem, e.g. 06_one_trace")
    run.add_argument("--out", type=Path, help="explicit output directory for export lessons")
    run.add_argument("--deck", type=Path, help="user-supplied GRDECL/EGRID for intake lesson")
    args = parser.parse_args(argv)
    if args.command == "list":
        for script in lessons():
            print(script.stem)
        return 0
    if args.command == "doctor":
        expected = upstream_revisions()
        for package, directory in (("cpg2unstructured", "cpg2python"),
                                   ("fault_projections", "fault-projections")):
            distribution = importlib.metadata.distribution(package.replace("_", "-"))
            print(f"{package}: {distribution.version}")
            sibling = EXAMPLES.parent.parent / directory
            if (sibling / ".git").exists():
                found = subprocess.run(
                    ["git", "-C", str(sibling), "rev-parse", "HEAD"],
                    check=True, capture_output=True, text=True,
                ).stdout.strip()
            else:
                direct_url = distribution.read_text("direct_url.json")
                info = json.loads(direct_url) if direct_url else {}
                found = info.get("vcs_info", {}).get("commit_id")
            if found and found != expected[package]:
                print(f"REVISION MISMATCH {package}: expected {expected[package]}, got {found}")
                return 1
            print(f"  source: {found or 'unverified (install from pinned source)'}")
        print(f"Python: {sys.version.split()[0]}; lessons: {len(lessons())}")
        return 0
    if args.command == "fetch":
        destination = fetch(args.url, args.out, args.sha256, max_bytes=args.max_bytes)
        provenance = destination.with_name(destination.name + ".source.json")
        provenance.write_text(
            json.dumps({
                "source": args.url, "sha256": args.sha256.lower(),
                "license": args.license_name, "attribution": args.attribution,
                "size_bytes": destination.stat().st_size,
            }, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"Verified {destination}; recorded source in {provenance}")
        return 0
    choices = {script.stem: script for script in lessons()}
    if args.lesson not in choices:
        parser.error(f"unknown lesson {args.lesson!r}; use 'list'")
    command = [sys.executable, str(choices[args.lesson])]
    if args.out is not None:
        command.extend(["--out", str(args.out)])
    if args.deck is not None:
        command.extend(["--deck", str(args.deck)])
    return subprocess.call(command)


if __name__ == "__main__":
    raise SystemExit(main())
