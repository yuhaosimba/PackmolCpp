#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile


def assert_pdb_basics(path: pathlib.Path) -> None:
    text = path.read_text()
    assert "HEADER" in text, "missing HEADER record"
    assert "END" in text, "missing END record"
    atom_lines = [line for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))]
    assert len(atom_lines) == 2, f"expected 2 atom lines, got {len(atom_lines)}"


def main() -> int:
    mode = sys.argv[1]
    probe = pathlib.Path(sys.argv[2])
    input_path = pathlib.Path(sys.argv[3])

    with tempfile.TemporaryDirectory() as tmpdir_name:
        tmpdir = pathlib.Path(tmpdir_name)

        if mode == "output":
            output_path = tmpdir / "generated.pdb"
            completed = subprocess.run(
                [str(probe), str(input_path), str(output_path)],
                capture_output=True,
                cwd=tmpdir,
            )
            stdout = completed.stdout.decode("utf-8", errors="replace")
            stderr = completed.stderr.decode("utf-8", errors="replace")
            if completed.returncode != 0:
                print(stdout, file=sys.stderr)
                print(stderr, file=sys.stderr)
                return 1
            if not output_path.exists():
                print("output probe did not create the requested pdb file", file=sys.stderr)
                return 1
            try:
                assert_pdb_basics(output_path)
            except AssertionError as exc:
                print(str(exc), file=sys.stderr)
                print(output_path.read_text(), file=sys.stderr)
                return 1
            return 0

        if mode == "checkpoint":
            completed = subprocess.run(
                [str(probe), str(input_path)],
                capture_output=True,
                cwd=tmpdir,
            )
            stdout = completed.stdout.decode("utf-8", errors="replace")
            stderr = completed.stderr.decode("utf-8", errors="replace")
            if completed.returncode != 0:
                print(stdout, file=sys.stderr)
                print(stderr, file=sys.stderr)
                return 1
            best = tmpdir / "output.pdb"
            forced = tmpdir / "output.pdb_FORCED"
            if not best.exists() or not forced.exists():
                print("checkpoint probe did not create both output files", file=sys.stderr)
                return 1
            if "ENDED WITHOUT PERFECT PACKING" not in stdout + stderr:
                print("checkpoint probe did not print the expected summary banner", file=sys.stderr)
                return 1
            try:
                assert_pdb_basics(best)
                assert_pdb_basics(forced)
            except AssertionError as exc:
                print(str(exc), file=sys.stderr)
                return 1
            return 0

    print(f"unsupported mode: {mode}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
