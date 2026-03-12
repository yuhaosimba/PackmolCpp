#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile


def main() -> int:
    probe = pathlib.Path(sys.argv[1])
    with tempfile.TemporaryDirectory() as tmpdir_name:
        tmpdir = pathlib.Path(tmpdir_name)
        completed = subprocess.run(
            [str(probe)],
            capture_output=True,
            text=True,
            cwd=tmpdir,
        )
        output = completed.stdout + completed.stderr
        log_path = tmpdir / "chkgrad.log"

        if completed.returncode != 0:
            print("comparegrad probe failed", file=sys.stderr)
            print(output, file=sys.stderr)
            return 1
        if "Comparing analytical and finite-difference" not in output:
            print("comparegrad probe did not announce the gradient comparison", file=sys.stderr)
            print(output, file=sys.stderr)
            return 1
        if "Done." not in output:
            print("comparegrad probe did not finish cleanly", file=sys.stderr)
            print(output, file=sys.stderr)
            return 1
        if not log_path.exists():
            print("comparegrad probe did not create chkgrad.log", file=sys.stderr)
            return 1

        log_text = log_path.read_text()
        if "Function Value =" not in log_text or "Maximum difference =" not in log_text:
            print("comparegrad log is missing expected summary lines", file=sys.stderr)
            print(log_text, file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
