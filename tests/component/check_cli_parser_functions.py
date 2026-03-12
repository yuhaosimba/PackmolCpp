#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import subprocess
import sys


def main() -> int:
    probe = pathlib.Path(sys.argv[1])
    completed = subprocess.run(
        [str(probe), "-i", "fixtures/input.inp", "-o", "fixtures/output.pdb"],
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        print("cli parser helper probe failed", file=sys.stderr)
        print(completed.stdout, file=sys.stderr)
        print(completed.stderr, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
