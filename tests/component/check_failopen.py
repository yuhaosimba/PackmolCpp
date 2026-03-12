#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import subprocess
import sys


def main() -> int:
    probe = pathlib.Path(sys.argv[1])
    completed = subprocess.run([str(probe)], capture_output=True, text=True)
    output = completed.stdout + completed.stderr
    if completed.returncode == 0:
        print("failopen probe unexpectedly exited with code 0", file=sys.stderr)
        return 1
    if "Could not open file" not in output:
        print("failopen probe did not emit the expected error message", file=sys.stderr)
        print(output, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
