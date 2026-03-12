#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import subprocess
import sys


def main() -> int:
    probe = pathlib.Path(sys.argv[1])
    completed = subprocess.run([str(probe)], capture_output=True, text=True)
    output = completed.stdout + completed.stderr

    if completed.returncode != 0:
        print("output messages probe failed", file=sys.stderr)
        print(output, file=sys.stderr)
        return 1

    required_snippets = [
        "PACKMOL - Packing optimization",
        "Version 21.2.1",
        "Packing solved for molecules of type",
        "Success!",
        "Final objective function value",
        "Maximum violation of target distance",
        "********************",
    ]
    for snippet in required_snippets:
        if snippet not in output:
            print(f"missing expected snippet: {snippet}", file=sys.stderr)
            print(output, file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
