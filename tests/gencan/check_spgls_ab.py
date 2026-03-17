#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import pathlib
import subprocess
import sys


MODE_ID = {
    "fortran": 0,
    "cpp": 1,
    "ab": 2,
}


def parse_output(text: str) -> dict[str, float | int]:
    values: dict[str, float | int] = {}
    for line in text.splitlines():
        tokens = line.strip().split()
        if tokens and tokens[0] == "mode" and len(tokens) >= 12:
            values["mode"] = int(tokens[1])
            values["inform"] = int(tokens[3])
            values["fcnt"] = int(tokens[5])
            values["f"] = float(tokens[7])
            values["xsum"] = float(tokens[9])
            values["dnorm2"] = float(tokens[11])
    return values


def run_probe(probe: pathlib.Path, mode: str) -> dict[str, float | int]:
    env = os.environ.copy()
    env["PACKMOL_GENCAN_IMPL"] = mode
    completed = subprocess.run(
        [str(probe)],
        capture_output=True,
        text=True,
        env=env,
    )
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        print(f"spgls probe failed for mode={mode}", file=sys.stderr)
        print(output, file=sys.stderr)
        raise RuntimeError("probe failed")
    parsed = parse_output(output)
    if parsed.get("mode") != MODE_ID[mode]:
        print(output, file=sys.stderr)
        raise RuntimeError(f"unexpected mode marker for mode={mode}")
    required = ("inform", "fcnt", "f", "xsum", "dnorm2")
    missing = [k for k in required if k not in parsed]
    if missing:
        print(output, file=sys.stderr)
        raise RuntimeError(f"missing keys for mode={mode}: {missing}")
    return parsed


def assert_close(name: str, lhs: float, rhs: float, atol: float = 1e-12, rtol: float = 1e-10) -> None:
    if not math.isclose(lhs, rhs, rel_tol=rtol, abs_tol=atol):
        raise RuntimeError(f"{name} mismatch: {lhs} vs {rhs}")


def main() -> int:
    probe = pathlib.Path(sys.argv[1])

    ref = run_probe(probe, "fortran")
    cpp = run_probe(probe, "cpp")
    ab = run_probe(probe, "ab")

    for key in ("inform", "fcnt"):
        if int(ref[key]) != int(cpp[key]) or int(ref[key]) != int(ab[key]):
            raise RuntimeError(f"{key} mismatch: fortran={ref[key]} cpp={cpp[key]} ab={ab[key]}")

    for key in ("f", "xsum", "dnorm2"):
        assert_close(key, float(ref[key]), float(cpp[key]))
        assert_close(f"{key}_ab", float(ref[key]), float(ab[key]))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
