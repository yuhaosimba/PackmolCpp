#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import pathlib
import shutil
import subprocess
import tempfile


def _parse_probe_output(text: str) -> dict[str, float | int]:
    values: dict[str, float | int] = {}
    for line in text.splitlines():
        tokens = line.strip().split()
        if not tokens:
            continue
        if tokens[0] == "mode" and len(tokens) >= 18:
            values["mode"] = int(tokens[1])
            values["inform"] = int(tokens[3])
            values["f"] = float(tokens[5])
            values["gpsupn"] = float(tokens[7])
            values["xsum"] = float(tokens[9])
            values["iter"] = int(tokens[11])
            values["fcnt"] = int(tokens[13])
            values["gcnt"] = int(tokens[15])
            values["cgcnt"] = int(tokens[17])
        elif tokens[0] == "gnorm2" and len(tokens) >= 2:
            values["gnorm2"] = float(tokens[1])
    return values


def _run_probe(
    probe: pathlib.Path,
    input_path: pathlib.Path,
    mode: str,
    workdir: pathlib.Path,
    handoff: bool,
    max_depth: int,
    cpp_replay: bool,
) -> dict[str, float | int]:
    env = os.environ.copy()
    env["PACKMOL_GENCAN_IMPL"] = mode
    env["PACKMOL_GENCAN_NUMERIC_CPP"] = "0"
    if handoff:
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_MAX_DEPTH"] = str(max_depth)
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY"] = "1" if cpp_replay else "0"

    completed = subprocess.run(
        [str(probe), str(input_path)],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"probe failed mode={mode} handoff={handoff}:\n{completed.stdout}\n{completed.stderr}"
        )
    parsed = _parse_probe_output(completed.stdout + completed.stderr)
    required = ("inform", "f", "gpsupn", "xsum", "gnorm2", "iter", "fcnt", "gcnt", "cgcnt")
    missing = [k for k in required if k not in parsed]
    if missing:
        raise RuntimeError(f"missing keys in probe output: {missing}")
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description="Diagnose TN post handoff drift against Fortran baseline.")
    parser.add_argument("probe", type=pathlib.Path)
    parser.add_argument("fixture", type=pathlib.Path)
    parser.add_argument("--max-depth", type=int, default=1)
    parser.add_argument("--strict", action="store_true", help="Exit non-zero when drift is detected.")
    parser.add_argument("--cpp-replay", action="store_true", help="Enable cpp replay in handoff canonicalization.")
    args = parser.parse_args()
    probe = args.probe.resolve()
    fixture = args.fixture.resolve()

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)
        structure_src = fixture.parent / "minimal_structure.pdb"
        if structure_src.exists():
            shutil.copy2(structure_src, tmpdir / structure_src.name)

        baseline = _run_probe(
            probe,
            staged_input,
            "fortran",
            tmpdir,
            handoff=False,
            max_depth=args.max_depth,
            cpp_replay=False,
        )
        handoff = _run_probe(
            probe,
            staged_input,
            "cpp",
            tmpdir,
            handoff=True,
            max_depth=args.max_depth,
            cpp_replay=args.cpp_replay,
        )

    float_keys = ("f", "gpsupn", "xsum", "gnorm2")
    int_keys = ("inform", "iter", "fcnt", "gcnt", "cgcnt")
    drift = False
    print("=== TN handoff drift ===")
    for key in int_keys:
        lhs = int(baseline[key])
        rhs = int(handoff[key])
        same = lhs == rhs
        drift = drift or (not same)
        print(f"{key:>6}: fortran={lhs:>6} handoff_cpp={rhs:>6} {'OK' if same else 'DIFF'}")
    for key in float_keys:
        lhs = float(baseline[key])
        rhs = float(handoff[key])
        same = math.isclose(lhs, rhs, rel_tol=1e-10, abs_tol=1e-12)
        drift = drift or (not same)
        print(f"{key:>6}: fortran={lhs:.16e} handoff_cpp={rhs:.16e} {'OK' if same else 'DIFF'}")

    if args.strict and drift:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
