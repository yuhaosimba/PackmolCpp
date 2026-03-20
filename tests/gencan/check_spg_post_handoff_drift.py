#!/usr/bin/env python3
from __future__ import annotations

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
    staged_input: pathlib.Path,
    workdir: pathlib.Path,
    mode: str,
    spg_handoff: bool,
) -> dict[str, float | int]:
    env = os.environ.copy()
    env["PACKMOL_GENCAN_IMPL"] = mode
    env["PACKMOL_GENCAN_NUMERIC_CPP"] = "0"
    if spg_handoff:
        env["PACKMOL_GENCAN_DEBUG"] = "1"
        env["PACKMOL_GENCAN_SPG_POST_CPP_HANDOFF"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_MAX_DEPTH"] = "1"

    completed = subprocess.run(
        [str(probe), str(staged_input)],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"probe failed mode={mode} spg_handoff={spg_handoff}:\n{completed.stdout}\n{completed.stderr}"
        )

    parsed = _parse_probe_output(completed.stdout + completed.stderr)
    required = ("inform", "f", "gpsupn", "xsum", "gnorm2", "iter", "fcnt", "gcnt", "cgcnt")
    missing = [key for key in required if key not in parsed]
    if missing:
        raise RuntimeError(f"missing keys in probe output: {missing}")
    return parsed


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    probe = repo / "build" / "tests" / "test_gencan_spg_post_ab_probe"
    fixture = repo / "tests" / "fixtures" / "tiny_packing.inp"

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)
        shutil.copy2(fixture.parent / "minimal_structure.pdb", tmpdir / "minimal_structure.pdb")

        baseline = _run_probe(probe, staged_input, tmpdir, mode="fortran", spg_handoff=False)
        handoff = _run_probe(probe, staged_input, tmpdir, mode="cpp", spg_handoff=True)

    float_keys = ("f", "gpsupn", "xsum", "gnorm2")
    int_keys = ("inform", "iter", "fcnt", "gcnt", "cgcnt")
    for key in int_keys:
        if int(baseline[key]) != int(handoff[key]):
            raise RuntimeError(f"drift on {key}: fortran={int(baseline[key])} cpp={int(handoff[key])}")
    for key in float_keys:
        if not math.isclose(float(baseline[key]), float(handoff[key]), rel_tol=1e-10, abs_tol=1e-12):
            raise RuntimeError(f"drift on {key}: fortran={float(baseline[key])} cpp={float(handoff[key])}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
