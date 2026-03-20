#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import re
import shutil
import subprocess
import tempfile
import math


SUMMARY_RE = re.compile(
    r"mode\s+(?P<mode>-?\d+)\s+"
    r"inform\s+(?P<inform>-?\d+)\s+"
    r"f\s+(?P<f>[-+0-9.Ee]+)\s+"
    r"gpsupn\s+(?P<gpsupn>[-+0-9.Ee]+)\s+"
    r"xsum\s+(?P<xsum>[-+0-9.Ee]+)\s+"
    r"iter\s+(?P<iter>-?\d+)\s+"
    r"fcnt\s+(?P<fcnt>-?\d+)\s+"
    r"gcnt\s+(?P<gcnt>-?\d+)\s+"
    r"cgcnt\s+(?P<cgcnt>-?\d+)"
)


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    probe = repo / "build" / "tests" / "test_gencan_ab_probe"
    fixture = repo / "tests" / "fixtures" / "tiny_packing.inp"
    structure = fixture.parent / "minimal_structure.pdb"

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)
        if structure.exists():
            shutil.copy2(structure, tmpdir / structure.name)

        base_env = os.environ.copy()
        base_env.update(
            {
                "PACKMOL_GENCAN_IMPL": "cpp",
                "PACKMOL_GENCAN_NUMERIC_CPP": "0",
                "PACKMOL_GENCAN_DISABLE_LEGACY_FORTRAN_FALLBACK": "1",
                "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF": "1",
                "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_MAX_DEPTH": "1",
                "PACKMOL_GENCAN_DEBUG": "1",
            }
        )

        safe_env = dict(base_env)
        safe_env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE"] = "0"
        safe_run = subprocess.run(
            [str(probe), str(staged_input)],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            env=safe_env,
            check=False,
        )
        safe_output = (safe_run.stdout + safe_run.stderr).replace("\x00", "")
        if safe_run.returncode != 0:
            raise RuntimeError(f"safe probe failed with code {safe_run.returncode}\n{safe_output}")

        unsafe_env = dict(base_env)
        unsafe_env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE"] = "1"
        unsafe_run = subprocess.run(
            [str(probe), str(staged_input)],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            env=unsafe_env,
            check=False,
        )
        unsafe_output = (unsafe_run.stdout + unsafe_run.stderr).replace("\x00", "")
        if unsafe_run.returncode != 0:
            raise RuntimeError(f"unsafe probe failed with code {unsafe_run.returncode}\n{unsafe_output}")

        if "[gencan-cpp-handoff-unsafe-downgraded]" not in unsafe_output:
            raise RuntimeError("unsafe handoff downgrade marker missing")

        safe_match = SUMMARY_RE.search(safe_output)
        unsafe_match = SUMMARY_RE.search(unsafe_output)
        if safe_match is None:
            raise RuntimeError(f"safe summary line missing\n{safe_output}")
        if unsafe_match is None:
            raise RuntimeError(f"unsafe summary line missing\n{unsafe_output}")

        mode = int(unsafe_match.group("mode"))
        inform = int(unsafe_match.group("inform"))
        iter_count = int(unsafe_match.group("iter"))
        fcnt = int(unsafe_match.group("fcnt"))
        gcnt = int(unsafe_match.group("gcnt"))
        cgcnt = int(unsafe_match.group("cgcnt"))
        if mode != 1:
            raise RuntimeError(f"unexpected mode={mode}")
        if inform < 0 or inform > 8:
            raise RuntimeError(f"unexpected inform={inform}")
        if min(iter_count, fcnt, gcnt, cgcnt) < 0:
            raise RuntimeError(
                f"negative counter(s): iter={iter_count} fcnt={fcnt} gcnt={gcnt} cgcnt={cgcnt}"
            )

        for key in ("inform", "iter", "fcnt", "gcnt", "cgcnt"):
            s = int(safe_match.group(key))
            u = int(unsafe_match.group(key))
            if s != u:
                raise RuntimeError(f"safe/unsafe drift in {key}: safe={s} unsafe={u}")

        for key in ("f", "gpsupn", "xsum"):
            s = float(safe_match.group(key))
            u = float(unsafe_match.group(key))
            if not math.isclose(s, u, rel_tol=1.0e-12, abs_tol=1.0e-12):
                raise RuntimeError(f"safe/unsafe drift in {key}: safe={s:.16e} unsafe={u:.16e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
