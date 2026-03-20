#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import tempfile


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    probe = repo / "build" / "tests" / "test_gencan_ab_probe"
    fixture = repo / "tests" / "fixtures" / "tiny_packing.inp"

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)
        shutil.copy2(fixture.parent / "minimal_structure.pdb", tmpdir / "minimal_structure.pdb")

        env = os.environ.copy()
        env.update({
            "PACKMOL_GENCAN_IMPL": "cpp",
            "PACKMOL_GENCAN_NUMERIC_CPP": "0",
            "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE": "1",
            "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY": "0",
            "PACKMOL_GENCAN_DEBUG": "1",
        })
        completed = subprocess.run(
            [str(probe), str(staged_input)],
            cwd=tmpdir,
            env=env,
            text=True,
            capture_output=True,
            check=False,
        )
        output = (completed.stdout + completed.stderr).replace("\x00", "")
        if completed.returncode != 0:
            raise RuntimeError(output)
        if "[gencan-cpp-handoff-canonicalize-cpp-replay-forced]" not in output:
            raise RuntimeError("missing forced cpp replay marker")
        if "[gencan-cpp-handoff-canonicalize-cpp-replay]" not in output:
            raise RuntimeError("missing cpp replay canonicalize marker")
        if "[gencan-cpp-handoff-canonicalize-fortran]" in output:
            raise RuntimeError("unexpected fortran canonicalize marker in forced replay mode")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
