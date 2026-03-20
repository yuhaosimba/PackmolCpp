#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import sys
import tempfile


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_spg_post_handoff.py <probe> <fixture>")

    probe = pathlib.Path(sys.argv[1]).resolve()
    fixture = pathlib.Path(sys.argv[2]).resolve()

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)

        structure_src = fixture.parent / "minimal_structure.pdb"
        if structure_src.exists():
            shutil.copy2(structure_src, tmpdir / structure_src.name)

        env = os.environ.copy()
        env.update({
            "PACKMOL_GENCAN_IMPL": "cpp",
            "PACKMOL_GENCAN_NUMERIC_CPP": "0",
            "PACKMOL_GENCAN_DEBUG": "1",
            "PACKMOL_GENCAN_SPG_POST_CPP_HANDOFF": "1",
            "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE": "1",
            "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY": "1",
            "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_MAX_DEPTH": "1",
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
        if "[gencan-cpp-handoff] reason=spg_post_nonterminal" not in output:
            raise RuntimeError("missing spg_post_nonterminal handoff marker")
        if "[gencan-cpp-handoff-cpp] reason=spg_post_nonterminal" not in output:
            raise RuntimeError("missing spg_post_nonterminal cpp handoff marker")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
