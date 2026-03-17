#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile


MODE_ID = {
    "fortran": 0,
    "cpp": 1,
    "ab": 2,
}

ALLOWED_FALLBACK_REASONS = {
    "tn_post_nonterminal",
}


def _run_probe(
    probe: pathlib.Path,
    input_path: pathlib.Path,
    mode: str,
    workdir: pathlib.Path,
) -> str:
    env = os.environ.copy()
    env["PACKMOL_GENCAN_IMPL"] = mode
    env["PACKMOL_GENCAN_DEBUG"] = "1"
    completed = subprocess.run(
        [str(probe), str(input_path)],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=env,
    )
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        print(f"probe failed in mode={mode}", file=sys.stderr)
        print(output, file=sys.stderr)
        raise RuntimeError("probe run failed")
    mode_marker = f"mode {MODE_ID[mode]} "
    if mode_marker not in output:
        print(output, file=sys.stderr)
        raise RuntimeError(f"unexpected mode marker for mode={mode}")
    return output


def main() -> int:
    probe = pathlib.Path(sys.argv[1])
    fixture = pathlib.Path(sys.argv[2])

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)
        structure_src = fixture.parent / "minimal_structure.pdb"
        if structure_src.exists():
            shutil.copy2(structure_src, tmpdir / structure_src.name)

        # Baseline still runs for parity workflow consistency.
        _run_probe(probe, staged_input, "fortran", tmpdir)
        cpp_output = _run_probe(probe, staged_input, "cpp", tmpdir)
        ab_output = _run_probe(probe, staged_input, "ab", tmpdir)

        pattern = re.compile(r"gencan-cpp-fallback\] reason=([a-zA-Z0-9_]+)")
        cpp_reasons = {match.group(1) for match in pattern.finditer(cpp_output)}
        ab_reasons = {match.group(1) for match in pattern.finditer(ab_output)}

        unexpected_cpp = sorted(reason for reason in cpp_reasons if reason not in ALLOWED_FALLBACK_REASONS)
        unexpected_ab = sorted(reason for reason in ab_reasons if reason not in ALLOWED_FALLBACK_REASONS)
        if unexpected_cpp:
            raise RuntimeError(f"unexpected fallback reason(s) in cpp mode: {unexpected_cpp}")
        if unexpected_ab:
            raise RuntimeError(f"unexpected fallback reason(s) in ab mode: {unexpected_ab}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
