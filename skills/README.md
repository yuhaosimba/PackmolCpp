# Repository Skills

This directory stores repository-local workflow skills.

- `fortran-baseline-testing`: maintain and extend the current Fortran baseline with `pixi`, `CMake`, `CTest`, fixtures, and inventory.
- `fortran-to-cpp-migration`: replace Fortran interfaces with C++ incrementally while keeping the baseline tests as the release gate.

These skills summarize the process already encoded in:

- `pixi.toml`
- `tests/CMakeLists.txt`
- `tests/inventory/`
- `tests/support/`
- `testing/`
