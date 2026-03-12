---
name: fortran-baseline-testing
description: Use when extending, repairing, or auditing the Packmol Fortran test baseline under pixi and CTest, especially when adding direct interface coverage or refreshing the procedure inventory
---

# Fortran Baseline Testing

## Overview

This repository treats the Fortran codebase as the reference implementation. The baseline is managed through `pixi`, built with `CMake`/`CTest`, and tracked by an interface inventory.

The rule is simple: lock behavior before changing it.

## When to Use

Use this skill when:

- adding tests for a Fortran procedure
- checking whether an interface is already covered
- extending fixtures or test shims
- fixing a failing baseline test
- refreshing the procedure inventory
- preparing a Fortran interface for later C++ replacement

Do not use this skill as the main guide for C++ implementation planning. Use `fortran-to-cpp-migration` for that.

## Canonical Entry Points

- Environment and task runner: `pixi.toml`
- Test registration: `tests/CMakeLists.txt`
- Inventory generator: `scripts/generate_inventory.py`
- Machine-readable inventory: `tests/inventory/procedures.json`
- Human checklist: `tests/inventory/CHECKLIST.md`
- Test support modules: `tests/support/`
- Small stable fixtures: `tests/fixtures/`
- Legacy regression layer: `testing/`

## Standard Workflow

1. Refresh inventory before claiming coverage changes.
2. Check `tests/inventory/CHECKLIST.md` to see whether the target interface already has a direct test.
3. Choose the smallest test layer that can verify the behavior:
   - `unit`: deterministic helpers and pure calculations
   - `component`: module state, file I/O, or multi-step setup
   - `gencan`: optimizer-internal routines
   - `regression`: existing CLI and end-to-end scripts in `testing/`
4. Reuse an existing fixture first. Only add a new fixture if current ones cannot express the behavior cleanly.
5. Add the test to `tests/CMakeLists.txt`.
6. Run the narrowest relevant test target first.
7. Run `pixi run test-all` before claiming success.

## Test Construction Rules

- Prefer one executable per focused behavior group.
- A test may cover multiple interfaces, but each interface should still have a direct and intentional call path.
- Use Fortran assertions for in-process checks.
- Use Python checker scripts when the assertion depends on:
  - subprocess exit status
  - stdout or stderr content
  - generated files
- Keep fixtures minimal and deterministic.
- Do not weaken assertions to “does not crash” unless the procedure is inherently a smoke-only boundary.

## Fixture Strategy

Prefer these existing fixture styles:

- `minimal_packmol.inp`: parser and I/O component tests
- `tiny_packing.inp`: runtime packing, `initial`, `movebad`, and related stateful flows
- empty toy problem setup: zero-atom or zero-gradient GENCAN and evaluator tests

If a new fixture is required:

- make it as small as possible
- keep it single-purpose
- place it under `tests/fixtures/`

## Shim Strategy

Use `tests/support/` for thin support code only:

- reset or initialize module state
- load minimal runtime problems
- normalize paths and test inputs

Do not move algorithmic behavior into support code.

## Verification

Default verification sequence:

```bash
pixi run inventory
pixi run build
ctest --test-dir build --output-on-failure -R '<new-test-name>'
pixi run test-all
```

If inventory changed, ensure both of these remain consistent:

- `tests/inventory/procedures.json`
- `tests/inventory/CHECKLIST.md`

## Common Mistakes

- Adding a test without refreshing inventory first.
- Treating regression scripts as a substitute for direct interface tests.
- Passing Fortran literals to routines that write back through arguments.
- Building large mixed-purpose harnesses that hide which interface failed.
- Encoding current interface counts into docs instead of regenerating inventory.

## Completion Standard

Work is complete only when:

- the new or changed interface is reflected in inventory
- the direct targeted test passes
- `pixi run test-all` passes

