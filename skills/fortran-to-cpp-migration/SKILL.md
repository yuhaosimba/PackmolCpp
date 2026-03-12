---
name: fortran-to-cpp-migration
description: Use when replacing Packmol Fortran procedures with C++ implementations, especially during mixed-language migration where interface-level tests must remain the release gate
---

# Fortran to C++ Migration

## Overview

Packmol migration should proceed interface by interface, not file by file and not by full rewrite. The Fortran implementation remains the behavioral reference until a C++ replacement passes the same baseline tests.

The rule is: no replacement without a stable baseline.

## When to Use

Use this skill when:

- choosing the order of C++ migration
- replacing one Fortran routine with a C++ implementation
- introducing mixed-language builds
- adding A/B checks between Fortran and C++ paths
- deciding whether to pause migration and add missing baseline tests first

Use `fortran-baseline-testing` when the task is mainly about baseline test creation or repair.

## Migration Order

Preferred order:

1. Pure helpers and deterministic utilities
2. Small component interfaces with controlled module state
3. Core evaluators and constraint kernels
4. Driver wrappers
5. High-coupling optimizer internals only after lower layers are stable

Do not start with whole-file rewrites of `gencan.f`, `getinp.f90`, or large runtime flows unless the relevant interface tests already exist.

## Required Inputs

Before replacing an interface, check:

- `tests/inventory/CHECKLIST.md`
- `tests/inventory/procedures.json`
- the existing direct test in `tests/unit/`, `tests/component/`, or `tests/gencan/`

If the target interface has no direct test, stop and add baseline coverage first.

## Replacement Workflow

1. Identify the exact interface boundary to replace.
2. Run the existing direct test for that interface.
3. Keep the public calling convention stable while introducing the C++ implementation.
4. Link the replacement in mixed-language mode.
5. Re-run the targeted direct test.
6. Re-run nearby component or regression tests that exercise the same path.
7. Run `pixi run test-all`.

## Test Discipline During Migration

- Reuse the existing fixtures whenever possible.
- Prefer behavior assertions over internal layout assertions.
- If both Fortran and C++ implementations are available, compare them through the same test inputs.
- Add new migration-only smoke tests only if existing direct tests cannot expose the language boundary safely.

## Mixed-Language Guidance

- Replace one procedure or one small group of tightly related procedures at a time.
- Avoid simultaneous rewrites across parser, evaluator, and optimizer layers.
- Keep ownership of file generation and process exits stable until the new layer is proven.
- Preserve current command entry points so `testing/` remains useful as an acceptance layer.

## High-Risk Areas

Treat these as last-stage migration targets:

- `gencan` main drivers and line-search internals
- parser flows with `stop`-based failure handling
- heavy module-state procedures that mutate many global arrays

For these areas, prefer:

- direct interface tests
- minimal harnesses
- explicit A/B comparison against Fortran behavior

## Verification

Minimum required verification for each replacement:

```bash
pixi run build
ctest --test-dir build --output-on-failure -R '<targeted-interface-test>'
pixi run test-all
```

If the change crosses parsing, output, or CLI behavior, also check the relevant regression target under `testing/`.

## Common Mistakes

- Replacing code before direct interface coverage exists.
- Migrating a whole source file because the current file structure is inconvenient.
- Using only end-to-end examples to judge equivalence.
- Changing interface shape and implementation at the same time.
- Ignoring Fortran-specific write-back arguments in mixed-language calls.

## Completion Standard

A migration step is complete only when:

- the C++ replacement passes the direct interface test
- related component or regression tests still pass
- `pixi run test-all` passes in the mixed-language build

