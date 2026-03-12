# Packmol Workflow Skills Design

## Summary

This repository keeps two local workflow skills to guide future testing and migration work:

1. `fortran-baseline-testing`
2. `fortran-to-cpp-migration`

The split is intentional. The testing skill defines how to maintain and extend the current Fortran baseline under `pixi` and `CTest`. The migration skill defines how to use that baseline during mixed-language migration so that C++ replacements are introduced with interface-level verification instead of full rewrites.

## Goals

- Keep the current `pixi`-managed test workflow discoverable inside the repository.
- Preserve the discipline that made the current baseline useful:
  - inventory first
  - test-first for new behavior
  - direct interface coverage
  - regression retention
- Give future work a stable migration path from Fortran to C++ without re-deriving the process from conversation history.

## Skill Layout

The skills live under `skills/`:

- `skills/fortran-baseline-testing/SKILL.md`
- `skills/fortran-to-cpp-migration/SKILL.md`

They are repository-local process references. They do not replace the code or test assets in `tests/`, `testing/`, `scripts/`, or `pixi.toml`; they explain how to use and extend them.

## Skill Boundaries

### fortran-baseline-testing

This skill covers:

- when to refresh inventory
- how to classify new interfaces
- how to add unit, component, optimizer-internal, and regression tests
- when to use direct executable tests versus Python subprocess checkers
- how to use `pixi run ...` tasks as the canonical workflow
- what must be verified before claiming new baseline coverage

This skill points to:

- `pixi.toml`
- `tests/CMakeLists.txt`
- `tests/support/`
- `tests/fixtures/`
- `tests/inventory/`
- `scripts/generate_inventory.py`

### fortran-to-cpp-migration

This skill covers:

- migration order
- how to replace one interface at a time
- how to keep mixed-language linking testable
- when to reuse current Fortran fixtures versus introducing dual-implementation tests
- when to stop and add baseline tests before replacing code

This skill points to:

- `tests/inventory/CHECKLIST.md`
- `tests/inventory/procedures.json`
- `tests/component/`
- `tests/gencan/`
- `testing/`

## Design Principles

- Prefer stable commands over narrative history.
- Prefer process constraints over informal advice.
- Avoid hardcoding counts that may drift; tell future work to regenerate inventory instead.
- Keep migration guidance interface-first, not file-first.
- Treat `pixi run test-all` as the default final verification step.

## Acceptance Criteria

- Both skills exist in the repository.
- Each skill has a clear trigger description and concrete workflow.
- The skills reference the actual repository entry points now in use.
- The skills are concise enough to scan during implementation work.

