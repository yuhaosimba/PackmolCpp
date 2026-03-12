# Remaining Fortran Tests Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Close the remaining interface-level Packmol Fortran test gaps so every discovered procedure has at least one direct executable test or subprocess contract test before the C++ migration.

**Architecture:** Keep the current `pixi -> CMake -> CTest` pipeline unchanged and add the missing tests in small batches. Prefer direct Fortran tests for deterministic interfaces, and use probe executables plus Python checkers for procedures that print, write files, or terminate the process.

**Tech Stack:** Pixi, CMake, CTest, Fortran, Python 3

---

### Task 1: Plan Baseline And Fixtures

**Files:**
- Create: `docs/plans/2026-03-12-remaining-fortran-tests.md`
- Modify: `tests/CMakeLists.txt`
- Modify: `tests/support/test_assertions.f90`
- Create: `tests/support/` helper modules only if required by later tasks

**Step 1: Write the failing test**

Add the missing test targets for the next uncovered procedures to `tests/CMakeLists.txt`.

**Step 2: Run test to verify it fails**

Run: `pixi run build`
Expected: configuration or build fails because the new test sources do not exist yet.

**Step 3: Write minimal implementation**

Create the missing test programs or subprocess checkers.

**Step 4: Run test to verify it passes**

Run: `ctest --test-dir build --output-on-failure -R '<new-test-regex>'`
Expected: all new tests pass.

**Step 5: Commit**

```bash
git add docs/plans/2026-03-12-remaining-fortran-tests.md tests
git commit -m "test: add next batch of Fortran baseline coverage"
```

### Task 2: Direct Core Evaluation Wrappers

**Files:**
- Modify: `tests/CMakeLists.txt`
- Create: `tests/component/test_cli_parser_functions_probe.f90`
- Create: `tests/component/check_cli_parser_functions.py`
- Create: `tests/component/test_core_evaluators.f90`
- Create: `tests/gencan/test_gencan_eval_wrappers.f90`

**Step 1: Write the failing test**

Register tests for `parse_command`, `get_filename`, `computef`, `computeg`, `evalal`, and `evalnal`.

**Step 2: Run test to verify it fails**

Run: `pixi run build`
Expected: missing source file errors for the new test targets.

**Step 3: Write minimal implementation**

Use the existing empty toy problem pattern to assert zero objective, zero gradient, and preserved reduced-space vectors. Use subprocess tests for CLI helper routines that rely on command arguments.

**Step 4: Run test to verify it passes**

Run: `ctest --test-dir build --output-on-failure -R 'cli_parser_functions|core_evaluators|gencan_eval_wrappers'`
Expected: all new tests pass.

**Step 5: Commit**

```bash
git add tests/CMakeLists.txt tests/component tests/gencan
git commit -m "test: cover core evaluators and CLI helper routines"
```

### Task 3: Core Geometry And Constraint Kernels

**Files:**
- Modify: `tests/CMakeLists.txt`
- Create: `tests/component/test_constraint_kernels.f90`
- Create: `tests/gencan/test_gencan_calchd.f90`
- Create: `tests/fixtures/` additions only if a geometry fixture is required

**Step 1: Write the failing test**

Register tests for `comprest`, `gparc`, `gwalls`, and `calchd`.

**Step 2: Run test to verify it fails**

Run: `pixi run build`
Expected: missing source file errors for the new test targets.

**Step 3: Write minimal implementation**

Build a tiny deterministic geometry fixture with one or two atoms and assert non-negative penalties, expected zero cases, and stable shrink/expand behavior around `calchd`.

**Step 4: Run test to verify it passes**

Run: `ctest --test-dir build --output-on-failure -R 'constraint_kernels|gencan_calchd'`
Expected: all new tests pass.

**Step 5: Commit**

```bash
git add tests/CMakeLists.txt tests/component tests/gencan tests/fixtures
git commit -m "test: cover constraint kernels"
```

### Task 4: Input And Output Mainline Paths

**Files:**
- Modify: `tests/CMakeLists.txt`
- Create: `tests/component/test_getinp_probe.f90`
- Create: `tests/component/check_getinp.py`
- Create: `tests/component/test_output_probe.f90`
- Create: `tests/component/check_output_files.py`
- Create: `tests/component/test_checkpoint_probe.f90`
- Create: `tests/component/check_checkpoint.py`

**Step 1: Write the failing test**

Register tests for `getinp`, `output`, `write_connect`, and `checkpoint`.

**Step 2: Run test to verify it fails**

Run: `pixi run build`
Expected: missing source file errors for the new test targets.

**Step 3: Write minimal implementation**

Use subprocess probes that execute in temporary directories and assert file creation, key text fragments, and non-destructive repeatability.

**Step 4: Run test to verify it passes**

Run: `ctest --test-dir build --output-on-failure -R 'getinp|output|checkpoint'`
Expected: all new tests pass.

**Step 5: Commit**

```bash
git add tests/CMakeLists.txt tests/component
git commit -m "test: cover input and output workflow procedures"
```

### Task 5: Heuristics And State Mutation Helpers

**Files:**
- Modify: `tests/CMakeLists.txt`
- Create: `tests/component/test_heuristics.f90`

**Step 1: Write the failing test**

Register tests for `initial`, `movebad`, `restmol`, and `swaptype`.

**Step 2: Run test to verify it fails**

Run: `pixi run build`
Expected: missing source file errors for the new test targets.

**Step 3: Write minimal implementation**

Construct the smallest mutable state fixture that lets these routines change positions or bookkeeping without requiring a full production input deck.

**Step 4: Run test to verify it passes**

Run: `ctest --test-dir build --output-on-failure -R 'heuristics'`
Expected: all new tests pass.

**Step 5: Commit**

```bash
git add tests/CMakeLists.txt tests/component
git commit -m "test: cover heuristic mutation routines"
```

### Task 6: Full Optimizer Driver Path

**Files:**
- Modify: `tests/CMakeLists.txt`
- Create: `tests/gencan/test_gencan_driver_probe.f90`
- Create: `tests/gencan/check_gencan_driver.py`
- Create: `tests/component/test_pgencan.f90`

**Step 1: Write the failing test**

Register tests for `pgencan`, `easygencan`, `gencan`, `spgls`, `cg`, and `tnls`.

**Step 2: Run test to verify it fails**

Run: `pixi run build`
Expected: missing source file errors for the new test targets.

**Step 3: Write minimal implementation**

Use a toy quadratic or empty feasible problem to assert clean exit codes, non-increasing objective behavior, and stable `inform/fcnt/gcnt/cgcnt` bookkeeping.

**Step 4: Run test to verify it passes**

Run: `ctest --test-dir build --output-on-failure -R 'gencan_driver|pgencan'`
Expected: all new tests pass.

**Step 5: Commit**

```bash
git add tests/CMakeLists.txt tests/gencan tests/component
git commit -m "test: cover optimizer driver routines"
```

### Task 7: Final Verification And Inventory Closure

**Files:**
- Modify: `tests/inventory/procedures.json`
- Modify: `tests/inventory/CHECKLIST.md`
- Modify: `README.md` only if the new test commands need documenting

**Step 1: Write the failing test**

Regenerate the inventory after the new tests exist and compare the remaining uncovered procedures against zero.

**Step 2: Run test to verify it fails**

Run: `python scripts/generate_inventory.py --root . --check`
Expected: fail only if the inventory file is stale.

**Step 3: Write minimal implementation**

Refresh the generated inventory assets and update any user-facing test documentation if needed.

**Step 4: Run test to verify it passes**

Run: `pixi run test-all`
Expected: `100% tests passed` and no uncovered procedures remain in the plan tracking.

**Step 5: Commit**

```bash
git add tests/inventory README.md
git commit -m "test: close remaining Fortran baseline gaps"
```
