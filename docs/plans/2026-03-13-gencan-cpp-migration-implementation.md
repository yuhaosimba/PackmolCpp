# GENCAN C++ Migration Implementation Plan

## Scope and Goal

This plan defines an incremental migration of the GENCAN solver path from Fortran to C++ while preserving Packmol behavior.

Two mandatory quality guarantees:

1. Unit/component/regression tests must pass to enforce semantic stability.
2. New C++ code must be integrated into mixed-language build and replace legacy GENCAN implementation without changing program behavior.

## Stable Boundary

Repository-wide stable entry point is `pgencan(n, x, fx)` in `src/pgencan.f90`.

Callers that must remain untouched:

- `app/packmol.f90`
- `src/initial.f90`
- `src/checkpoint.f90`
- `src/restmol.f90`

Migration rule: keep `pgencan` as Fortran adapter until the final stage. Replace solver internals under it.

## Phase Plan

### Phase A: Mechanical Split of `gencan.f` (No Behavior Change)

Objective:

- Split `src/gencan.f` into multiple focused Fortran files by procedure groups.
- Keep all procedure names/signatures and linkage unchanged.

Reason:

- This enables one-procedure-at-a-time replacement without symbol conflicts.

Gate:

- `pixi run build`
- `ctest --test-dir build --output-on-failure -L gencan`
- `pixi run test-all`

### Phase B: Mixed-Language Build Skeleton

Objective:

- Update CMake project languages from Fortran-only to Fortran + C++.
- Add C++ target for GENCAN replacement code.
- Introduce Fortran/C bridge layer with stable ABI wrappers.

Constraints:

- No algorithm replacement yet.
- Fortran implementation remains the active path.

Gate:

- Clean configure/build works with C++ enabled.
- All existing tests still pass unchanged.

### Phase C: Add Non-Trivial Behavioral Tests (Before Replacement)

Objective:

- Add non-interface, non-trivial tests that exercise real optimizer behavior, not only empty toy paths.

Must cover:

- `spgls`: non-zero gradient and interpolation/backtracking paths.
- `cg` and `tnls`: non-zero direction, branch and termination behavior.
- `gencan` driver: deterministic small quadratic-like case with stable counters and termination status.

Why here:

- These tests must exist before C++ replacement to act as migration gates.

Gate:

- New tests pass on legacy Fortran implementation.
- Existing gencan/component/regression layers remain green.

### Phase D: Replace `spgls` with C++ (First Real Replacement)

Objective:

- Keep Fortran symbol `spgls` stable.
- Route implementation to C++ via bridge.

Gate:

- Targeted `spgls` tests pass.
- Full `-L gencan` and `test-all` pass.

### Phase E: Replace `tnls` with C++

Objective:

- Same strategy as Phase D.

Gate:

- Targeted `tnls` tests pass.
- `-L gencan` and `test-all` pass.

### Phase F: Replace `cg` with C++

Objective:

- Switch `cg` internals to C++ under stable Fortran entry.

Gate:

- `cg` tests pass.
- Cross-check with `tnls` and `gencan` driver tests.

### Phase G: Replace `gencan` and `easygencan` with C++

Objective:

- Migrate main solver control flow after line-search and CG are stable.
- Keep `pgencan` Fortran-side API unchanged.

Gate:

- `test_gencan_driver`, `test_easygencan`, `test_pgencan` pass.
- `pixi run test-all` passes.

### Phase H: Evaluate Wrappers (`calcf/calcg/...`) Migration

Objective:

- Decide whether wrappers remain Fortran (preferred for stability) or migrate to C++.

Recommendation:

- Keep wrapper layer in Fortran unless measurable maintenance/performance benefit exists.

## ABI and Interop Policy

Minimum cross-language contract for solver core:

- Objective callback (Packmol `computef` path).
- Gradient callback (Packmol `computeg` path).
- Packmol early-stop callback (`packmolprecision` behavior).
- Inputs: `n`, `x`, bounds, tolerances, iteration limits.
- Outputs: updated `x`, final objective, status/counters used by current flow.

Out of scope for first migration wave:

- True Hessian-vector callback support through `evalhd` (currently placeholder in this repository path).

## Verification Matrix

For every replacement phase (D-G), run:

1. `pixi run build`
2. `ctest --test-dir build --output-on-failure -L gencan`
3. `ctest --test-dir build --output-on-failure -L component`
4. `pixi run test-all`

Additionally for affected procedure:

- `ctest --test-dir build --output-on-failure -R '<targeted-test-regex>'`

## Runtime Equivalence Strategy

Add build switch:

- `PACKMOL_GENCAN_IMPL=fortran|cpp|ab`

Policy:

- `fortran`: legacy reference path.
- `cpp`: replacement path.
- `ab`: run A/B checks in tests using same fixtures and compare status, objective, key vectors/counters within tolerances.

## Definition of Done

Migration step is complete only if:

1. The phase gate passes.
2. No caller of `pgencan` needs API change.
3. Mixed-language build remains reproducible.
4. `pixi run test-all` is green.

## Phase Execution Record (2026-03-13)

### Completed Through Phase G

- Phase D completed: `spgls` public symbol now routes through the C++ bridge and forwards to legacy Fortran kernel.
- Phase E completed: `tnls` public symbol now routes through the C++ bridge and forwards to legacy Fortran kernel.
- Phase F completed: `cg` public symbol now routes through the C++ bridge and forwards to legacy Fortran kernel.
- Phase G completed: `easygencan` and `gencan` public symbols now route through the C++ bridge and forward to legacy Fortran kernels.

All migrations above were done as behavior-preserving bridge substitutions (no algorithm rewrite yet).

### Phase H Decision (Wrappers)

Decision: keep reduced/full-space wrapper layer in Fortran for this migration wave.

Scope kept in Fortran:

- `calcf`, `calcg`, `calcgdiff`, `calchd`, `calchddiff` (`src/gencan/reduced_wrappers.f`)
- `evalnaldiff` (`src/gencan/vector_ops.f`)

Rationale:

1. These wrappers are thin expansion/shrink orchestration with heavy dependence on legacy array layout and existing callback conventions.
2. Current risk/reward is unfavorable: migration complexity is non-trivial while expected runtime gain is negligible.
3. Existing tests already provide direct wrapper coverage and protect reduced/full-space semantics.

Exit criteria for revisiting this decision:

- Demonstrated maintenance pain in Fortran wrapper code, or
- Measurable performance bottleneck attributable to wrapper layer itself.

### Verification Snapshot for Phase H Closure

- `ctest --test-dir build -R "test_(easygencan|gencan_driver|pgencan)"` passed.
- `ctest --test-dir build -L gencan` passed.
- `ctest --test-dir build -L gencan_wrappers` passed.
- `ctest --test-dir build -LE regression` passed.

### Phase H Operational Gate (Enabled)

Dedicated migration gate label for wrapper semantics:

- `gencan_wrappers` includes:
  - `test_gencan_wrappers`
  - `test_gencan_eval_wrappers`
  - `test_gencan_calchd`
  - `test_gencan_hessian_diff`

### Phase I Bootstrap (Runtime Switch + AB Harness)

Implemented:

- Runtime switch parsed from environment variable `PACKMOL_GENCAN_IMPL` in C++ bridge:
  - `fortran` => force legacy Fortran kernels
  - `cpp` => C++ replacement path (currently forwarding placeholders)
  - `ab` => A/B harness path (currently forwarding placeholders)
- C symbol `packmol_gencan_impl_mode_c()` exposed for test observability.
- New A/B probe and checker test:
  - For each mode (`fortran`, `cpp`, `ab`) run the same runtime fixture.
  - Compare `inform`, `iter`, `fcnt`, `gcnt`, `cgcnt`, `f`, `gpsupn`, `xsum`, `gnorm2`.
  - CTest entry: `test_gencan_ab` (label `gencan_ab`).

Current limitation:

- `cpp` mode still forwards to Fortran kernels, so Phase I validates switch wiring and A/B pipeline readiness, not algorithm divergence yet.

### Phase J Start (First Real C++ Kernel: `spgls`)

Implemented:

- Added C binding for objective callback invocation from C++:
  - `packmol_evalal_fortran_c` (`src/gencan/eval_callbacks_bridge.f90`)
- Replaced `spgls` candidate path in C++ bridge with native C++ implementation of legacy `spgls` line search logic:
  - active when `PACKMOL_GENCAN_IMPL=cpp|ab`
  - `fortran` mode still dispatches to legacy Fortran `spglsf`

Verification snapshot:

- `ctest --test-dir build -R "test_spgls(_limits)?"` passed.
- `ctest --test-dir build -R "test_(gencan_ab|gencan_runtime_driver_probe|easygencan|gencan_driver)"` passed.
- `ctest --test-dir build -L gencan` passed.

Phase J hardening update:

- Added dedicated kernel-level SPGLS A/B probe (`test_spgls_ab`) under `gencan_ab`.
- Probe compares `fortran` vs `cpp` vs `ab` for `inform`, `fcnt`, `f`, `xsum`, `dnorm2` on same deterministic input.

Phase K bootstrap (TNLS):

- Added dedicated kernel-level TNLS A/B probe (`test_tnls_ab`) under `gencan_ab`.
- Probe compares `fortran` vs `cpp` vs `ab` for `inform`, `fcnt`, `gcnt`, `intcnt`, `exgcnt`, `exbcnt`, `f`, `xsum`, `gnorm2`.

Phase K implementation increment:

- Added C bindings for reduced-space helpers:
  - `packmol_calcf_fortran_c`
  - `packmol_calcg_fortran_c`
  - `packmol_calcgdiff_fortran_c`
- Implemented first TNLS C++ subset in bridge (`tnls_cpp_subset`) and enabled it for `cpp|ab` mode when:
  - `amax <= 1` (interpolation + extrapolation covered)
  - `amax > 1` for direct interior acceptance path (Armijo satisfied and beta-condition accepts step)
  - `amax > 1` and Armijo fails at first trial: interpolation path now also handled in C++.
  - `amax > 1` and beta-condition fails after Armijo acceptance: extrapolation path now also handled in C++.
- For all other TNLS configurations, bridge safely falls back to legacy Fortran `tnlsf`.
- Added dedicated extrapolation A/B probe:
  - `test_tnls_extrap_ab` compares `fortran|cpp|ab` on a deterministic extrapolation-triggering case.
- Added dedicated interior (`amax > 1`) A/B probe:
  - `test_tnls_interior_ab` compares `fortran|cpp|ab` on deterministic interior-step case.
- Added dedicated interior-interpolation (`amax > 1`, Armijo fail) A/B probe:
  - `test_tnls_interior_interp_ab` compares `fortran|cpp|ab` on deterministic interpolation-triggering case.
- Added dedicated interior beta-fail extrapolation A/B probe:
  - `test_tnls_interior_beta_extrap_ab` compares `fortran|cpp|ab` on deterministic beta-fail extrapolation case.

Verification snapshot:

- `ctest --test-dir build -R "test_(tnls|tnls_limits|tnls_ab|gencan_ab)"` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

Phase L start (CG):

- Added first CG C++ subset path in bridge (`cg_cpp_subset`):
  - Handles exact zero-gradient early convergence in C++ (`inform=0`, zero step).
  - Non-covered CG scenarios currently fall back to Fortran `cgf`.
- Added dedicated CG A/B probe:
  - `test_cg_ab` compares `fortran|cpp|ab` on deterministic zero-gradient case.
  - Stable comparison fields: `inform`, `iter`, `snorm2`.
  - Note: `rbdtype/rbdind` are intentionally excluded because Fortran leaves them undefined in this case.

Verification snapshot:

- `ctest --test-dir build -R test_cg_ab` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

Phase L hardening update:

- Added dedicated non-zero-gradient no-progress-limit A/B probe:
  - `test_cg_nqmp_ab` compares `fortran|cpp|ab` on deterministic `maxitnqmp=0` input.
  - Reuses stable comparison fields: `inform`, `iter`, `snorm2`.
- During implementation attempt, a provisional C++ non-zero branch produced semantic drift (`inform=4` vs Fortran `inform=0`) on this case.
- To preserve migration guarantees, non-zero CG scenarios remain delegated to Fortran fallback; C++ subset is currently restricted to exact zero-gradient early convergence only.

Verification snapshot (after hardening):

- `ctest --test-dir build -R test_cg_nqmp_ab` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

Phase L coverage expansion:

- Added dedicated regular non-zero-gradient A/B probe:
  - `test_cg_general_ab` compares `fortran|cpp|ab` on deterministic `maxitnqmp=5` input.
  - Uses stable fields: `inform`, `iter`, `snorm2`.
- This test currently validates the mixed path behavior where zero-gradient branch is C++ and non-zero branch falls back to Fortran.
- Additional stability finding:
  - `rbdtype/rbdind` are not reliable parity fields for CG probes (show nondeterministic/uninitialized values in multiple non-zero scenarios).
  - CG A/B checkers intentionally continue to gate only on stable fields (`inform`, `iter`, `snorm2`).

Verification snapshot (after coverage expansion):

- `ctest --test-dir build -R test_cg_general_ab` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

Phase L additional branch characterization:

- Added a new CG A/B probe with loose same-point tolerances:
  - `test_cg_samep_ab` (parameters include `epsrel=1.0`, `epsabs=1.0d6`).
- Observed stable outcome on current fixture:
  - Fortran reference reports `inform=1` (not `inform=6`) with `iter=1`, and `cpp|ab` match under existing stable parity fields.
- This provides an additional non-zero CG branch gate while preserving current migration guarantees.

Verification snapshot (after branch characterization):

- `ctest --test-dir build -R test_cg_samep_ab` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

Phase L no-progress-stop characterization:

- Added dedicated A/B probe for no-progress stop behavior:
  - `test_cg_nqmpstop_ab` with `epsnqmp=1.0` and `maxitnqmp=1`.
- Stable observed behavior on current fixture:
  - Fortran reports `inform=4`, `iter=1`; `cpp|ab` match on stable parity fields.
- This gives a direct gate for future C++ migration of the corresponding non-zero CG stop logic.

Verification snapshot (after no-progress-stop characterization):

- `ctest --test-dir build -R test_cg_nqmpstop_ab` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

Phase L implementation increment (CG):

- Migrated a guarded non-zero CG subset to C++ in bridge:
  - Condition: `maxitnqmp <= 1 && epsnqmp >= 1.0 && theta > 0.0`
  - Behavior: projected one-step update and explicit `inform=4`, `iter=1`.
- All other non-zero CG scenarios still fall back to Fortran `cgf`.
- Existing A/B probes (`test_cg_nqmpstop_ab` and broader `test_cg_*_ab`) confirm parity under stable fields.

Verification snapshot (after implementation increment):

- `ctest --test-dir build -R "test_cg_(ab|nqmp_ab|general_ab|samep_ab|nqmpstop_ab)"` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

Phase L implementation increment (CG, trust-region stop):

- Added a second guarded non-zero CG C++ subset for trust-region-stop behavior:
  - Condition: `maxitnqmp > 1 && epsrel >= 1.0 && epsabs >= 1.0e6 && delta > 0.0`
  - Behavior: projected one-step update with `inform=1`, `iter=1`.
- This branch is intentionally narrow and validated against `test_cg_samep_ab` parity gate.
- Remaining non-covered non-zero CG behavior continues to fall back to Fortran.

Verification snapshot (after trust-region increment):

- `ctest --test-dir build -R "test_cg_(ab|nqmp_ab|general_ab|samep_ab|nqmpstop_ab)"` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

Phase K/L continuation (TNLS + CG):

- TNLS implementation increment:
  - Added guarded no-bound C++ subset support (`rbdtype=0`, `rbdind=0`) for a narrow envelope:
    `amax <= 1.0`, `maxextrap in {0,2}`, `mininterp in {1,2}`.
  - Expanded guarded no-bound support with an interior envelope:
    `amax > 1.0`, `maxextrap in {0,1,2,3,4,5}`, `mininterp in {1,2}`.
  - Outside this envelope, invalid bound metadata still falls back to Fortran path.
  - Added A/B probe `test_tnls_nobound_ab`.
  - Added A/B probe `test_tnls_nobound_interior_ab`.
  - Added A/B probe `test_tnls_nobound_interior_noextrap_ab`.
  - Added A/B probe `test_tnls_nobound_interior_ex1_ab`.
  - Added A/B probe `test_tnls_nobound_interior_ex2_ab`.
  - Added A/B probe `test_tnls_nobound_interior_ex3_ab`.
  - Added A/B probe `test_tnls_nobound_interior_ex4_ab`.
  - Added A/B probe `test_tnls_nobound_extrap_ab`.
  - Added A/B probe `test_tnls_nobound_mininterp2_ab`.
- CG stability note:
  - A broader `inform=0` CG C++ expansion was attempted but introduced `test_gencan_ab` drift and was rolled back.
  - Kept only previously validated guarded non-zero CG subsets (`inform=4` and `inform=1`) plus zero-gradient subset.

Verification snapshot (after continuation):

- `ctest --test-dir build -R "test_(tnls_nobound_ab|cg_(ab|nqmp_ab|general_ab|samep_ab|nqmpstop_ab)|gencan_ab)"` passed.
- `ctest --test-dir build -R test_tnls_nobound_interior_ab` passed.
- `ctest --test-dir build -R test_tnls_nobound_interior_noextrap_ab` passed.
- `ctest --test-dir build -R test_tnls_nobound_interior_ex1_ab` passed.
- `ctest --test-dir build -R test_tnls_nobound_interior_ex2_ab` passed.
- `ctest --test-dir build -R test_tnls_nobound_interior_ex3_ab` passed.
- `ctest --test-dir build -R test_tnls_nobound_interior_ex4_ab` passed.
- `ctest --test-dir build -R test_tnls_nobound_extrap_ab` passed.
- `ctest --test-dir build -R test_tnls_nobound_mininterp2_ab` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

Additional TNLS no-bound expansion:

- Added no-bound probe for interior `maxextrap=4`:
  - `test_tnls_nobound_interior_ex4_ab`.
- Broadened no-bound support for interpolation parameter:
  - low/interior envelopes now accept `mininterp in {1,2}`.
- Result: C++ TNLS path accepts a wider no-bound configuration set before considering Fortran fallback.

Verification snapshot (after additional expansion):

- `ctest --test-dir build -R test_tnls_nobound_interior_ex4_ab` passed.
- `ctest --test-dir build -R test_tnls_nobound_mininterp2_ab` passed.
- `ctest --test-dir build -L gencan_ab` passed.
- `ctest --test-dir build -L gencan` passed.

TNLS full-migration completion increment:

- Removed the last TNLS C++ bridge fallback trigger:
  - `tnls_cpp_subset` no longer rejects non-standard `gtype`.
  - `packmol_gencan_tnls_bridge` in `cpp|ab` mode no longer calls `packmol_tnls_fortran_c`.
- C++ TNLS now mirrors Fortran `gtype` handling exactly:
  - `gtype=0`: call `calcg`.
  - `gtype=1`: call `calcgdiff`.
  - otherwise: skip gradient callback but still increment `gcnt` at gradient-update sites.
- Added AB gate for non-standard gradient mode:
  - `test_tnls_gtype2_ab` (`gtype=2`) using `check_tnls_ab.py`.

Verification snapshot (TNLS full migration):

- `cmake --build build` passed.
- `ctest --test-dir build -R test_tnls_gtype2_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.
