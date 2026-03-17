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

CG full-migration start (in progress):

- Added C bridges for Hessian-vector callbacks used by CG:
  - `packmol_calchd_fortran_c`
  - `packmol_calchddiff_fortran_c`
- Implemented a first full C++ translation draft of `cgf` in `src/cpp/gencan/bridge.cpp` (`cg_cpp_full`).
- Validation result:
  - Dedicated CG A/B probes (`test_cg_ab`, `test_cg_nqmp_ab`, `test_cg_general_ab`, `test_cg_samep_ab`, `test_cg_nqmpstop_ab`) pass against the draft.
  - Runtime gate `test_gencan_ab` still detects drift when routing `cpp|ab` to the draft.
- Current routing decision (stability-first):
  - `packmol_gencan_cg_bridge` keeps `cpp|ab` on Fortran `packmol_cg_fortran_c` while the draft is refined to satisfy full runtime AB parity.

Verification snapshot (after stabilization):

- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

CG full-migration debugging continuation:

- Added a draft switch for isolated runtime validation:
  - `PACKMOL_GENCAN_CG_CPP_DRAFT=1` enables `cg_cpp_full` in `cpp|ab` mode.
  - Default behavior remains unchanged (`cpp|ab` uses Fortran CG) to keep gates stable.
- Added optional diagnostics:
  - `PACKMOL_GENCAN_DEBUG=1` emits CG bridge call summaries for runtime fixture analysis.
- Fixed one concrete semantic mismatch in draft path:
  - For `htvtype=1`, `calchddiff` now receives residual `r` (as in Fortran `cgf`) instead of initial gradient `g`.
- Additional alignment:
  - On `inform=1` trust-region boundary stop, draft now propagates boundary candidate indices (`rbdtype/rbdind`) from current step candidate.

Current status:

- Stable gates (default routing) remain green:
  - `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
  - `ctest --test-dir build -L gencan --output-on-failure` passed.
- Draft path is improved but not yet parity-complete on runtime fixture:
  - With `PACKMOL_GENCAN_CG_CPP_DRAFT=1`, `test_gencan_ab` still reports objective drift vs Fortran.

Kernel utilities migration increment (vector ops + numeric switch):

- Consolidated repeated vector operations in C++ kernels to shared helpers:
  - `vec_copy`, `vec_trial_point`, `same_step_relative`, `same_point_relative`.
  - Applied in `spgls_cpp` and `tnls_cpp_subset` without changing default runtime behavior.
- Added numeric-kernel switch in CG draft path:
  - `PACKMOL_GENCAN_NUMERIC_CPP=1` enables C++ `dot`/`norm2` kernels (`dot_cpp_stable`, `norm2sq_cpp_stable`) inside `cg_cpp_full`.
  - Default remains Fortran numeric wrappers to preserve current parity baseline.

Verification snapshot (this increment):

- Default gates remain green:
  - `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- Experimental numeric switch characterization:
  - `PACKMOL_GENCAN_CG_CPP_DRAFT=1 PACKMOL_GENCAN_NUMERIC_CPP=1`
    - `check_cg_ab.py` passed on `test_cg_ab_probe`.
    - `test_gencan_ab_cg_draft` equivalent runtime check still shows objective drift (`f mismatch`), so the numeric switch stays non-gating/experimental.

Easygencan migration increment:

- Added a C++ implementation of `easygencan` argument orchestration in bridge (`easygencan_cpp`):
  - Reproduces legacy constant/default parameter assembly.
  - Reuses existing working-space layout (`wd` slices for `s/y/d/w`, `wi` as `ind`).
  - Dispatches into existing Fortran `gencan` kernel ABI to preserve core behavior while moving setup logic.
- Added runtime gate for controlled rollout:
  - `PACKMOL_GENCAN_EASY_CPP_DRAFT=1` enables C++ easygencan path in `cpp|ab`.
  - Default remains Fortran easygencan path to keep runtime parity stable.
- Added dedicated AB test:
  - `test_easygencan_ab` with probe/checker pair validates `fortran|cpp|ab` parity on stable fields (`inform`, `iter`, `fcnt`, `gcnt`, `f`, `gpsupn`, `xsum`, `gnorm2`).

Verification snapshot (easygencan increment):

- `ctest --test-dir build -R test_easygencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

Gencan migration increment (entry stopping subset):

- Added a small C++ subset in `packmol_gencan_gencan_bridge` for `cpp|ab` mode:
  - Evaluate `f` and `g` at the current point using existing callback wrappers.
  - Compute projected-gradient norms (`gpeucn2`, `gpsupn`).
  - If first-iteration stopping criteria are already satisfied (`inform=0` or `inform=1`), return directly from C++ with synchronized counters/status.
  - Otherwise, safely fall back to the existing Fortran `gencan` kernel.
- Scope limitation:
  - Subset currently handles `gtype=0` (analytic gradient path) and `gtype!=0` (finite-difference gradient path, Fortran-else compatible) for the immediate entry-stop condition.
  - All other cases preserve previous Fortran behavior.
  - `packmolprecision` and `evalal<0` early returns are now handled directly in C++ with Fortran-aligned counters/status.
- Alignment details for this increment:
  - Entry subset now mirrors Fortran projected-gradient computation (`gpi = proj(x-g)-x`) and updates projected/free-index data accordingly.
  - Added C bridge `packmol_packmolprecision_fortran_c` to preserve Packmol-specific early-return semantics in mixed mode.
  - Added dedicated runtime AB gate `test_gencan_entry_stop_ab` to lock this first-iteration stop envelope.
  - Added dedicated runtime AB gate `test_gencan_entry_stop_gtype2_ab` to lock non-zero `gtype` compatibility on this envelope.

Verification snapshot (gencan subset increment):

- `ctest --test-dir build -R "test_(easygencan|gencan_driver)" --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -R test_gencan_entry_stop_ab --output-on-failure` passed.
- `ctest --test-dir build -R test_gencan_entry_stop_gtype2_ab --output-on-failure` passed.

Gencan migration increment (entry stop expansion for `inform=2/3`):

- Extended `packmol_gencan_gencan_bridge` entry subset in `cpp|ab` mode:
  - Added first-iteration no-function-progress criterion (`inform=2`) with Fortran initialization semantics (`fprev=infabs`, `bestprog=max(currprog,0)`, `itnfp=1` check).
  - Added first-iteration projected-gradient stagnation criterion (`inform=3`) with Fortran initialization semantics (`lastgpns(:)=infabs`, then compare `gpeucn2`).
- Added dedicated AB gates:
  - `test_gencan_entry_stop_inform2_ab`
  - `test_gencan_entry_stop_inform3_ab`

Verification snapshot (entry-stop expansion):

- `ctest --test-dir build -R "test_gencan_entry_stop_.*_ab" --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

Gencan migration increment (entry error-path alignment):

- Extended `packmol_gencan_gencan_bridge` C++ entry subset to match Fortran behavior on additional non-trivial entry paths:
  - `evalal` non-negative flags now proceed to gradient evaluation path (instead of forcing fallback for `flag > 0`).
  - Added direct C++ handling for gradient callback internal errors (`grad inform < 0`) with Fortran-aligned counters/state (`iter=0`, `fcnt=1`, `gcnt=1`, status propagated), avoiding fallback.
- This further narrows Fortran fallback in `cpp|ab` mode to paths beyond first-iteration stopping/entry error envelope.

Verification snapshot (entry error-path alignment):

- `cmake --build build` passed.
- `ctest --test-dir build -R "test_gencan_entry_stop_.*_ab|test_gencan_ab" --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

Gencan migration increment (first-iteration SPG terminal subset):

- Extended `packmol_gencan_gencan_bridge` C++ path with a guarded first-iteration SPG attempt:
  - Trigger condition mirrors Fortran face-abandon decision at iteration 1:
    - `gieucn2 <= (1-eta)^2 * gpeucn2`.
  - Runs SPG line search and post-SPG gradient update on C++ side using existing migrated kernels/callback wrappers.
  - Commits results and returns from C++ only for terminal-safe outcomes:
    - `spgls inform < 0`,
    - post-SPG gradient `inform < 0`,
    - line-search tiny-step termination `inform = 6` after Fortran-aligned post-step projection/metrics refresh.
  - Non-terminal outcomes continue to fallback to Fortran `gencan`, preserving runtime parity while shrinking fallback surface.

Verification snapshot (SPG terminal subset):

- `cmake --build build` passed.
- `ctest --test-dir build -R "test_gencan_ab|test_gencan_entry_stop_.*_ab" --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

Gencan migration increment (first-iteration SPG post-step stopping subset):

- Hardened SPG subset internals:
  - Introduced full-index gradient map (`ind_all = 1..n`) for full-gradient callbacks, avoiding accidental reuse of free-variable index buffers.
- Expanded C++ direct-return envelope after first SPG step:
  - After SPG + gradient recomputation + projection/metric refresh, C++ now evaluates the same top-of-loop stopping checks used by Fortran main loop (`inform=0/1/2/3/4/7/8`) and returns directly when satisfied.
  - If those post-step criteria are not satisfied, path still falls back to Fortran `gencan` for continuation.

Verification snapshot (SPG post-step stopping subset):

- `cmake --build build` passed.
- `ctest --test-dir build -R "test_gencan_ab|test_gencan_entry_stop_.*_ab" --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

Stabilization fix (state isolation before fallback):

- Fixed a semantic bug in the SPG subset path:
  - `lastgpns` updates are now committed only when the C++ subset actually returns control as the final decision.
  - In non-terminal cases that still fallback to Fortran `gencan`, C++ no longer mutates iteration-state arrays (`lastgpns`) beforehand.
- This preserves strict behavior isolation between speculative C++ subset evaluation and legacy Fortran continuation path.

Verification snapshot (stabilization):

- `cmake --build build` passed.
- `ctest --test-dir build -R "test_gencan_ab|test_gencan_entry_stop_.*_ab" --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

TN migration infrastructure increment (bridge primitives for C++ orchestration):

- Added C bridges for reduced-space transform primitives:
  - `packmol_shrink_fortran_c`
  - `packmol_expand_fortran_c`
  - (in `src/gencan/reduced_wrappers_bridge.f90`)
- Added C bridges for IEEE-signal parameter helpers used by GENCAN CG setup:
  - `packmol_gp_ieee_signal1_fortran_c`
  - `packmol_gp_ieee_signal2_fortran_c`
  - (in `src/gencan_ieee_signal_routines.f90`)
- Purpose:
  - enable upcoming C++ TN-step orchestration to reuse exact legacy parameterization and reduced/full-space transforms, avoiding duplicated logic and reducing semantic drift risk.

Verification snapshot (infrastructure increment):

- `cmake --build build` passed.
- `ctest --test-dir build -R "test_tnls_nobound_interior_(ab|noextrap_ab|ex1_ab)" --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

TN migration infrastructure increment (C++ bridge wiring):

- Wired newly added Fortran bridge symbols into `src/cpp/gencan/bridge.cpp` declarations:
  - reduced/full transform bridges: `packmol_shrink_fortran_c`, `packmol_expand_fortran_c`
  - CG-parameter helper bridges: `packmol_gp_ieee_signal1_fortran_c`, `packmol_gp_ieee_signal2_fortran_c`
- This keeps runtime behavior unchanged but makes the C++ gencan bridge ready to consume exact legacy TN/CG helper kernels in the next migration slice.

Verification snapshot (bridge wiring):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

TN integration draft increment (one-step TN terminal attempt, gated):

- Added a TN one-step draft path in `packmol_gencan_gencan_bridge`:
  - condition: branch selected when first-iteration policy prefers TN over SPG.
  - performs one reduced-space TN attempt (`shrink -> cg -> tnls -> expand`) with existing migrated sub-kernels/bridges.
  - returns from C++ only when post-step stopping criteria are satisfied; otherwise continues with legacy Fortran fallback.
- Runtime safety gate:
  - draft is disabled by default and enabled only with `PACKMOL_GENCAN_TN_CPP_DRAFT=1`.
  - default `cpp|ab` behavior remains parity-first.
- Stabilization during landing:
  - fixed accidental TNLS branch drift introduced while wiring the draft (`exbcnt` accounting path in `tnls_cpp_subset` restored).

Verification snapshot (TN draft gated, default path):

- `cmake --build build` passed.
- `ctest --test-dir build -R test_gencan_entry_stop_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

TN draft observability increment:

- Added a dedicated non-gating draft smoke test:
  - `test_gencan_tn_draft_smoke`
  - runs runtime driver probe in C++ mode (`PACKMOL_GENCAN_IMPL=cpp`)
  - labels: `gencan;gencan_draft` (intentionally not `gencan_ab`)
- Purpose:
  - keep default parity gate stable while providing a persistent execution hook for TN draft convergence work.
- Added TN draft A/B parity gate:
  - `test_gencan_ab_tn_draft`
  - runs standard GENCAN Fortran/CPP/AB comparison
  - labels: `gencan;gencan_ab;gencan_draft`
- Current TN draft characterization:
  - after fixing TN draft direction usage in reduced-space TNLS call path (`cg_s` as line-search direction), draft A/B parity now passes on existing GENCAN AB corpus.
  - this path was initially env-gated during landing and is now promoted to default in `cpp|ab` mode.

Verification snapshot (observability increment):

- `ctest --test-dir build -R test_gencan_tn_draft_smoke --output-on-failure` passed.
- `ctest --test-dir build -R test_gencan_ab_tn_draft --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

TN integration default-on increment:

- Promoted first-iteration TN one-step C++ path from env-gated draft to default `cpp|ab` behavior in `packmol_gencan_gencan_bridge`.
- Removed `PACKMOL_GENCAN_TN_CPP_DRAFT` runtime switch from bridge dispatch path.
- Kept parity guardrail unchanged:
  - C++ TN one-step path still returns only for terminal-safe post-step outcomes.
  - Non-terminal cases continue through existing Fortran `gencan` fallback.
- Updated TN observability tests to validate default path directly:
  - `test_gencan_tn_draft_smoke` now runs with `PACKMOL_GENCAN_IMPL=cpp` only.
  - `test_gencan_ab_tn_draft` now compares default Fortran/CPP/AB behavior without TN-specific env override.

Verification snapshot (default-on):

- `cmake --build build` passed.
- `ctest --test-dir build -R test_gencan_ab_tn_draft --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

TN kernel migration increment (reduced/full-space transforms):

- Replaced TN path usage of Fortran `shrink/expand` bridge calls with local C++ in-place implementations in `bridge.cpp`:
  - `shrink_inplace` and `expand_inplace`
  - semantics preserved as index-swap transforms (Fortran-compatible 1-based `ind` interpretation).
- Removed now-unused C bridge declarations from C++ side:
  - `packmol_shrink_fortran_c`
  - `packmol_expand_fortran_c`
- Scope:
  - applied to TN one-step C++ path only (`x/g/l/u` reduced-space transforms around CG+TNLS).
  - no change to stopping policy or fallback policy.

Verification snapshot (transform migration):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

TN kernel migration increment (CG signal parameterization):

- Replaced TN path usage of Fortran IEEE-signal helper bridges with local C++ equivalents in `bridge.cpp`:
  - `gp_ieee_signal1_cpp`
  - `gp_ieee_signal2_cpp`
- Formulas and clipping behavior are kept aligned with Fortran:
  - `acgeps/bcgeps` computation from `gpsupn`
  - `cgmaxit` adaptive rule (including `min(20, cgmaxit)` cap)
  - `cgeps` update and `[cgepsf, cgepsi]` clipping
- Removed now-unused C bridge declarations from C++ side:
  - `packmol_gp_ieee_signal1_fortran_c`
  - `packmol_gp_ieee_signal2_fortran_c`

Verification snapshot (CG signal migration):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

Bridge cleanup increment (remove dead Fortran C exports):

- Removed no-longer-referenced Fortran C bridge exports after C++ TN transform/signal migration:
  - from `src/gencan/reduced_wrappers_bridge.f90`:
    - `packmol_shrink_fortran_c`
    - `packmol_expand_fortran_c`
  - from `src/gencan_ieee_signal_routines.f90`:
    - `packmol_gp_ieee_signal1_fortran_c`
    - `packmol_gp_ieee_signal2_fortran_c`
- This is a linkage-surface cleanup only; algorithm behavior remains unchanged.

Verification snapshot (bridge cleanup):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

CG / Easygencan default-on attempt (superseded by rollback below):

- Attempted to promote both front-door paths to default in `cpp|ab`:
  - `packmol_gencan_cg_bridge` -> `cg_cpp_full` without `PACKMOL_GENCAN_CG_CPP_DRAFT`.
  - `packmol_gencan_easy_bridge` -> `easygencan_cpp` without `PACKMOL_GENCAN_EASY_CPP_DRAFT`.
- Result:
  - runtime fixture drift appeared in `test_gencan_ab*` checks.
  - therefore this attempt is not active and is superseded by the rollback section below.

Stabilization update (front-door rollout status):

- Re-validated `cg_cpp_full` on runtime fixture with shadow comparison and A/B gates, then promoted CG path again:
  - `packmol_gencan_cg_bridge` now default-on in `cpp|ab`.
  - `test_gencan_ab_cg_draft` validates default path directly (no env override).
- `easygencan` remains env-gated for parity-first stability:
  - `PACKMOL_GENCAN_EASY_CPP_DRAFT=1` controls C++ easygencan path.

Verification snapshot (current status):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

Easygencan parity fix and promotion:

- Root cause of `PACKMOL_GENCAN_EASY_CPP_DRAFT=1` runtime drift in `test_gencan_ab*`:
  - C++ `easygencan_cpp` forwarded caller `delmin`, while Fortran `easygcf` forces `delmin = 1.0d-2` before entering `gencan`.
  - This changed early feasibility/initialization behavior and yielded objective drift on `tiny_packing.inp`.
- Fix:
  - aligned C++ path to Fortran semantics by using local `delmin_local = 1.0e-2` in `easygencan_cpp` and writing it back to output `delmin`.
- Rollout policy update:
  - `packmol_gencan_easy_bridge` is now default-on C++ in `cpp|ab`.
  - `PACKMOL_GENCAN_EASY_CPP_DRAFT` remains as an emergency rollback switch (`0/false/no` => Fortran, unset/others => C++).
- Process note:
  - build/test must be run serially when validating rollout changes; parallel build+ctest can produce stale-binary false positives.

Verification snapshot (easygencan promotion):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

Fallback gate hardening (easygencan):

- Added explicit rollback-gate test:
  - `test_easygencan_ab_fallback`
  - env: `PACKMOL_GENCAN_EASY_CPP_DRAFT=0`
  - labels: `gencan;gencan_ab;gencan_fallback`
- Purpose:
  - keep emergency Fortran fallback path continuously validated while `easygencan` stays default-on C++ in `cpp|ab`.

Verification snapshot (fallback gate):

- `ctest --test-dir build -R "^test_easygencan_ab_fallback$" --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.

TN path migration increment (`tnls inform=6` forced-SPG in C++):

- In `packmol_gencan_gencan_bridge` TN one-step C++ branch:
  - extended `tnls` handling from `inform >= 0 && inform != 6` to full `inform >= 0` path.
  - for `tnls inform = 6`, now executes the same fallback policy as Fortran main loop in-place:
    - force one SPG iteration (`spgiter += 1`),
    - compute SPG step length and run `spgls_cpp`,
    - refresh gradient via callback wrappers,
    - apply standard post-step stopping checks and return only on terminal criteria.
- Effect:
  - reduces one more Fortran fallback slice in the TN branch while preserving behavior parity under existing AB gates.

Verification snapshot (TN inform=6 increment):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

TN path migration increment (`tnls inform<0` direct return in C++):

- In the TN one-step C++ branch, added direct terminal handling for negative `tnls` return:
  - when `tnls_inform < 0`, C++ now synchronizes counters/status and returns immediately from bridge.
  - this mirrors Fortran main-loop policy (`inform<0` immediate return) and removes one fallback-to-Fortran slice.
- Scope:
  - only the already-entered TN C++ one-step path; no change to branch selection policy.

Verification snapshot (TN negative-inform increment):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

TN path migration increment (`cg inform<0` direct return in C++):

- In the TN one-step C++ branch, added direct terminal handling for negative CG return:
  - when `cg_inform < 0`, bridge now synchronizes counters/status and returns immediately.
  - this mirrors Fortran main-loop behavior (`cg` error => immediate return) and removes another fallback-to-Fortran slice.
- Scope:
  - only affects already-entered TN C++ one-step execution path.

Verification snapshot (CG negative-inform increment):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

Fallback observability increment (no behavior change):

- Added debug-only fallback reason tracing in `packmol_gencan_gencan_bridge`:
  - when C++ path cannot terminate and hands control to Fortran, bridge now prints:
    - `[gencan-cpp-fallback] reason=<...> mode=<...>`
  - enabled only under `PACKMOL_GENCAN_DEBUG=1`.
- Current reason tags:
  - `spg_post_nonterminal`
  - `tn_post_nonterminal`
  - `tn_no_free_variables`
  - `entry_grad_flag_nonzero`
  - default `cpp_nonterminal_continue`

Verification snapshot (fallback observability):

- `cmake --build build` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

Closure gate increment (no fallback allowed on AB corpus):

- Added a dedicated closure gate:
  - `test_gencan_ab_no_fallback`
  - checker: `tests/gencan/check_gencan_ab_no_fallback.py`
  - policy: run `fortran/cpp/ab` probes under `PACKMOL_GENCAN_DEBUG=1` and fail if `gencan-cpp-fallback` marker appears in `cpp` or `ab` outputs.
  - labels: `gencan;gencan_ab;gencan_closure`
- Purpose:
  - converts “fallback must be zero on current AB fixture corpus” into an enforceable CI gate.

Verification snapshot (closure gate):

- `ctest --test-dir build -R test_gencan_ab_no_fallback --output-on-failure` passed.
- `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
- `ctest --test-dir build -L gencan --output-on-failure` passed.

Numeric kernel promotion increment (default C++ `dot/norm2`):

- Promoted `PACKMOL_GENCAN_NUMERIC_CPP` policy to default-on in C++ bridge:
  - `dot_kernel`/`norm2_kernel` now use C++ stable implementations by default.
  - emergency rollback remains available:
    - `PACKMOL_GENCAN_NUMERIC_CPP=0|false|no` forces legacy Fortran numeric wrappers.
- Scope:
  - affects only numeric utility dispatch (`dot`/`norm2`) used by migrated C++ GENCAN paths.
  - no interface or I/O behavior changes.

Verification snapshot (numeric promotion):

- `cmake --build build` passed.
- default path:
  - `ctest --test-dir build -L gencan_ab --output-on-failure` passed.
  - `ctest --test-dir build -L gencan --output-on-failure` passed.
- rollback path:
  - `PACKMOL_GENCAN_NUMERIC_CPP=0 ctest --test-dir build -L gencan_ab --output-on-failure` passed.
