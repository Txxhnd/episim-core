# EPISIM-CORE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver `episim-core` v0.1.0 as a dependency-light Rust library that provides validated deterministic SIR and SEIR simulations using a shared classical RK4 solver.

**Architecture:** The public API exposes concrete `simulate_sir` and `simulate_seir` functions. Both call a shared internal simulation driver and fixed-size-array RK4 integrator; model equations, input validation, trajectory recording, diagnostics, and metric calculation remain separate modules. No CLI, GUI, network service, persistence, stochastic model, or public plugin trait is included.

**Tech Stack:** Rust 2024, MSRV 1.85, Cargo, `f64`, standard library, built-in test framework, Criterion as a development-only benchmark dependency.

**Spec:** `TECHNICAL_SPEC.md` in the repository root (implements Part II of `Product_Vision_&_MVP_Product_Requirements_Document.md`; product direction is in Part I of the same file).

> **Frozen at v0.1.0 (November 14, 2026).** This document records what release 0.1.0 set out to do and how it was built. It is no longer updated. For current behaviour, see `docs/`, the API documentation and `CHANGELOG.md`; design decisions after 0.1.0 are recorded under `design/`.

**Status:** Complete. Tasks 1 to 13 are implemented. Tasks 1 to 10 were reviewed as a whole by an independent reviewer, and Tasks 11 and 12 came out of that review. Tasks 11 and 12 were then reviewed in turn, which led to Task 13. Task 13 itself has not been independently reviewed; its guard was checked by the randomised stability evidence test that improvement plan 002 added. Where the implementation departed from the task text, the Execution Notes section near the end records what was done and why. Work after this plan is summarised in "After This Plan" and tracked in `plans/README.md`.

## Global Constraints

- Pure Rust library; no CLI, GUI, HTTP API, database, or WebAssembly binding in v0.1.0.
- Rust edition is `2024`; declared `rust-version` is `1.85`.
- Core simulation functionality has zero mandatory runtime dependencies.
- `#![forbid(unsafe_code)]` is required.
- Numerical precision is `f64`; population counts may be fractional expected counts.
- Initial models are deterministic SIR and SEIR with a closed population, homogeneous mixing, constant parameters, and permanent immunity.
- Transmission is frequency-dependent: infection flow is `beta * (S / N) * I`.
- Classical RK4 with a fixed timestep is the only solver in the MVP.
- Final partial timestep must end exactly at the requested duration within the documented floating-point time handling.
- Default maximum number of integration steps is 1,000,000; explicit overrides must be supported.
- Initial and accepted compartment values must be finite; input compartments must be non-negative.
- Beta must be finite and non-negative; gamma and (for SEIR) sigma must be finite and strictly positive.
- A timestep for which the model's fastest decay rate multiplied by the step taken exceeds 2 must be rejected before integration with `TimeStepTooLarge { max_step_days }`. The fastest rate is `max(beta, gamma)` for SIR and `max(beta, (sigma + gamma) / 2 + sqrt(((sigma - gamma) / 2)^2 + sigma * beta))` for SEIR (specification section 6.3).
- Materially negative numerical states, non-finite intermediate values, or population conservation drift above the specified tolerance must return typed errors; do not silently clip values.
- Relative population drift target is `1e-9`; tolerated negative roundoff uses `64 * f64::EPSILON * max(N0, 1)` and must be counted in diagnostics.
- Results contain the initial observation, each accepted step, and the final observation.
- Peak infectious metrics use recorded observations; ties select the earliest time. If every recorded infectious value is zero, peak infectious time is day zero.
- Cumulative infections since simulation start are evaluated at the final simulation time `T` as `S(0) - S(T)`; cumulative incidence fraction divides that value by `N(0)`.
- `R0` is the model-derived theoretical value `beta / gamma`; do not expose an empirical `Rt` metric.
- Metrics named cumulative infections and cumulative incidence fraction refer to the final simulation time, not every individual observation.
- Public functions return `Result`; normal invalid user input must not panic.
- Use `Vec::try_reserve`/`try_reserve_exact` for trajectory reservation and map reservation failure to a typed error.
- Every task must follow TDD, run its targeted tests, and end in a focused Git commit.

## Review Focus

1. **Finite compartment values whose sum overflows:** reject an initial state when its computed population total is non-finite; test SIR and SEIR separately.
2. **Unrepresentably small timestep:** detect when adding the next step cannot advance the floating-point time; return `TimeStepTooSmall` rather than looping indefinitely.
3. **Derived metric overflow:** a finite beta divided by a very small positive gamma can overflow; return `MetricOverflow` instead of returning an infinite `R0`.
4. **Final-step rounding and step-count boundaries:** ensure the last observation is exactly the requested duration, and enforce the maximum-step limit without off-by-one errors.
5. **Non-finite intermediate RK4 stages:** detect overflow/NaN in `k1` through `k4`, not only in the final state; return a numerical error with the affected stage/context where practical.

---

## File Structure

The plan assumes a new repository. If a repository already exists, inspect its conventions first and preserve compatible structure rather than overwriting it.

- Create: `Cargo.toml` — crate metadata, edition, MSRV, benchmark dev-dependency.
- Create: `src/lib.rs` — crate attributes and stable public re-exports.
- Create: `src/error.rs` — public typed errors and `Display`/`Error` implementations.
- Create: `src/settings.rs` — validated simulation duration, timestep, and maximum-step configuration.
- Create: `src/validation.rs` — shared finite/range/population validation helpers.
- Create: `src/result.rs` — public observations, simulation results, and diagnostics.
- Create: `src/metrics.rs` — model-derived epidemiological metrics.
- Create: `src/models/mod.rs` — private model module declarations and model-specific interfaces.
- Create: `src/models/sir.rs` — SIR public types, equations, and wrapper function.
- Create: `src/models/seir.rs` — SEIR public types, equations, and wrapper function.
- Create: `src/solver/mod.rs` — private solver exports.
- Create: `src/solver/rk4.rs` — fixed-size-array RK4 stepper.
- Create: `src/simulation/mod.rs` — private simulation module exports.
- Create: `src/simulation/driver.rs` — shared integration loop, step counting, trajectory reservation, and diagnostics.
- Create: `tests/crate_smoke_test.rs` — compile-only library smoke test.
- Create: `tests/validation_tests.rs` — invalid inputs and configuration boundaries.
- Create: `tests/rk4_tests.rs` — solver correctness and intermediate-stage failures.
- Create: `tests/sir_tests.rs` — SIR API and analytical scenarios.
- Create: `tests/seir_tests.rs` — SEIR API and analytical scenarios.
- Create: `tests/numerical_tests.rs` — conservation, positivity, and diagnostic behavior.
- Create: `tests/metrics_tests.rs` — metrics tested against hand-constructed trajectories.
- Create: `tests/convergence_tests.rs` — timestep refinement against independent references.
- Create: `examples/basic_sir.rs` and `examples/basic_seir.rs` — downstream-consumable examples.
- Create: `benches/simulation_bench.rs` — Criterion benchmarks for 10k, 100k, and 1m steps.
- Create: `README.md` — quick start, assumptions, scope, and limitations.
- Create: `docs/architecture.md`, `docs/mathematical_models.md`, `docs/numerical_methods.md`, and `docs/api_design.md` — technical documentation.

## API Contracts Used Throughout the Plan

### Simulation settings

Use validated construction to prevent invalid settings from entering the driver:

```rust
pub struct SimulationSettings {
    duration_days: f64,
    step_days: f64,
    max_steps: usize,
}

impl SimulationSettings {
    pub fn new(duration_days: f64, step_days: f64) -> Result<Self, SimulationError>;
    pub fn with_max_steps(self, max_steps: usize) -> Result<Self, SimulationError>;
    pub fn duration_days(&self) -> f64;
    pub fn step_days(&self) -> f64;
    pub fn max_steps(&self) -> usize;
}
```

`new` sets `max_steps` to `1_000_000`. `with_max_steps` rejects zero. This is the concrete implementation of the PRD's default and override behavior.

### Public model functions

```rust
pub fn simulate_sir(
    parameters: SirParameters,
    initial_state: SirState,
    settings: SimulationSettings,
) -> Result<SimulationResult<SirState>, SimulationError>;

pub fn simulate_seir(
    parameters: SeirParameters,
    initial_state: SeirState,
    settings: SimulationSettings,
) -> Result<SimulationResult<SeirState>, SimulationError>;
```

### Public results

```rust
pub struct Observation<S> {
    pub time_days: f64,
    pub state: S,
}

pub struct SimulationResult<S> {
    pub trajectory: Vec<Observation<S>>,
    pub metrics: EpidemicMetrics,
    pub diagnostics: SimulationDiagnostics,
}

pub struct EpidemicMetrics {
    pub basic_reproduction_number: f64,
    pub peak_infectious: f64,
    pub peak_infectious_time_days: f64,
    pub cumulative_infections_since_start: f64,
    pub cumulative_incidence_fraction: f64,
}

pub struct SimulationDiagnostics {
    pub steps_executed: usize,
    pub maximum_absolute_population_drift: f64,
    pub maximum_relative_population_drift: f64,
    pub tolerated_negative_value_count: usize,
}
```

### Internal solver contract

```rust
pub(crate) fn rk4_step<const N: usize, F>(
    state: [f64; N],
    step: f64,
    derivative: F,
) -> Result<[f64; N], SimulationError>
where
    F: Fn([f64; N]) -> [f64; N];
```

The derivative callback is evaluated at all four RK4 stages. Every stage derivative and intermediate/final state must be checked for finiteness. The solver itself does not clip negative states; the driver applies model-aware validation and diagnostics. Use `NonFiniteIntermediate { stage }` for a non-finite RK4 intermediate stage or derivative, and `NonFiniteNumericalState { compartment }` for an invalid accepted model state. The driver must check initial and post-step population sums using the fixed initial population supplied by the model wrapper.

---

## Task 1: Initialize the Rust Crate and Quality Gates

**Files:**
- Create: `Cargo.toml`
- Create: `src/lib.rs`
- Create: `.gitignore`

**Interfaces:**
- Produces a compilable library crate named `episim-core`.
- Exposes no simulation API yet.

- [x] **Step 1: Write the crate smoke test**

Create `tests/crate_smoke_test.rs` with `use episim_core as _;` and a minimal `#[test] fn crate_compiles_as_library() {}`. This is a compile smoke test and must not require a placeholder public API.

- [x] **Step 2: Run the test to verify the scaffold is missing**

Run: `cargo test --test crate_smoke_test`
Expected: FAIL because the crate manifest/source has not yet been created.

- [x] **Step 3: Create the minimum crate scaffold**

Set package name `episim-core`, version `0.1.0`, edition `2024`, and `rust-version = "1.85"`. Add `#![forbid(unsafe_code)]`, `#![warn(missing_docs)]`, and `#![deny(missing_debug_implementations)]` to `src/lib.rs`. Add standard Rust build outputs to `.gitignore`.

- [x] **Step 4: Run the smoke test and baseline quality checks**

Run: `cargo test --test crate_smoke_test`
Expected: PASS.

Run: `cargo fmt --all -- --check && cargo clippy --all-targets -- -D warnings`
Expected: PASS with no warnings.

- [x] **Step 5: Commit**

```bash
git add Cargo.toml src/lib.rs .gitignore tests/crate_smoke_test.rs
git commit -m "chore: initialize episim-core crate"
```

## Task 2: Define Typed Errors and Validated Simulation Settings

**Files:**
- Create: `src/error.rs`
- Create: `src/settings.rs`
- Modify: `src/lib.rs`
- Create: `tests/validation_tests.rs`

**Interfaces:**
- Produces `SimulationError` implementing `Debug`, `Display`, and `std::error::Error`.
- Produces `SimulationSettings::new`, `with_max_steps`, and getters exactly as listed in API Contracts.
- Default maximum steps: `1_000_000`.

- [x] **Step 1: Write failing settings tests**

Add tests named:
- `settings_accepts_positive_finite_duration_and_step`
- `settings_rejects_zero_negative_nan_and_infinite_duration`
- `settings_rejects_zero_negative_nan_and_infinite_step`
- `settings_defaults_to_one_million_max_steps`
- `settings_rejects_zero_max_steps`
- `settings_accepts_explicit_positive_max_steps`

Assert valid settings preserve the supplied values and invalid settings return the correct error variant.

- [x] **Step 2: Run targeted tests to verify failure**

Run: `cargo test --test validation_tests settings_`
Expected: FAIL because settings and error types do not exist.

- [x] **Step 3: Implement settings and error types**

Implement private settings fields, validated constructors, getters, and these error variants: `InvalidParameter { name }`, `InvalidInitialState { compartment }`, `InvalidPopulation`, `InvalidDuration`, `InvalidTimeStep`, `InvalidMaximumSteps`, `StepLimitExceeded`, `TimeStepTooSmall`, `NonFiniteIntermediate { stage }`, `NonFiniteNumericalState { compartment }`, `NegativeCompartment { compartment }`, `PopulationConservationViolation`, `MetricOverflow`, and `AllocationFailure`. Implement `Display` manually to avoid a mandatory error-handling dependency. Add contextual numeric fields to error variants where they materially help diagnosis, without exposing internal implementation types.

- [x] **Step 4: Run settings tests**

Run: `cargo test --test validation_tests settings_`
Expected: PASS.

- [x] **Step 5: Commit**

```bash
git add src/error.rs src/settings.rs src/lib.rs tests/validation_tests.rs
git commit -m "feat: add simulation settings and typed errors"
```

## Task 3: Define Public Model Types and Input Validation

**Files:**
- Create: `src/validation.rs`
- Create: `src/models/mod.rs`
- Create: `src/models/sir.rs`
- Create: `src/models/seir.rs`
- Modify: `src/lib.rs`
- Modify: `tests/validation_tests.rs`

**Interfaces:**
- `SirParameters { beta_per_day: f64, gamma_per_day: f64 }`.
- `SirState { susceptible: f64, infected: f64, recovered: f64 }`.
- `SeirParameters { beta_per_day: f64, sigma_per_day: f64, gamma_per_day: f64 }`.
- `SeirState { susceptible: f64, exposed: f64, infected: f64, recovered: f64 }`.
- Model state and parameter structs derive `Debug`, `Clone`, `Copy`, and `PartialEq`.
- Internal validation helpers return typed errors and calculate a finite positive initial population.

- [x] **Step 1: Write failing model validation tests**

Add tests:
- `sir_accepts_zero_beta_and_positive_gamma`
- `sir_rejects_negative_or_non_finite_beta`
- `sir_rejects_non_positive_or_non_finite_gamma`
- `seir_accepts_zero_beta_and_positive_sigma_gamma`
- `seir_rejects_invalid_sigma`
- `seir_rejects_invalid_gamma`
- `sir_rejects_negative_or_non_finite_compartments`
- `seir_rejects_negative_or_non_finite_compartments`
- `sir_rejects_zero_or_overflowed_total_population`
- `seir_rejects_zero_or_overflowed_total_population`

The overflow tests must use individually finite compartment values whose sum becomes non-finite.

- [x] **Step 2: Run targeted tests to verify failure**

Run: `cargo test --test validation_tests sir_`
Run: `cargo test --test validation_tests seir_`
Expected: FAIL because model types and validators are not implemented.

- [x] **Step 3: Implement model types and validators**

Implement the public structs and internal validation functions. Validate each compartment before summation; use a finite check on the summed population and require `N0 > 0`. Do not place simulation logic in this task. Export the public structs from the crate root only after their documentation is complete.

- [x] **Step 4: Run validation tests**

Run: `cargo test --test validation_tests`
Expected: PASS.

- [x] **Step 5: Commit**

```bash
git add src/validation.rs src/models src/lib.rs tests/validation_tests.rs
git commit -m "feat: add epidemiological model types and validation"
```

## Task 4: Implement the Fixed-Size RK4 Solver

**Files:**
- Create: `src/solver/mod.rs`
- Create: `src/solver/rk4.rs`
- Create: `tests/rk4_tests.rs`

**Interfaces:**
- Internal `rk4_step<const N: usize, F>(state: [f64; N], step: f64, derivative: F) -> Result<[f64; N], SimulationError>`.
- `F: Fn([f64; N]) -> [f64; N]`.
- No heap allocation in one RK4 step.

- [x] **Step 1: Write failing RK4 reference tests**

Add tests:
- `rk4_integrates_scalar_exponential_decay_with_fourth_order_accuracy`
- `rk4_integrates_constant_derivative_exactly_within_roundoff`
- `rk4_supports_three_and_four_component_arrays`
- `rk4_rejects_non_finite_step`
- `rk4_rejects_non_finite_intermediate_stage`
- `rk4_rejects_non_finite_final_state`

For exponential decay, compare against `exp(-t)` at a fixed endpoint and verify an error-reduction ratio consistent with fourth-order convergence over selected sufficiently small steps.

- [x] **Step 2: Run tests to verify failure**

Run: `cargo test --test rk4_tests`
Expected: FAIL because `rk4_step` is not implemented.

- [x] **Step 3: Implement `rk4_step`**

Calculate `k1`, `k2`, `k3`, and `k4` using fixed-size arrays. Check every intermediate state, derivative array, and final value for finiteness. Return a typed numerical error rather than propagating NaN or infinity.

- [x] **Step 4: Run RK4 tests**

Run: `cargo test --test rk4_tests`
Expected: PASS.

- [x] **Step 5: Commit**

```bash
git add src/solver tests/rk4_tests.rs
git commit -m "feat: implement fixed-step rk4 integrator"
```

## Task 5: Implement the Shared Simulation Driver and Trajectory Diagnostics

**Files:**
- Create: `src/simulation/mod.rs`
- Create: `src/simulation/driver.rs`
- Create: `src/result.rs`
- Create: `tests/numerical_tests.rs`
- Modify: `src/lib.rs`

**Interfaces:**
- `Observation<S> { time_days: f64, state: S }`, deriving `Debug` and `Clone` (and `PartialEq` when supported by `S`).
- `SimulationDiagnostics` with step count and population drift fields, deriving `Debug`, `Clone`, and `PartialEq`.
- Internal `SimulationOutput<S> { trajectory: Vec<Observation<S>>, diagnostics: SimulationDiagnostics }`.
- Internal `run_simulation<const N: usize, S, D, C>(initial: [f64; N], population: f64, settings: SimulationSettings, derivative: D, to_state: C, compartment_names: [&'static str; N]) -> Result<SimulationOutput<S>, SimulationError>` where `D: Fn([f64; N]) -> [f64; N]` and `C: Fn([f64; N]) -> S`.
- The driver does not calculate epidemiological metrics; Task 6 adds metrics and the public `SimulationResult<S>` wrapper.

- [x] **Step 1: Write failing driver and diagnostics tests**

Add tests:
- `driver_records_initial_state`
- `driver_records_final_state_at_requested_duration`
- `driver_uses_short_final_partial_step`
- `driver_counts_executed_steps`
- `driver_rejects_step_limit_exceeded_before_integration`
- `driver_rejects_time_step_that_cannot_advance_clock`
- `driver_checks_population_drift_at_every_accepted_step`
- `driver_rejects_materially_negative_compartment`
- `driver_counts_tolerated_negative_roundoff_without_clipping`
- `driver_maps_trajectory_reservation_failure_to_typed_error` (use an injected reservation abstraction or a narrowly scoped test seam; do not attempt to exhaust system memory in a test).

- [x] **Step 2: Run tests to verify failure**

Run: `cargo test --test numerical_tests`
Expected: FAIL because result and driver types are not implemented.

- [x] **Step 3: Implement result structures and driver**

Reserve `required_steps + 1` observations using checked arithmetic and `try_reserve_exact`. Calculate the required number of steps without an unchecked float-to-integer cast: reject a non-finite duration/step ratio or any ratio exceeding `max_steps` before converting it, then use checked arithmetic for the observation capacity. Enforce the configured step limit before integration. Advance time by `min(step_days, duration_days - current_time)` and detect a step that cannot advance the floating-point clock. Check finiteness, compartment negativity tolerance, and population drift after each accepted step. Do not alter tolerated negative values; count them.

- [x] **Step 4: Run numerical diagnostics tests**

Run: `cargo test --test numerical_tests`
Expected: PASS.

- [x] **Step 5: Commit**

```bash
git add src/result.rs src/simulation src/lib.rs tests/numerical_tests.rs
git commit -m "feat: add shared simulation driver and diagnostics"
```

## Task 6: Implement Result Types and Epidemiological Metrics

**Files:**
- Create: `src/metrics.rs`
- Modify: `src/result.rs`
- Modify: `src/lib.rs`
- Create: `tests/metrics_tests.rs`

**Interfaces:**
- `EpidemicMetrics` and `SimulationResult<S>` exactly as defined in API Contracts.
- Internal `calculate_metrics<S, FS, FI>(trajectory: &[Observation<S>], initial_population: f64, beta_per_day: f64, gamma_per_day: f64, susceptible_of: FS, infectious_of: FI) -> Result<EpidemicMetrics, SimulationError>` where `FS: Fn(&S) -> f64` and `FI: Fn(&S) -> f64`.

- [x] **Step 1: Write failing metric tests**

Add tests:
- `metrics_calculate_basic_reproduction_number`
- `metrics_calculate_peak_infectious_and_earliest_peak_time`
- `metrics_calculate_cumulative_infections_since_start`
- `metrics_calculate_cumulative_incidence_fraction`
- `metrics_zero_transmission_returns_zero_r0`
- `metrics_reject_non_finite_reproduction_number`
- `metrics_handle_zero_peak_infection`
- `metrics_reject_non_finite_trajectory_values`

Use small hand-constructed trajectories to test metric calculations independently of numerical integration.

- [x] **Step 2: Run tests to verify failure**

Run: `cargo test --test metrics_tests`
Expected: FAIL because metrics and public result types are not implemented.

- [x] **Step 3: Implement result and metric structures**

Implement `EpidemicMetrics` and `SimulationResult<S>`. Calculate `R0 = beta / gamma` and reject non-finite results. Calculate peak infectious value from recorded observations; if there is a tie, retain the earliest observation. Calculate cumulative infections at the final simulation time as `S(0) - S(T)` and cumulative incidence fraction as that value divided by `N0`. Validate that metric inputs and outputs are finite. Do not add an `Rt` field.

- [x] **Step 4: Run metric tests**

Run: `cargo test --test metrics_tests`
Expected: PASS.

- [x] **Step 5: Commit**

```bash
git add src/metrics.rs src/result.rs src/lib.rs tests/metrics_tests.rs
git commit -m "feat: add epidemiological result metrics"
```

## Task 7: Implement the Public SIR Simulation API

**Files:**
- Modify: `src/models/sir.rs`
- Modify: `src/lib.rs`
- Create: `tests/sir_tests.rs`

**Interfaces:**
- `simulate_sir(parameters: SirParameters, initial_state: SirState, settings: SimulationSettings) -> Result<SimulationResult<SirState>, SimulationError>`.
- SIR compartment array order is `[susceptible, infected, recovered]`.
- Infectious metric accessor returns `state.infected`; susceptible accessor returns `state.susceptible`.

- [x] **Step 1: Write failing SIR integration tests**

Add tests:
- `sir_returns_initial_and_final_observations`
- `sir_zero_transmission_keeps_susceptible_constant`
- `sir_zero_initial_infection_has_no_new_infections`
- `sir_beta_zero_matches_analytical_recovery_solution`
- `sir_preserves_total_population_within_tolerance`
- `sir_returns_error_for_invalid_parameters_and_state`

Use beta = 0 and gamma > 0 as an analytical reference scenario. Verify `S(t) = S0`, `I(t) = I0 * exp(-gamma * t)`, and `R(t) = R0_initial + I0 - I(t)` within an explicitly chosen numerical tolerance.

- [x] **Step 2: Run tests to verify failure**

Run: `cargo test --test sir_tests`
Expected: FAIL because `simulate_sir` is not implemented.

- [x] **Step 3: Implement the SIR wrapper and derivative**

Convert `SirState` to `[f64; 3]`, compute the frequency-dependent infection flow and recovery flow, call `run_simulation` with compartment names `["susceptible", "infected", "recovered"]`, then calculate metrics from the returned trajectory and assemble `SimulationResult<SirState>`. Re-export `simulate_sir`, `SirParameters`, and `SirState` from the crate root.

- [x] **Step 4: Run SIR tests**

Run: `cargo test --test sir_tests`
Expected: PASS.

- [x] **Step 5: Commit**

```bash
git add src/models/sir.rs src/lib.rs tests/sir_tests.rs
git commit -m "feat: implement deterministic sir simulation"
```

## Task 8: Implement the Public SEIR Simulation API

**Files:**
- Modify: `src/models/seir.rs`
- Modify: `src/lib.rs`
- Create: `tests/seir_tests.rs`

**Interfaces:**
- `simulate_seir(parameters: SeirParameters, initial_state: SeirState, settings: SimulationSettings) -> Result<SimulationResult<SeirState>, SimulationError>`.
- SEIR compartment array order is `[susceptible, exposed, infected, recovered]`.
- Infectious metric accessor returns `state.infected`; susceptible accessor returns `state.susceptible`.

- [x] **Step 1: Write failing SEIR integration tests**

Add tests:
- `seir_returns_initial_and_final_observations`
- `seir_zero_transmission_keeps_susceptible_constant`
- `seir_zero_exposed_and_infected_has_no_new_infections`
- `seir_beta_zero_matches_analytical_exposed_and_infected_solution`
- `seir_preserves_total_population_within_tolerance`
- `seir_returns_error_for_invalid_parameters_and_state`

For beta = 0 and sigma != gamma, compare exposed and infected trajectories against the analytical solution in the technical specification.

- [x] **Step 2: Run tests to verify failure**

Run: `cargo test --test seir_tests`
Expected: FAIL because `simulate_seir` is not implemented.

- [x] **Step 3: Implement the SEIR wrapper and derivative**

Convert `SeirState` to `[f64; 4]`, implement infection, exposed progression, and recovery flows, call `run_simulation` with compartment names `["susceptible", "exposed", "infected", "recovered"]`, then calculate metrics and assemble `SimulationResult<SeirState>`. Re-export `simulate_seir`, `SeirParameters`, and `SeirState` from the crate root.

- [x] **Step 4: Run SEIR tests**

Run: `cargo test --test seir_tests`
Expected: PASS.

- [x] **Step 5: Commit**

```bash
git add src/models/seir.rs src/lib.rs tests/seir_tests.rs
git commit -m "feat: implement deterministic seir simulation"
```

## Task 9: Add Cross-Model Numerical and Convergence Verification

**Files:**
- Create: `tests/convergence_tests.rs`
- Modify: `tests/numerical_tests.rs`
- Modify: `tests/sir_tests.rs`
- Modify: `tests/seir_tests.rs`

**Interfaces:**
- No public API changes.
- Tests consume `simulate_sir` and `simulate_seir` only.

- [x] **Step 1: Write convergence and edge-case tests**

Add tests:
- `sir_timestep_refinement_reduces_reference_error`
- `seir_timestep_refinement_reduces_reference_error`
- `sir_population_drift_stays_within_release_threshold`
- `seir_population_drift_stays_within_release_threshold`
- `sir_final_time_equals_requested_duration`
- `seir_final_time_equals_requested_duration`
- `metrics_and_diagnostics_remain_finite_for_reference_scenarios`

Use independently calculated analytical reference scenarios where possible. For nontrivial scenarios, compare against a much finer-step reference and document that it is a numerical reference, not an exact solution. Assert a reduction factor of at least 8 for the designated smooth convergence fixtures when halving the timestep, while avoiding the floating-point error floor.

- [x] **Step 2: Run the new tests to verify failure where appropriate**

Run: `cargo test --test convergence_tests`
Expected: Initial convergence assertions may expose tolerance or solver defects; do not weaken thresholds without documenting numerical evidence.

- [x] **Step 3: Correct implementation defects only**

Fix solver, time advancement, model derivative, or diagnostic defects revealed by the tests. Do not clip invalid compartment values or change the stated model equations to make tests pass.

- [x] **Step 4: Run all functional and numerical tests**

Run: `cargo test --all-targets`
Expected: PASS.

- [x] **Step 5: Commit**

```bash
git add tests
git commit -m "test: verify numerical accuracy and convergence"
```

## Task 10: Add Examples, Documentation, and Release Quality Gates

**Files:**
- Create: `examples/basic_sir.rs`
- Create: `examples/basic_seir.rs`
- Create: `README.md`
- Create: `docs/architecture.md`
- Create: `docs/mathematical_models.md`
- Create: `docs/numerical_methods.md`
- Create: `docs/api_design.md`
- Modify: `Cargo.toml`
- Create: `benches/simulation_bench.rs`

**Interfaces:**
- Examples use only the public API.
- Criterion is a development-only dependency.
- No production API changes.

- [x] **Step 1: Write executable example checks**

Create both examples and ensure each prints a small summary of model configuration, peak infectious population, peak time, cumulative infections, and numerical diagnostics.

- [x] **Step 2: Run examples to verify the public workflow**

Run: `cargo check --examples`
Expected: FAIL until examples and all required public re-exports are correctly connected; then PASS.

- [x] **Step 3: Write documentation**

Document equations, units, model assumptions, frequency-dependent transmission, RK4 stages, fixed-step limitations, negative-state behavior, result metric definitions, and a downstream quick start. Clearly state that the library is a simulation engine, not a calibrated real-world forecasting product.

- [x] **Step 4: Add benchmark workloads**

Add Criterion as a dev-dependency and configure `[[bench]] name = "simulation_bench"` with `harness = false`. Benchmark SIR and SEIR at 10,000, 100,000, and 1,000,000 steps using release mode. Record compiler, CPU, operating system, and workload parameters with the baseline; do not invent absolute runtime targets.

- [x] **Step 5: Run full quality and packaging checks**

Run:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets
cargo doc --no-deps --all-features
cargo check --examples
cargo package --list
cargo package
```

Expected: all commands pass. Resolve warnings and packaging omissions before release.

- [x] **Step 6: Commit**

```bash
git add Cargo.toml Cargo.lock README.md docs examples benches
git commit -m "docs: prepare episim-core v0.1.0 release"
```

## Task 11: Reject Timesteps Beyond the RK4 Stability Limit

> **Superseded in part by Task 13.** A review found that checking each rate separately against 2.785 does not prevent instability in SEIR. Task 13 replaced the rule, the limit, and the shape of the error variant. The text below is kept as the record of what Task 11 did.

Added after the final review of Tasks 1 to 10. The review found that a step too large for a rate can return `Ok` with a wrong trajectory: with beta = 0, gamma = 3, and a step of one day, the infected compartment grows by 37.5% per step, stays positive, and conserves the population. The same happens for susceptibles when `beta * step_days` is about 3 and nearly everyone is infected. See Section 6.3 of the technical specification.

**Files:**
- Modify: `src/error.rs`
- Modify: `src/validation.rs`
- Modify: `src/models/sir.rs`
- Modify: `src/models/seir.rs`
- Modify: `tests/sir_tests.rs`
- Modify: `tests/seir_tests.rs`
- Modify: `README.md`, `docs/numerical_methods.md`, `docs/api_design.md`

**Interfaces:**
- Adds the error variant `TimeStepTooLarge { parameter: &'static str }`.
- Internal `stable_step(step_days: f64, rates: &[(&'static str, f64)]) -> Result<(), SimulationError>` in `src/validation.rs`, using the limit `2.785`. It returns the first rate, in the order given, whose product with `step_days` is greater than the limit.
- `simulate_sir` checks `beta_per_day` and `gamma_per_day`; `simulate_seir` checks `beta_per_day`, `sigma_per_day`, and `gamma_per_day`. The check runs after parameter and initial-state validation and before integration.
- No other public API changes.

- [x] **Step 1: Write failing stability tests**

Add tests:
- `sir_rejects_step_too_large_for_stable_integration`
- `seir_rejects_step_too_large_for_stable_integration`

For each rate parameter, assert that a product of 2.786 with the step returns `TimeStepTooLarge` naming that parameter, and that the same rate with half the step runs. Pin the boundary: a product of exactly 2.785 must not return `TimeStepTooLarge`. In SIR with beta = 0 that run succeeds and the infected decay. In SEIR with beta = 0 and sigma at the limit the run passes the guard and is then stopped with `NegativeCompartment { compartment: "infected" }`, which is an acceptable outcome because it is an error and not a wrong trajectory.

- [x] **Step 2: Run tests to verify failure**

Run: `cargo test --no-fail-fast --test sir_tests --test seir_tests too_large`
Expected: FAIL because the variant does not exist; once the variant alone is added, FAIL because the unstable configurations return `Ok`.

- [x] **Step 3: Implement the guard**

Add the variant and its `Display` text, add `stable_step`, and call it from both public functions. Update the `# Errors` documentation of both functions, and describe the guard and its limits in `README.md`, `docs/numerical_methods.md`, and `docs/api_design.md`. State clearly that it is a stability guard and does not detect inaccuracy.

- [x] **Step 4: Run all quality checks**

Run: `cargo test`
Expected: PASS.

Run: `cargo fmt --all -- --check && cargo clippy --all-targets --all-features -- -D warnings`
Expected: PASS with no warnings.

- [x] **Step 5: Commit**

```bash
git add src/error.rs src/validation.rs src/models tests/sir_tests.rs tests/seir_tests.rs README.md docs
git commit -m "feat: reject timesteps beyond the rk4 stability limit"
```

## Task 12: Release Hardening from the Final Review

The final review of Tasks 1 to 10 left a list of minor findings. This task closes the ones that need no owner decision. Choosing a licence stays with the project owner.

**Files:**
- Modify: `src/error.rs`
- Modify: `src/metrics.rs`
- Modify: `src/models/sir.rs`
- Modify: `src/models/seir.rs`
- Modify: `src/simulation/driver.rs`
- Modify: `src/solver/rk4.rs`
- Modify: `tests/sir_tests.rs`
- Modify: `tests/seir_tests.rs`
- Modify: `tests/numerical_tests.rs`
- Create: `CHANGELOG.md`
- Modify: `Cargo.toml`
- Modify: `docs/api_design.md`
- Modify: `implementation_plan.md`

**Interfaces:**
- `SimulationError` becomes `#[non_exhaustive]`, matching the documented promise that it may gain variants.
- Internal `basic_reproduction_number(beta_per_day: f64, gamma_per_day: f64) -> Result<f64, SimulationError>` in `src/metrics.rs`, used by `calculate_metrics` and called by both public functions before integration.
- No other API changes.

- [x] **Step 1: Write the tests**

Add tests:
- `sir_reports_metric_overflow_before_integrating` and `seir_reports_metric_overflow_before_integrating`: with an overflowing `beta / gamma` and settings that would exceed the step limit, the error is `MetricOverflow`, not `StepLimitExceeded`.
- `sir_surfaces_step_limit_clock_and_population_errors` and `seir_surfaces_step_limit_clock_and_population_errors`: `StepLimitExceeded`, `TimeStepTooSmall` (64-bit targets), and an overflowed-sum `InvalidPopulation` reach the caller through the public functions.
- `driver_rejects_non_finite_accepted_state`: `check_state` returns `NonFiniteNumericalState` naming the compartment.
- `driver_applies_negative_tolerance_at_its_boundary`: a value of exactly `-64 * EPSILON * max(N0, 1)` is accepted and counted, a value just below it is rejected, for a population above 1 and one below 1.
- `driver_applies_drift_threshold_near_its_boundary`: a relative drift of `0.9e-9` is accepted and `1.1e-9` is rejected.

Extend `rk4_rejects_non_finite_intermediate_stage` with a case that fails at `k3`. In `metrics_and_diagnostics_remain_finite_for_reference_scenarios`, replace the assertion `peak_infectious >= 10.0`, which cannot fail, with one that requires growth when beta is positive and an unchanged peak at day zero when beta is zero.

- [x] **Step 2: Run the tests**

Run: `cargo test`
Expected: the two `*_reports_metric_overflow_before_integrating` tests FAIL with `StepLimitExceeded`. The other new tests describe behaviour that already exists and pass; prove each can fail by breaking the code it covers, then restore the code.

- [x] **Step 3: Implement**

Extract `basic_reproduction_number` and call it in both public functions after the stability check and before `run_simulation`. Add `#[non_exhaustive]` to `SimulationError`. Add `CHANGELOG.md` with a 0.1.0 entry and include it in the package.

- [x] **Step 4: Run all quality checks**

Run:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
cargo doc --no-deps --all-features
cargo package --list
```

Expected: all pass.

- [x] **Step 5: Commit**

```bash
git add src tests docs CHANGELOG.md Cargo.toml implementation_plan.md
git commit -m "test: harden error paths and boundaries for release"
```

## Task 13: Correct the Stability Guard After Review

An independent review of Tasks 11 and 12 found that the Task 11 guard does not do its job for SEIR. Exposed and infected feed each other, so together they decay faster than sigma or gamma alone, and a run could pass the per-rate check and still return `Ok` with the infected growing where the model says they must decline (beta = 1.0, sigma = 0.25, gamma = 2.78, step 1: reported peak 370, correct peak 100). The review also noted that a limit of 2.785 sits where RK4 no longer damps at all.

**Files:**
- Modify: `src/error.rs`
- Modify: `src/validation.rs`
- Modify: `src/models/sir.rs`
- Modify: `src/models/seir.rs`
- Modify: `src/simulation/driver.rs`
- Modify: `tests/sir_tests.rs`
- Modify: `tests/seir_tests.rs`
- Modify: `README.md`, `CHANGELOG.md`, `docs/numerical_methods.md`, `docs/api_design.md`
- Modify: `TECHNICAL_SPEC.md`, `Product_Vision_&_MVP_Product_Requirements_Document.md`, `implementation_plan.md`

**Interfaces:**
- The error variant becomes `TimeStepTooLarge { max_step_days: f64 }`: the largest step accepted for the same parameters. A combined rate has no single parameter to name.
- Internal `stable_step(step_days: f64, fastest_rate_per_day: f64) -> Result<(), SimulationError>`, rejecting a product above `2.0`.
- `simulate_sir` passes `max(beta, gamma)`; `simulate_seir` passes `max(beta, (sigma + gamma) / 2 + sqrt(((sigma - gamma) / 2)^2 + sigma * beta))`. Both pass `min(step_days, duration_days)` as the step.

- [x] **Step 1: Write failing tests**

Replace the two Task 11 tests and add:
- `sir_rejects_step_too_large_for_stable_integration` and `seir_rejects_step_too_large_for_stable_integration`: each rate at 4 per day with a step of one day returns `TimeStepTooLarge { max_step_days: 0.5 }`, and the same rate runs with a step of 0.5. In SIR, a product of exactly 2 is accepted and one step leaves one third of the infected; 2.01 is rejected.
- `seir_rejects_step_too_large_for_the_coupled_exposed_and_infected_rate`: with all three rates at 1.5 and a step of one day, every single product is 1.5 but the run is rejected with `max_step_days` of 2/3; with all rates at 1 the coupled product is exactly 2 and the guard does not fire.
- `seir_rejects_subcritical_outbreak_that_a_coarse_step_would_inflate`: the review's input is rejected, and with a step of 0.5 the infected never exceed their initial value.
- `sir_stability_check_uses_the_step_actually_taken` and `seir_stability_check_uses_the_step_actually_taken`: a duration of one day with a nominal step of seven days is judged on the one-day step.

Tighten `driver_applies_drift_threshold_near_its_boundary` to 0.99e-9 and 1.01e-9.

- [x] **Step 2: Run tests to verify failure**

Run: `cargo test --no-fail-fast`
Expected: the six stability tests FAIL against the Task 11 guard.

- [x] **Step 3: Implement**

Change the variant, `stable_step`, and both call sites. Correct every statement about the guard in the README, changelog, docs, specification section 6.3, the product requirements, and this plan.

- [x] **Step 4: Verify**

Run the full quality checks from Task 12 Step 4. Then confirm the guard empirically: in a scratch crate outside the repository, run random configurations of both models against runs at 1/64 of the step and count accepted runs in which a compartment moves in the wrong direction or the final incidence is off by more than 0.1.
Expected: all checks pass, and the count is zero. Measured: 226,590 accepted SIR runs and 107,633 accepted SEIR runs out of 250,000 configurations each, none flagged. Before the fix the same measure found such runs from a product of about 2.2 (SEIR) and 2.4 (SIR).

- [x] **Step 5: Commit**

```bash
git add src tests docs README.md CHANGELOG.md TECHNICAL_SPEC.md "Product_Vision_&_MVP_Product_Requirements_Document.md" implementation_plan.md
git commit -m "fix: guard the fastest coupled rate and lower the stability limit"
```

---

## Execution Notes

These record where the implementation of Tasks 1 to 11 departed from the task text above. The task text is kept as written.

- **Document locations.** The specification, product requirements, and this plan live in the repository root, not under `docs/`.
- **Tests of internal code (Tasks 3 to 6).** The plan placed them in `tests/validation_tests.rs`, `tests/rk4_tests.rs`, `tests/numerical_tests.rs`, and `tests/metrics_tests.rs`. Integration tests cannot see `pub(crate)` items, and the specification keeps the solver, driver, validators, and metric calculation private. Those tests are unit tests inside `src/models/sir.rs`, `src/models/seir.rs`, `src/solver/rk4.rs`, `src/simulation/driver.rs`, and `src/metrics.rs`, with the test names the plan lists. Run them with `cargo test --lib`. `tests/rk4_tests.rs` and `tests/metrics_tests.rs` do not exist; `tests/numerical_tests.rs` was created in Task 9 with public-API tests.
- **Error variants (Task 2).** No extra numeric context fields were added. `stage` and `parameter` are `&'static str`.
- **Solver errors (Task 4).** A non-finite step returns `InvalidTimeStep`. A non-finite stage input or slope returns `NonFiniteIntermediate` with stage `"k1"` to `"k4"`; a non-finite combined result uses stage `"final"`, because the solver has no compartment names.
- **Derives (Tasks 5 and 6).** `Observation`, `SimulationDiagnostics`, and `EpidemicMetrics` also derive `Copy`, as specification section 5.1 asks.
- **Time handling (Task 5).** Observation times are `index * step_days`, with the last pinned to the duration, not a running sum. The step count is `duration / step` reduced by four ulps and then rounded up, so a duration that is a whole number of steps does not gain a step of rounding-error length (2.7 / 0.3 gives 9 steps). `TimeStepTooSmall` is detected before integration as `duration + step == duration`, after the step-limit check. The allocation-failure test calls the reservation function with a capacity of `usize::MAX`.
- **Metrics (Task 6).** An empty trajectory returns `MetricOverflow`; it cannot occur through the public API.
- **Model tests (Tasks 7 and 8).** `sir_first_step_matches_model_equations` and `seir_first_step_matches_model_equations` were added after the final review, because no listed test would have caught a wrong transmission term.
- **Convergence (Task 9).** The fixtures use the nonlinear models against a run at 1/64 day. Measured error-reduction factors: 13.9 and 14.9 for SIR, 15.6 and 15.8 for SEIR.
- **Release work (Task 10).** Criterion 0.7 is used with default features off except `cargo_bench_support`. `Cargo.toml` has an `include` list so only the crate and its documentation are packaged. The benchmark baseline is in `docs/architecture.md`; its three workloads use different timesteps so that each covers 120 to 200 simulated days. The Rust examples in `README.md` and `docs/api_design.md` run as doctests. No licence had been chosen at the time; improvement plan 001 later added the MIT licence.
- **Self-review numbering.** In the checklist below, the model equations are Tasks 7 and 8 and the non-finite `R0` test is Task 6.

---

## After This Plan

The task text and the API contracts above are kept as written. Later changes came through improvement plans 001 to 012. The advisor wrote them in `plans/`, executor agents carried them out, and the advisor reviewed each one. `plans/` is in the working copy and is not committed. Up to release 0.1.0 these plans changed the contracts as follows:

| Plan | Change to the contracts above |
|---|---|
| 001 | MIT licence, package metadata, GitHub Actions CI on stable and on Rust 1.85.0 |
| 002 | `tests/stability_evidence_tests.rs`, the ignored measurement behind the stability limit |
| 003 | `SimulationSettings` gains `record_every`, with `with_record_every` and `record_every()`, and the error `InvalidRecordInterval`. The trajectory holds every N-th step plus the initial and final states. |
| 004 | `estimate_sir_step_error`, `estimate_seir_step_error` and `StepErrorEstimate`, in `src/accuracy.rs` |
| 005 | `EpidemicMetrics` gains `herd_immunity_threshold` and `final_size_fraction`. `SirParameters` and `SeirParameters` gain `effective_reproduction_number`, which is model-derived and is not `Rt`. |
| 006 | Examples `outbreak_curve` and `r0_sweep` |
| 007 | Design spike for SIS, SIRS and SEIRS (`design/007-waning-immunity-models.md`); no code change |
| 008 | `Observation` gains `cumulative_infections`, integrated as `dC/dt = beta * S * I / N`. The driver treats the last array entry as that accumulator, outside the population checks, and the cumulative metrics come from it. |
| 009 | `EpidemicMetrics` becomes `#[non_exhaustive]`; `herd_immunity_threshold` and `final_size_fraction` become `Option<f64>`; the peak time ignores roundoff at a level peak |
| 010 | `simulate_sis`, `SisParameters`, `SisState`, `estimate_sis_step_error`, and the metric `endemic_prevalence` |
| 011 | `simulate_sirs`, `SirsParameters` (reusing `SirState`), `estimate_sirs_step_error` |
| 012 | `simulate_seirs`, `SeirsParameters` (reusing `SeirState`), `estimate_seirs_step_error` |

The test suite grew from the 80 tests at the end of Task 13 to 144 passed and 2 ignored at release 0.1.0, tagged `v0.1.0` on November 13, 2026.

---

## Self-Review Checklist

### Spec coverage

- SIR and SEIR model equations: Tasks 6 and 7.
- Shared RK4 solver: Task 4.
- Input and configuration validation: Tasks 2 and 3.
- Fixed timestep and final partial step: Task 5.
- Maximum-step limit and allocation protection: Task 5.
- Time-series output: Task 5.
- Population conservation and negativity diagnostics: Tasks 5 and 9.
- Epidemiological metrics: Task 6.
- Analytical reference scenarios: Tasks 7, 8, and 9.
- Numerical convergence: Task 9.
- API documentation and examples: Task 10.
- Performance benchmarks: Task 10.
- Package verification: Task 10.
- Timestep stability guard: Task 11.
- Error-path, boundary, and release-note hardening: Task 12.
- Corrected stability guard (coupled SEIR rate, limit of 2): Task 13.

### Type consistency

- `SimulationSettings` is passed by value to the public simulation functions.
- Public result types are generic over the concrete public state type.
- SIR compartment order is `[S, I, R]`; SEIR compartment order is `[S, E, I, R]`.
- The infectious compartment index is 1 for SIR and 2 for SEIR; susceptible is index 0 for both.
- Both public functions return `Result<SimulationResult<...>, SimulationError>`.
- Metric names and units are consistent across the API and tests.

### Review focus coverage

- Overflowed initial population: Task 3 tests.
- Time-step unable to advance clock: Task 5 tests.
- Non-finite `R0`: Task 8 tests.
- Final partial step and maximum-step boundary: Task 5 and Task 9 tests.
- Non-finite RK4 intermediate stage: Task 4 tests.

### Technical notes to resolve during implementation

- Keep allocation failure tests deterministic through a test seam; never attempt to exhaust real system memory.
- Use checked arithmetic for `required_steps + 1` before reserving trajectory storage.
- Ensure the step-count calculation and actual integration loop use the same endpoint convention to avoid off-by-one differences.
- The specified convergence ratio is a conservative release criterion, not a proof of fourth-order behavior; retain a solver-level exponential-decay test for the expected RK4 order.
- A tolerated tiny negative roundoff value is counted but not clipped. Document that downstream consumers may observe such a value within the stated tolerance.
- If the existing repository uses different directory conventions, adapt paths while keeping the module boundaries and test ownership intact.

---

## Execution Handoff

This plan was executed in order, with one focused test cycle and commit per task, after the product and technical specification documents were in the repository and the plan had been reviewed. It is complete; new work is planned in `plans/`.
