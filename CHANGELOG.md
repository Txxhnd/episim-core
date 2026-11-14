# Changelog

All notable changes to this project are recorded in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- `TransmissionSchedule::multiplier_at` agrees with the multiplier the run
  applies at a recorded time on a decimal step grid: a change on day 0.9 now
  counts at the recorded time `3 * 0.3 = 0.8999999999999999`, within the
  same four-ulp tolerance as the step-grid check (issue #4).
- A schedule's multiplier is found by binary search, so a schedule with a
  change at every step no longer costs a scan of all changes per step
  (40,000 changes: 2.4 s before, 27 ms after; issue #5).
- The SEIR stability limit computes the coupled rate on rates scaled by a
  power of two, so it no longer overflows or underflows when every rate is
  scaled by a very large or very small factor. Ordinary inputs give the same
  limit to the last bit (issue #6).
- `final_size_fraction` solves for the newly infected fraction directly,
  with a series for `exp(-x) - 1 + x` near the threshold, so tiny seeds near
  `R0 = 1` are no longer lost to rounding (issue #7).

### Added

- `repository` in the package metadata.

## [0.2.0] - 2026-11-14

Transmission schedules for every model, a stricter SIRS stability limit, and
fixes to three metrics and the step-error estimate. Under Cargo's rules for
0.x versions this release may change the public API; the one breaking change
is marked below. Tagged `v0.2.0`; not published to crates.io.

### Added

- `TransmissionSchedule`, a piecewise-constant multiplier on the
  transmission rate, attached with
  `SimulationSettings::with_transmission_schedule`, for interventions such
  as lockdowns. Works with every model and with the step-error estimates.
  The stability limit uses the largest transmission rate in force.
- `SimulationError::InvalidSchedule` and
  `SimulationError::ScheduleChangeOffStepGrid { day }`.
- A GitHub Actions workflow, "Stability evidence", that runs the release-mode
  stability evidence test every Monday and on demand.

### Changed

- **Breaking:** `SimulationSettings` is `Clone` but no longer `Copy`, because
  it can hold a transmission schedule. Clone settings to reuse them after a
  call.
- Under a transmission schedule that changes the multiplier,
  `final_size_fraction` is `None`, and `endemic_prevalence` uses the
  multiplier in force at the end of the run.
- The SIRS stability limit is `beta + gamma + omega`, as for SEIRS, instead
  of `max(beta + gamma, omega)`. The old limit accepted steps that ended a
  quarter of the population away from the exact solution when `omega` was
  close to `beta + gamma`, and 1,215 of the evidence's accepted runs failed
  with `NegativeCompartment`. Steps it accepted may now return
  `TimeStepTooLarge`; at typical rates the largest step barely changes
  (4.9 days instead of 5 at `beta = 0.3`, `gamma = 0.1`, `omega = 0.01`).

### Fixed

- The step-error estimates split every step the run actually takes, including
  a final step shortened to the duration. They used to halve the nominal
  step, so a run whose only step was shorter than the nominal one reported an
  error of 0 (issue #1).
- SEIRS `endemic_prevalence` is computed from rate ratios, so it no longer
  overflows to 0 or loses precision when every rate is scaled by a very large
  or very small factor (issue #2).
- `final_size_fraction` is solved in fractions of the population, so it no
  longer overflows for populations near `f64::MAX` (issue #3).

- `endemic_prevalence` is 0 for SIS, SIRS and SEIRS runs that start with
  nobody infected (for SEIRS, nobody exposed or infected). It used to report
  the endemic level, `1 - 1/R0` for SIS, although such a run never leaves
  the disease-free state.

## [0.1.0] - 2026-11-13

The first version. Tagged `v0.1.0`; not published to crates.io.

### Added

- `simulate_sir` and `simulate_seir`: deterministic SIR and SEIR models with
  frequency-dependent transmission, a closed population, and constant
  parameters.
- `simulate_sis` and `estimate_sis_step_error`: the SIS model, in which
  recovery gives no immunity. Its `EpidemicMetrics` report
  `endemic_prevalence` instead of a herd-immunity threshold and final size.
- `simulate_sirs` and `estimate_sirs_step_error`: SIR with immunity that
  wanes at rate `omega_per_day`, using `SirState`. It reports
  `endemic_prevalence` and a herd-immunity threshold, and no final size.
- `simulate_seirs` and `estimate_seirs_step_error`: SEIR with immunity that
  wanes, using `SeirState`, with the same metrics as SIRS.
- A fixed-step classical fourth-order Runge-Kutta solver shared by all
  models. The final step is shortened so a run ends exactly at the requested
  duration.
- `SimulationSettings` with a validated duration and timestep, and a step
  limit that defaults to 1,000,000 and can be changed with `with_max_steps`.
- `SimulationSettings::with_record_every`, to store only every N-th step of
  a run while still integrating and checking every step.
- `estimate_sir_step_error` and `estimate_seir_step_error`, which estimate
  the error the timestep introduces by comparing a run with the same run at
  half the step.
- Model-derived threshold quantities: `effective_reproduction_number` on
  `SirParameters` and `SeirParameters` (`R0 * S / N` in any state), and the
  `herd_immunity_threshold` and `final_size_fraction` metrics. Both are
  `Option<f64>`, because models with waning immunity will not report them.
  `EpidemicMetrics` is `#[non_exhaustive]`.
- `Observation::cumulative_infections`: new infections since day zero,
  integrated with the compartments as `dC/dt = beta * S * I / N`.
  `cumulative_infections_since_start` is its value at the final time.
- `SimulationResult` with the full trajectory, `EpidemicMetrics` (basic
  reproduction number, peak infectious population and its time, cumulative
  infections, cumulative incidence fraction), and `SimulationDiagnostics`
  (steps executed, population drift, tolerated negative values).
- `SimulationError`, a typed, non-exhaustive error for invalid parameters,
  invalid initial states, unusable settings, numerical failures, resource
  limits, and non-finite metrics.
- Input checks before integration: parameter ranges, a finite and positive
  initial population, a timestep small enough for stable integration (the
  model's fastest decay rate times the step at most 2), and a finite
  `beta / gamma`.
- Checks after every step: finite values, no materially negative
  compartment, and relative population drift of at most `1e-9`. Values are
  never clipped.
- Examples `basic_sir` and `basic_seir`; teaching examples `outbreak_curve`
  (one outbreak as CSV, one row per day) and `r0_sweep` (peak, timing and
  final size across R0 values); Criterion benchmarks; and
  documentation of the models, the numerical method, the architecture, and
  the API.
