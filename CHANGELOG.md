# Changelog

All notable changes to this project are recorded in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - unreleased

The first version. It has not been published yet.

### Added

- `simulate_sir` and `simulate_seir`: deterministic SIR and SEIR models with
  frequency-dependent transmission, a closed population, and constant
  parameters.
- `simulate_sis` and `estimate_sis_step_error`: the SIS model, in which
  recovery gives no immunity. Its `EpidemicMetrics` report
  `endemic_prevalence` instead of a herd-immunity threshold and final size.
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
