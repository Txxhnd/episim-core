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
- A fixed-step classical fourth-order Runge-Kutta solver shared by both
  models. The final step is shortened so a run ends exactly at the requested
  duration.
- `SimulationSettings` with a validated duration and timestep, and a step
  limit that defaults to 1,000,000 and can be changed with `with_max_steps`.
- `SimulationResult` with the full trajectory, `EpidemicMetrics` (basic
  reproduction number, peak infectious population and its time, cumulative
  infections, cumulative incidence fraction), and `SimulationDiagnostics`
  (steps executed, population drift, tolerated negative values).
- `SimulationError`, a typed, non-exhaustive error for invalid parameters,
  invalid initial states, unusable settings, numerical failures, resource
  limits, and non-finite metrics.
- Input checks before integration: parameter ranges, a finite and positive
  initial population, a timestep within the solver's stability limit (any
  rate times the step at most 2.785), and a finite `beta / gamma`.
- Checks after every step: finite values, no materially negative
  compartment, and relative population drift of at most `1e-9`. Values are
  never clipped.
- Examples `basic_sir` and `basic_seir`, Criterion benchmarks, and
  documentation of the models, the numerical method, the architecture, and
  the API.
