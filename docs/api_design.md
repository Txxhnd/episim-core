# API design

The public API is small and concrete: two functions, the plain data types
they take and return, and one error type. Everything else in the crate is
private.

## Public items

| Item | Kind | Purpose |
|---|---|---|
| `simulate_sir` | function | Run an SIR simulation. |
| `simulate_seir` | function | Run an SEIR simulation. |
| `estimate_sir_step_error`, `estimate_seir_step_error` | function | Estimate the error the timestep introduces, by comparing with a run at half the step. |
| `SirParameters`, `SeirParameters` | struct | Rate parameters, per day. |
| `SirState`, `SeirState` | struct | Compartment values; used for the initial state and for each observation. |
| `SimulationSettings` | struct | Duration, timestep, and step limit. |
| `SimulationResult<S>` | struct | Trajectory, metrics, and diagnostics. |
| `Observation<S>` | struct | A time in days, the cumulative infections since day zero, and the state at that time. |
| `EpidemicMetrics` | struct | Model-derived summary metrics. |
| `SirParameters::effective_reproduction_number`, `SeirParameters::effective_reproduction_number` | method | Model-derived effective reproduction number `R0 * S / N` in a given state. |
| `SimulationDiagnostics` | struct | Step count, population drift, tolerated negative values. |
| `StepErrorEstimate` | struct | Estimated largest absolute and relative error of a run. |
| `SimulationError` | enum | Every failure the library can report. |

```rust,ignore
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

## Conventions

- **Concrete functions, no public model trait.** Each model has its own
  function and types. The solver and driver they share are internal. A
  public extension trait is deferred until there is a real need for one.
- **Arguments by value.** Parameter, state, and settings types are small and
  `Copy`.
- **Units in names.** Rates end in `_per_day` and times in `_days`.
- **Plain data for parameters and states.** Their fields are public and are
  validated when a simulation starts.
- **Validated settings.** `SimulationSettings` has private fields.
  `SimulationSettings::new(duration_days, step_days)` rejects values that
  are not finite and positive, and sets the step limit to 1,000,000.
  `with_max_steps` changes the limit and rejects zero.
  `with_record_every` records only every N-th step and rejects zero.
- **No panics on bad input.** The simulation functions and the settings
  constructors return `Result`.
- **No hidden state.** There are no globals and no randomness; the same
  inputs give the same result in the same environment.
- **No silent correction.** The library does not clip, rescale, or otherwise
  repair a numerical state. It reports the problem.

## Errors

`SimulationError` implements `Debug`, `Display`, `Clone`, `PartialEq`, and
`std::error::Error`. Its variants fall into five groups.

| Group | Variants | When |
|---|---|---|
| Invalid configuration | `InvalidParameter { name }`, `InvalidDuration`, `InvalidTimeStep`, `InvalidMaximumSteps`, `InvalidRecordInterval`, `TimeStepTooLarge { max_step_days }` | A rate, duration, timestep, step limit, or recording interval is outside its valid range, or the timestep is too large for stable integration (the model's fastest decay rate times the step is above 2; `max_step_days` is the largest step that would be accepted). `name` is the field name, for example `"gamma_per_day"`. |
| Invalid initial population | `InvalidInitialState { compartment }`, `InvalidPopulation` | A compartment is negative or not finite, or the total is zero or not finite. |
| Numerical instability | `NonFiniteIntermediate { stage }`, `NonFiniteNumericalState { compartment }`, `NegativeCompartment { compartment }`, `PopulationConservationViolation` | The integration produced a value the library will not accept. |
| Resource limits | `StepLimitExceeded`, `TimeStepTooSmall`, `AllocationFailure` | The run needs more steps than allowed, a step too small for `f64` time, or more memory than could be reserved. |
| Invalid metrics | `MetricOverflow` | A metric is not finite, for example `beta / gamma` with an extremely small `gamma`. That case is reported before integration starts. |

`NonFiniteNumericalState` is the error for a non-finite value in an accepted
state. The solver checks its own output first, so in the current version a
non-finite value is reported as `NonFiniteIntermediate`.

Handling an error by group:

```rust
use episim_core::{SimulationError, SimulationSettings, SirParameters, SirState, simulate_sir};

let parameters = SirParameters { beta_per_day: 0.3, gamma_per_day: 0.0 };
let initial_state = SirState { susceptible: 990.0, infected: 10.0, recovered: 0.0 };
let settings = SimulationSettings::new(160.0, 0.25).expect("valid settings");

match simulate_sir(parameters, initial_state, settings) {
    Ok(result) => println!("{} observations", result.trajectory.len()),
    Err(SimulationError::InvalidParameter { name }) => println!("fix parameter {name}"),
    Err(other) => println!("simulation failed: {other}"),
}
```

## Stability

This is a 0.1 release. The public API may still change before 1.0.
`SimulationError` is marked `#[non_exhaustive]`: it may gain variants, so
the compiler requires a wildcard arm when you match on it.
`EpidemicMetrics` is `#[non_exhaustive]` too: it may gain fields, so code
outside the crate reads its fields but cannot build one with a struct
literal. Metrics that only some models define, such as
`herd_immunity_threshold` and `final_size_fraction`, are `Option<f64>`.

Not part of the API in this version: estimation of `Rt` from case data,
serialization, and user-defined models.
