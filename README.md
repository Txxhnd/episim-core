# episim-core

A small Rust library for deterministic epidemiological simulation. It
integrates the SIR and SEIR compartmental models with a fixed-step classical
fourth-order Runge-Kutta (RK4) solver and returns the full trajectory, a few
model-derived metrics, and numerical diagnostics.

- Pure Rust, no runtime dependencies, no `unsafe` code.
- The simulation functions and the settings constructors return a `Result`
  with a typed error; invalid input does not panic.
- Numerical problems the library detects are reported, never patched over:
  compartment values are not clipped, and population conservation is checked
  after every step. A fixed-step solver cannot detect every inaccurate
  result; see the limitations below.

`episim-core` is a simulation engine, not a forecasting product. Its output
describes what the stated model does under the stated assumptions. It is not
calibrated against, or an estimate of, any real-world outbreak.

## Quick start

The crate is not published on crates.io yet, so depend on it by path or git:

```toml
[dependencies]
episim-core = { path = "../episim-core" }
```

```rust
use episim_core::{SimulationError, SimulationSettings, SirParameters, SirState, simulate_sir};

fn main() -> Result<(), SimulationError> {
    // Rates are per day; compartments are population counts.
    let parameters = SirParameters {
        beta_per_day: 0.3,
        gamma_per_day: 0.1,
    };
    let initial_state = SirState {
        susceptible: 990.0,
        infected: 10.0,
        recovered: 0.0,
    };
    // Simulate 160 days with a fixed step of a quarter day.
    let settings = SimulationSettings::new(160.0, 0.25)?;

    let result = simulate_sir(parameters, initial_state, settings)?;

    for observation in result.trajectory.iter().step_by(80) {
        let state = observation.state;
        println!(
            "day {:>3}: S = {:6.1}, I = {:6.1}, R = {:6.1}",
            observation.time_days, state.susceptible, state.infected, state.recovered
        );
    }
    println!(
        "peak of {:.1} infectious on day {}",
        result.metrics.peak_infectious, result.metrics.peak_infectious_time_days
    );
    Ok(())
}
```

SEIR works the same way through `simulate_seir`, `SeirParameters`, and
`SeirState`. Two runnable examples print a fuller summary:

```sh
cargo run --example basic_sir
cargo run --example basic_seir
```

## What a simulation returns

`SimulationResult` has three public fields.

| Field | Contents |
|---|---|
| `trajectory` | One `Observation` (time in days plus the model state) for day zero and for every recorded integration step: every step by default, or every N-th step with `SimulationSettings::with_record_every`. The last one is always the final state, exactly at the requested duration. |
| `metrics` | Basic reproduction number `R0 = beta / gamma` and the herd-immunity threshold `1 - 1/R0`, peak infectious population and the day it is first reached, new infections since the start `S(0) - S(T)` and that number as a fraction of the initial population, and the fraction the model infects by the time the epidemic is over (its final size). |
| `diagnostics` | Steps executed, the largest absolute and relative population drift, and how many recorded values were negative within roundoff tolerance. |

Peak metrics are read from the recorded observations, so their time
resolution is the recording interval, which is the timestep by default.

The effective reproduction number at any point of a run, `R0 * S / N`, comes
from the parameters:

```rust
use episim_core::{SimulationError, SimulationSettings, SirParameters, SirState, simulate_sir};

let parameters = SirParameters { beta_per_day: 0.3, gamma_per_day: 0.1 };
let initial_state = SirState { susceptible: 990.0, infected: 10.0, recovered: 0.0 };
let result = simulate_sir(parameters, initial_state, SimulationSettings::new(160.0, 0.25)?)?;

for observation in result.trajectory.iter().step_by(40) {
    let effective_r = parameters.effective_reproduction_number(&observation.state)?;
    println!("day {:>3}: effective R = {effective_r:.2}", observation.time_days);
}
# Ok::<(), SimulationError>(())
```

All of these are properties of the model. They are not estimates from case
data: in particular, the model's effective R is not the `Rt` that
public-health agencies estimate from surveillance.

## Model assumptions

Both models assume:

1. A closed population: no births, deaths, or migration.
2. Homogeneous mixing.
3. Parameters that stay constant for the whole simulation.
4. Permanent immunity after recovery.
5. Transmission only from the infected compartment.
6. Deterministic dynamics, with fractional "expected" population counts.

Transmission is frequency-dependent: new infections occur at the rate
`beta * S * I / N`, where `N` is the total of the initial state and stays
fixed. In SEIR the exposed compartment is a latent, non-infectious stage; it
is not the clinical incubation period.

## Scope and limitations

Version 0.1.0 contains the SIR and SEIR models and one solver. It has no
command-line interface, serialization, stochastic or agent-based models, age
structure, interventions, time-varying parameters, parameter fitting, or
estimation of `Rt` from case data.

The solver uses a fixed step and has no error control. Choose a step that is
small compared with `1 / beta`, `1 / sigma`, and `1 / gamma`. A run is
rejected with `TimeStepTooLarge` if the model's fastest decay rate times the
step exceeds 2, because RK4 damps too weakly beyond that and can turn decay
into growth; the error reports the largest step that would be accepted. A run
also stops with an error if a value becomes non-finite, a compartment goes
materially negative, or the population drifts. Below that limit a coarse
step is still inaccurate, and a plain run does not tell you so. To find out
how much error the step introduces, ask for an estimate:

```rust
use episim_core::{
    SimulationError, SimulationSettings, SirParameters, SirState, estimate_sir_step_error,
};

let parameters = SirParameters {
    beta_per_day: 0.3,
    gamma_per_day: 0.1,
};
let initial_state = SirState {
    susceptible: 990.0,
    infected: 10.0,
    recovered: 0.0,
};
let settings = SimulationSettings::new(160.0, 0.25)?;

let estimate = estimate_sir_step_error(parameters, initial_state, settings)?;
// With a quarter-day step no compartment is off by more than 0.001 people.
assert!(estimate.max_absolute_error < 1e-3);
# Ok::<(), SimulationError>(())
```

`estimate_seir_step_error` does the same for SEIR. The estimate reruns the
simulation at half the step and compares the two, so it costs about three
times as much as one run. It measures the error from the step size only,
not how well the model matches reality.
Every recorded observation is kept in memory (32 bytes each for SIR, 40 for
SEIR). To use a small step without storing every step, record every N-th
one:

```rust
use episim_core::{SimulationError, SimulationSettings};

// 16,000 steps of 0.01 days, recorded once per simulated day.
let settings = SimulationSettings::new(160.0, 0.01)?.with_record_every(100)?;
assert_eq!(settings.record_every(), 100);
# Ok::<(), SimulationError>(())
```

Every step is still integrated and checked. A run is limited to 1,000,000
steps unless you raise the limit with `SimulationSettings::with_max_steps`.

## Documentation

- [Mathematical models](docs/mathematical_models.md): equations, units,
  assumptions, and metric definitions.
- [Numerical methods](docs/numerical_methods.md): the RK4 scheme, time
  handling, and the checks applied to every step.
- [Architecture](docs/architecture.md): module layout, data flow, and the
  performance baseline.
- [API design](docs/api_design.md): public types, conventions, and errors.

API reference: `cargo doc --no-deps --open`.

## Development

Requires Rust 1.85 or newer.

```sh
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets
cargo test --doc
cargo doc --no-deps --all-features
cargo bench
```

## License

Licensed under the [MIT License](LICENSE).
