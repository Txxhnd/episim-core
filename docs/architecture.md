# Architecture

`episim-core` is a single library crate with no runtime dependencies. The
public API exposes one function per model. Internally all models share the
same validation helpers, RK4 solver, simulation driver, and metric
calculation.

## Data flow

```text
simulate_sir / _seir / _sis / _sirs / _seirs   src/models/
  -> simulate (shared pipeline)       src/models/mod.rs
  1. validate parameters              src/validation.rs
  2. validate the initial state and
     compute the fixed population N
     apply the stability limit and check R0
  3. run_simulation                   src/simulation/driver.rs
       - count steps, enforce the step limit
       - reserve the trajectory
       - for each step:
           look up the step's transmission multiplier
           rk4_step                   src/solver/rk4.rs
           check the new state
           record an observation
  4. calculate_metrics                src/metrics.rs
  5. assemble SimulationResult        src/result.rs
```

## Modules

| Module | Visibility | Responsibility |
|---|---|---|
| `lib.rs` | public | Crate attributes and re-exports. |
| `error.rs` | type is public | `SimulationError` and its `Display` text. |
| `settings.rs` | type is public | `SimulationSettings`: validated duration, step, step limit, recording interval, and transmission schedule, with the per-step multiplier lookups. |
| `schedule.rs` | type is public | `TransmissionSchedule`: a validated piecewise-constant multiplier on `beta`. |
| `validation.rs` | private | Range checks for rates and for the initial population. |
| `models/mod.rs` | private | The private `Model` trait and the shared `simulate` and `estimate` pipeline. |
| `models/sir.rs`, `models/seir.rs`, `models/sis.rs`, `models/sirs.rs`, `models/seirs.rs` | types and functions are public | Parameter and state types, each model's `Model` implementation (equations, stability rate, model-specific metrics), and the public `simulate_*` and `estimate_*_step_error` functions, which call the shared pipeline. |
| `solver/rk4.rs` | private | One RK4 step over a fixed-size array, with finiteness checks. |
| `simulation/driver.rs` | private | Step counting, time handling, trajectory storage, and per-state checks. |
| `metrics.rs` | type is public | `EpidemicMetrics` and its calculation from a trajectory. |
| `accuracy.rs` | type is public | `StepErrorEstimate` and the comparison of a run with its half-step run. |
| `result.rs` | types are public | `Observation`, `SimulationDiagnostics`, `SimulationResult`. |

The modules themselves are private. Public types are re-exported from the
crate root, so users never name a module path.

## Separation of concerns

- **Models know equations, not integration.** A model supplies a derivative
  function over a fixed-size array, a conversion from that array back to its
  state type, and compartment names for error messages.
- **The solver knows nothing about epidemiology.** `rk4_step` is generic
  over the array length and the derivative function. It does not judge
  negative values; it only rejects non-finite ones.
- **The driver owns the loop.** It decides the time grid, calls the solver,
  applies the negativity and conservation checks, and collects diagnostics.
  It does not compute epidemiological metrics.
- **Metrics are computed afterwards** from the recorded trajectory, through
  an accessor that picks the infectious value out of a state, and each
  model then fills in the metrics only it reports.

The internal generics (`const N: usize` for the number of compartments, and
closures for the derivative and conversions) are resolved at compile time.
None of them appears in the public API.

## Design decisions

| Decision | Reason |
|---|---|
| Concrete public functions, no public model trait | A few fixed models do not justify a plugin interface; internals can still be shared. |
| A private `Model` trait behind one shared pipeline | Validation, the stability limit, the `R0` check, integration and metrics run in one place, so a change to the pipeline, such as a transmission schedule, is made once rather than in five copies. The trait is not exported. |
| `f64` everywhere | Compartments are expected counts, not integers. |
| Population `N` fixed from the initial state | Recomputing it each step would feed numerical drift back into the transmission term. |
| Fixed-size arrays in the solver | No heap allocation inside a step. |
| Every step stored by default | Simple, and what the metrics and callers need; memory is linear in the number of recorded steps. `with_record_every` stores fewer. |
| Step limit with a default of 1,000,000 | Guards against accidentally huge runs. |
| Fallible trajectory reservation | An oversized run returns an error and does not abort the process. |
| Errors, never clipping | A clipped state would hide a step size that is too large. |

## Tests

| Location | Covers |
|---|---|
| Unit tests in `src/` | Internal code that integration tests cannot reach: validators, `rk4_step`, the driver, and metric calculation. |
| `tests/validation_tests.rs` | `SimulationSettings` construction. |
| `tests/sir_tests.rs`, `tests/seir_tests.rs`, `tests/sis_tests.rs`, `tests/sirs_tests.rs`, `tests/seirs_tests.rs` | The public functions, including closed-form reference cases. |
| `tests/numerical_tests.rs` | Population drift and finite metrics across reference scenarios. |
| `tests/convergence_tests.rs` | Error reduction under timestep refinement. |
| Doctests | The examples in the API docs, `README.md`, and `docs/api_design.md`. |

## Cost

For a fixed model, run time and memory are both linear in the number of
steps. Each observation takes 40 bytes for SIR and 48 bytes for SEIR (the
time, the cumulative infections, and the compartments), so a 1,000,000-step
run stores about 40 MB or 48 MB when every step is recorded.

## Performance baseline

Measured with `cargo bench` (Criterion 0.7, 10 samples per case) on
2026-11-13, for release 0.1.0. These numbers are a baseline for comparison on the same machine,
not a performance target or guarantee.

| | |
|---|---|
| CPU | Intel Core i5-5350U, 1.80 GHz |
| Memory | 8 GB |
| Operating system | macOS 12.7.6 |
| Compiler | rustc 1.99.1 |
| Build | Cargo `bench` profile (release, `opt-level = 3`), no custom flags |

Workloads: SIR with `beta = 0.3`, `gamma = 0.1`, initial state
`(990, 10, 0)`; SEIR with `beta = 0.5`, `sigma = 0.2`, `gamma = 0.1`,
initial state `(980, 10, 10, 0)`. Each run includes validation, the fallible
trajectory reservation, integration, and metric calculation.

| Steps | Timestep (days) | Duration (days) | SIR median | SEIR median |
|---:|---:|---:|---:|---:|
| 10,000 | 1/64 | 156.25 | 0.82 ms | 0.81 ms |
| 100,000 | 1/512 | 195.31 | 9.3 ms | 9.4 ms |
| 1,000,000 | 1/8192 | 122.07 | 90 ms | 88 ms |

That is roughly 80 to 95 nanoseconds per step. Criterion's interval for the
1,000,000-step SIR case was 86 to 96 ms. On 2026-11-12, before observations
carried the cumulative infections, the same workloads took 0.69, 9.2 and
88 ms (SIR) and 0.87, 11.0 and 108 ms (SEIR). That run's interval for the
1,000,000-step SIR case was 78 to 105 ms, so treat differences of less than
about 20% as noise on this machine; the extra integrated value made no
measurable difference. Memory use was not measured; the trajectory accounts
for 40 bytes (SIR) or 48 bytes (SEIR) per recorded step. SIS, SIRS and SEIRS
are not benchmarked.

The timesteps differ between rows so that each run covers a similar span of
an active outbreak. An earlier run that kept the timestep at 1/64 day and
stretched the 1,000,000-step case to 15,625 simulated days took 560 ms for
SIR and 788 ms for SEIR. Long after the outbreak ends, the infectious
compartment decays to a subnormal `f64` value and stays there, and
arithmetic on subnormal values is much slower on this CPU. Very long
horizons can therefore cost several times more per step than the table
suggests.
