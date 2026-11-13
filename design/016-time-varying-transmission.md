# Design 016: Time-varying transmission

Design study for version 0.2.0, on branch `spike/016-time-varying-transmission`
from `1681ff0`. The prototype on that branch is not for merging; this
document, and the build plans written from it, move on.

The prototype is in these files:

| File | What it is |
|---|---|
| `src/solver/rk4.rs` | `rk4_step_timed`: RK4 whose derivative also receives the stage's position in the step (0, 1/2 or 1). `rk4_step` now delegates to it; all 156 tests pass unchanged. |
| `src/spike.rs` | `#[doc(hidden)]` `integrate`, a bare fixed-step loop whose derivative receives the step's start time, the step length and the stage fraction. No checks. |
| `tests/time_varying_spike_tests.rs` | Measurements M1 to M4 below. M3 is ignored; run it with `--release -- --include-ignored`. |

All commands run from the repository root of that branch:

```sh
cargo test --release --test time_varying_spike_tests -- --include-ignored --nocapture
```

## Summary

1. **Describe changes as a piecewise-constant multiplier on `beta`**, for
   example "contacts at 40% from day 20 to day 50". That covers lockdowns,
   school closures and reopenings. It keeps `R0 = beta / gamma` meaningful as
   the baseline, and makes the largest `beta` in force exact and easy to find.
2. **Read the multiplier at the start of each step.** The solver does not
   change: the driver gives the derivative the step's start time.
3. **Require every change to fall on the step grid.** On the grid the run
   keeps fourth-order accuracy. Off it, every way of reading the multiplier
   loses it, with errors of 4 to 19 people.
4. **Stability limit:** the model's existing fastest-rate formula, with
   `beta` replaced by the largest `beta` in force during the run. Measured:
   no flagged run in 250,000 SIR and SEIR configurations with random
   schedules. Using only the `beta` at day 0 flags 97 and 57.
5. **Metrics:** `R0` and the herd-immunity threshold stay the baseline
   values. The final size is `None` once the schedule changes anything. The
   endemic prevalence uses the multiplier in force at the end of the run.
6. **Share the five models' pipeline first**, as a private refactor, so the
   schedule is handled in one place.
7. **Leave smooth seasonality out of 0.2.0.** It needs a different way of
   reading the multiplier (measured below), and a schedule of changes covers
   the interventions the product document names.

## Q1 What can vary, and how a user says so

| Option | For | Against |
|---|---|---|
| **A. Piecewise-constant multiplier `m(t)` on `beta`** (recommended) | Matches how interventions are described ("cut contacts by 60%"). `R0` keeps its meaning as the value without interventions. The largest `beta` is `beta * max m`, known before the run. Fully validatable. | Cannot express gradual change except as many small steps. |
| B. Piecewise-constant absolute `beta` values | Same numerics as A. | `R0` and `beta_per_day` in the parameters become ambiguous: which one is "the" `beta`? |
| C. Smooth seasonal curve, `m(t) = 1 + a sin(2 pi t / P + phi)` | Classic for SIRS and influenza. | Needs stage-time lookup (Q2). A different family of tests; not named in the product document. |
| D. A user closure `Fn(f64) -> f64` | Most flexible. | Cannot be validated (finite, non-negative everywhere), its maximum is unknowable so the stability limit cannot be applied, and it cannot derive `Debug` or `PartialEq`, which the crate requires of public types. Rejected. |

Recommended shape, all names provisional:

```rust
/// Piecewise-constant multiplier on beta; 1 before the first change.
pub struct TransmissionSchedule { /* private: Vec<(f64, f64)> */ }

impl TransmissionSchedule {
    /// Changes as (day, multiplier). Days finite, positive and strictly
    /// increasing; multipliers finite and non-negative.
    pub fn new(changes: &[(f64, f64)]) -> Result<Self, SimulationError>;
    /// The multiplier in force on `day`.
    pub fn multiplier_at(&self, day: f64) -> f64;
}
```

## Q2 How the solver sees time

Today the derivative sees only the state. The prototype adds
`rk4_step_timed`, whose derivative also learns where in the step each stage
is. That allows three ways of reading `m(t)`:

- **stage time**, `start + fraction * h`;
- **step start**, `start` for every stage;
- **inward**, `start + fraction * h / 2`, which stays inside the step.

M1 and M2 run SIR (`beta = 0.3`, `gamma = 0.1`, `S = 990`, `I = 10`, 160 days)
at steps of 1, 1/2 and 1/4 day. They compare every compartment and the
cumulative infections at every step against a fine run: 1/256 day for M1,
and 1/1024 or 1/768 day for M2, so the reference has the changes on its own
grid. "Ratio" is how much the error falls when the step is halved; RK4
should give about 16.

**M1, smooth seasonal multiplier** `1 + 0.3 sin(2 pi t / 60)`:

| Lookup | Error at h = 1, 1/2, 1/4 (people) | Ratios |
|---|---|---|
| stage time | 1.1e-2, 7.7e-4, 5.0e-5 | 14.8, 15.4 |
| step start | 10.0, 5.0, 2.5 | 2.0, 2.0 |
| inward | 4.9, 2.5, 1.2 | 2.0, 2.0 |

**M2, lockdown** (multiplier 0.4 from day 20 to day 50):

| Lookup | Changes on the step grid (days 20, 50) | Changes off it (days 20 1/3, 50 1/3) |
|---|---|---|
| step start | 4.2e-3, 2.8e-4, 1.8e-5; ratios 14.9, 15.5 | 19.2, 4.9, 4.9; ratios 3.9, 1.0 |
| inward | identical to step start | 14.5, 4.9, 3.7; ratios 3.0, 1.3 |
| stage time | 5.0, 2.5, 1.2; ratios 2.0, 2.0 | 4.9, 2.4, 1.2; ratios 2.0, 2.0 |

What this says:

- **A schedule of changes needs step-start lookup and changes on the grid.**
  Then every step lies inside one segment, so RK4 sees a constant `beta` and
  keeps fourth order. Stage-time lookup is about 1,000 times worse: the k4
  stage of the step that ends on a change already reads the new value.
- **Changes off the grid lose fourth order however the multiplier is
  read.** Either reject them or split the step at the change (Q3).
- **Smooth seasonality is the opposite case:** it needs stage-time lookup,
  and step-start lookup turns RK4 into a first-order method.
- With step-start lookup the solver needs no change at all. The driver
  passes the step's start time to the derivative, `Fn(f64, [f64; N])`, and
  calls `rk4_step` with a closure that captures it. `rk4_step_timed` is only
  needed if smooth seasonality is ever added.

## Q3 Changes that do not fall on the step grid

| Option | For | Against |
|---|---|---|
| **Reject** (recommended): a change day must equal `k * step` within the four-ulp tolerance the step count already uses, or the run fails with a new error | Simple. Every step lies in one segment, so the M2 on-grid accuracy holds. The observation grid stays regular. | Users must choose a step that divides their change days. Days are usually whole numbers and steps are usually 1, 1/2, 1/4 or 0.1, so this rarely bites. |
| Split the step at each change | Any change day works. | The driver's step count, step limit, recording interval and observation times all become irregular, and the half-step error estimate must split at the same days. |

If the schedule lives in `SimulationSettings` (Q6), the check runs when the
schedule is added, so the error appears at configuration time.

## Q4 Stability limit

Each step lies inside one segment, so within a step the model is an ordinary
constant-parameter model with `beta * m`. The existing limit (fastest rate
times step at most 2) should therefore hold if the fastest rate is computed
with the largest `beta * m` among the segments that start before the end of
the run.

M3 checks this with the existing harness design. It draws 250,000
configurations from seed 42 with the library harness's draws, plus a
three-segment schedule: changes on whole days inside the run, multipliers in
`[0, 2)`. Each configuration runs at one day and at 1/64 day. A run is
flagged if, over some one-day step, a compartment or the cumulative
infections move against the exact flow from the same state (64 substeps),
or any final value misses the fine run by more than 0.1 N. Runs that go
materially negative count as loud failures, not as accepted.

Each configuration is judged once; each candidate limit then accepts the
runs whose rate times step is at most 2. "Smallest flagged product" is over
all flagged runs, accepted or not. A value above 2 is the candidate's margin.

| Model | Limit computed from | Accepted | Flagged | Smallest flagged rate x step |
|---|---|---|---|---|
| SIR | **largest `beta * m` in force** | 161,036 | **0** | 2.075 |
| SIR | `beta * m` at day 0 only | 197,466 | 97 | 0.022 |
| SEIR | **largest `beta * m` in force** | 79,102 | **0** | 2.281 |
| SEIR | `beta * m` at day 0 only | 110,716 | 57 | 0.467 |

M3 took 194 s in release mode. The largest `beta` in force gives no flagged
run, with margins of 2.075 and 2.281, against 2.035 for SEIR without a
schedule. Using only the `beta` at day 0 lets through runs whose later
segment raises `beta` past the limit; a flagged product of 0.022 means a run
accepted on day 0's `beta` went wrong after a large increase. SIS, SIRS and
SEIRS were not measured with schedules. The same argument applies (each step
lies in one segment with a constant `beta`), and the build plan adds
schedule rows for all five models to the library harness.

## Q5 Metrics

M4: SIR as above, run for 400 days at 1/16 day.

| Run | Incidence by day 400 | Peak (day) | Susceptible left |
|---|---|---|---|
| no schedule | 0.931 | 303.8 (day 26.6) | 58.8 |
| contacts at 40% from day 20 to 50 | 0.817 | 239.9 (day 20) | 172.9 |

The final-size formula with `R0 = 3` gives 0.931 in both cases. It assumes a
constant `beta`; under the lockdown it is wrong by 0.11.

| Metric | Recommendation |
|---|---|
| `basic_reproduction_number` | Unchanged: `beta / gamma`, the value without the schedule. The multiplier is 1 before the first change, so this is also `R0` at day 0. |
| `herd_immunity_threshold` | Unchanged, from that `R0`. The doc says it is the threshold without interventions. |
| `final_size_fraction` | `None` once the schedule has any change: there is no closed form. A schedule with no changes behaves exactly as today. |
| `endemic_prevalence` (SIS, SIRS, SEIRS) | Computed with `m_last * R0`, the multiplier in force after the last change, because that is the equilibrium the run settles to. |
| peak, cumulative infections, incidence | Unchanged; they are read from the trajectory. |
| effective R | `R_eff(t) = m(t) * R0 * S / N`. Keep the existing method (the baseline) and add `TransmissionSchedule::multiplier_at`, so `m(t) * R_eff` is one multiplication. Its doc must say so. |

## Q6 Where the schedule goes

| Option | For | Against |
|---|---|---|
| **`SimulationSettings::with_transmission_schedule(schedule)`** (recommended) | No `simulate_*` or `estimate_*` signature changes. The settings know the step, so the grid check (Q3) runs when the schedule is added. Half-step runs keep the schedule, which stays on the halved grid. | `SimulationSettings` stops being `Copy` because it holds a `Vec`. That breaks 5 places in the library (each step-error estimate uses the settings twice) and any user code that reuses settings after a call. `.clone()` fixes each. |
| A new argument on all ten public functions | Explicit. | Breaks every call. Users without a schedule must pass one anyway. |
| A field on each of the five parameter types | Close to `beta`. | Parameters stop being `Copy`, every struct literal breaks, and five types carry the same field. |

A fixed-capacity array (say 16 changes) would keep `Copy`, at the cost of an
arbitrary limit. Not recommended.

## Q7 The shared model pipeline

Each `simulate_*` function is 63 to 81 lines, 368 in all, and each
`estimate_*` is 14 to 21 lines. Of each `simulate_*`, only the derivative,
the fastest-rate formula, the state conversion and three to eight metric
lines differ. Adding the schedule by hand means about ten more repeated lines
per model:

- reading the multiplier;
- using `beta * max m` in the limit;
- the `final_size_fraction` and `endemic_prevalence` rules.

Recommendation: before the schedule, move the pipeline into one private
function driven by a private per-model trait. The trait has the compartment
array and names, validation, the fastest rate for a given `beta`, the
derivative for a given `beta`, the state conversion and the model-specific
metrics. The public API stays concrete: five `simulate_*` and five
`estimate_*` functions, now a few lines each. The technical specification's
decision against a *public* model trait is unaffected. Plan 014's tests are
the check that the move preserves behaviour.

## Owner decisions

On 2026-11-14 the owner accepted the recommended answer to all eight
questions below:

1. A piecewise-constant multiplier on `beta` only; no smooth seasonality in 0.2.0.
2. A multiplier, 1 before the first change.
3. Changes off the step grid are rejected with a new error naming the day.
4. The stability limit uses the largest `beta * m` in force.
5. `R0` and the herd-immunity threshold stay the baseline values; the final
   size is `None` once the schedule changes anything; the endemic prevalence
   uses the multiplier in force at the end of the run.
6. `TransmissionSchedule::multiplier_at` is added; the existing effective R
   method stays the baseline.
7. The schedule lives in `SimulationSettings`, which stops being `Copy`.
8. The model pipeline is shared first, in its own plan, with no behaviour change.

One change was made while writing the build plans. Proposed plans 2 and 3
are merged into one plan (018). The driver passing time to the derivative
only matters once a schedule exists, and in that plan the driver passes the
multiplier itself rather than the time.

## Open questions for the owner

1. **What may vary?** A piecewise-constant multiplier on `beta` only, or also
   smooth seasonality? *Recommended: the multiplier only, for 0.2.0.*
2. **Multiplier or absolute values?** *Recommended: a multiplier, 1 before
   the first change.*
3. **Changes off the step grid:** reject, or split the step? *Recommended:
   reject, with a new error naming the day.*
4. **Stability limit:** the largest `beta * m` in force, in the model's
   existing formula. *Recommended: measured with no flags (Q4).*
5. **Metrics:** baseline `R0` and herd-immunity threshold; final size `None`
   once the schedule changes anything; endemic prevalence from the last
   multiplier. *Recommended: as listed.*
6. **Effective R:** add `TransmissionSchedule::multiplier_at` and keep the
   existing method, or add a method that takes the schedule and a day?
   *Recommended: `multiplier_at`.*
7. **Where the schedule lives:** in `SimulationSettings`, which stops being
   `Copy`? *Recommended: yes.*
8. **Share the pipeline first?** *Recommended: yes, as its own plan, with no
   behaviour change.*

## Proposed build plans

In this order. Each leaves `cargo test` green and updates `docs/`,
`CHANGELOG.md` and the stability evidence. The frozen documents are not
edited.

1. **Shared model pipeline** (effort M). A private trait and one
   `simulate`/`estimate` function. No behaviour change: every test passes
   unchanged, the stability evidence prints the same counts, and the example
   outputs are byte-identical.
2. **The driver passes the step start time** (effort S). The derivative
   becomes `Fn(f64, [f64; N])`; the models ignore the time. No behaviour
   change.
3. **`TransmissionSchedule`** (effort M):
   - the type, and `SimulationSettings::with_transmission_schedule`;
   - the grid check and the new error;
   - the stability limit with the largest `beta`;
   - the metric rules;
   - tests: on-grid accuracy against a fine run, schedule-free runs
     unchanged bit for bit, the lockdown example, every error path;
   - schedule rows in `tests/stability_evidence_tests.rs`.
4. **Teaching example** (effort S, optional): `lockdown_curve`, CSV of the
   same outbreak with and without a lockdown.

The prototype's `rk4_step_timed` does not move on unless question 1 adds
seasonality.
