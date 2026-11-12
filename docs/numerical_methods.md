# Numerical methods

`episim-core` has one integrator: the classical fourth-order Runge-Kutta
method (RK4) with a fixed timestep. This page describes the scheme, how
simulation time is handled, the checks applied to every step, and the
limits of the approach.

## The RK4 step

For `dy/dt = f(y)` and a step of size `h`:

```text
k1 = f(y)
k2 = f(y + h/2 * k1)
k3 = f(y + h/2 * k2)
k4 = f(y + h   * k3)

y_next = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
```

The step works on fixed-size arrays (three values for SIR, four for SEIR)
and allocates nothing. Each stage input, each slope `k1` to `k4`, and
`y_next` is checked for NaN and infinity. A failure is reported as
`NonFiniteIntermediate { stage }`, where `stage` is `"k1"`, `"k2"`, `"k3"`,
`"k4"`, or `"final"`.

The global error of RK4 is `O(h^4)` for smooth problems: halving the step
should divide the error by about 16 until floating-point roundoff takes
over.

## Time handling

The number of steps is `duration / step` rounded up, and the run uses the
requested step for all of them except possibly the last, which is shortened
to `duration - time` so the simulation never runs past the requested
duration.

Observation times are computed as `index * step`, not by adding the step
repeatedly, so rounding error does not accumulate. The final observation's
time is set to the duration itself, so `time_days == duration_days` holds
exactly.

Two floating-point details are handled explicitly:

- **Step count rounding.** A duration that is a whole number of steps often
  divides to just above that number: `2.7 / 0.3` evaluates to
  `9.000000000000002`. Rounding that up would add a tenth step of
  rounding-error length and two nearly identical final observations. The
  quotient is therefore reduced by four ulps before it is rounded up, so
  `2.7 / 0.3` gives 9 steps. The consequence is that a genuine final
  interval shorter than about `1e-15` times the duration does not get a
  step of its own; the last full step absorbs it.
- **Steps too small for the clock.** If `duration + step == duration` in
  `f64`, the step cannot advance simulation time near the end of the run and
  the run is rejected with `TimeStepTooSmall`. Such a configuration needs
  more than about `2^53` steps, so it is only reached when the step limit
  has been raised far above its default; otherwise `StepLimitExceeded` is
  reported first.

## Resource limits

A run may take at most `max_steps` steps: 1,000,000 by default, changeable
with `SimulationSettings::with_max_steps`. The limit is compared with the
step count described above, before any integration or allocation, and
exceeding it returns `StepLimitExceeded`.

Memory for the recorded trajectory (the initial state plus one observation
per recorded step) is reserved up front with a fallible reservation. By
default every step is recorded; `SimulationSettings::with_record_every`
records every N-th step and always the final state. If it cannot be reserved, the run returns
`AllocationFailure`.

## Checks on every accepted state

The initial state and the state after every step go through the same
checks, in this order, whether or not the step is recorded.

1. **Finite values.** A NaN or infinite compartment is an error.
2. **Negative values.** RK4 does not guarantee non-negative compartments.
   With tolerance `64 * f64::EPSILON * max(N0, 1)`:
   - a value below `-tolerance` returns `NegativeCompartment`;
   - a value in `[-tolerance, 0)` is accepted, left unchanged in the
     trajectory, and counted in `tolerated_negative_value_count`.

   Values are never clipped to zero. Code that consumes a trajectory may
   therefore see a negative value of roundoff size.
3. **Population conservation.** With `N(t)` the sum of the compartments, a
   relative drift `|N(t) - N0| / N0` above `1e-9` returns
   `PopulationConservationViolation`. The largest absolute and relative
   drift seen are reported in the diagnostics.

## Limitations of a fixed step

- There is no adaptive step control and no error estimate. Accuracy depends
  entirely on the step you choose.
- RK4 is an explicit method. If the step is large compared with the fastest
  time scale (`1 / beta`, `1 / sigma`, or `1 / gamma`), the solution becomes
  inaccurate and, past a limit, unstable.
- **Steps that are too large are rejected up front.** For a compartment
  that should decay at rate `r`, RK4 multiplies it each step by
  `1 - z + z^2/2 - z^3/6 + z^4/24` with `z = r * step`. That factor is 1/3 at
  `z = 2`, rises steeply after it, and exceeds 1 once `z` is above about
  2.785. Past that point the compartment grows when it should shrink, while
  staying positive and conserving the population, so the per-step checks
  cannot see it: with `beta = 0`, `gamma = 3`, and a step of one day the
  infected compartment would grow by 37.5% per step.

  A run is therefore rejected with `TimeStepTooLarge { max_step_days }`
  before integration if the model's fastest decay rate times the step
  exceeds 2. `max_step_days` is the largest step accepted for the same
  parameters, `2 / fastest rate`. The fastest rate is:

  | Model | Fastest decay rate |
  |---|---|
  | SIR | `max(beta, gamma)` |
  | SEIR | `max(beta, (sigma + gamma) / 2 + sqrt(((sigma - gamma) / 2)^2 + sigma * beta))` |

  Susceptibles empty at up to `beta` per day. In SEIR the exposed and
  infected compartments feed each other, so together they decay faster than
  `sigma` or `gamma` alone; the second SEIR term is the fast rate of that
  pair in a fully susceptible population. With `beta = 0` it reduces to
  `max(sigma, gamma)`, and with `beta = sigma = gamma` it is twice the rate.
  The step compared with the limit is the step actually taken: a duration
  shorter than the timestep is integrated in one step of that duration.

  The limit of 2 was chosen from measurements, not derived. In 250,000
  random configurations per model, compared with runs at 1/64 of the step,
  no accepted run moved a compartment in the wrong direction or missed the
  final incidence by more than 0.1 (226,590 accepted SIR runs and 107,633
  accepted SEIR runs). With a limit at the stability boundary of 2.785 such
  runs do occur, from a product of about 2.2. This is evidence, not a proof
  that every accepted run is qualitatively right. The measurement is kept
  as an ignored test; repeat it with
  `cargo test --release --test stability_evidence_tests -- --ignored --nocapture`.
- **Inaccuracy is not detected.** The limit is a stability guard, not an
  accuracy target. At the limit a run still has large errors: with
  `gamma * step = 2` the infected compartment falls to 33% per step where
  the exact value is 14%. Keep every rate times the step well below 1.
- Positivity is not enforced by the method.
- Peak metrics are taken from recorded observations, so the peak time is
  resolved only to the nearest recorded step.

A practical check is to rerun with half the step and compare.
`estimate_sir_step_error` and `estimate_seir_step_error` do exactly that.

## Estimating the step error

Halving the step divides RK4's error by about 16. The difference between a
run and the same run at half the step is therefore about 15/16 of the first
run's error, and

```text
estimated error = 16/15 * largest difference between the two runs
```

taken over every compartment at every recorded time. The half-step run
records the same times as the original: its recording interval and its step
limit are both doubled.

| Field | Meaning |
|---|---|
| `max_absolute_error` | Estimated largest error of any compartment at any recorded time, in individuals. |
| `max_relative_error` | The same, divided by the initial population. |

Limits of the estimate:

- It costs about three times one run.
- It assumes the step is small enough for the sixteenfold rule to hold. In
  the reference scenarios (SIR with `beta = 0.5`, `gamma = 0.1`; SEIR with
  `beta = 0.5`, `sigma = 0.2`, `gamma = 0.1`; 40 days) it was within 4% of
  the error measured against a run at 1/256 day, for steps from 4 days, the
  largest the stability limit accepts there, down to 1/8 day.
- It covers the recorded times only. With a recording interval above one, an
  error peak between recorded steps is not seen.
- It is an estimate of discretisation error. It does not judge the model.

## Verification

The test suite checks the solver in three ways.

- **Order of accuracy.** For `dy/dt = -y` the error ratio between
  successive step halvings must lie between 14 and 18.
- **Closed-form references.** With `beta = 0` the SIR and SEIR trajectories
  are compared with the exact exponential solutions in
  [mathematical models](mathematical_models.md).
- **Refinement on the nonlinear models.** Runs with steps of 1, 0.5, and
  0.25 days are compared with a reference run at 1/64 day. The reference is
  numerical, not exact; its own error is about `16^4` times smaller than
  that of the finest step tested. The release criterion is an
  error-reduction factor of at least 8 per halving. The measured factors
  were 13.9 and 14.9 for SIR, and 15.6 and 15.8 for SEIR.

Population drift is checked across whole trajectories for populations of
`1e3` and `1e9`. The `1e-9` threshold is far above the drift observed in
these scenarios, which is around `1e-15`.
