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

The number of steps is `ceil(duration / step)`, and the run uses the
requested step for all of them except possibly the last, which is shortened
to `duration - time` so the simulation never runs past the requested
duration.

Observation times are computed as `index * step`, not by adding the step
repeatedly, so rounding error does not accumulate. The final observation's
time is set to the duration itself, so `time_days == duration_days` holds
exactly.

Two floating-point details are handled explicitly:

- **Step count rounding.** `duration / step` can round up to a count whose
  second-to-last grid point already equals the duration. For example,
  `609.6 / 0.6` evaluates to 1017, while `1016 * 0.6` is already `609.6`.
  The count is reduced by one in that case; otherwise the run would end with
  a zero-length step and two observations at the same time.
- **Steps too small for the clock.** If `duration + step == duration` in
  `f64`, the step cannot advance simulation time near the end of the run and
  the run is rejected with `TimeStepTooSmall`. Such a configuration needs
  more than about `2^53` steps, so it is only reached when the step limit
  has been raised far above its default; otherwise `StepLimitExceeded` is
  reported first.

## Resource limits

A run may take at most `max_steps` steps: 1,000,000 by default, changeable
with `SimulationSettings::with_max_steps`. The limit is checked before any
integration or allocation, and exceeding it returns `StepLimitExceeded`.

Memory for the whole trajectory (`steps + 1` observations) is reserved up
front with a fallible reservation. If it cannot be reserved, the run returns
`AllocationFailure`.

## Checks on every recorded state

The initial state and the state after every step go through the same
checks, in this order.

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
  time scale (`1 / beta`, `1 / sigma`, or `1 / gamma`), the solution can
  become inaccurate or unstable. Instability shows up as one of the errors
  above, not as silently wrong output. An inaccurate but stable run is not
  detected, so compare against a smaller step if accuracy matters.
- Positivity is not enforced by the method.
- Peak metrics are taken from recorded observations, so the peak time is
  resolved only to the nearest step.

A practical check is to rerun with half the step and compare.

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
