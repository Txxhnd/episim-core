use crate::error::SimulationError;
use crate::result::{Observation, SimulationDiagnostics};
use crate::settings::SimulationSettings;
use crate::solver::rk4::rk4_step;

/// Largest accepted value of `|N(t) - N(0)| / N(0)`.
const MAX_RELATIVE_POPULATION_DRIFT: f64 = 1e-9;

/// Trajectory and diagnostics of a completed integration, before metrics.
#[derive(Debug)]
pub(crate) struct SimulationOutput<S> {
    pub(crate) trajectory: Vec<Observation<S>>,
    pub(crate) diagnostics: SimulationDiagnostics,
}

/// Integrates `derivative` from `initial` over the configured duration,
/// recording the initial state, the state after every `record_every`-th
/// accepted step, and the final state.
///
/// The derivative receives the transmission multiplier for the step, from
/// the settings' schedule (1 without one), and the state. The last entry of
/// the state is the cumulative infection count `C`, which the derivative
/// integrates as `dC/dt = beta * S * I / N`. The entries
/// before it are the compartments, named by `compartment_names`. `C` is
/// recorded in each [`Observation`] and is not part of the population.
///
/// `population` is the fixed initial population the caller validated; the
/// compartments of every accepted state are checked against it, recorded or
/// not.
pub(crate) fn run_simulation<const N: usize, S, D, C>(
    initial: [f64; N],
    population: f64,
    settings: SimulationSettings,
    derivative: D,
    to_state: C,
    compartment_names: &[&'static str],
) -> Result<SimulationOutput<S>, SimulationError>
where
    D: Fn(f64, [f64; N]) -> [f64; N],
    C: Fn([f64; N]) -> S,
{
    assert_eq!(compartment_names.len() + 1, N, "one name per compartment");
    let compartments = N - 1;
    let observe = |time_days, state: [f64; N]| Observation {
        time_days,
        cumulative_infections: state[compartments],
        state: to_state(state),
    };
    let duration = settings.duration_days();
    let step = settings.step_days();
    let required_steps = required_steps(&settings)?;
    let record_every = settings.record_every();
    let substeps = settings.substeps();
    // The initial state, every `record_every`-th step, and the final step.
    let capacity = required_steps
        .div_ceil(record_every)
        .checked_add(1)
        .ok_or(SimulationError::AllocationFailure)?;
    let mut trajectory = reserve_trajectory(capacity)?;
    let mut diagnostics = SimulationDiagnostics {
        steps_executed: 0,
        maximum_absolute_population_drift: 0.0,
        maximum_relative_population_drift: 0.0,
        tolerated_negative_value_count: 0,
    };

    let mut state = initial;
    let mut time = 0.0;
    check_state(
        &state[..compartments],
        compartment_names,
        population,
        &mut diagnostics,
    )?;
    trajectory.push(observe(time, state));
    for index in 1..=required_steps {
        // The schedule's changes lie on the step grid, so one multiplier holds
        // for the whole step and RK4 keeps its order.
        let multiplier = settings.multiplier_for_step(index - 1);
        // A step-error estimate's refined run splits every step, including a
        // shortened final one, into equal substeps; each is checked.
        let substep = step.min(duration - time) / substeps as f64;
        for substep_index in 0..substeps {
            state = rk4_step(state, substep, |values| derivative(multiplier, values))?;
            if substep_index + 1 < substeps {
                check_state(
                    &state[..compartments],
                    compartment_names,
                    population,
                    &mut diagnostics,
                )?;
            }
        }
        // Times come from the step index, not from repeated addition, so
        // rounding does not accumulate; the last one is the duration itself.
        time = if index == required_steps {
            duration
        } else {
            index as f64 * step
        };
        check_state(
            &state[..compartments],
            compartment_names,
            population,
            &mut diagnostics,
        )?;
        if index % record_every == 0 || index == required_steps {
            trajectory.push(observe(time, state));
        }
        diagnostics.steps_executed = index * substeps;
    }
    Ok(SimulationOutput {
        trajectory,
        diagnostics,
    })
}

/// Number of integration steps needed to reach the duration, rejecting
/// configurations that exceed the step limit or cannot advance the clock.
fn required_steps(settings: &SimulationSettings) -> Result<usize, SimulationError> {
    let duration = settings.duration_days();
    let step = settings.step_days();
    // A whole number of steps often divides to just above that number
    // (2.7 / 0.3 is 9.000000000000002). Shaving a few ulps off before rounding
    // up keeps such a run from gaining a final step of rounding-error length.
    // A duration far smaller than the step underflows to zero; it still needs
    // one step.
    // shortcut: the slack grows with the quotient and nears a whole step
    // around 2^50 steps; a run that long would take years to compute, so
    // revisit this only if steps become several orders of magnitude cheaper.
    let steps = (duration / step * (1.0 - 4.0 * f64::EPSILON))
        .ceil()
        .max(1.0);
    if !steps.is_finite() || steps > settings.max_steps() as f64 {
        return Err(SimulationError::StepLimitExceeded);
    }
    // If one step cannot move the clock near the end of the run, the time
    // grid collapses onto repeated values.
    if duration + step == duration {
        return Err(SimulationError::TimeStepTooSmall);
    }
    // In range: `steps` is at least 1 and at most `max_steps`.
    Ok(steps as usize)
}

fn reserve_trajectory<S>(capacity: usize) -> Result<Vec<Observation<S>>, SimulationError> {
    let mut trajectory = Vec::new();
    trajectory
        .try_reserve_exact(capacity)
        .map_err(|_| SimulationError::AllocationFailure)?;
    Ok(trajectory)
}

/// Rejects a non-finite or materially negative compartment and excessive
/// population drift; counts tolerated negative roundoff without changing it.
fn check_state(
    compartments: &[f64],
    compartment_names: &[&'static str],
    population: f64,
    diagnostics: &mut SimulationDiagnostics,
) -> Result<(), SimulationError> {
    let negative_tolerance = 64.0 * f64::EPSILON * population.max(1.0);
    for (&value, &compartment) in compartments.iter().zip(compartment_names) {
        if !value.is_finite() {
            return Err(SimulationError::NonFiniteNumericalState { compartment });
        }
        if value < -negative_tolerance {
            return Err(SimulationError::NegativeCompartment { compartment });
        }
        if value < 0.0 {
            diagnostics.tolerated_negative_value_count += 1;
        }
    }

    let absolute_drift = (compartments.iter().sum::<f64>() - population).abs();
    let relative_drift = absolute_drift / population;
    if !relative_drift.is_finite() || relative_drift > MAX_RELATIVE_POPULATION_DRIFT {
        return Err(SimulationError::PopulationConservationViolation);
    }
    diagnostics.maximum_absolute_population_drift = diagnostics
        .maximum_absolute_population_drift
        .max(absolute_drift);
    diagnostics.maximum_relative_population_drift = diagnostics
        .maximum_relative_population_drift
        .max(relative_drift);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    type Output = Result<SimulationOutput<[f64; 2]>, SimulationError>;

    fn settings(duration_days: f64, step_days: f64) -> SimulationSettings {
        SimulationSettings::new(duration_days, step_days).unwrap()
    }

    /// Runs a two-compartment system whose population is the initial sum,
    /// with an accumulator that stays at zero.
    fn run<D>(initial: [f64; 2], settings: SimulationSettings, derivative: D) -> Output
    where
        D: Fn([f64; 2]) -> [f64; 2],
    {
        let population = initial.iter().sum();
        run_simulation(
            [initial[0], initial[1], 0.0],
            population,
            settings,
            |_, [a, b, _]| {
                let [da, db] = derivative([a, b]);
                [da, db, 0.0]
            },
            |[a, b, _]| [a, b],
            &["a", "b"],
        )
    }

    /// Moves one individual per day from `a` to `b`, conserving the total.
    fn unit_flow(_: [f64; 2]) -> [f64; 2] {
        [-1.0, 1.0]
    }

    fn times(output: &SimulationOutput<[f64; 2]>) -> Vec<f64> {
        output.trajectory.iter().map(|o| o.time_days).collect()
    }

    #[test]
    fn driver_records_initial_state() {
        let output = run([100.0, 0.0], settings(10.0, 1.0), unit_flow).unwrap();
        assert_eq!(
            output.trajectory[0],
            Observation {
                time_days: 0.0,
                cumulative_infections: 0.0,
                state: [100.0, 0.0]
            }
        );
    }

    #[test]
    fn driver_records_the_accumulator_outside_the_population() {
        // One individual a day moves from a to b while the accumulator grows
        // by five a day: it is recorded, but not counted in the population.
        let output = run_simulation(
            [100.0, 0.0, 0.0],
            100.0,
            settings(10.0, 1.0),
            |_, _| [-1.0, 1.0, 5.0],
            |[a, b, _]| [a, b],
            &["a", "b"],
        )
        .unwrap();
        let last = output.trajectory[10];
        assert_eq!(last.state, [90.0, 10.0]);
        assert_eq!(last.cumulative_infections, 50.0);
        assert_eq!(output.diagnostics.maximum_absolute_population_drift, 0.0);
    }

    #[test]
    fn driver_records_final_state_at_requested_duration() {
        // 0.1 accumulates to just under 1.0 after ten additions. The other
        // durations are whole multiples of the step whose quotient rounds to
        // just above the whole number (2.7 / 0.3 is 9.000000000000002).
        for (duration, step, steps) in [
            (1.0, 0.1, 10),
            (609.6, 0.6, 1016),
            (2.7, 0.3, 9),
            (4.2, 0.7, 6),
        ] {
            let output = run([1000.0, 0.0], settings(duration, step), unit_flow).unwrap();
            let times = times(&output);
            assert_eq!(times.len(), steps + 1, "duration {duration}");
            assert_eq!(times[steps], duration);
            assert!(
                times.windows(2).all(|pair| pair[0] < pair[1]),
                "times must strictly increase for duration {duration}"
            );
        }
    }

    #[test]
    fn driver_uses_short_final_partial_step() {
        let output = run([100.0, 0.0], settings(10.0, 3.0), unit_flow).unwrap();
        assert_eq!(times(&output), [0.0, 3.0, 6.0, 9.0, 10.0]);
        // A full final step would have moved 12 individuals instead of 10.
        assert_eq!(output.trajectory[4].state, [90.0, 10.0]);
    }

    #[test]
    fn driver_counts_executed_steps() {
        let output = run([100.0, 0.0], settings(10.0, 3.0), unit_flow).unwrap();
        assert_eq!(output.diagnostics.steps_executed, 4);
    }

    #[test]
    fn driver_rejects_step_limit_exceeded_before_integration() {
        let calls = Cell::new(0);
        let counted = |state| {
            calls.set(calls.get() + 1);
            unit_flow(state)
        };
        // Ten steps are required: a limit of nine is exceeded, a limit of ten is not.
        let nine = settings(10.0, 1.0).with_max_steps(9).unwrap();
        assert_eq!(
            run([100.0, 0.0], nine, counted).unwrap_err(),
            SimulationError::StepLimitExceeded
        );
        assert_eq!(calls.get(), 0, "no derivative evaluation before rejection");

        let ten = settings(10.0, 1.0).with_max_steps(10).unwrap();
        let output = run([100.0, 0.0], ten, unit_flow).unwrap();
        assert_eq!(output.diagnostics.steps_executed, 10);

        // The limit applies to the steps actually taken: 609.6 / 0.6 needs
        // 1016 steps even though the quotient rounds to just above 1016.
        let exact = settings(609.6, 0.6).with_max_steps(1016).unwrap();
        let output = run([1000.0, 0.0], exact, unit_flow).unwrap();
        assert_eq!(output.diagnostics.steps_executed, 1016);
        let one_short = settings(609.6, 0.6).with_max_steps(1015).unwrap();
        assert_eq!(
            run([1000.0, 0.0], one_short, unit_flow).unwrap_err(),
            SimulationError::StepLimitExceeded
        );
    }

    #[test]
    fn driver_counts_whole_steps_for_decimal_multiples_of_the_step() {
        for step in [0.1, 0.15, 0.3, 0.6, 0.7] {
            for multiple in 1..=5000_usize {
                // The decimal number a user would type, not `multiple * step` in binary.
                let duration: f64 = format!("{:.2}", multiple as f64 * step).parse().unwrap();
                assert_eq!(
                    required_steps(&settings(duration, step)),
                    Ok(multiple),
                    "duration {duration}, step {step}"
                );
            }
        }
    }

    // On 32-bit targets every such configuration already exceeds the largest step limit.
    #[cfg(target_pointer_width = "64")]
    #[test]
    fn driver_rejects_time_step_that_cannot_advance_clock() {
        // Near 1e17 adjacent f64 values are 16 apart, so adding 1.0 changes nothing.
        let stalled = settings(1e17, 1.0).with_max_steps(usize::MAX).unwrap();
        assert_eq!(
            run([100.0, 0.0], stalled, unit_flow).unwrap_err(),
            SimulationError::TimeStepTooSmall
        );
    }

    #[test]
    fn driver_checks_population_drift_at_every_accepted_step() {
        // Conserves the total for two steps, then starts creating individuals.
        let leaks_later = |[a, _]: [f64; 2]| [-1.0, if a < 97.5 { 2.0 } else { 1.0 }];
        assert_eq!(
            run([100.0, 0.0], settings(10.0, 1.0), leaks_later).unwrap_err(),
            SimulationError::PopulationConservationViolation
        );

        // The initial state is checked against the supplied population as well.
        let mismatched = run_simulation(
            [100.0, 0.0, 0.0],
            200.0,
            settings(10.0, 1.0),
            |_, _| [-1.0, 1.0, 0.0],
            |[a, b, _]| [a, b],
            &["a", "b"],
        );
        assert_eq!(
            mismatched.unwrap_err(),
            SimulationError::PopulationConservationViolation
        );

        let conserving = run([100.0, 0.0], settings(10.0, 1.0), unit_flow).unwrap();
        assert_eq!(
            conserving.diagnostics.maximum_absolute_population_drift,
            0.0
        );
        assert_eq!(
            conserving.diagnostics.maximum_relative_population_drift,
            0.0
        );
    }

    #[test]
    fn driver_rejects_materially_negative_compartment() {
        assert_eq!(
            run([0.5, 100.0], settings(10.0, 1.0), unit_flow).unwrap_err(),
            SimulationError::NegativeCompartment { compartment: "a" }
        );
    }

    #[test]
    fn driver_counts_tolerated_negative_roundoff_without_clipping() {
        // The tolerance for a population of 1000 is 64 * EPSILON * 1000, about 1.4e-11.
        let tiny_flow = |_| [-1e-12, 1e-12];
        let output = run([0.0, 1000.0], settings(2.0, 1.0), tiny_flow).unwrap();
        assert_eq!(output.diagnostics.tolerated_negative_value_count, 2);
        assert_eq!(output.trajectory[1].state[0], -1e-12);
        assert_eq!(output.trajectory[2].state[0], -2e-12);
    }

    /// Runs the per-state checks alone and returns the tolerated-negative count.
    fn check(state: [f64; 2], population: f64) -> Result<usize, SimulationError> {
        let mut diagnostics = SimulationDiagnostics {
            steps_executed: 0,
            maximum_absolute_population_drift: 0.0,
            maximum_relative_population_drift: 0.0,
            tolerated_negative_value_count: 0,
        };
        check_state(&state, &["a", "b"], population, &mut diagnostics)?;
        Ok(diagnostics.tolerated_negative_value_count)
    }

    #[test]
    fn driver_rejects_non_finite_accepted_state() {
        // The solver rejects non-finite results itself, so this check is
        // only reachable directly.
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                check([100.0, value], 100.0),
                Err(SimulationError::NonFiniteNumericalState { compartment: "b" }),
                "value {value}"
            );
        }
    }

    #[test]
    fn driver_applies_negative_tolerance_at_its_boundary() {
        // The tolerance is 64 * EPSILON * max(N0, 1): it scales with a
        // population above 1 and stops shrinking below 1.
        let cases = [
            (1000.0, 64.0 * f64::EPSILON * 1000.0),
            (0.5, 64.0 * f64::EPSILON),
        ];
        for (population, tolerance) in cases {
            assert_eq!(
                check([population + tolerance, -tolerance], population),
                Ok(1),
                "population {population}"
            );
            let beyond = tolerance * 1.000001;
            assert_eq!(
                check([population + beyond, -beyond], population),
                Err(SimulationError::NegativeCompartment { compartment: "b" }),
                "population {population}"
            );
        }
    }

    #[test]
    fn driver_applies_drift_threshold_near_its_boundary() {
        // In a population of 1000 a relative drift of 1e-9 is 1e-6 individuals.
        for accepted in [1000.0 + 0.99e-6, 1000.0 - 0.99e-6] {
            assert_eq!(check([accepted, 0.0], 1000.0), Ok(0), "total {accepted}");
        }
        for rejected in [1000.0 + 1.01e-6, 1000.0 - 1.01e-6] {
            assert_eq!(
                check([rejected, 0.0], 1000.0),
                Err(SimulationError::PopulationConservationViolation),
                "total {rejected}"
            );
        }
    }

    #[test]
    fn driver_records_every_nth_step_and_always_the_final_state() {
        let every = |steps| settings(10.0, 1.0).with_record_every(steps).unwrap();

        // Ten steps, every third recorded, plus the final one.
        let output = run([100.0, 0.0], every(3), unit_flow).unwrap();
        assert_eq!(times(&output), [0.0, 3.0, 6.0, 9.0, 10.0]);
        assert_eq!(output.trajectory[4].state, [90.0, 10.0]);
        assert_eq!(output.diagnostics.steps_executed, 10);

        // The final step is not recorded twice when the interval divides the count.
        let output = run([100.0, 0.0], every(5), unit_flow).unwrap();
        assert_eq!(times(&output), [0.0, 5.0, 10.0]);

        // An interval longer than the run leaves the two end points.
        let output = run([100.0, 0.0], every(1000), unit_flow).unwrap();
        assert_eq!(times(&output), [0.0, 10.0]);
    }

    #[test]
    fn driver_checks_and_counts_steps_that_are_not_recorded() {
        let sparse = |duration| settings(duration, 1.0).with_record_every(1000).unwrap();

        // The leak starts at the third step, which is never recorded.
        let leaks_later = |[a, _]: [f64; 2]| [-1.0, if a < 97.5 { 2.0 } else { 1.0 }];
        assert_eq!(
            run([100.0, 0.0], sparse(10.0), leaks_later).unwrap_err(),
            SimulationError::PopulationConservationViolation
        );

        // Both steps produce a tolerated negative value; only the last is recorded.
        let tiny_flow = |_| [-1e-12, 1e-12];
        let output = run([0.0, 1000.0], sparse(2.0), tiny_flow).unwrap();
        assert_eq!(output.trajectory.len(), 2);
        assert_eq!(output.diagnostics.tolerated_negative_value_count, 2);
    }

    #[test]
    fn driver_halved_settings_split_every_step_actually_taken() {
        // Issue #1: the half-step run must split each step the run takes,
        // including a final step shortened to the duration and a single step
        // shorter than the nominal one, and record the same times.
        for (duration, step, steps) in [(10.0, 3.0, 4), (1.0, 2.0, 1)] {
            let settings = settings(duration, step);
            let run = run([100.0, 0.0], settings.clone(), unit_flow).unwrap();
            let half = run_simulation_with(&settings.halved());
            assert_eq!(times(&half), times(&run), "duration {duration}");
            assert_eq!(run.diagnostics.steps_executed, steps);
            assert_eq!(
                half.diagnostics.steps_executed,
                2 * steps,
                "duration {duration}"
            );
        }
    }

    /// Runs the unit flow from `[100, 0]` with the given settings.
    fn run_simulation_with(settings: &SimulationSettings) -> SimulationOutput<[f64; 2]> {
        run([100.0, 0.0], settings.clone(), unit_flow).unwrap()
    }

    #[test]
    fn driver_maps_trajectory_reservation_failure_to_typed_error() {
        // A capacity of usize::MAX observations overflows the size computation,
        // so the reservation fails without touching any memory.
        assert_eq!(
            reserve_trajectory::<[f64; 2]>(usize::MAX).unwrap_err(),
            SimulationError::AllocationFailure
        );
    }
}
