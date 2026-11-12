use crate::result::Observation;

/// The error of classical RK4 falls by about `2^4 = 16` when the step is
/// halved, so the difference between a run and the same run at half the step
/// is about 15/16 of the first run's error.
const RICHARDSON_FACTOR: f64 = 16.0 / 15.0;

/// Estimated error that the timestep introduces into a simulation.
///
/// The estimate compares the run with the same run at half the timestep. It
/// measures discretisation error only: how far the numerical trajectory is
/// from the exact solution of the model equations. It says nothing about how
/// well the model describes a real outbreak.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepErrorEstimate {
    /// Estimated largest absolute error, in individuals, of any compartment
    /// at any recorded time.
    pub max_absolute_error: f64,
    /// `max_absolute_error` divided by the initial population.
    pub max_relative_error: f64,
}

/// Compares a run with its half-step run, which must hold observations at
/// the same times (see `SimulationSettings::halved`).
pub(crate) fn estimate_step_error<S, const N: usize>(
    run: &[Observation<S>],
    half_step_run: &[Observation<S>],
    population: f64,
    compartments: impl Fn(&S) -> [f64; N],
) -> StepErrorEstimate {
    debug_assert_eq!(run.len(), half_step_run.len());
    let mut max_difference = 0.0_f64;
    for (coarse, fine) in run.iter().zip(half_step_run) {
        debug_assert_eq!(coarse.time_days, fine.time_days);
        let (coarse, fine) = (compartments(&coarse.state), compartments(&fine.state));
        for index in 0..N {
            max_difference = max_difference.max((coarse[index] - fine[index]).abs());
        }
    }
    let max_absolute_error = RICHARDSON_FACTOR * max_difference;
    StepErrorEstimate {
        max_absolute_error,
        max_relative_error: max_absolute_error / population,
    }
}
