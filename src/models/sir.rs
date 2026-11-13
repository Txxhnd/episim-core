use crate::accuracy::{StepErrorEstimate, estimate_step_error};
use crate::error::SimulationError;
use crate::metrics::{
    basic_reproduction_number, calculate_metrics, effective_reproduction_number,
    final_size_fraction, herd_immunity_threshold,
};
use crate::result::SimulationResult;
use crate::settings::SimulationSettings;
use crate::simulation::driver::run_simulation;
use crate::validation;

/// Rate parameters of the SIR model, constant for the whole simulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SirParameters {
    /// Transmission rate in inverse days; must be finite and non-negative.
    pub beta_per_day: f64,
    /// Recovery rate in inverse days; must be finite and strictly positive.
    pub gamma_per_day: f64,
}

/// Population counts in each SIR compartment.
///
/// Counts are expected values and may be fractional. Every compartment must
/// be finite and non-negative, and the total must be finite and strictly
/// positive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SirState {
    /// Individuals who can become infected.
    pub susceptible: f64,
    /// Individuals who are currently infectious.
    pub infected: f64,
    /// Individuals who have recovered and are assumed immune.
    pub recovered: f64,
}

impl SirParameters {
    pub(crate) fn validate(&self) -> Result<(), SimulationError> {
        validation::non_negative_rate("beta_per_day", self.beta_per_day)?;
        validation::positive_rate("gamma_per_day", self.gamma_per_day)
    }

    /// Model-derived effective reproduction number in `state`: `R0 * S / N`,
    /// where `N` is the total of `state`.
    ///
    /// It is the average number of people one infected person infects while a
    /// fraction `S / N` of the population is still susceptible. Because
    /// `dI/dt = gamma * I * (effective R - 1)`, the infected grow while it is
    /// above 1 and decline while it is below 1. For a state returned by
    /// [`simulate_sir`], `N` equals the initial population to within the
    /// population drift tolerance of `1e-9`.
    ///
    /// This is a property of the model. It is not an estimate of the effective
    /// reproduction number `Rt` that public-health agencies infer from case
    /// data.
    ///
    /// # Errors
    ///
    /// Returns [`SimulationError::InvalidParameter`] for an invalid rate,
    /// [`SimulationError::InvalidPopulation`] if the total of `state` is not
    /// finite and positive, and [`SimulationError::MetricOverflow`] if the
    /// result is not finite.
    ///
    /// # Examples
    ///
    /// ```
    /// use episim_core::{SimulationSettings, SirParameters, SirState, simulate_sir};
    ///
    /// let parameters = SirParameters { beta_per_day: 0.3, gamma_per_day: 0.1 };
    /// let initial_state = SirState { susceptible: 990.0, infected: 10.0, recovered: 0.0 };
    /// let result = simulate_sir(parameters, initial_state, SimulationSettings::new(160.0, 0.25)?)?;
    ///
    /// // R0 is 3 and 99% of the population starts susceptible.
    /// let start = parameters.effective_reproduction_number(&result.trajectory[0].state)?;
    /// assert!((start - 2.97).abs() < 1e-9);
    /// // By day 160 too few susceptibles remain to sustain transmission.
    /// let end = parameters.effective_reproduction_number(&result.trajectory.last().unwrap().state)?;
    /// assert!(end < 1.0);
    /// # Ok::<(), episim_core::SimulationError>(())
    /// ```
    pub fn effective_reproduction_number(&self, state: &SirState) -> Result<f64, SimulationError> {
        self.validate()?;
        let r0 = basic_reproduction_number(self.beta_per_day, self.gamma_per_day)?;
        effective_reproduction_number(
            r0,
            state.susceptible,
            state.susceptible + state.infected + state.recovered,
        )
    }
}

impl SirState {
    /// Validates the state and returns the total population it defines.
    pub(crate) fn initial_population(&self) -> Result<f64, SimulationError> {
        validation::initial_population(&[
            ("susceptible", self.susceptible),
            ("infected", self.infected),
            ("recovered", self.recovered),
        ])
    }
}

/// Simulates the deterministic SIR model with a fixed-step classical RK4 integrator.
///
/// The model assumes a closed, homogeneously mixing population with constant
/// parameters and permanent immunity. Transmission is frequency-dependent,
/// with `N` fixed at the initial total population:
///
/// ```text
/// dS/dt = -beta * S * I / N
/// dI/dt =  beta * S * I / N - gamma * I
/// dR/dt =  gamma * I
/// ```
///
/// The returned trajectory starts at day zero, holds one observation per
/// integration step, and ends exactly at the requested duration.
///
/// # Errors
///
/// Returns a [`SimulationError`] if a parameter or the initial state is
/// invalid, if the timestep is too large for stable integration with these
/// rates, if the settings require more steps than their limit or a
/// timestep too small to advance the clock, if the integration produces a
/// non-finite value, a materially negative compartment, or excessive
/// population drift, if a metric is not finite, or if the trajectory cannot
/// be allocated.
///
/// # Examples
///
/// ```
/// use episim_core::{SimulationSettings, SirParameters, SirState, simulate_sir};
///
/// let parameters = SirParameters { beta_per_day: 0.3, gamma_per_day: 0.1 };
/// let initial_state = SirState { susceptible: 990.0, infected: 10.0, recovered: 0.0 };
/// let settings = SimulationSettings::new(160.0, 0.25)?;
///
/// let result = simulate_sir(parameters, initial_state, settings)?;
///
/// assert_eq!(result.trajectory.last().unwrap().time_days, 160.0);
/// assert!(result.metrics.peak_infectious > initial_state.infected);
/// # Ok::<(), episim_core::SimulationError>(())
/// ```
pub fn simulate_sir(
    parameters: SirParameters,
    initial_state: SirState,
    settings: SimulationSettings,
) -> Result<SimulationResult<SirState>, SimulationError> {
    parameters.validate()?;
    let population = initial_state.initial_population()?;
    let SirParameters {
        beta_per_day,
        gamma_per_day,
    } = parameters;
    // Susceptibles empty at up to beta per day and infected at gamma per
    // day; no state of the model decays faster than the larger of the two.
    // A duration shorter than the step is integrated in one shorter step.
    validation::stable_step(
        settings.step_days().min(settings.duration_days()),
        beta_per_day.max(gamma_per_day),
    )?;
    // Checked now so that an unusable R0 is reported before any integration.
    basic_reproduction_number(beta_per_day, gamma_per_day)?;

    let output = run_simulation(
        [
            initial_state.susceptible,
            initial_state.infected,
            initial_state.recovered,
            0.0,
        ],
        population,
        settings,
        |[susceptible, infected, _, _]| {
            // Dividing first avoids multiplying two large population values.
            let infection_flow = beta_per_day * (susceptible / population) * infected;
            let recovery_flow = gamma_per_day * infected;
            [
                -infection_flow,
                infection_flow - recovery_flow,
                recovery_flow,
                infection_flow,
            ]
        },
        |[susceptible, infected, recovered, _]| SirState {
            susceptible,
            infected,
            recovered,
        },
        &["susceptible", "infected", "recovered"],
    )?;
    let mut metrics = calculate_metrics(
        &output.trajectory,
        population,
        beta_per_day,
        gamma_per_day,
        |state| state.infected,
    )?;
    let r0 = metrics.basic_reproduction_number;
    metrics.herd_immunity_threshold = Some(herd_immunity_threshold(r0));
    metrics.final_size_fraction = Some(final_size_fraction(
        initial_state.susceptible,
        initial_state.infected,
        population,
        r0,
    ));
    Ok(SimulationResult {
        trajectory: output.trajectory,
        metrics,
        diagnostics: output.diagnostics,
    })
}

/// Estimates the error that the timestep introduces into [`simulate_sir`] for
/// the same arguments.
///
/// The simulation is run twice, once as given and once at half the timestep,
/// and the two trajectories are compared at every recorded time. This costs
/// about three times as much as one run. The half-step run is allowed twice
/// the configured step limit.
///
/// The estimate assumes the timestep is small enough for the solver's error
/// to fall sixteenfold when the step is halved. In the reference scenarios
/// it was within 4% of the true error for every step the stability limit
/// accepts.
///
/// # Errors
///
/// Returns any [`SimulationError`] that [`simulate_sir`] returns for either run.
///
/// # Examples
///
/// ```
/// use episim_core::{SimulationSettings, SirParameters, SirState, estimate_sir_step_error};
///
/// let parameters = SirParameters { beta_per_day: 0.3, gamma_per_day: 0.1 };
/// let initial_state = SirState { susceptible: 990.0, infected: 10.0, recovered: 0.0 };
///
/// let coarse = estimate_sir_step_error(parameters, initial_state, SimulationSettings::new(160.0, 1.0)?)?;
/// let fine = estimate_sir_step_error(parameters, initial_state, SimulationSettings::new(160.0, 0.25)?)?;
///
/// // A quarter of the step gives roughly 4^4 = 256 times less error.
/// assert!(fine.max_absolute_error < coarse.max_absolute_error / 100.0);
/// assert!(fine.max_relative_error < 1e-6);
/// # Ok::<(), episim_core::SimulationError>(())
/// ```
pub fn estimate_sir_step_error(
    parameters: SirParameters,
    initial_state: SirState,
    settings: SimulationSettings,
) -> Result<StepErrorEstimate, SimulationError> {
    let run = simulate_sir(parameters, initial_state, settings)?;
    let half_step_run = simulate_sir(parameters, initial_state, settings.halved())?;
    Ok(estimate_step_error(
        &run.trajectory,
        &half_step_run.trajectory,
        initial_state.initial_population()?,
        |state| [state.susceptible, state.infected, state.recovered],
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const INVALID_VALUES: [f64; 4] = [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

    const VALID_STATE: SirState = SirState {
        susceptible: 990.0,
        infected: 10.0,
        recovered: 0.0,
    };

    #[test]
    fn sir_accepts_zero_beta_and_positive_gamma() {
        let parameters = SirParameters {
            beta_per_day: 0.0,
            gamma_per_day: 0.1,
        };
        assert_eq!(parameters.validate(), Ok(()));
        assert_eq!(VALID_STATE.initial_population(), Ok(1000.0));
    }

    #[test]
    fn sir_rejects_negative_or_non_finite_beta() {
        for beta_per_day in INVALID_VALUES {
            let parameters = SirParameters {
                beta_per_day,
                gamma_per_day: 0.1,
            };
            assert_eq!(
                parameters.validate(),
                Err(SimulationError::InvalidParameter {
                    name: "beta_per_day"
                }),
                "beta {beta_per_day}"
            );
        }
    }

    #[test]
    fn sir_rejects_non_positive_or_non_finite_gamma() {
        for gamma_per_day in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let parameters = SirParameters {
                beta_per_day: 0.3,
                gamma_per_day,
            };
            assert_eq!(
                parameters.validate(),
                Err(SimulationError::InvalidParameter {
                    name: "gamma_per_day"
                }),
                "gamma {gamma_per_day}"
            );
        }
    }

    #[test]
    fn sir_rejects_negative_or_non_finite_compartments() {
        for invalid in INVALID_VALUES {
            let cases = [
                (
                    "susceptible",
                    SirState {
                        susceptible: invalid,
                        ..VALID_STATE
                    },
                ),
                (
                    "infected",
                    SirState {
                        infected: invalid,
                        ..VALID_STATE
                    },
                ),
                (
                    "recovered",
                    SirState {
                        recovered: invalid,
                        ..VALID_STATE
                    },
                ),
            ];
            for (compartment, state) in cases {
                assert_eq!(
                    state.initial_population(),
                    Err(SimulationError::InvalidInitialState { compartment }),
                    "{compartment} = {invalid}"
                );
            }
        }
    }

    #[test]
    fn sir_rejects_zero_or_overflowed_total_population() {
        let zero = SirState {
            susceptible: 0.0,
            infected: 0.0,
            recovered: 0.0,
        };
        assert_eq!(
            zero.initial_population(),
            Err(SimulationError::InvalidPopulation)
        );

        // Each compartment is finite; only the sum overflows.
        let overflowed = SirState {
            susceptible: f64::MAX,
            infected: f64::MAX,
            recovered: 0.0,
        };
        assert_eq!(
            overflowed.initial_population(),
            Err(SimulationError::InvalidPopulation)
        );
    }

    #[test]
    fn sir_half_step_run_records_the_same_times() {
        let parameters = SirParameters {
            beta_per_day: 0.3,
            gamma_per_day: 0.1,
        };
        // Whole and fractional step counts, including quotients that round
        // across a whole number, with several recording intervals.
        let grids = [
            (10.0, 3.0),
            (7.3, 0.25),
            (1.0, 0.1),
            (609.6, 0.6),
            (2.7, 0.3),
            (40.0, 1.0),
        ];
        for (duration, step) in grids {
            for record_every in [1, 2, 3, 7, 1000] {
                let settings = SimulationSettings::new(duration, step)
                    .unwrap()
                    .with_record_every(record_every)
                    .unwrap();
                let run = simulate_sir(parameters, VALID_STATE, settings).unwrap();
                let half = simulate_sir(parameters, VALID_STATE, settings.halved()).unwrap();
                let times = |result: &SimulationResult<SirState>| -> Vec<f64> {
                    result.trajectory.iter().map(|o| o.time_days).collect()
                };
                assert_eq!(
                    times(&run),
                    times(&half),
                    "duration {duration}, step {step}, every {record_every}"
                );
                assert!(
                    half.diagnostics.steps_executed >= 2 * run.diagnostics.steps_executed - 1,
                    "duration {duration}, step {step}"
                );
            }
        }
    }
}
