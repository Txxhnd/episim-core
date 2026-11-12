use crate::error::SimulationError;
use crate::metrics::calculate_metrics;
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
/// invalid, if the timestep is too large for stable integration with one of
/// the rates, if the settings require more steps than their limit or a
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
    validation::stable_step(
        settings.step_days(),
        &[
            ("beta_per_day", beta_per_day),
            ("gamma_per_day", gamma_per_day),
        ],
    )?;

    let output = run_simulation(
        [
            initial_state.susceptible,
            initial_state.infected,
            initial_state.recovered,
        ],
        population,
        settings,
        |[susceptible, infected, _]| {
            // Dividing first avoids multiplying two large population values.
            let infection_flow = beta_per_day * (susceptible / population) * infected;
            let recovery_flow = gamma_per_day * infected;
            [
                -infection_flow,
                infection_flow - recovery_flow,
                recovery_flow,
            ]
        },
        |[susceptible, infected, recovered]| SirState {
            susceptible,
            infected,
            recovered,
        },
        ["susceptible", "infected", "recovered"],
    )?;
    let metrics = calculate_metrics(
        &output.trajectory,
        population,
        beta_per_day,
        gamma_per_day,
        |state| state.susceptible,
        |state| state.infected,
    )?;
    Ok(SimulationResult {
        trajectory: output.trajectory,
        metrics,
        diagnostics: output.diagnostics,
    })
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
}
