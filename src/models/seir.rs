use crate::error::SimulationError;
use crate::metrics::calculate_metrics;
use crate::result::SimulationResult;
use crate::settings::SimulationSettings;
use crate::simulation::driver::run_simulation;
use crate::validation;

/// Rate parameters of the SEIR model, constant for the whole simulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeirParameters {
    /// Transmission rate in inverse days; must be finite and non-negative.
    pub beta_per_day: f64,
    /// Rate of progression from exposed to infectious in inverse days; must
    /// be finite and strictly positive.
    pub sigma_per_day: f64,
    /// Recovery rate in inverse days; must be finite and strictly positive.
    pub gamma_per_day: f64,
}

/// Population counts in each SEIR compartment.
///
/// Counts are expected values and may be fractional. Every compartment must
/// be finite and non-negative, and the total must be finite and strictly
/// positive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeirState {
    /// Individuals who can become infected.
    pub susceptible: f64,
    /// Infected individuals who are not yet infectious.
    pub exposed: f64,
    /// Individuals who are currently infectious.
    pub infected: f64,
    /// Individuals who have recovered and are assumed immune.
    pub recovered: f64,
}

impl SeirParameters {
    pub(crate) fn validate(&self) -> Result<(), SimulationError> {
        validation::non_negative_rate("beta_per_day", self.beta_per_day)?;
        validation::positive_rate("sigma_per_day", self.sigma_per_day)?;
        validation::positive_rate("gamma_per_day", self.gamma_per_day)
    }
}

impl SeirState {
    /// Validates the state and returns the total population it defines.
    pub(crate) fn initial_population(&self) -> Result<f64, SimulationError> {
        validation::initial_population(&[
            ("susceptible", self.susceptible),
            ("exposed", self.exposed),
            ("infected", self.infected),
            ("recovered", self.recovered),
        ])
    }
}

/// Simulates the deterministic SEIR model with a fixed-step classical RK4 integrator.
///
/// The model assumes a closed, homogeneously mixing population with constant
/// parameters and permanent immunity. Only the infected compartment
/// transmits; the exposed compartment is a latent, non-infectious stage.
/// Transmission is frequency-dependent, with `N` fixed at the initial total
/// population:
///
/// ```text
/// dS/dt = -beta * S * I / N
/// dE/dt =  beta * S * I / N - sigma * E
/// dI/dt =  sigma * E - gamma * I
/// dR/dt =  gamma * I
/// ```
///
/// The returned trajectory starts at day zero, holds one observation per
/// integration step, and ends exactly at the requested duration.
///
/// # Errors
///
/// Returns a [`SimulationError`] if a parameter or the initial state is
/// invalid, if the settings require more steps than their limit or a timestep
/// too small to advance the clock, if the integration produces a non-finite
/// value, a materially negative compartment, or excessive population drift,
/// if a metric is not finite, or if the trajectory cannot be allocated.
///
/// # Examples
///
/// ```
/// use episim_core::{SeirParameters, SeirState, SimulationSettings, simulate_seir};
///
/// let parameters = SeirParameters {
///     beta_per_day: 0.5,
///     sigma_per_day: 0.2,
///     gamma_per_day: 0.1,
/// };
/// let initial_state = SeirState {
///     susceptible: 980.0,
///     exposed: 10.0,
///     infected: 10.0,
///     recovered: 0.0,
/// };
/// let settings = SimulationSettings::new(160.0, 0.25)?;
///
/// let result = simulate_seir(parameters, initial_state, settings)?;
///
/// assert_eq!(result.trajectory.last().unwrap().time_days, 160.0);
/// assert!(result.metrics.peak_infectious > initial_state.infected);
/// # Ok::<(), episim_core::SimulationError>(())
/// ```
pub fn simulate_seir(
    parameters: SeirParameters,
    initial_state: SeirState,
    settings: SimulationSettings,
) -> Result<SimulationResult<SeirState>, SimulationError> {
    parameters.validate()?;
    let population = initial_state.initial_population()?;
    let SeirParameters {
        beta_per_day,
        sigma_per_day,
        gamma_per_day,
    } = parameters;

    let output = run_simulation(
        [
            initial_state.susceptible,
            initial_state.exposed,
            initial_state.infected,
            initial_state.recovered,
        ],
        population,
        settings,
        |[susceptible, exposed, infected, _]| {
            // Dividing first avoids multiplying two large population values.
            let infection_flow = beta_per_day * (susceptible / population) * infected;
            let progression_flow = sigma_per_day * exposed;
            let recovery_flow = gamma_per_day * infected;
            [
                -infection_flow,
                infection_flow - progression_flow,
                progression_flow - recovery_flow,
                recovery_flow,
            ]
        },
        |[susceptible, exposed, infected, recovered]| SeirState {
            susceptible,
            exposed,
            infected,
            recovered,
        },
        ["susceptible", "exposed", "infected", "recovered"],
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

    const INVALID_RATES: [f64; 5] = [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

    const VALID_PARAMETERS: SeirParameters = SeirParameters {
        beta_per_day: 0.3,
        sigma_per_day: 0.2,
        gamma_per_day: 0.1,
    };

    const VALID_STATE: SeirState = SeirState {
        susceptible: 980.0,
        exposed: 10.0,
        infected: 10.0,
        recovered: 0.0,
    };

    #[test]
    fn seir_accepts_zero_beta_and_positive_sigma_gamma() {
        let parameters = SeirParameters {
            beta_per_day: 0.0,
            ..VALID_PARAMETERS
        };
        assert_eq!(parameters.validate(), Ok(()));
        assert_eq!(VALID_STATE.initial_population(), Ok(1000.0));
    }

    #[test]
    fn seir_rejects_invalid_sigma() {
        for sigma_per_day in INVALID_RATES {
            let parameters = SeirParameters {
                sigma_per_day,
                ..VALID_PARAMETERS
            };
            assert_eq!(
                parameters.validate(),
                Err(SimulationError::InvalidParameter {
                    name: "sigma_per_day"
                }),
                "sigma {sigma_per_day}"
            );
        }
    }

    #[test]
    fn seir_rejects_invalid_gamma() {
        for gamma_per_day in INVALID_RATES {
            let parameters = SeirParameters {
                gamma_per_day,
                ..VALID_PARAMETERS
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
    fn seir_rejects_negative_or_non_finite_compartments() {
        for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let cases = [
                (
                    "susceptible",
                    SeirState {
                        susceptible: invalid,
                        ..VALID_STATE
                    },
                ),
                (
                    "exposed",
                    SeirState {
                        exposed: invalid,
                        ..VALID_STATE
                    },
                ),
                (
                    "infected",
                    SeirState {
                        infected: invalid,
                        ..VALID_STATE
                    },
                ),
                (
                    "recovered",
                    SeirState {
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
    fn seir_rejects_zero_or_overflowed_total_population() {
        let zero = SeirState {
            susceptible: 0.0,
            exposed: 0.0,
            infected: 0.0,
            recovered: 0.0,
        };
        assert_eq!(
            zero.initial_population(),
            Err(SimulationError::InvalidPopulation)
        );

        // Each compartment is finite; only the sum overflows.
        let overflowed = SeirState {
            susceptible: f64::MAX,
            exposed: 0.0,
            infected: 0.0,
            recovered: f64::MAX,
        };
        assert_eq!(
            overflowed.initial_population(),
            Err(SimulationError::InvalidPopulation)
        );
    }
}
