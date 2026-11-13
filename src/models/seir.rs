use crate::accuracy::StepErrorEstimate;
use crate::error::SimulationError;
use crate::metrics::{
    EpidemicMetrics, basic_reproduction_number, effective_reproduction_number, final_size_fraction,
    herd_immunity_threshold,
};
use crate::models::{Model, estimate, simulate};
use crate::result::SimulationResult;
use crate::settings::SimulationSettings;
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

/// Population counts in each SEIR compartment, used by both
/// [`simulate_seir`] and [`simulate_seirs`](crate::simulate_seirs).
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
    /// Individuals who have recovered and are immune: for good in SEIR, until
    /// their immunity wanes in SEIRS.
    pub recovered: f64,
}

impl SeirParameters {
    pub(crate) fn validate(&self) -> Result<(), SimulationError> {
        validation::non_negative_rate("beta_per_day", self.beta_per_day)?;
        validation::positive_rate("sigma_per_day", self.sigma_per_day)?;
        validation::positive_rate("gamma_per_day", self.gamma_per_day)
    }

    /// Model-derived effective reproduction number in `state`: `R0 * S / N`,
    /// where `N` is the total of `state`.
    ///
    /// It is the average number of people one infected person infects while a
    /// fraction `S / N` of the population is still susceptible. Because
    /// `d(E + I)/dt = gamma * I * (effective R - 1)`, the exposed and infected
    /// together grow while it is above 1 and decline while it is below 1. For
    /// a state returned by [`simulate_seir`], `N` equals the initial population
    /// to within the population drift tolerance of `1e-9`.
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
    pub fn effective_reproduction_number(&self, state: &SeirState) -> Result<f64, SimulationError> {
        self.validate()?;
        let r0 = basic_reproduction_number(self.beta_per_day, self.gamma_per_day)?;
        effective_reproduction_number(
            r0,
            state.susceptible,
            state.susceptible + state.exposed + state.infected + state.recovered,
        )
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
    simulate(&parameters, initial_state, settings)
}

/// Estimates the error that the timestep introduces into [`simulate_seir`] for
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
/// Returns any [`SimulationError`] that [`simulate_seir`] returns for either run.
///
/// # Examples
///
/// ```
/// use episim_core::{SeirParameters, SeirState, SimulationSettings, estimate_seir_step_error};
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
///
/// let coarse = estimate_seir_step_error(parameters, initial_state, SimulationSettings::new(160.0, 1.0)?)?;
/// let fine = estimate_seir_step_error(parameters, initial_state, SimulationSettings::new(160.0, 0.25)?)?;
///
/// // A quarter of the step gives roughly 4^4 = 256 times less error.
/// assert!(fine.max_absolute_error < coarse.max_absolute_error / 100.0);
/// assert!(fine.max_relative_error < 1e-6);
/// # Ok::<(), episim_core::SimulationError>(())
/// ```
pub fn estimate_seir_step_error(
    parameters: SeirParameters,
    initial_state: SeirState,
    settings: SimulationSettings,
) -> Result<StepErrorEstimate, SimulationError> {
    estimate(&parameters, initial_state, settings)
}

impl Model<5> for SeirParameters {
    type State = SeirState;
    const COMPARTMENTS: &'static [&'static str] =
        &["susceptible", "exposed", "infected", "recovered"];

    fn validate(&self) -> Result<(), SimulationError> {
        SeirParameters::validate(self)
    }

    fn population(state: &SeirState) -> Result<f64, SimulationError> {
        state.initial_population()
    }

    fn beta(&self) -> f64 {
        self.beta_per_day
    }

    fn gamma(&self) -> f64 {
        self.gamma_per_day
    }

    fn fastest_rate(&self) -> f64 {
        let (beta_per_day, sigma_per_day, gamma_per_day) =
            (self.beta_per_day, self.sigma_per_day, self.gamma_per_day);
        // Exposed and infected feed each other (E -> I by progression, I -> E by
        // transmission), so together they decay faster than sigma or gamma alone.
        // This is the fast eigenvalue of that pair in a fully susceptible
        // population, where the coupling is strongest.
        let coupled_rate = (sigma_per_day + gamma_per_day) / 2.0
            + (((sigma_per_day - gamma_per_day) / 2.0).powi(2) + sigma_per_day * beta_per_day)
                .sqrt();
        // Susceptibles empty at up to beta per day.
        beta_per_day.max(coupled_rate)
    }

    fn to_array(state: &SeirState) -> [f64; 5] {
        [
            state.susceptible,
            state.exposed,
            state.infected,
            state.recovered,
            0.0,
        ]
    }

    fn from_array([susceptible, exposed, infected, recovered, _]: [f64; 5]) -> SeirState {
        SeirState {
            susceptible,
            exposed,
            infected,
            recovered,
        }
    }

    fn derivative(
        &self,
        population: f64,
        [susceptible, exposed, infected, _, _]: [f64; 5],
    ) -> [f64; 5] {
        // Dividing first avoids multiplying two large population values.
        let infection_flow = self.beta_per_day * (susceptible / population) * infected;
        let progression_flow = self.sigma_per_day * exposed;
        let recovery_flow = self.gamma_per_day * infected;
        [
            -infection_flow,
            infection_flow - progression_flow,
            progression_flow - recovery_flow,
            recovery_flow,
            infection_flow,
        ]
    }

    fn infectious(state: &SeirState) -> f64 {
        state.infected
    }

    fn model_metrics(&self, metrics: &mut EpidemicMetrics, initial: &SeirState, population: f64) {
        let r0 = metrics.basic_reproduction_number;
        metrics.herd_immunity_threshold = Some(herd_immunity_threshold(r0));
        metrics.final_size_fraction = Some(final_size_fraction(
            initial.susceptible,
            initial.exposed + initial.infected,
            population,
            r0,
        ));
    }
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
