use crate::error::SimulationError;
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
