use crate::error::SimulationError;
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
