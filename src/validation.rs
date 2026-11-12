use crate::error::SimulationError;

/// Accepts a rate that is finite and greater than or equal to zero.
pub(crate) fn non_negative_rate(name: &'static str, value: f64) -> Result<(), SimulationError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(SimulationError::InvalidParameter { name })
    }
}

/// Accepts a rate that is finite and strictly greater than zero.
pub(crate) fn positive_rate(name: &'static str, value: f64) -> Result<(), SimulationError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(SimulationError::InvalidParameter { name })
    }
}

/// Validates every compartment, then returns the finite, strictly positive total.
pub(crate) fn initial_population(
    compartments: &[(&'static str, f64)],
) -> Result<f64, SimulationError> {
    for &(compartment, value) in compartments {
        if !(value.is_finite() && value >= 0.0) {
            return Err(SimulationError::InvalidInitialState { compartment });
        }
    }
    // Finite compartments can still sum to infinity, so the total is checked separately.
    let population: f64 = compartments.iter().map(|&(_, value)| value).sum();
    if population.is_finite() && population > 0.0 {
        Ok(population)
    } else {
        Err(SimulationError::InvalidPopulation)
    }
}
