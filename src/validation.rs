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

/// Largest `rate * step` for which classical RK4 still damps a decaying
/// quantity. Beyond it the per-step factor `1 - z + z^2/2 - z^3/6 + z^4/24`
/// exceeds 1 (its crossing point is near 2.7853), so decay turns into growth.
const RK4_STABILITY_LIMIT: f64 = 2.785;

/// Rejects a timestep that is too large for any of the given rates.
///
/// Each rate bounds how fast one compartment can empty, so this guards
/// against instability; it does not make a coarse step accurate.
pub(crate) fn stable_step(
    step_days: f64,
    rates: &[(&'static str, f64)],
) -> Result<(), SimulationError> {
    for &(parameter, rate) in rates {
        if rate * step_days > RK4_STABILITY_LIMIT {
            return Err(SimulationError::TimeStepTooLarge { parameter });
        }
    }
    Ok(())
}
