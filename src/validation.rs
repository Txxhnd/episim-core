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

/// Largest accepted product of the model's fastest decay rate and the step.
///
/// Classical RK4 multiplies a quantity decaying at rate `r` by
/// `1 - z + z^2/2 - z^3/6 + z^4/24` per step, with `z = r * step`. That factor
/// is 1/3 at `z = 2`, climbs steeply afterwards, and passes 1 near 2.785,
/// where decay turns into growth without any negative value or population
/// drift. In randomised comparisons against fine-step runs, no run with
/// `z <= 2` moved a compartment in the wrong direction or missed the final
/// incidence by more than 0.1; such runs appear from about 2.2.
const MAX_RATE_TIMES_STEP: f64 = 2.0;

/// Rejects a timestep that is too large for `fastest_rate_per_day`, the
/// fastest rate at which the model's compartments can decay.
///
/// This guards against instability; it does not make a coarse step accurate.
pub(crate) fn stable_step(
    step_days: f64,
    fastest_rate_per_day: f64,
) -> Result<(), SimulationError> {
    if fastest_rate_per_day * step_days > MAX_RATE_TIMES_STEP {
        return Err(SimulationError::TimeStepTooLarge {
            max_step_days: MAX_RATE_TIMES_STEP / fastest_rate_per_day,
        });
    }
    Ok(())
}
