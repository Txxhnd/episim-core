use crate::error::SimulationError;
use crate::result::Observation;

/// Model-derived summary metrics of a completed simulation.
///
/// These are properties of the simulated model under its assumptions, not
/// estimates derived from observed surveillance data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EpidemicMetrics {
    /// Theoretical basic reproduction number `R0 = beta / gamma`
    /// (dimensionless) for a fully susceptible population.
    pub basic_reproduction_number: f64,
    /// Largest recorded infectious population, in individuals. Its precision
    /// is limited by the timestep because only recorded observations count.
    pub peak_infectious: f64,
    /// Time in days of the earliest recorded observation at the infectious
    /// peak; day zero if the infectious population never exceeds its initial
    /// value.
    pub peak_infectious_time_days: f64,
    /// New infections between the start and the final simulation time `T`,
    /// in individuals: `S(0) - S(T)`. Infections already present at the start
    /// are not counted.
    pub cumulative_infections_since_start: f64,
    /// Fraction of the initial population newly infected by the final
    /// simulation time: `(S(0) - S(T)) / N(0)`.
    pub cumulative_incidence_fraction: f64,
}

/// Theoretical `R0 = beta / gamma`, or [`SimulationError::MetricOverflow`] if
/// the quotient is not finite.
pub(crate) fn basic_reproduction_number(
    beta_per_day: f64,
    gamma_per_day: f64,
) -> Result<f64, SimulationError> {
    let value = beta_per_day / gamma_per_day;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SimulationError::MetricOverflow)
    }
}

/// Derives the metrics from a recorded trajectory, returning
/// [`SimulationError::MetricOverflow`] if any input or result is not finite.
pub(crate) fn calculate_metrics<S, FS, FI>(
    trajectory: &[Observation<S>],
    initial_population: f64,
    beta_per_day: f64,
    gamma_per_day: f64,
    susceptible_of: FS,
    infectious_of: FI,
) -> Result<EpidemicMetrics, SimulationError>
where
    FS: Fn(&S) -> f64,
    FI: Fn(&S) -> f64,
{
    let finite = |value: f64| {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(SimulationError::MetricOverflow)
        }
    };
    let (Some(first), Some(last)) = (trajectory.first(), trajectory.last()) else {
        return Err(SimulationError::MetricOverflow);
    };

    let mut peak_infectious = finite(infectious_of(&first.state))?;
    let mut peak_infectious_time_days = first.time_days;
    for observation in trajectory {
        let infectious = finite(infectious_of(&observation.state))?;
        // Strictly greater, so a tie keeps the earliest observation.
        if infectious > peak_infectious {
            peak_infectious = infectious;
            peak_infectious_time_days = observation.time_days;
        }
    }

    let cumulative_infections_since_start =
        finite(susceptible_of(&first.state) - susceptible_of(&last.state))?;
    Ok(EpidemicMetrics {
        basic_reproduction_number: basic_reproduction_number(beta_per_day, gamma_per_day)?,
        peak_infectious,
        peak_infectious_time_days,
        cumulative_infections_since_start,
        cumulative_incidence_fraction: finite(
            cumulative_infections_since_start / initial_population,
        )?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (time in days, susceptible, infectious); the infectious peak of 80 is
    /// reached on day 1 and held on day 2.
    const EPIDEMIC: [(f64, f64, f64); 4] = [
        (0.0, 990.0, 10.0),
        (1.0, 900.0, 80.0),
        (2.0, 700.0, 80.0),
        (3.0, 650.0, 40.0),
    ];

    /// Metrics of a hand-built trajectory in a population of 1000.
    fn metrics(
        points: &[(f64, f64, f64)],
        beta_per_day: f64,
        gamma_per_day: f64,
    ) -> Result<EpidemicMetrics, SimulationError> {
        let trajectory: Vec<_> = points
            .iter()
            .map(|&(time_days, susceptible, infectious)| Observation {
                time_days,
                state: [susceptible, infectious],
            })
            .collect();
        calculate_metrics(
            &trajectory,
            1000.0,
            beta_per_day,
            gamma_per_day,
            |state| state[0],
            |state| state[1],
        )
    }

    #[test]
    fn metrics_calculate_basic_reproduction_number() {
        let metrics = metrics(&EPIDEMIC, 0.5, 0.25).unwrap();
        assert_eq!(metrics.basic_reproduction_number, 2.0);
    }

    #[test]
    fn metrics_calculate_peak_infectious_and_earliest_peak_time() {
        let metrics = metrics(&EPIDEMIC, 0.5, 0.25).unwrap();
        assert_eq!(metrics.peak_infectious, 80.0);
        assert_eq!(metrics.peak_infectious_time_days, 1.0);
    }

    #[test]
    fn metrics_calculate_cumulative_infections_since_start() {
        let metrics = metrics(&EPIDEMIC, 0.5, 0.25).unwrap();
        assert_eq!(metrics.cumulative_infections_since_start, 340.0);
    }

    #[test]
    fn metrics_calculate_cumulative_incidence_fraction() {
        let metrics = metrics(&EPIDEMIC, 0.5, 0.25).unwrap();
        assert_eq!(metrics.cumulative_incidence_fraction, 0.34);
    }

    #[test]
    fn metrics_zero_transmission_returns_zero_r0() {
        let metrics = metrics(&EPIDEMIC, 0.0, 0.25).unwrap();
        assert_eq!(metrics.basic_reproduction_number, 0.0);
    }

    #[test]
    fn metrics_reject_non_finite_reproduction_number() {
        // Both rates are valid inputs, but their quotient overflows to infinity.
        assert_eq!(
            metrics(&EPIDEMIC, 10.0, f64::MIN_POSITIVE),
            Err(SimulationError::MetricOverflow)
        );
    }

    #[test]
    fn metrics_handle_zero_peak_infection() {
        let no_infection = [(0.0, 1000.0, 0.0), (1.0, 1000.0, 0.0), (2.0, 1000.0, 0.0)];
        let metrics = metrics(&no_infection, 0.5, 0.25).unwrap();
        assert_eq!(metrics.peak_infectious, 0.0);
        assert_eq!(metrics.peak_infectious_time_days, 0.0);
        assert_eq!(metrics.cumulative_infections_since_start, 0.0);
    }

    #[test]
    fn metrics_reject_non_finite_trajectory_values() {
        let nan_infectious = [
            (0.0, 990.0, 10.0),
            (1.0, 900.0, f64::NAN),
            (2.0, 700.0, 80.0),
        ];
        let infinite_susceptible = [(0.0, 990.0, 10.0), (1.0, f64::INFINITY, 80.0)];
        // No observations at all leaves nothing to derive a finite metric from.
        for points in [&nan_infectious[..], &infinite_susceptible[..], &[]] {
            assert_eq!(
                metrics(points, 0.5, 0.25),
                Err(SimulationError::MetricOverflow)
            );
        }
    }
}
