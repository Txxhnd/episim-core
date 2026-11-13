use crate::error::SimulationError;
use crate::result::Observation;

/// Model-derived summary metrics of a completed simulation.
///
/// These are properties of the simulated model under its assumptions, not
/// estimates derived from observed surveillance data.
///
/// A metric that only some models define is an `Option`, `None` for the
/// models where it does not apply. The struct is `#[non_exhaustive]`, so
/// metrics can be added without breaking code that reads it; it cannot be
/// built with a struct literal outside this crate.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct EpidemicMetrics {
    /// Theoretical basic reproduction number `R0 = beta / gamma`
    /// (dimensionless) for a fully susceptible population.
    pub basic_reproduction_number: f64,
    /// Fraction of the population that must be immune for the infected to
    /// decline: `1 - 1/R0`, or 0 when `R0 <= 1`. `Some` for SIR and SEIR.
    pub herd_immunity_threshold: Option<f64>,
    /// Largest recorded infectious population, in individuals. Its precision
    /// is limited by the timestep because only recorded observations count.
    pub peak_infectious: f64,
    /// Time in days of the earliest recorded observation within
    /// `64 * f64::EPSILON * max(N0, 1)` individuals of the infectious peak;
    /// day zero if the infectious population never exceeds its initial value
    /// by more than that. The tolerance keeps roundoff from choosing the time
    /// when the infectious population levels off at its peak.
    pub peak_infectious_time_days: f64,
    /// New infections between the start and the final simulation time `T`,
    /// in individuals: the cumulative infection count `C(T)` of the last
    /// observation, which for SIR and SEIR equals `S(0) - S(T)`. Infections
    /// already present at the start are not counted.
    pub cumulative_infections_since_start: f64,
    /// New infections by the final simulation time per head of the initial
    /// population: `C(T) / N(0)`.
    pub cumulative_incidence_fraction: f64,
    /// Fraction of the initial population newly infected by the time the
    /// epidemic is over, `(S(0) - S(inf)) / N(0)`, from the final size
    /// equation. `cumulative_incidence_fraction` approaches it as the
    /// simulated duration grows. `Some` for SIR and SEIR.
    pub final_size_fraction: Option<f64>,
    /// Fraction of the population infected once the model settles at its
    /// endemic equilibrium, `I* / N`, or 0 when `R0 <= 1`. It follows from the
    /// parameters alone, like `basic_reproduction_number`. `Some` for models
    /// in which people become susceptible again: SIS.
    pub endemic_prevalence: Option<f64>,
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

/// Herd-immunity threshold `1 - 1/R0`, or 0 when `R0 <= 1`.
pub(crate) fn herd_immunity_threshold(r0: f64) -> f64 {
    if r0 > 1.0 { 1.0 - 1.0 / r0 } else { 0.0 }
}

/// Model-derived effective reproduction number `R0 * S / N`, or an error if
/// `N` is not finite and positive or the result is not finite.
pub(crate) fn effective_reproduction_number(
    r0: f64,
    susceptible: f64,
    population: f64,
) -> Result<f64, SimulationError> {
    if !(population.is_finite() && population > 0.0) {
        return Err(SimulationError::InvalidPopulation);
    }
    let value = r0 * (susceptible / population);
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SimulationError::MetricOverflow)
    }
}

/// Fraction of the population `N` newly infected by the end of the epidemic,
/// `(S(0) - S(inf)) / N`, where `S(inf)` solves the final size equation
///
/// ```text
/// S(inf) = S(0) * exp(-R0 * (S(0) + seed - S(inf)) / N)
/// ```
///
/// and `seed` is the number infected, plus for SEIR exposed, at the start.
/// In both models `dS/dR = -R0 * S / N`, and everyone infected or exposed
/// eventually recovers, which gives the equation.
pub(crate) fn final_size_fraction(susceptible: f64, seed: f64, population: f64, r0: f64) -> f64 {
    // With no seed the equation has a second root that no epidemic reaches.
    if seed == 0.0 || r0 == 0.0 || susceptible == 0.0 {
        return 0.0;
    }
    let gap = |remaining: f64| {
        remaining - susceptible * (-r0 * (susceptible + seed - remaining) / population).exp()
    };
    // gap(0) < 0 <= gap(S(0)) and gap is concave, so exactly one root lies
    // between; bisection narrows it to adjacent floating-point values.
    let (mut low, mut high) = (0.0_f64, susceptible);
    loop {
        let middle = low + (high - low) / 2.0;
        if middle <= low || middle >= high {
            break;
        }
        if gap(middle) < 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    (susceptible - high) / population
}

/// Derives the metrics every model shares from a recorded trajectory,
/// returning [`SimulationError::MetricOverflow`] if any input or result is not
/// finite. The model-specific metrics are `None`; each model fills in its own.
pub(crate) fn calculate_metrics<S, FI>(
    trajectory: &[Observation<S>],
    initial_population: f64,
    beta_per_day: f64,
    gamma_per_day: f64,
    infectious_of: FI,
) -> Result<EpidemicMetrics, SimulationError>
where
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
    for observation in trajectory {
        peak_infectious = peak_infectious.max(finite(infectious_of(&observation.state))?);
    }
    // Where the infectious population levels off at its peak, roundoff alone
    // decides which observation is largest; take the first one that is
    // within roundoff of the peak. The tolerance is the driver's.
    let tolerance = 64.0 * f64::EPSILON * initial_population.max(1.0);
    let peak_infectious_time_days = trajectory
        .iter()
        .find(|observation| infectious_of(&observation.state) >= peak_infectious - tolerance)
        .map_or(first.time_days, |observation| observation.time_days);

    let cumulative_infections_since_start =
        finite(last.cumulative_infections - first.cumulative_infections)?;
    let r0 = basic_reproduction_number(beta_per_day, gamma_per_day)?;
    Ok(EpidemicMetrics {
        basic_reproduction_number: r0,
        herd_immunity_threshold: None,
        peak_infectious,
        peak_infectious_time_days,
        cumulative_infections_since_start,
        cumulative_incidence_fraction: finite(
            cumulative_infections_since_start / initial_population,
        )?,
        final_size_fraction: None,
        endemic_prevalence: None,
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

    /// Metrics of a hand-built trajectory in a population of 1000, with, as in
    /// SIR, the drop in susceptibles as the cumulative infections.
    fn metrics(
        points: &[(f64, f64, f64)],
        beta_per_day: f64,
        gamma_per_day: f64,
    ) -> Result<EpidemicMetrics, SimulationError> {
        let trajectory: Vec<_> = points
            .iter()
            .map(|&(time_days, susceptible, infectious)| Observation {
                time_days,
                cumulative_infections: points[0].1 - susceptible,
                state: [susceptible, infectious],
            })
            .collect();
        calculate_metrics(
            &trajectory,
            1000.0,
            beta_per_day,
            gamma_per_day,
            |state: &[f64; 2]| state[1],
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
    fn metrics_peak_time_ignores_roundoff_at_a_level_peak() {
        // The tolerance for a population of 1000 is 64 * EPSILON * 1000,
        // about 1.4e-11 individuals.
        let level = 666.666_666_666_666_7;
        let roundoff = [
            (0.0, 990.0, 10.0),
            (1.0, 500.0, level),
            (2.0, 400.0, level + 1e-12),
            (3.0, 350.0, level + 2e-12),
        ];
        let level_peak = metrics(&roundoff, 0.3, 0.1).unwrap();
        assert_eq!(level_peak.peak_infectious, level + 2e-12);
        assert_eq!(level_peak.peak_infectious_time_days, 1.0);

        // A rise larger than the tolerance is a real peak.
        let rising = [
            (0.0, 990.0, 10.0),
            (1.0, 500.0, level),
            (2.0, 400.0, level + 1e-10),
        ];
        assert_eq!(
            metrics(&rising, 0.3, 0.1)
                .unwrap()
                .peak_infectious_time_days,
            2.0
        );
    }

    #[test]
    fn metrics_calculate_cumulative_infections_since_start() {
        let metrics = metrics(&EPIDEMIC, 0.5, 0.25).unwrap();
        assert_eq!(metrics.cumulative_infections_since_start, 340.0);
    }

    #[test]
    fn metrics_count_cumulative_infections_from_the_accumulator() {
        // Susceptibles that stay level while infections accumulate: people
        // are infected and replaced, as when immunity wanes.
        let trajectory =
            [(0.0, 0.0), (1.0, 30.0), (2.0, 70.0)].map(|(time_days, infections)| Observation {
                time_days,
                cumulative_infections: infections,
                state: [500.0, 100.0],
            });
        let metrics = calculate_metrics(&trajectory, 1000.0, 0.5, 0.25, |state| state[1]).unwrap();
        assert_eq!(metrics.cumulative_infections_since_start, 70.0);
        assert_eq!(metrics.cumulative_incidence_fraction, 0.07);
    }

    #[test]
    fn metrics_calculate_cumulative_incidence_fraction() {
        let metrics = metrics(&EPIDEMIC, 0.5, 0.25).unwrap();
        assert_eq!(metrics.cumulative_incidence_fraction, 0.34);
    }

    #[test]
    fn metrics_leave_model_specific_metrics_to_the_model() {
        let metrics = metrics(&EPIDEMIC, 0.5, 0.25).unwrap();
        assert_eq!(metrics.herd_immunity_threshold, None);
        assert_eq!(metrics.final_size_fraction, None);
        assert_eq!(metrics.endemic_prevalence, None);
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

    #[test]
    fn metrics_herd_immunity_threshold_is_one_minus_inverse_r0() {
        assert_eq!(herd_immunity_threshold(2.0), 0.5);
        assert_eq!(herd_immunity_threshold(4.0), 0.75);
        // Without an epidemic threshold to cross, no immunity is needed.
        for r0 in [1.0, 0.8, 0.0] {
            assert_eq!(herd_immunity_threshold(r0), 0.0, "R0 {r0}");
        }
    }

    #[test]
    fn metrics_final_size_solves_the_final_size_equation() {
        // (susceptible, seed, recovered, R0) at the start.
        let cases = [
            (999.0, 1.0, 0.0, 2.0),
            (985.0, 10.0, 5.0, 3.0),
            (990.0, 10.0, 0.0, 0.5),
            (500.0, 1.0, 499.0, 4.0),
            (999_999_990.0, 10.0, 0.0, 1.5),
        ];
        for (susceptible, seed, recovered, r0) in cases {
            let population = susceptible + seed + recovered;
            let fraction = final_size_fraction(susceptible, seed, population, r0);
            // S(inf) = S(0) * exp(-R0 * (S(0) + seed - S(inf)) / N)
            let remaining = susceptible - fraction * population;
            let predicted =
                susceptible * (-r0 * (susceptible + seed - remaining) / population).exp();
            assert!(
                (remaining - predicted).abs() <= 1e-9 * population,
                "case {susceptible}, {seed}, {recovered}, {r0}: {remaining} vs {predicted}"
            );
            assert!(fraction > 0.0 && fraction <= susceptible / population);
        }
        // The textbook value: with R0 = 2 and a vanishing seed, z = 1 - exp(-2z)
        // gives z = 0.796812...
        let fraction = final_size_fraction(1e9 - 1.0, 1.0, 1e9, 2.0);
        assert!((fraction - 0.796_812).abs() < 1e-6, "{fraction}");
    }

    #[test]
    fn metrics_final_size_is_zero_without_seed_transmission_or_susceptibles() {
        // Nobody to start an epidemic, no transmission, or nobody to infect.
        assert_eq!(final_size_fraction(990.0, 0.0, 1000.0, 3.0), 0.0);
        assert_eq!(final_size_fraction(990.0, 10.0, 1000.0, 0.0), 0.0);
        assert_eq!(final_size_fraction(0.0, 10.0, 1000.0, 3.0), 0.0);
    }
}
