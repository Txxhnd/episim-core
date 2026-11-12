use episim_core::{
    EpidemicMetrics, SeirParameters, SeirState, SimulationDiagnostics, SimulationSettings,
    SirParameters, SirState, simulate_seir, simulate_sir,
};

const DRIFT_THRESHOLD: f64 = 1e-9;

/// (transmission rate, population, timestep): a typical outbreak, the same
/// outbreak in a population of a billion with a coarse step, and pure recovery.
const SCENARIOS: [(f64, f64, f64); 3] = [(0.3, 1e3, 0.25), (0.5, 1e9, 1.0), (0.0, 1e3, 0.5)];

fn settings(step_days: f64) -> SimulationSettings {
    SimulationSettings::new(200.0, step_days).unwrap()
}

#[test]
fn sir_population_drift_stays_within_release_threshold() {
    for (beta_per_day, population, step) in SCENARIOS {
        let parameters = SirParameters {
            beta_per_day,
            gamma_per_day: 0.1,
        };
        let initial = SirState {
            susceptible: population - 10.0,
            infected: 10.0,
            recovered: 0.0,
        };
        let result = simulate_sir(parameters, initial, settings(step)).unwrap();

        for observation in &result.trajectory {
            let state = observation.state;
            let total = state.susceptible + state.infected + state.recovered;
            assert!(
                ((total - population) / population).abs() <= DRIFT_THRESHOLD,
                "population {population}: total {total} at day {}",
                observation.time_days
            );
        }
        assert!(result.diagnostics.maximum_relative_population_drift <= DRIFT_THRESHOLD);
    }
}

#[test]
fn seir_population_drift_stays_within_release_threshold() {
    for (beta_per_day, population, step) in SCENARIOS {
        let parameters = SeirParameters {
            beta_per_day,
            sigma_per_day: 0.2,
            gamma_per_day: 0.1,
        };
        let initial = SeirState {
            susceptible: population - 20.0,
            exposed: 10.0,
            infected: 10.0,
            recovered: 0.0,
        };
        let result = simulate_seir(parameters, initial, settings(step)).unwrap();

        for observation in &result.trajectory {
            let state = observation.state;
            let total = state.susceptible + state.exposed + state.infected + state.recovered;
            assert!(
                ((total - population) / population).abs() <= DRIFT_THRESHOLD,
                "population {population}: total {total} at day {}",
                observation.time_days
            );
        }
        assert!(result.diagnostics.maximum_relative_population_drift <= DRIFT_THRESHOLD);
    }
}

fn assert_finite_and_plausible(
    metrics: EpidemicMetrics,
    diagnostics: SimulationDiagnostics,
    observations: usize,
    infectious_must_grow: bool,
) {
    for value in [
        metrics.basic_reproduction_number,
        metrics.peak_infectious,
        metrics.peak_infectious_time_days,
        metrics.cumulative_infections_since_start,
        metrics.cumulative_incidence_fraction,
        diagnostics.maximum_absolute_population_drift,
        diagnostics.maximum_relative_population_drift,
    ] {
        assert!(value.is_finite(), "{metrics:?} {diagnostics:?}");
    }
    // Every scenario starts with 10 infected.
    if infectious_must_grow {
        assert!(metrics.peak_infectious > 10.0, "{metrics:?}");
        assert!(metrics.peak_infectious_time_days > 0.0, "{metrics:?}");
    } else {
        assert_eq!(metrics.peak_infectious, 10.0);
        assert_eq!(metrics.peak_infectious_time_days, 0.0);
    }
    assert!(
        (0.0..=1.0).contains(&metrics.cumulative_incidence_fraction),
        "{metrics:?}"
    );
    assert_eq!(diagnostics.steps_executed, observations - 1);
}

#[test]
fn metrics_and_diagnostics_remain_finite_for_reference_scenarios() {
    for (beta_per_day, population, step) in SCENARIOS {
        let sir = simulate_sir(
            SirParameters {
                beta_per_day,
                gamma_per_day: 0.1,
            },
            SirState {
                susceptible: population - 10.0,
                infected: 10.0,
                recovered: 0.0,
            },
            settings(step),
        )
        .unwrap();
        // Without transmission the SIR infected can only recover.
        assert_finite_and_plausible(
            sir.metrics,
            sir.diagnostics,
            sir.trajectory.len(),
            beta_per_day > 0.0,
        );

        let seir = simulate_seir(
            SeirParameters {
                beta_per_day,
                sigma_per_day: 0.2,
                gamma_per_day: 0.1,
            },
            SeirState {
                susceptible: population - 20.0,
                exposed: 10.0,
                infected: 10.0,
                recovered: 0.0,
            },
            settings(step),
        )
        .unwrap();
        // The initially exposed become infectious even without transmission.
        assert_finite_and_plausible(seir.metrics, seir.diagnostics, seir.trajectory.len(), true);
    }
}
