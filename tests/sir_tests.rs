use episim_core::{
    Observation, SimulationError, SimulationSettings, SirParameters, SirState, simulate_sir,
};

const PARAMETERS: SirParameters = SirParameters {
    beta_per_day: 0.3,
    gamma_per_day: 0.1,
};

const INITIAL: SirState = SirState {
    susceptible: 990.0,
    infected: 10.0,
    recovered: 0.0,
};

fn settings(duration_days: f64, step_days: f64) -> SimulationSettings {
    SimulationSettings::new(duration_days, step_days).unwrap()
}

#[test]
fn sir_returns_initial_and_final_observations() {
    let result = simulate_sir(PARAMETERS, INITIAL, settings(10.0, 3.0)).unwrap();

    assert_eq!(result.trajectory.len(), 5);
    assert_eq!(
        result.trajectory[0],
        Observation {
            time_days: 0.0,
            state: INITIAL
        }
    );
    let last = result.trajectory[4];
    assert_eq!(last.time_days, 10.0);
    // With beta > gamma the epidemic grows: susceptibles fall, recoveries accumulate.
    assert!(last.state.susceptible < INITIAL.susceptible);
    assert!(last.state.recovered > INITIAL.recovered);
    assert_eq!(result.diagnostics.steps_executed, 4);
    assert!((result.metrics.basic_reproduction_number - 3.0).abs() < 1e-12);
}

#[test]
fn sir_zero_transmission_keeps_susceptible_constant() {
    let parameters = SirParameters {
        beta_per_day: 0.0,
        ..PARAMETERS
    };
    let result = simulate_sir(parameters, INITIAL, settings(30.0, 0.5)).unwrap();

    for observation in &result.trajectory {
        assert_eq!(observation.state.susceptible, INITIAL.susceptible);
    }
    assert_eq!(result.metrics.basic_reproduction_number, 0.0);
    assert_eq!(result.metrics.cumulative_infections_since_start, 0.0);
}

#[test]
fn sir_zero_initial_infection_has_no_new_infections() {
    let initial = SirState {
        susceptible: 1000.0,
        infected: 0.0,
        recovered: 0.0,
    };
    let result = simulate_sir(PARAMETERS, initial, settings(30.0, 0.5)).unwrap();

    for observation in &result.trajectory {
        assert_eq!(observation.state, initial);
    }
    assert_eq!(result.metrics.cumulative_infections_since_start, 0.0);
    assert_eq!(result.metrics.peak_infectious, 0.0);
    assert_eq!(result.metrics.peak_infectious_time_days, 0.0);
}

#[test]
fn sir_beta_zero_matches_analytical_recovery_solution() {
    let parameters = SirParameters {
        beta_per_day: 0.0,
        gamma_per_day: 0.1,
    };
    let initial = SirState {
        susceptible: 900.0,
        infected: 100.0,
        recovered: 5.0,
    };
    let result = simulate_sir(parameters, initial, settings(30.0, 0.5)).unwrap();

    // RK4's global error for pure decay with gamma * h = 0.05 peaks near
    // 2e-6 individuals here; 1e-5 leaves a fivefold margin.
    let tolerance = 1e-5;
    for observation in &result.trajectory {
        let infected = initial.infected * (-parameters.gamma_per_day * observation.time_days).exp();
        let recovered = initial.recovered + initial.infected - infected;
        let state = observation.state;
        assert_eq!(state.susceptible, initial.susceptible);
        assert!(
            (state.infected - infected).abs() < tolerance,
            "infected at day {}: {} vs {infected}",
            observation.time_days,
            state.infected
        );
        assert!(
            (state.recovered - recovered).abs() < tolerance,
            "recovered at day {}: {} vs {recovered}",
            observation.time_days,
            state.recovered
        );
    }
}

#[test]
fn sir_preserves_total_population_within_tolerance() {
    let result = simulate_sir(PARAMETERS, INITIAL, settings(160.0, 0.25)).unwrap();

    for observation in &result.trajectory {
        let state = observation.state;
        let total = state.susceptible + state.infected + state.recovered;
        assert!(
            ((total - 1000.0) / 1000.0).abs() <= 1e-9,
            "total {total} at day {}",
            observation.time_days
        );
    }
    assert!(result.diagnostics.maximum_relative_population_drift <= 1e-9);
}

#[test]
fn sir_returns_error_for_invalid_parameters_and_state() {
    let run = |parameters, initial| simulate_sir(parameters, initial, settings(10.0, 1.0));

    let negative_beta = SirParameters {
        beta_per_day: -0.3,
        ..PARAMETERS
    };
    assert_eq!(
        run(negative_beta, INITIAL).unwrap_err(),
        SimulationError::InvalidParameter {
            name: "beta_per_day"
        }
    );

    let zero_gamma = SirParameters {
        gamma_per_day: 0.0,
        ..PARAMETERS
    };
    assert_eq!(
        run(zero_gamma, INITIAL).unwrap_err(),
        SimulationError::InvalidParameter {
            name: "gamma_per_day"
        }
    );

    let negative_infected = SirState {
        infected: -10.0,
        ..INITIAL
    };
    assert_eq!(
        run(PARAMETERS, negative_infected).unwrap_err(),
        SimulationError::InvalidInitialState {
            compartment: "infected"
        }
    );

    let empty = SirState {
        susceptible: 0.0,
        infected: 0.0,
        recovered: 0.0,
    };
    assert_eq!(
        run(PARAMETERS, empty).unwrap_err(),
        SimulationError::InvalidPopulation
    );
}

#[test]
fn sir_final_time_equals_requested_duration() {
    // Durations that are not a whole number of steps, including two whose
    // step count is sensitive to floating-point rounding.
    for (duration, step) in [(10.0, 3.0), (7.3, 0.25), (1.0, 0.1), (609.6, 0.6)] {
        let result = simulate_sir(PARAMETERS, INITIAL, settings(duration, step)).unwrap();
        let times: Vec<f64> = result.trajectory.iter().map(|o| o.time_days).collect();

        assert_eq!(times[0], 0.0);
        assert_eq!(*times.last().unwrap(), duration);
        assert!(
            times.windows(2).all(|pair| pair[0] < pair[1]),
            "times must strictly increase for duration {duration}"
        );
        assert_eq!(result.diagnostics.steps_executed, times.len() - 1);
    }
}
