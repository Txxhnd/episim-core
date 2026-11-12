use episim_core::{
    Observation, SeirParameters, SeirState, SimulationError, SimulationSettings, simulate_seir,
};

const PARAMETERS: SeirParameters = SeirParameters {
    beta_per_day: 0.5,
    sigma_per_day: 0.2,
    gamma_per_day: 0.1,
};

const INITIAL: SeirState = SeirState {
    susceptible: 980.0,
    exposed: 10.0,
    infected: 10.0,
    recovered: 0.0,
};

fn settings(duration_days: f64, step_days: f64) -> SimulationSettings {
    SimulationSettings::new(duration_days, step_days).unwrap()
}

#[test]
fn seir_returns_initial_and_final_observations() {
    let result = simulate_seir(PARAMETERS, INITIAL, settings(10.0, 3.0)).unwrap();

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
    assert_eq!(result.metrics.basic_reproduction_number, 5.0);
}

#[test]
fn seir_zero_transmission_keeps_susceptible_constant() {
    let parameters = SeirParameters {
        beta_per_day: 0.0,
        ..PARAMETERS
    };
    let result = simulate_seir(parameters, INITIAL, settings(30.0, 0.5)).unwrap();

    for observation in &result.trajectory {
        assert_eq!(observation.state.susceptible, INITIAL.susceptible);
    }
    assert_eq!(result.metrics.basic_reproduction_number, 0.0);
    assert_eq!(result.metrics.cumulative_infections_since_start, 0.0);
}

#[test]
fn seir_zero_exposed_and_infected_has_no_new_infections() {
    let initial = SeirState {
        susceptible: 1000.0,
        exposed: 0.0,
        infected: 0.0,
        recovered: 0.0,
    };
    let result = simulate_seir(PARAMETERS, initial, settings(30.0, 0.5)).unwrap();

    for observation in &result.trajectory {
        assert_eq!(observation.state, initial);
    }
    assert_eq!(result.metrics.cumulative_infections_since_start, 0.0);
    assert_eq!(result.metrics.peak_infectious, 0.0);
    assert_eq!(result.metrics.peak_infectious_time_days, 0.0);
}

#[test]
fn seir_beta_zero_matches_analytical_exposed_and_infected_solution() {
    let parameters = SeirParameters {
        beta_per_day: 0.0,
        sigma_per_day: 0.2,
        gamma_per_day: 0.1,
    };
    let initial = SeirState {
        susceptible: 850.0,
        exposed: 50.0,
        infected: 100.0,
        recovered: 0.0,
    };
    let result = simulate_seir(parameters, initial, settings(30.0, 0.25)).unwrap();

    let (sigma, gamma) = (parameters.sigma_per_day, parameters.gamma_per_day);
    // RK4's global error with sigma * h = 0.05 stays near 1e-6 individuals
    // here; 1e-5 leaves roughly a tenfold margin.
    let tolerance = 1e-5;
    for observation in &result.trajectory {
        let time = observation.time_days;
        let exposed = initial.exposed * (-sigma * time).exp();
        let infected = initial.infected * (-gamma * time).exp()
            + sigma * initial.exposed / (gamma - sigma)
                * ((-sigma * time).exp() - (-gamma * time).exp());
        let state = observation.state;
        assert_eq!(state.susceptible, initial.susceptible);
        assert!(
            (state.exposed - exposed).abs() < tolerance,
            "exposed at day {time}: {} vs {exposed}",
            state.exposed
        );
        assert!(
            (state.infected - infected).abs() < tolerance,
            "infected at day {time}: {} vs {infected}",
            state.infected
        );
    }
}

#[test]
fn seir_preserves_total_population_within_tolerance() {
    let result = simulate_seir(PARAMETERS, INITIAL, settings(160.0, 0.25)).unwrap();

    for observation in &result.trajectory {
        let state = observation.state;
        let total = state.susceptible + state.exposed + state.infected + state.recovered;
        assert!(
            ((total - 1000.0) / 1000.0).abs() <= 1e-9,
            "total {total} at day {}",
            observation.time_days
        );
    }
    assert!(result.diagnostics.maximum_relative_population_drift <= 1e-9);
}

#[test]
fn seir_returns_error_for_invalid_parameters_and_state() {
    let run = |parameters, initial| simulate_seir(parameters, initial, settings(10.0, 1.0));

    let negative_beta = SeirParameters {
        beta_per_day: -0.5,
        ..PARAMETERS
    };
    assert_eq!(
        run(negative_beta, INITIAL).unwrap_err(),
        SimulationError::InvalidParameter {
            name: "beta_per_day"
        }
    );

    let zero_sigma = SeirParameters {
        sigma_per_day: 0.0,
        ..PARAMETERS
    };
    assert_eq!(
        run(zero_sigma, INITIAL).unwrap_err(),
        SimulationError::InvalidParameter {
            name: "sigma_per_day"
        }
    );

    let zero_gamma = SeirParameters {
        gamma_per_day: 0.0,
        ..PARAMETERS
    };
    assert_eq!(
        run(zero_gamma, INITIAL).unwrap_err(),
        SimulationError::InvalidParameter {
            name: "gamma_per_day"
        }
    );

    let negative_exposed = SeirState {
        exposed: -10.0,
        ..INITIAL
    };
    assert_eq!(
        run(PARAMETERS, negative_exposed).unwrap_err(),
        SimulationError::InvalidInitialState {
            compartment: "exposed"
        }
    );

    let empty = SeirState {
        susceptible: 0.0,
        exposed: 0.0,
        infected: 0.0,
        recovered: 0.0,
    };
    assert_eq!(
        run(PARAMETERS, empty).unwrap_err(),
        SimulationError::InvalidPopulation
    );
}

#[test]
fn seir_final_time_equals_requested_duration() {
    // Durations that are not a whole number of steps, including two whose
    // step count is sensitive to floating-point rounding.
    for (duration, step) in [(10.0, 3.0), (7.3, 0.25), (1.0, 0.1), (609.6, 0.6)] {
        let result = simulate_seir(PARAMETERS, INITIAL, settings(duration, step)).unwrap();
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

#[test]
fn seir_first_step_matches_model_equations() {
    // Over one very short step each compartment changes at the rate the
    // equations give for the initial state: with S = 980, E = 10, I = 10,
    // N = 1000 the infection flow is 0.5 * 980 * 10 / 1000 = 4.9, the
    // progression flow is 0.2 * 10 = 2, and the recovery flow is 0.1 * 10 = 1.
    let step = 1e-3;
    let result = simulate_seir(PARAMETERS, INITIAL, settings(step, step)).unwrap();
    let next = result.trajectory[1].state;

    let rates = [
        ("susceptible", next.susceptible - INITIAL.susceptible, -4.9),
        ("exposed", next.exposed - INITIAL.exposed, 4.9 - 2.0),
        ("infected", next.infected - INITIAL.infected, 2.0 - 1.0),
        ("recovered", next.recovered - INITIAL.recovered, 1.0),
    ];
    for (compartment, change, expected) in rates {
        // The rates themselves drift by well under 0.1% during the step.
        assert!(
            (change / step - expected).abs() < 2e-3,
            "{compartment}: {} per day, expected {expected}",
            change / step
        );
    }
}

#[test]
fn seir_rejects_step_too_large_for_stable_integration() {
    // RK4 turns decay into growth once a rate times the step exceeds about
    // 2.785, without any negative value or drift to give it away.
    let run = |beta_per_day, sigma_per_day, gamma_per_day, step| {
        let parameters = SeirParameters {
            beta_per_day,
            sigma_per_day,
            gamma_per_day,
        };
        simulate_seir(parameters, INITIAL, settings(10.0, step))
    };
    let too_large = |parameter| SimulationError::TimeStepTooLarge { parameter };

    assert_eq!(
        run(2.786, 0.2, 0.1, 1.0).unwrap_err(),
        too_large("beta_per_day")
    );
    assert_eq!(
        run(0.0, 2.786, 0.1, 1.0).unwrap_err(),
        too_large("sigma_per_day")
    );
    assert_eq!(
        run(0.0, 0.2, 2.786, 1.0).unwrap_err(),
        too_large("gamma_per_day")
    );

    // Exactly at the limit this guard does not fire. The run is still too
    // coarse to succeed: the per-step checks stop it when a compartment goes
    // negative, so it fails with an error and not with a wrong trajectory.
    assert_eq!(
        run(0.0, 2.785, 0.1, 1.0).unwrap_err(),
        SimulationError::NegativeCompartment {
            compartment: "infected"
        }
    );
    // The same rate is fine with a smaller step.
    assert!(run(0.0, 2.786, 0.1, 0.5).is_ok());
}

#[test]
fn seir_reports_metric_overflow_before_integrating() {
    // beta / gamma overflows to infinity, and these settings would also need
    // ten billion steps. The parameter problem is reported first, without
    // running or sizing the simulation.
    let parameters = SeirParameters {
        beta_per_day: 10.0,
        sigma_per_day: 0.2,
        gamma_per_day: f64::MIN_POSITIVE,
    };
    assert_eq!(
        simulate_seir(parameters, INITIAL, settings(1e9, 0.1)).unwrap_err(),
        SimulationError::MetricOverflow
    );
}

#[test]
fn seir_surfaces_step_limit_clock_and_population_errors() {
    let limited = settings(10.0, 1.0).with_max_steps(9).unwrap();
    assert_eq!(
        simulate_seir(PARAMETERS, INITIAL, limited).unwrap_err(),
        SimulationError::StepLimitExceeded
    );

    // Near 1e17 adjacent f64 values are 16 apart, so a step of 1 cannot
    // advance the clock. Only a 64-bit step limit can be raised this far.
    #[cfg(target_pointer_width = "64")]
    {
        let stalled = settings(1e17, 1.0).with_max_steps(usize::MAX).unwrap();
        assert_eq!(
            simulate_seir(PARAMETERS, INITIAL, stalled).unwrap_err(),
            SimulationError::TimeStepTooSmall
        );
    }

    // Each compartment is finite; only their sum overflows.
    let overflowed = SeirState {
        susceptible: f64::MAX,
        exposed: 0.0,
        infected: f64::MAX,
        recovered: 0.0,
    };
    assert_eq!(
        simulate_seir(PARAMETERS, overflowed, settings(10.0, 1.0)).unwrap_err(),
        SimulationError::InvalidPopulation
    );
}
