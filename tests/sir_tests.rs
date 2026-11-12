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

#[test]
fn sir_first_step_matches_model_equations() {
    // Over one very short step each compartment changes at the rate the
    // equations give for the initial state: with S = 990, I = 10, N = 1000
    // the infection flow is 0.3 * 990 * 10 / 1000 = 2.97 and the recovery
    // flow is 0.1 * 10 = 1.
    let step = 1e-3;
    let result = simulate_sir(PARAMETERS, INITIAL, settings(step, step)).unwrap();
    let next = result.trajectory[1].state;

    let rates = [
        ("susceptible", next.susceptible - INITIAL.susceptible, -2.97),
        ("infected", next.infected - INITIAL.infected, 2.97 - 1.0),
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
fn sir_rejects_step_too_large_for_stable_integration() {
    // The fastest rate times the step may be at most 2. Beyond that RK4
    // damps too weakly, or not at all, and nothing else would notice.
    let run = |beta_per_day, gamma_per_day, duration, step| {
        let parameters = SirParameters {
            beta_per_day,
            gamma_per_day,
        };
        simulate_sir(parameters, INITIAL, settings(duration, step))
    };
    // A rate of 4 per day allows a step of at most 2 / 4 = 0.5 days.
    let too_large = SimulationError::TimeStepTooLarge { max_step_days: 0.5 };

    assert_eq!(run(0.0, 4.0, 10.0, 1.0).unwrap_err(), too_large);
    assert_eq!(run(4.0, 0.1, 10.0, 1.0).unwrap_err(), too_large);
    // The step named in the error is accepted, for either rate.
    assert!(run(0.0, 4.0, 10.0, 0.5).is_ok());
    assert!(run(4.0, 0.1, 10.0, 0.5).is_ok());

    // Exactly at the limit the run is accepted, and pure recovery still
    // removes two thirds of the infected per step (RK4's factor at 2 is 1/3).
    let at_limit = run(0.0, 2.0, 10.0, 1.0).unwrap();
    let after_one_step = at_limit.trajectory[1].state.infected;
    assert!((after_one_step - INITIAL.infected / 3.0).abs() < 1e-9);
    assert!(matches!(
        run(0.0, 2.01, 10.0, 1.0),
        Err(SimulationError::TimeStepTooLarge { .. })
    ));
}

#[test]
fn sir_stability_check_uses_the_step_actually_taken() {
    // A duration shorter than the step is integrated in one step of that
    // duration, so the nominal step of 7 days is never used.
    let parameters = SirParameters {
        beta_per_day: 0.3,
        gamma_per_day: 0.5,
    };
    let result = simulate_sir(parameters, INITIAL, settings(1.0, 7.0)).unwrap();
    assert_eq!(result.diagnostics.steps_executed, 1);

    // The single step of 1 day is still too large for a rate of 4 per day.
    let fast = SirParameters {
        beta_per_day: 0.0,
        gamma_per_day: 4.0,
    };
    assert_eq!(
        simulate_sir(fast, INITIAL, settings(1.0, 7.0)).unwrap_err(),
        SimulationError::TimeStepTooLarge { max_step_days: 0.5 }
    );
}

#[test]
fn sir_reports_metric_overflow_before_integrating() {
    // beta / gamma overflows to infinity, and these settings would also need
    // ten billion steps. The parameter problem is reported first, without
    // running or sizing the simulation.
    let parameters = SirParameters {
        beta_per_day: 10.0,
        gamma_per_day: f64::MIN_POSITIVE,
    };
    assert_eq!(
        simulate_sir(parameters, INITIAL, settings(1e9, 0.1)).unwrap_err(),
        SimulationError::MetricOverflow
    );
}

#[test]
fn sir_surfaces_step_limit_clock_and_population_errors() {
    let limited = settings(10.0, 1.0).with_max_steps(9).unwrap();
    assert_eq!(
        simulate_sir(PARAMETERS, INITIAL, limited).unwrap_err(),
        SimulationError::StepLimitExceeded
    );

    // Near 1e17 adjacent f64 values are 16 apart, so a step of 1 cannot
    // advance the clock. Only a 64-bit step limit can be raised this far.
    #[cfg(target_pointer_width = "64")]
    {
        let stalled = settings(1e17, 1.0).with_max_steps(usize::MAX).unwrap();
        assert_eq!(
            simulate_sir(PARAMETERS, INITIAL, stalled).unwrap_err(),
            SimulationError::TimeStepTooSmall
        );
    }

    // Each compartment is finite; only their sum overflows.
    let overflowed = SirState {
        susceptible: f64::MAX,
        infected: f64::MAX,
        recovered: 0.0,
    };
    assert_eq!(
        simulate_sir(PARAMETERS, overflowed, settings(10.0, 1.0)).unwrap_err(),
        SimulationError::InvalidPopulation
    );
}

#[test]
fn sir_records_every_nth_step_without_changing_the_integration() {
    let every_step = settings(160.0, 0.25);
    let daily = every_step.with_record_every(4).unwrap();
    let full = simulate_sir(PARAMETERS, INITIAL, every_step).unwrap();
    let thinned = simulate_sir(PARAMETERS, INITIAL, daily).unwrap();

    // 640 steps of a quarter day, recorded once per day.
    assert_eq!(thinned.trajectory.len(), 161);
    assert_eq!(thinned.diagnostics.steps_executed, 640);
    for (day, observation) in thinned.trajectory.iter().enumerate() {
        assert_eq!(observation.time_days, day as f64);
        // Each recorded observation is exactly the one the full run has.
        assert_eq!(*observation, full.trajectory[day * 4]);
    }
    assert_eq!(thinned.diagnostics, full.diagnostics);
    assert_eq!(
        thinned.metrics.cumulative_infections_since_start,
        full.metrics.cumulative_infections_since_start
    );
    // The peak is read from fewer observations, so it can only be lower or equal.
    assert!(thinned.metrics.peak_infectious <= full.metrics.peak_infectious);
}

/// Largest compartment difference between a run and a reference run whose
/// step is `stride` times smaller, at the times of the first run.
fn error_against_reference(
    run: &[Observation<SirState>],
    reference: &[Observation<SirState>],
    stride: usize,
) -> f64 {
    run.iter()
        .enumerate()
        .flat_map(|(index, observation)| {
            let (state, exact) = (observation.state, reference[index * stride].state);
            [
                (state.susceptible - exact.susceptible).abs(),
                (state.infected - exact.infected).abs(),
                (state.recovered - exact.recovered).abs(),
            ]
        })
        .fold(0.0, f64::max)
}

#[test]
fn sir_step_error_estimate_matches_the_error_against_a_fine_reference() {
    let parameters = SirParameters {
        beta_per_day: 0.5,
        gamma_per_day: 0.1,
    };
    let run = |step| simulate_sir(parameters, INITIAL, settings(40.0, step)).unwrap();
    let estimate =
        episim_core::estimate_sir_step_error(parameters, INITIAL, settings(40.0, 1.0)).unwrap();

    // A step 256 times smaller stands in for the exact solution.
    let actual = error_against_reference(&run(1.0).trajectory, &run(1.0 / 256.0).trajectory, 256);
    let ratio = estimate.max_absolute_error / actual;
    assert!(
        (0.9..1.1).contains(&ratio),
        "estimated {}, actual {actual}",
        estimate.max_absolute_error
    );
    // INITIAL holds a population of 1000.
    assert_eq!(
        estimate.max_relative_error,
        estimate.max_absolute_error / 1000.0
    );
}

#[test]
fn sir_step_error_estimate_falls_about_sixteenfold_when_the_step_is_halved() {
    let parameters = SirParameters {
        beta_per_day: 0.5,
        gamma_per_day: 0.1,
    };
    let estimate = |step| {
        episim_core::estimate_sir_step_error(parameters, INITIAL, settings(40.0, step))
            .unwrap()
            .max_absolute_error
    };
    for step in [1.0, 0.5] {
        let ratio = estimate(step) / estimate(step / 2.0);
        assert!((12.0..18.0).contains(&ratio), "step {step}: ratio {ratio}");
    }
}

#[test]
fn sir_step_error_estimate_follows_the_recording_interval_and_reports_errors() {
    let every_step = settings(40.0, 0.5);
    let sparse = every_step.with_record_every(8).unwrap();
    let all = episim_core::estimate_sir_step_error(PARAMETERS, INITIAL, every_step).unwrap();
    let some = episim_core::estimate_sir_step_error(PARAMETERS, INITIAL, sparse).unwrap();
    // The largest error over fewer recorded times cannot be larger.
    assert!(some.max_absolute_error > 0.0);
    assert!(some.max_absolute_error <= all.max_absolute_error);

    // A run the simulation rejects is rejected here with the same error.
    let limited = settings(10.0, 1.0).with_max_steps(9).unwrap();
    assert_eq!(
        episim_core::estimate_sir_step_error(PARAMETERS, INITIAL, limited).unwrap_err(),
        SimulationError::StepLimitExceeded
    );
}
