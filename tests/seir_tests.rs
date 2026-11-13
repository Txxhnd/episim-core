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
            cumulative_infections: 0.0,
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
fn seir_cumulative_infections_equal_the_drop_in_susceptibles() {
    // Everyone who leaves S is infected, so C(t) = S(0) - S(t) exactly; the
    // two are integrated separately and may differ only by roundoff.
    for (beta_per_day, sigma_per_day, gamma_per_day) in [(0.5, 0.2, 0.1), (1.2, 0.5, 0.25)] {
        let parameters = SeirParameters {
            beta_per_day,
            sigma_per_day,
            gamma_per_day,
        };
        let result = simulate_seir(parameters, INITIAL, settings(200.0, 0.25)).unwrap();
        assert_eq!(result.trajectory[0].cumulative_infections, 0.0);
        let mut largest_gap = 0.0_f64;
        for pair in result.trajectory.windows(2) {
            assert!(pair[1].cumulative_infections >= pair[0].cumulative_infections);
        }
        for observation in &result.trajectory {
            let drop = INITIAL.susceptible - observation.state.susceptible;
            largest_gap = largest_gap.max((observation.cumulative_infections - drop).abs());
        }
        assert!(largest_gap < 1e-9, "beta {beta_per_day}: {largest_gap}");
        let last = result.trajectory.last().unwrap();
        assert_eq!(
            result.metrics.cumulative_infections_since_start,
            last.cumulative_infections
        );
    }
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
    // The fastest rate times the step may be at most 2. Beyond that RK4
    // damps too weakly, or not at all, and nothing else would notice.
    let run = |beta_per_day, sigma_per_day, gamma_per_day, step| {
        let parameters = SeirParameters {
            beta_per_day,
            sigma_per_day,
            gamma_per_day,
        };
        simulate_seir(parameters, INITIAL, settings(10.0, step))
    };
    let too_large = |max_step_days| SimulationError::TimeStepTooLarge { max_step_days };

    // Without transmission the fastest rate is the larger of sigma and gamma.
    assert_eq!(run(0.0, 4.0, 0.1, 1.0).unwrap_err(), too_large(0.5));
    assert_eq!(run(0.0, 0.2, 4.0, 1.0).unwrap_err(), too_large(0.5));
    // Susceptibles can empty at up to beta per day.
    assert_eq!(run(4.0, 0.2, 0.1, 1.0).unwrap_err(), too_large(0.5));
    // The step named in the error is accepted.
    assert!(run(0.0, 4.0, 0.1, 0.5).is_ok());
    assert!(run(0.0, 0.2, 4.0, 0.5).is_ok());
}

#[test]
fn seir_rejects_step_too_large_for_the_coupled_exposed_and_infected_rate() {
    // Exposed and infected feed each other, so together they decay faster
    // than sigma or gamma alone. With sigma = gamma = g the fastest rate is
    // g + sqrt(g * beta).
    let run = |rate, step| {
        let parameters = SeirParameters {
            beta_per_day: rate,
            sigma_per_day: rate,
            gamma_per_day: rate,
        };
        simulate_seir(parameters, INITIAL, settings(10.0, step))
    };

    // Every rate times the step is 1.5, but the coupled rate is 3 per day.
    assert_eq!(
        run(1.5, 1.0).unwrap_err(),
        SimulationError::TimeStepTooLarge {
            max_step_days: 2.0 / 3.0
        }
    );
    assert!(run(1.5, 0.5).is_ok());
    // With all rates 1 the coupled rate is exactly 2: at the limit, not over it.
    assert!(!matches!(
        run(1.0, 1.0),
        Err(SimulationError::TimeStepTooLarge { .. })
    ));
}

#[test]
fn seir_rejects_subcritical_outbreak_that_a_coarse_step_would_inflate() {
    // R0 = beta / gamma = 0.36, so the infected can only decline. With a
    // step of one day every single rate looks acceptable, yet the coupled
    // rate is about 2.88 per day; unguarded, RK4 reports a peak of 370
    // infected and 80% of the population infected.
    let parameters = SeirParameters {
        beta_per_day: 1.0,
        sigma_per_day: 0.25,
        gamma_per_day: 2.78,
    };
    let initial = SeirState {
        susceptible: 800.0,
        exposed: 100.0,
        infected: 100.0,
        recovered: 0.0,
    };
    assert!(matches!(
        simulate_seir(parameters, initial, settings(100.0, 1.0)),
        Err(SimulationError::TimeStepTooLarge { .. })
    ));

    // With a step the guard accepts, the infected never exceed their start.
    let result = simulate_seir(parameters, initial, settings(100.0, 0.5)).unwrap();
    assert_eq!(result.metrics.peak_infectious, 100.0);
    assert_eq!(result.metrics.peak_infectious_time_days, 0.0);
}

#[test]
fn seir_stability_check_uses_the_step_actually_taken() {
    // A duration shorter than the step is integrated in one step of that
    // duration, so the nominal step of 7 days is never used.
    let result = simulate_seir(PARAMETERS, INITIAL, settings(1.0, 7.0)).unwrap();
    assert_eq!(result.diagnostics.steps_executed, 1);
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

#[test]
fn seir_records_every_nth_step_without_changing_the_integration() {
    let every_step = settings(160.0, 0.25);
    let daily = every_step.with_record_every(4).unwrap();
    let full = simulate_seir(PARAMETERS, INITIAL, every_step).unwrap();
    let thinned = simulate_seir(PARAMETERS, INITIAL, daily).unwrap();

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
    run: &[Observation<SeirState>],
    reference: &[Observation<SeirState>],
    stride: usize,
) -> f64 {
    run.iter()
        .enumerate()
        .flat_map(|(index, observation)| {
            let (state, exact) = (observation.state, reference[index * stride].state);
            [
                (state.susceptible - exact.susceptible).abs(),
                (state.exposed - exact.exposed).abs(),
                (state.infected - exact.infected).abs(),
                (state.recovered - exact.recovered).abs(),
            ]
        })
        .fold(0.0, f64::max)
}

#[test]
fn seir_step_error_estimate_matches_the_error_against_a_fine_reference() {
    let parameters = PARAMETERS;
    let run = |step| simulate_seir(parameters, INITIAL, settings(40.0, step)).unwrap();
    let estimate =
        episim_core::estimate_seir_step_error(parameters, INITIAL, settings(40.0, 1.0)).unwrap();

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
fn seir_step_error_estimate_falls_about_sixteenfold_when_the_step_is_halved() {
    let parameters = PARAMETERS;
    let estimate = |step| {
        episim_core::estimate_seir_step_error(parameters, INITIAL, settings(40.0, step))
            .unwrap()
            .max_absolute_error
    };
    for step in [1.0, 0.5] {
        let ratio = estimate(step) / estimate(step / 2.0);
        assert!((12.0..18.0).contains(&ratio), "step {step}: ratio {ratio}");
    }
}

#[test]
fn seir_step_error_estimate_follows_the_recording_interval_and_reports_errors() {
    let every_step = settings(40.0, 0.5);
    let sparse = every_step.with_record_every(8).unwrap();
    let all = episim_core::estimate_seir_step_error(PARAMETERS, INITIAL, every_step).unwrap();
    let some = episim_core::estimate_seir_step_error(PARAMETERS, INITIAL, sparse).unwrap();
    // The largest error over fewer recorded times cannot be larger.
    assert!(some.max_absolute_error > 0.0);
    assert!(some.max_absolute_error <= all.max_absolute_error);

    // A run the simulation rejects is rejected here with the same error.
    let limited = settings(10.0, 1.0).with_max_steps(9).unwrap();
    assert_eq!(
        episim_core::estimate_seir_step_error(PARAMETERS, INITIAL, limited).unwrap_err(),
        SimulationError::StepLimitExceeded
    );
}

#[test]
fn seir_final_size_matches_a_long_simulation() {
    let result = simulate_seir(PARAMETERS, INITIAL, settings(400.0, 0.1)).unwrap();
    let metrics = result.metrics;
    let final_size = metrics.final_size_fraction.unwrap();
    // By day 400 the epidemic is over, so the simulated incidence has reached
    // the limit the final size equation predicts.
    assert!(
        (metrics.cumulative_incidence_fraction - final_size).abs() < 1e-6,
        "simulated {} vs final size {final_size}",
        metrics.cumulative_incidence_fraction,
    );
    // R0 = 5, so four fifths of the population must be immune.
    assert!((metrics.herd_immunity_threshold.unwrap() - 0.8).abs() < 1e-12);
}

#[test]
fn seir_effective_reproduction_number_falls_as_susceptibles_are_used_up() {
    let result = simulate_seir(PARAMETERS, INITIAL, settings(160.0, 0.25)).unwrap();
    let values: Vec<f64> = result
        .trajectory
        .iter()
        .map(|o| PARAMETERS.effective_reproduction_number(&o.state).unwrap())
        .collect();

    // R0 = 5 and 98% of the population starts susceptible.
    assert!((values[0] - 4.9).abs() < 1e-9);
    // Susceptibles only ever decrease; the slack covers roundoff in N.
    assert!(values.windows(2).all(|pair| pair[1] <= pair[0] + 1e-12));
    assert!(*values.last().unwrap() < 1.0);
}

#[test]
fn seir_effective_reproduction_number_rejects_invalid_input() {
    let no_progression = SeirParameters {
        sigma_per_day: 0.0,
        ..PARAMETERS
    };
    assert_eq!(
        no_progression
            .effective_reproduction_number(&INITIAL)
            .unwrap_err(),
        SimulationError::InvalidParameter {
            name: "sigma_per_day"
        }
    );

    let empty = SeirState {
        susceptible: 0.0,
        exposed: 0.0,
        infected: 0.0,
        recovered: 0.0,
    };
    assert_eq!(
        PARAMETERS
            .effective_reproduction_number(&empty)
            .unwrap_err(),
        SimulationError::InvalidPopulation
    );
}
