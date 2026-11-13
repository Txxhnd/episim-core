use episim_core::{
    SimulationError, SimulationSettings, SisParameters, SisState, estimate_sis_step_error,
    simulate_sis,
};

const PARAMETERS: SisParameters = SisParameters {
    beta_per_day: 0.3,
    gamma_per_day: 0.1,
};

const INITIAL: SisState = SisState {
    susceptible: 990.0,
    infected: 10.0,
};

fn settings(duration_days: f64, step_days: f64) -> SimulationSettings {
    SimulationSettings::new(duration_days, step_days).unwrap()
}

/// The exact SIS solution is logistic. With `r = beta - gamma`,
/// `K = N (1 - gamma / beta)` and `A = K / I(0) - 1`, it returns `I(t)` and
/// the exact cumulative infections `C(t) = I(t) - I(0) + gamma * integral of I`.
/// Both formulas also hold for `R0 < 1`, where `K < 0`.
fn logistic(parameters: SisParameters, initial: SisState, time: f64) -> (f64, f64) {
    let SisParameters {
        beta_per_day: beta,
        gamma_per_day: gamma,
    } = parameters;
    let population = initial.susceptible + initial.infected;
    let rate = beta - gamma;
    let capacity = population * (1.0 - gamma / beta);
    let a = capacity / initial.infected - 1.0;
    let infected = capacity / (1.0 + a * (-rate * time).exp());
    let integral =
        capacity * time + capacity / rate * ((1.0 + a * (-rate * time).exp()) / (1.0 + a)).ln();
    (infected, infected - initial.infected + gamma * integral)
}

#[test]
fn sis_matches_the_logistic_closed_form() {
    // (beta, gamma, initially infected): an endemic case, one barely above
    // the threshold, and one below it.
    for (beta_per_day, gamma_per_day, infected) in
        [(0.3, 0.1, 10.0), (1.5, 1.4, 10.0), (0.1, 0.3, 500.0)]
    {
        let parameters = SisParameters {
            beta_per_day,
            gamma_per_day,
        };
        let initial = SisState {
            susceptible: 1000.0 - infected,
            infected,
        };
        let result = simulate_sis(parameters, initial, settings(60.0, 0.25)).unwrap();
        for observation in &result.trajectory {
            let (infected, cumulative) = logistic(parameters, initial, observation.time_days);
            // The advisor measured at most 2.1e-5 people for I and 4.1e-5 for C.
            assert!(
                (observation.state.infected - infected).abs() < 2e-4,
                "beta {beta_per_day}, day {}: I {} vs {infected}",
                observation.time_days,
                observation.state.infected
            );
            assert!(
                (observation.cumulative_infections - cumulative).abs() < 2e-4,
                "beta {beta_per_day}, day {}: C {} vs {cumulative}",
                observation.time_days,
                observation.cumulative_infections
            );
        }
    }
}

#[test]
fn sis_settles_at_the_endemic_equilibrium() {
    // (beta, gamma, expected I* / N): R0 = 3 gives 1 - 1/3; R0 = 1/3 dies out.
    for (beta_per_day, gamma_per_day, prevalence) in [(0.3, 0.1, 2.0 / 3.0), (0.1, 0.3, 0.0)] {
        let parameters = SisParameters {
            beta_per_day,
            gamma_per_day,
        };
        let result = simulate_sis(parameters, INITIAL, settings(400.0, 0.25)).unwrap();
        let metrics = result.metrics;
        assert!(
            (metrics.endemic_prevalence.unwrap() - prevalence).abs() < 1e-12,
            "{metrics:?}"
        );
        let last = result.trajectory.last().unwrap().state;
        assert!(
            (last.infected / 1000.0 - prevalence).abs() < 1e-6,
            "beta {beta_per_day}: {last:?}"
        );
        // Nobody is ever immune, and nothing is final.
        assert_eq!(metrics.herd_immunity_threshold, None);
        assert_eq!(metrics.final_size_fraction, None);
    }
}

#[test]
fn sis_conserves_the_population() {
    let result = simulate_sis(PARAMETERS, INITIAL, settings(200.0, 0.25)).unwrap();
    for observation in &result.trajectory {
        let total = observation.state.susceptible + observation.state.infected;
        assert!(
            (total - 1000.0).abs() <= 1e-6,
            "day {}",
            observation.time_days
        );
    }
    assert!(result.diagnostics.maximum_relative_population_drift <= 1e-9);
}

#[test]
fn sis_counts_reinfections() {
    let result = simulate_sis(PARAMETERS, INITIAL, settings(200.0, 0.25)).unwrap();
    assert_eq!(result.trajectory[0].cumulative_infections, 0.0);
    for pair in result.trajectory.windows(2) {
        assert!(pair[1].cumulative_infections >= pair[0].cumulative_infections);
    }
    // Close to the closed form's 12.59 infections per head: far more than
    // one, because people recover and are infected again.
    let (_, exact) = logistic(PARAMETERS, INITIAL, 200.0);
    let metrics = result.metrics;
    assert!((metrics.cumulative_infections_since_start - exact).abs() < 1e-3);
    assert!(
        (metrics.cumulative_incidence_fraction - 12.59).abs() < 0.01,
        "{metrics:?}"
    );
}

#[test]
fn sis_peak_time_is_when_the_infected_level_off() {
    // I rises to I* and then stays there, moved only by roundoff. The peak
    // time is the first observation within roundoff of the peak, not
    // wherever roundoff happens to put the largest value.
    let result = simulate_sis(PARAMETERS, INITIAL, settings(200.0, 0.25)).unwrap();
    let metrics = result.metrics;
    let tolerance = 64.0 * f64::EPSILON * 1000.0;
    let first_at_peak = result
        .trajectory
        .iter()
        .find(|o| o.state.infected >= metrics.peak_infectious - tolerance)
        .unwrap();
    assert_eq!(metrics.peak_infectious_time_days, first_at_peak.time_days);
    // Near the plateau the logistic curve is I* - I(t) = I* * A * exp(-r t),
    // with A = I* / I(0) - 1 and r = beta - gamma, so it comes within the
    // tolerance at t = ln(I* * A / tolerance) / r, about day 178.7. Without
    // the tolerance the peak time is day 191, wherever roundoff put the
    // largest value.
    let equilibrium = 1000.0 * (1.0 - 1.0 / 3.0);
    let a = equilibrium / INITIAL.infected - 1.0;
    let reached = (equilibrium * a / tolerance).ln() / 0.2;
    assert!(
        (metrics.peak_infectious_time_days - reached).abs() < 2.0,
        "{metrics:?}, plateau reached at day {reached}"
    );
}

#[test]
fn sis_effective_reproduction_number_sets_the_direction_of_the_infected() {
    // Above the equilibrium the infected decline, below it they grow, and in
    // both cases effective R tends to 1.
    let above = SisState {
        susceptible: 100.0,
        infected: 900.0,
    };
    for initial in [INITIAL, above] {
        let result = simulate_sis(PARAMETERS, initial, settings(200.0, 0.25)).unwrap();
        for pair in result.trajectory.windows(2) {
            let r = PARAMETERS
                .effective_reproduction_number(&pair[0].state)
                .unwrap();
            let change = pair[1].state.infected - pair[0].state.infected;
            if (r - 1.0).abs() > 1e-6 {
                assert_eq!(change > 0.0, r > 1.0, "R {r}, change {change}");
            }
        }
        let last = result.trajectory.last().unwrap().state;
        let r = PARAMETERS.effective_reproduction_number(&last).unwrap();
        assert!((r - 1.0).abs() < 1e-6, "{r}");
    }
}

#[test]
fn sis_rejects_invalid_parameters_and_state() {
    let run = |parameters, initial| simulate_sis(parameters, initial, settings(10.0, 0.5));
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        let beta = SisParameters {
            beta_per_day: invalid,
            ..PARAMETERS
        };
        assert_eq!(
            run(beta, INITIAL).unwrap_err(),
            SimulationError::InvalidParameter {
                name: "beta_per_day"
            }
        );
        let state = SisState {
            infected: invalid,
            ..INITIAL
        };
        assert_eq!(
            run(PARAMETERS, state).unwrap_err(),
            SimulationError::InvalidInitialState {
                compartment: "infected"
            }
        );
    }
    let no_recovery = SisParameters {
        gamma_per_day: 0.0,
        ..PARAMETERS
    };
    assert_eq!(
        run(no_recovery, INITIAL).unwrap_err(),
        SimulationError::InvalidParameter {
            name: "gamma_per_day"
        }
    );
    let empty = SisState {
        susceptible: 0.0,
        infected: 0.0,
    };
    assert_eq!(
        run(PARAMETERS, empty).unwrap_err(),
        SimulationError::InvalidPopulation
    );
}

#[test]
fn sis_rejects_step_too_large_for_stable_integration() {
    let run = |beta_per_day, gamma_per_day, step| {
        let parameters = SisParameters {
            beta_per_day,
            gamma_per_day,
        };
        simulate_sis(parameters, INITIAL, settings(10.0, step))
    };
    // A fastest rate of 4 per day allows a step of at most 2 / 4 = 0.5 days:
    // beta + gamma = 3 + 1, and with no transmission 1.25 * gamma = 1.25 * 3.2.
    let too_large = SimulationError::TimeStepTooLarge { max_step_days: 0.5 };
    for (beta, gamma) in [(3.0, 1.0), (0.0, 3.2)] {
        assert_eq!(
            run(beta, gamma, 0.51).unwrap_err(),
            too_large,
            "beta {beta}"
        );
        assert!(run(beta, gamma, 0.5).is_ok(), "beta {beta}");
    }
}

#[test]
fn sis_step_error_estimate_matches_the_error_against_the_closed_form() {
    let parameters = SisParameters {
        beta_per_day: 0.5,
        gamma_per_day: 0.1,
    };
    let estimate =
        |step| estimate_sis_step_error(parameters, INITIAL, settings(40.0, step)).unwrap();
    let run = simulate_sis(parameters, INITIAL, settings(40.0, 1.0)).unwrap();
    let actual = run
        .trajectory
        .iter()
        .flat_map(|o| {
            let (infected, cumulative) = logistic(parameters, INITIAL, o.time_days);
            [
                (o.state.infected - infected).abs(),
                (o.cumulative_infections - cumulative).abs(),
            ]
        })
        .fold(0.0, f64::max);
    let coarse = estimate(1.0);
    let ratio = coarse.max_absolute_error / actual;
    assert!(
        (0.9..1.1).contains(&ratio),
        "estimated {coarse:?}, actual {actual}"
    );
    // Halving the step divides the error by about 2^4 = 16.
    let halving = coarse.max_absolute_error / estimate(0.5).max_absolute_error;
    assert!((12.0..20.0).contains(&halving), "{halving}");
}

#[test]
fn sis_without_initial_infection_never_becomes_endemic() {
    // R0 = 3, but with nobody infected the infection can never start, so
    // the run stays where it began and the endemic prevalence is 0.
    let initial = SisState {
        susceptible: 1000.0,
        infected: 0.0,
    };
    let result = simulate_sis(PARAMETERS, initial, settings(400.0, 0.25)).unwrap();
    for observation in &result.trajectory {
        assert_eq!(observation.state, initial);
        assert_eq!(observation.cumulative_infections, 0.0);
    }
    assert_eq!(result.metrics.endemic_prevalence, Some(0.0));
}

#[test]
fn sis_stability_check_uses_the_step_actually_taken() {
    // A duration shorter than the step is integrated in one step of that
    // duration, so the nominal step of 7 days is never used.
    let result = simulate_sis(PARAMETERS, INITIAL, settings(1.0, 7.0)).unwrap();
    assert_eq!(result.diagnostics.steps_executed, 1);

    // The single step of 1 day is still too large for 1.25 * gamma = 4 per day.
    let fast = SisParameters {
        beta_per_day: 0.0,
        gamma_per_day: 3.2,
    };
    assert_eq!(
        simulate_sis(fast, INITIAL, settings(1.0, 7.0)).unwrap_err(),
        SimulationError::TimeStepTooLarge { max_step_days: 0.5 }
    );
}

#[test]
fn sis_reports_metric_overflow_before_integrating() {
    // beta / gamma overflows to infinity, and these settings would also need
    // ten billion steps. The parameter problem is reported first.
    let parameters = SisParameters {
        beta_per_day: 10.0,
        gamma_per_day: f64::MIN_POSITIVE,
    };
    assert_eq!(
        simulate_sis(parameters, INITIAL, settings(1e9, 0.1)).unwrap_err(),
        SimulationError::MetricOverflow
    );
}

#[test]
fn sis_effective_reproduction_number_rejects_invalid_input() {
    let no_recovery = SisParameters {
        gamma_per_day: 0.0,
        ..PARAMETERS
    };
    assert_eq!(
        no_recovery
            .effective_reproduction_number(&INITIAL)
            .unwrap_err(),
        SimulationError::InvalidParameter {
            name: "gamma_per_day"
        }
    );
    let empty = SisState {
        susceptible: 0.0,
        infected: 0.0,
    };
    assert_eq!(
        PARAMETERS
            .effective_reproduction_number(&empty)
            .unwrap_err(),
        SimulationError::InvalidPopulation
    );
    let overflowing = SisParameters {
        beta_per_day: 10.0,
        gamma_per_day: f64::MIN_POSITIVE,
    };
    assert_eq!(
        overflowing
            .effective_reproduction_number(&INITIAL)
            .unwrap_err(),
        SimulationError::MetricOverflow
    );
}
