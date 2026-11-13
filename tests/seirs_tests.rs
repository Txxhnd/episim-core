use episim_core::{
    Observation, SeirParameters, SeirState, SeirsParameters, SimulationError, SimulationSettings,
    estimate_seirs_step_error, simulate_seir, simulate_seirs,
};

const PARAMETERS: SeirsParameters = SeirsParameters {
    beta_per_day: 0.5,
    sigma_per_day: 0.2,
    gamma_per_day: 0.1,
    omega_per_day: 0.01,
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

fn parameters(beta: f64, sigma: f64, gamma: f64, omega: f64) -> SeirsParameters {
    SeirsParameters {
        beta_per_day: beta,
        sigma_per_day: sigma,
        gamma_per_day: gamma,
        omega_per_day: omega,
    }
}

#[test]
fn seirs_without_waning_is_seir() {
    // With omega = 0 every waning term is exactly zero, so the runs agree
    // to the last bit, cumulative infections included.
    let seirs = simulate_seirs(
        parameters(0.5, 0.2, 0.1, 0.0),
        INITIAL,
        settings(200.0, 0.25),
    )
    .unwrap();
    let seir = simulate_seir(
        SeirParameters {
            beta_per_day: 0.5,
            sigma_per_day: 0.2,
            gamma_per_day: 0.1,
        },
        INITIAL,
        settings(200.0, 0.25),
    )
    .unwrap();
    assert_eq!(seirs.trajectory, seir.trajectory);
    assert_eq!(seirs.diagnostics, seir.diagnostics);
    assert_eq!(
        seirs.metrics.herd_immunity_threshold,
        seir.metrics.herd_immunity_threshold
    );
    // Without waning nobody stays infected for good.
    assert_eq!(seirs.metrics.endemic_prevalence, Some(0.0));
    assert_eq!(seirs.metrics.final_size_fraction, None);
}

#[test]
fn seirs_settles_at_the_endemic_equilibrium() {
    // At equilibrium S* = N / R0 and I* / N = (1 - 1/R0) / (1 + gamma/sigma + gamma/omega),
    // or S* = N and I* = 0 when R0 <= 1. The last set has R0 = 1/2.
    for (beta, sigma, gamma, omega, infected) in [
        (0.5, 0.2, 0.1, 0.01, 69.565_217_391),
        (0.8, 0.5, 0.25, 0.05, 105.769_230_769),
        (1.8, 1.0, 0.3, 1.2, 537.634_408_602),
        (0.1, 0.2, 0.2, 0.05, 0.0),
    ] {
        let result = simulate_seirs(
            parameters(beta, sigma, gamma, omega),
            INITIAL,
            settings(5000.0, 0.25),
        )
        .unwrap();
        let metrics = result.metrics;
        assert!(
            (metrics.endemic_prevalence.unwrap() * 1000.0 - infected).abs() < 1e-6,
            "{metrics:?}"
        );
        let last = result.trajectory.last().unwrap().state;
        let susceptible = 1000.0 * (gamma / beta).min(1.0);
        assert!((last.infected - infected).abs() < 1e-6 * 1000.0, "{last:?}");
        assert!(
            (last.susceptible - susceptible).abs() < 1e-6 * 1000.0,
            "{last:?}"
        );
        assert_eq!(metrics.final_size_fraction, None);
    }
}

#[test]
fn seirs_conserves_the_population() {
    let result = simulate_seirs(
        parameters(1.8, 1.0, 0.3, 1.2),
        INITIAL,
        settings(200.0, 0.25),
    )
    .unwrap();
    for observation in &result.trajectory {
        let state = observation.state;
        let total = state.susceptible + state.exposed + state.infected + state.recovered;
        assert!(
            (total - 1000.0).abs() <= 1e-6,
            "day {}",
            observation.time_days
        );
    }
    assert!(result.diagnostics.maximum_relative_population_drift <= 1e-9);
}

#[test]
fn seirs_counts_reinfections() {
    let result = simulate_seirs(PARAMETERS, INITIAL, settings(1000.0, 0.25)).unwrap();
    assert_eq!(result.trajectory[0].cumulative_infections, 0.0);
    // C - (S(0) - S) grows at omega * R: the people who lost their immunity
    // and were counted again when reinfected.
    let reinfections = |o: &Observation<SeirState>| {
        o.cumulative_infections - (INITIAL.susceptible - o.state.susceptible)
    };
    for pair in result.trajectory.windows(2) {
        assert!(pair[1].cumulative_infections >= pair[0].cumulative_infections);
        assert!(reinfections(&pair[1]) >= reinfections(&pair[0]) - 1e-9);
    }
    assert!(
        result.metrics.cumulative_incidence_fraction > 4.0,
        "{:?}",
        result.metrics
    );
    // R0 = 5; the threshold keeps its meaning at any given time.
    assert_eq!(result.metrics.herd_immunity_threshold, Some(0.8));
}

#[test]
fn seirs_effective_reproduction_number_sets_the_direction_of_the_exposed_and_infected() {
    // E + I rise and fall in damped waves around the equilibrium, rising
    // while effective R is above 1 and falling while it is below. A short
    // step keeps effective R nearly constant over each step.
    let result = simulate_seirs(PARAMETERS, INITIAL, settings(400.0, 0.05)).unwrap();
    let mut crossings = 0;
    for pair in result.trajectory.windows(2) {
        let r = PARAMETERS
            .effective_reproduction_number(&pair[0].state)
            .unwrap();
        let next = PARAMETERS
            .effective_reproduction_number(&pair[1].state)
            .unwrap();
        crossings += usize::from((r - 1.0) * (next - 1.0) < 0.0);
        if (r - 1.0).abs() > 1e-2 {
            let (before, after) = (pair[0].state, pair[1].state);
            let change = after.exposed + after.infected - before.exposed - before.infected;
            assert_eq!(change > 0.0, r > 1.0, "day {}: R {r}", pair[0].time_days);
        }
    }
    // Effective R falls below 1 and rises above it again: waning immunity
    // refills the susceptibles.
    assert!(crossings >= 3, "{crossings}");
}

#[test]
fn seirs_rejects_invalid_parameters() {
    let run = |parameters| simulate_seirs(parameters, INITIAL, settings(10.0, 0.5));
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            run(parameters(0.5, 0.2, 0.1, invalid)).unwrap_err(),
            SimulationError::InvalidParameter {
                name: "omega_per_day"
            }
        );
    }
    assert_eq!(
        run(parameters(0.5, 0.0, 0.1, 0.01)).unwrap_err(),
        SimulationError::InvalidParameter {
            name: "sigma_per_day"
        }
    );
    assert_eq!(
        run(parameters(0.5, 0.2, 0.0, 0.01)).unwrap_err(),
        SimulationError::InvalidParameter {
            name: "gamma_per_day"
        }
    );
    let negative = SeirState {
        exposed: -1.0,
        ..INITIAL
    };
    assert_eq!(
        simulate_seirs(PARAMETERS, negative, settings(10.0, 0.5)).unwrap_err(),
        SimulationError::InvalidInitialState {
            compartment: "exposed"
        }
    );
}

#[test]
fn seirs_rejects_step_too_large_for_stable_integration() {
    // The fastest rate is the sum of the four, here 4 per day, which allows
    // a step of at most 2 / 4 = 0.5 days.
    let run = |step| {
        simulate_seirs(
            parameters(1.0, 1.5, 1.0, 0.5),
            INITIAL,
            settings(10.0, step),
        )
    };
    assert_eq!(
        run(0.51).unwrap_err(),
        SimulationError::TimeStepTooLarge { max_step_days: 0.5 }
    );
    assert!(run(0.5).is_ok());
}

#[test]
fn seirs_step_error_estimate_matches_the_error_against_a_fine_reference() {
    let parameters = parameters(0.5, 0.2, 0.1, 0.05);
    let run = |step| simulate_seirs(parameters, INITIAL, settings(40.0, step)).unwrap();
    let estimate =
        |step| estimate_seirs_step_error(parameters, INITIAL, settings(40.0, step)).unwrap();
    // A step 256 times smaller stands in for the exact solution.
    let (coarse, reference) = (run(1.0), run(1.0 / 256.0));
    let actual = coarse
        .trajectory
        .iter()
        .enumerate()
        .flat_map(|(index, observation)| {
            let exact = reference.trajectory[index * 256];
            let (state, exact_state) = (observation.state, exact.state);
            [
                (state.susceptible - exact_state.susceptible).abs(),
                (state.exposed - exact_state.exposed).abs(),
                (state.infected - exact_state.infected).abs(),
                (state.recovered - exact_state.recovered).abs(),
                (observation.cumulative_infections - exact.cumulative_infections).abs(),
            ]
        })
        .fold(0.0, f64::max);
    let estimated = estimate(1.0).max_absolute_error;
    assert!(
        (0.9..1.1).contains(&(estimated / actual)),
        "estimated {estimated}, actual {actual}"
    );
    // Halving the step divides the error by about 2^4 = 16.
    let halving = estimated / estimate(0.5).max_absolute_error;
    assert!((12.0..20.0).contains(&halving), "{halving}");
}

#[test]
fn seirs_without_initial_infection_never_becomes_endemic() {
    // R0 = 5, but with nobody exposed or infected the infection can never
    // start, and the endemic prevalence is 0.
    let initial = SeirState {
        susceptible: 900.0,
        exposed: 0.0,
        infected: 0.0,
        recovered: 100.0,
    };
    let result = simulate_seirs(PARAMETERS, initial, settings(400.0, 0.25)).unwrap();
    for observation in &result.trajectory {
        assert_eq!(observation.state.exposed, 0.0);
        assert_eq!(observation.state.infected, 0.0);
        assert_eq!(observation.cumulative_infections, 0.0);
    }
    assert_eq!(result.metrics.endemic_prevalence, Some(0.0));
}

#[test]
fn seirs_stability_check_uses_the_step_actually_taken() {
    // A duration shorter than the step is integrated in one step of that
    // duration, so the nominal step of 7 days is never used.
    let result = simulate_seirs(PARAMETERS, INITIAL, settings(1.0, 7.0)).unwrap();
    assert_eq!(result.diagnostics.steps_executed, 1);

    // The single step of 1 day is still too large for rates summing to 4 per day.
    assert_eq!(
        simulate_seirs(parameters(1.0, 1.5, 1.0, 0.5), INITIAL, settings(1.0, 7.0)).unwrap_err(),
        SimulationError::TimeStepTooLarge { max_step_days: 0.5 }
    );
}

#[test]
fn seirs_reports_metric_overflow_before_integrating() {
    // beta / gamma overflows to infinity, and these settings would also need
    // ten billion steps. The parameter problem is reported first.
    assert_eq!(
        simulate_seirs(
            parameters(10.0, 0.1, f64::MIN_POSITIVE, 0.0),
            INITIAL,
            settings(1e9, 0.1)
        )
        .unwrap_err(),
        SimulationError::MetricOverflow
    );
}

#[test]
fn seirs_effective_reproduction_number_rejects_invalid_input() {
    assert_eq!(
        parameters(0.5, 0.0, 0.1, 0.01)
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
    assert_eq!(
        parameters(10.0, 0.1, f64::MIN_POSITIVE, 0.0)
            .effective_reproduction_number(&INITIAL)
            .unwrap_err(),
        SimulationError::MetricOverflow
    );
}
