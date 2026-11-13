use episim_core::{
    SimulationError, SimulationSettings, SirParameters, SirState, SirsParameters,
    estimate_sirs_step_error, simulate_sir, simulate_sirs,
};

const PARAMETERS: SirsParameters = SirsParameters {
    beta_per_day: 0.3,
    gamma_per_day: 0.1,
    omega_per_day: 0.01,
};

const INITIAL: SirState = SirState {
    susceptible: 990.0,
    infected: 10.0,
    recovered: 0.0,
};

fn settings(duration_days: f64, step_days: f64) -> SimulationSettings {
    SimulationSettings::new(duration_days, step_days).unwrap()
}

fn parameters(beta_per_day: f64, gamma_per_day: f64, omega_per_day: f64) -> SirsParameters {
    SirsParameters {
        beta_per_day,
        gamma_per_day,
        omega_per_day,
    }
}

#[test]
fn sirs_without_waning_is_sir() {
    // With omega = 0 every waning term is exactly zero, so the runs agree
    // to the last bit, cumulative infections included.
    let sirs = simulate_sirs(parameters(0.3, 0.1, 0.0), INITIAL, settings(160.0, 0.25)).unwrap();
    let sir = simulate_sir(
        SirParameters {
            beta_per_day: 0.3,
            gamma_per_day: 0.1,
        },
        INITIAL,
        settings(160.0, 0.25),
    )
    .unwrap();
    assert_eq!(sirs.trajectory, sir.trajectory);
    assert_eq!(sirs.diagnostics, sir.diagnostics);
    assert_eq!(
        sirs.metrics.herd_immunity_threshold,
        sir.metrics.herd_immunity_threshold
    );
    // Without waning nobody stays infected for good.
    assert_eq!(sirs.metrics.endemic_prevalence, Some(0.0));
    assert_eq!(sirs.metrics.final_size_fraction, None);
}

#[test]
fn sirs_settles_at_the_endemic_equilibrium() {
    // At equilibrium S* = N / R0 and I* / N = (1 - 1/R0) * omega / (omega + gamma),
    // or S* = N and I* = 0 when R0 <= 1. The last set has R0 = 1/2.
    for (beta, gamma, omega, infected) in [
        (0.3, 0.1, 0.01, 60.606_060_606),
        (0.5, 0.2, 0.05, 120.0),
        (1.8, 0.3, 1.2, 666.666_666_667),
        (0.1, 0.2, 0.05, 0.0),
    ] {
        let result = simulate_sirs(
            parameters(beta, gamma, omega),
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
fn sirs_conserves_the_population() {
    let result = simulate_sirs(parameters(1.8, 0.3, 1.2), INITIAL, settings(200.0, 0.25)).unwrap();
    for observation in &result.trajectory {
        let state = observation.state;
        let total = state.susceptible + state.infected + state.recovered;
        assert!(
            (total - 1000.0).abs() <= 1e-6,
            "day {}",
            observation.time_days
        );
    }
    assert!(result.diagnostics.maximum_relative_population_drift <= 1e-9);
}

#[test]
fn sirs_counts_reinfections() {
    let result = simulate_sirs(PARAMETERS, INITIAL, settings(1000.0, 0.25)).unwrap();
    assert_eq!(result.trajectory[0].cumulative_infections, 0.0);
    // C - (S(0) - S) grows at omega * R: the people who lost their immunity
    // and were counted again when reinfected.
    let reinfections = |o: &episim_core::Observation<SirState>| {
        o.cumulative_infections - (INITIAL.susceptible - o.state.susceptible)
    };
    for pair in result.trajectory.windows(2) {
        assert!(pair[1].cumulative_infections >= pair[0].cumulative_infections);
        assert!(reinfections(&pair[1]) >= reinfections(&pair[0]) - 1e-9);
    }
    // At equilibrium gamma * I* = 6.06 infections a day, so over 1000 days
    // each person is infected several times.
    assert!(
        result.metrics.cumulative_incidence_fraction > 4.0,
        "{:?}",
        result.metrics
    );
    // R0 = 3; the threshold keeps its meaning at any given time.
    assert_eq!(result.metrics.herd_immunity_threshold, Some(2.0 / 3.0));
}

#[test]
fn sirs_effective_reproduction_number_sets_the_direction_of_the_infected() {
    // The infected rise and fall in damped waves around the equilibrium,
    // always rising while effective R is above 1 and falling while it is
    // below. A short step keeps effective R nearly constant over each step.
    let result = simulate_sirs(PARAMETERS, INITIAL, settings(400.0, 0.05)).unwrap();
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
            let change = pair[1].state.infected - pair[0].state.infected;
            assert_eq!(change > 0.0, r > 1.0, "day {}: R {r}", pair[0].time_days);
        }
    }
    // Effective R falls below 1 and rises above it again: waning immunity
    // refills the susceptibles.
    assert!(crossings >= 3, "{crossings}");
}

#[test]
fn sirs_rejects_invalid_parameters() {
    let run = |parameters| simulate_sirs(parameters, INITIAL, settings(10.0, 0.5));
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            run(parameters(0.3, 0.1, invalid)).unwrap_err(),
            SimulationError::InvalidParameter {
                name: "omega_per_day"
            }
        );
        assert_eq!(
            run(parameters(invalid, 0.1, 0.01)).unwrap_err(),
            SimulationError::InvalidParameter {
                name: "beta_per_day"
            }
        );
    }
    assert_eq!(
        run(parameters(0.3, 0.0, 0.01)).unwrap_err(),
        SimulationError::InvalidParameter {
            name: "gamma_per_day"
        }
    );
    let negative = SirState {
        recovered: -1.0,
        ..INITIAL
    };
    assert_eq!(
        simulate_sirs(PARAMETERS, negative, settings(10.0, 0.5)).unwrap_err(),
        SimulationError::InvalidInitialState {
            compartment: "recovered"
        }
    );
}

#[test]
fn sirs_rejects_step_too_large_for_stable_integration() {
    let run = |beta, gamma, omega, step| {
        simulate_sirs(
            parameters(beta, gamma, omega),
            INITIAL,
            settings(10.0, step),
        )
    };
    // The fastest rate is beta + gamma + omega. Each case sums to 4 per day,
    // allowing a step of at most 2 / 4 = 0.5 days, although neither
    // beta + gamma nor omega alone reaches 4.
    let too_large = SimulationError::TimeStepTooLarge { max_step_days: 0.5 };
    for (beta, gamma, omega) in [(2.5, 1.0, 0.5), (0.2, 0.1, 3.7)] {
        assert_eq!(run(beta, gamma, omega, 0.51).unwrap_err(), too_large);
        assert!(run(beta, gamma, omega, 0.5).is_ok(), "omega {omega}");
    }
}

#[test]
fn sirs_step_error_estimate_matches_the_error_against_a_fine_reference() {
    let parameters = parameters(0.5, 0.1, 0.05);
    let run = |step| simulate_sirs(parameters, INITIAL, settings(40.0, step)).unwrap();
    let estimate =
        |step| estimate_sirs_step_error(parameters, INITIAL, settings(40.0, step)).unwrap();
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
fn sirs_without_initial_infection_never_becomes_endemic() {
    // R0 = 3, and immunity wanes, but with nobody infected the infection can
    // never start: people only move from R back to S, and the endemic
    // prevalence is 0.
    let initial = SirState {
        susceptible: 900.0,
        infected: 0.0,
        recovered: 100.0,
    };
    let result = simulate_sirs(PARAMETERS, initial, settings(400.0, 0.25)).unwrap();
    for observation in &result.trajectory {
        assert_eq!(observation.state.infected, 0.0);
        assert_eq!(observation.cumulative_infections, 0.0);
    }
    assert_eq!(result.metrics.endemic_prevalence, Some(0.0));
}

#[test]
fn sirs_stability_check_uses_the_step_actually_taken() {
    // A duration shorter than the step is integrated in one step of that
    // duration, so the nominal step of 7 days is never used.
    let result = simulate_sirs(PARAMETERS, INITIAL, settings(1.0, 7.0)).unwrap();
    assert_eq!(result.diagnostics.steps_executed, 1);

    // The single step of 1 day is still too large for rates summing to 4 per day.
    assert_eq!(
        simulate_sirs(parameters(0.2, 0.1, 3.7), INITIAL, settings(1.0, 7.0)).unwrap_err(),
        SimulationError::TimeStepTooLarge { max_step_days: 0.5 }
    );
}

#[test]
fn sirs_reports_metric_overflow_before_integrating() {
    // beta / gamma overflows to infinity, and these settings would also need
    // ten billion steps. The parameter problem is reported first.
    assert_eq!(
        simulate_sirs(
            parameters(10.0, f64::MIN_POSITIVE, 0.0),
            INITIAL,
            settings(1e9, 0.1)
        )
        .unwrap_err(),
        SimulationError::MetricOverflow
    );
}

#[test]
fn sirs_effective_reproduction_number_rejects_invalid_input() {
    assert_eq!(
        parameters(0.3, 0.1, -1.0)
            .effective_reproduction_number(&INITIAL)
            .unwrap_err(),
        SimulationError::InvalidParameter {
            name: "omega_per_day"
        }
    );
    let empty = SirState {
        susceptible: 0.0,
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
        parameters(10.0, f64::MIN_POSITIVE, 0.0)
            .effective_reproduction_number(&INITIAL)
            .unwrap_err(),
        SimulationError::MetricOverflow
    );
}

#[test]
fn sirs_rejects_a_step_that_leaves_a_quarter_of_the_population_wrong() {
    // Found by the stability evidence: with these rates and most people
    // infected, one step of a day ended 26% of the population away from a
    // fine run (recovered fell by 93 where it should rise by 166). The old
    // limit, max(beta + gamma, omega) = 1.87, accepted it; the sum, 3.73,
    // does not.
    let parameters = parameters(0.3929114889456895, 1.4664779429151176, 1.8734350589317061);
    let initial = SirState {
        susceptible: 46.97145797103499,
        infected: 856.3605939210153,
        recovered: 96.66794810794975,
    };
    assert!(matches!(
        simulate_sirs(parameters, initial, settings(1.0, 1.0)),
        Err(SimulationError::TimeStepTooLarge { .. })
    ));
}
