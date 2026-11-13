//! Evidence for the timestep stability limit.
//!
//! Every `simulate_*` function rejects a run when the model's fastest decay
//! rate times the step exceeds 2. That limit was chosen from
//! measurement, not derived, so this file keeps the measurement runnable.
//! It draws random configurations, runs each one at a step of one day and
//! again at 1/64 day, and counts accepted runs whose coarse result is
//! qualitatively wrong.
//!
//! The SIS rows reuse the SIR draws, so adding them changed no other row. The
//! exact SIS solution is logistic: `I` moves monotonically towards its
//! equilibrium `I*` and never crosses it, and the cumulative infections never
//! fall. A run is flagged if it breaks either property, or if its incidence
//! or final infected fraction misses the fine run by more than 0.1.
//!
//! The SIRS rows also reuse the SIR draws, and take `omega` from a second
//! generator. No compartment of SIRS is monotone everywhere, so a run is
//! flagged if, over some coarse step, a compartment or the cumulative
//! infections move against the exact flow from the same starting state,
//! where that flow is monotone over the step. The exact flow is a run of 64
//! substeps from that state. A run is also flagged if its incidence or any
//! final compartment, as a fraction of `N`, misses the fine run by more than
//! 0.1.
//!
//! The measurement is ignored by default because it needs a release build to
//! finish in about a minute. `docs/numerical_methods.md` quotes its result.
//!
//! ```sh
//! cargo test --release --test stability_evidence_tests -- --ignored --nocapture
//! ```
//!
//! A smaller sample is not kept in the normal suite: with 500 configurations
//! it could not detect either a loosened limit or a missing coupled SEIR rate,
//! while this one detects both.

use episim_core::{
    Observation, SeirParameters, SeirState, SimulationSettings, SirParameters, SirState,
    SirsParameters, SisParameters, SisState, simulate_seir, simulate_sir, simulate_sirs,
    simulate_sis,
};

const POPULATION: f64 = 1000.0;
/// Changes smaller than this are rounding, not a wrong direction.
const TOLERANCE: f64 = 1e-9 * POPULATION;
/// Largest accepted difference in the cumulative incidence fraction.
const MAX_INCIDENCE_ERROR: f64 = 0.1;

/// A small deterministic generator, so the evidence is reproducible without
/// adding a dependency.
struct Generator(u64);

impl Generator {
    /// Next value in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[derive(Debug, Default)]
struct Evidence {
    accepted_sir: u32,
    accepted_seir: u32,
    accepted_sis: u32,
    accepted_sirs: u32,
    /// One line per accepted run that is qualitatively wrong.
    flagged: Vec<String>,
}

/// Whether some entry moved from `before` to `after` against the exact flow
/// from `before`, over a step on which that flow is monotone.
fn against_exact_flow<const N: usize>(
    before: [f64; N],
    after: [f64; N],
    exact: &[[f64; N]],
) -> bool {
    (0..N).any(|entry| {
        let rising = exact.windows(2).all(|pair| pair[1][entry] > pair[0][entry]);
        let falling = exact.windows(2).all(|pair| pair[1][entry] < pair[0][entry]);
        let change = after[entry] - before[entry];
        (rising && change < -TOLERANCE) || (falling && change > TOLERANCE)
    })
}

/// SIR compartments and cumulative infections of an observation.
fn sir_entries(observation: &Observation<SirState>) -> [f64; 4] {
    let state = observation.state;
    [
        state.susceptible,
        state.infected,
        state.recovered,
        observation.cumulative_infections,
    ]
}

/// Runs `configurations` random cases of each model and collects the evidence.
fn gather(configurations: u32) -> Evidence {
    let mut generator = Generator(42);
    // A second generator, so that adding waning immunity changed no other draw.
    let mut omega_generator = Generator(7);
    let mut evidence = Evidence::default();
    for _ in 0..configurations {
        // Rates up to 2.1 per day with a one-day step reach well past the
        // limit, so both the accepting and the rejecting side are exercised.
        let beta = generator.unit() * 2.1;
        let sigma = (generator.unit() * 2.1).max(1e-3);
        let gamma = (generator.unit() * 2.1).max(1e-3);
        let mut duration = 2.0 + (generator.unit() * 40.0).floor();
        // Some durations end with a shortened final step.
        if generator.unit() < 0.3 {
            duration += generator.unit();
        }
        let (a, b, c, d) = (
            generator.unit(),
            generator.unit(),
            generator.unit() * 0.2,
            generator.unit(),
        );
        let omega = (omega_generator.unit() * 2.1).max(1e-3);
        let coarse = SimulationSettings::new(duration, 1.0).unwrap();
        let fine = SimulationSettings::new(duration, 1.0 / 64.0).unwrap();

        let parameters = SirParameters {
            beta_per_day: beta,
            gamma_per_day: gamma,
        };
        let total = a + b + c;
        let initial = SirState {
            susceptible: POPULATION * a / total,
            infected: POPULATION * b / total,
            recovered: POPULATION * c / total,
        };
        if let (Ok(run), Ok(reference)) = (
            simulate_sir(parameters, initial, coarse),
            simulate_sir(parameters, initial, fine),
        ) {
            evidence.accepted_sir += 1;
            let wrong_direction = run.trajectory.windows(2).any(|pair| {
                let (before, after) = (pair[0].state, pair[1].state);
                after.recovered < before.recovered - TOLERANCE
                    || after.susceptible > before.susceptible + TOLERANCE
                    // Below the epidemic threshold the infected must decline.
                    || (beta * before.susceptible / POPULATION < gamma
                        && after.infected > before.infected + TOLERANCE)
            });
            let incidence_error = (run.metrics.cumulative_incidence_fraction
                - reference.metrics.cumulative_incidence_fraction)
                .abs();
            if wrong_direction || incidence_error > MAX_INCIDENCE_ERROR {
                evidence.flagged.push(format!(
                    "SIR {parameters:?} {initial:?} duration {duration}: \
                     wrong direction {wrong_direction}, incidence error {incidence_error}"
                ));
            }
        }

        let parameters = SeirParameters {
            beta_per_day: beta,
            sigma_per_day: sigma,
            gamma_per_day: gamma,
        };
        let total = a + b + c + d;
        let initial = SeirState {
            susceptible: POPULATION * a / total,
            exposed: POPULATION * d / total,
            infected: POPULATION * b / total,
            recovered: POPULATION * c / total,
        };
        if let (Ok(run), Ok(reference)) = (
            simulate_seir(parameters, initial, coarse),
            simulate_seir(parameters, initial, fine),
        ) {
            evidence.accepted_seir += 1;
            let wrong_direction = run.trajectory.windows(2).any(|pair| {
                let (before, after) = (pair[0].state, pair[1].state);
                after.recovered < before.recovered - TOLERANCE
                    || after.susceptible > before.susceptible + TOLERANCE
                    // Below the epidemic threshold the infected and exposed
                    // together must decline.
                    || (beta * before.susceptible / POPULATION < gamma
                        && after.exposed + after.infected
                            > before.exposed + before.infected + TOLERANCE)
            });
            let incidence_error = (run.metrics.cumulative_incidence_fraction
                - reference.metrics.cumulative_incidence_fraction)
                .abs();
            if wrong_direction || incidence_error > MAX_INCIDENCE_ERROR {
                evidence.flagged.push(format!(
                    "SEIR {parameters:?} {initial:?} duration {duration}: \
                     wrong direction {wrong_direction}, incidence error {incidence_error}"
                ));
            }
        }

        let parameters = SisParameters {
            beta_per_day: beta,
            gamma_per_day: gamma,
        };
        let initial = SisState {
            susceptible: POPULATION * a / (a + b),
            infected: POPULATION * b / (a + b),
        };
        if let (Ok(run), Ok(reference)) = (
            simulate_sis(parameters, initial, coarse),
            simulate_sis(parameters, initial, fine),
        ) {
            evidence.accepted_sis += 1;
            let equilibrium = run.metrics.endemic_prevalence.unwrap() * POPULATION;
            let wrong_direction = run.trajectory.windows(2).any(|pair| {
                let (before, after) = (pair[0].state.infected, pair[1].state.infected);
                // Away from the equilibrium, or across it.
                (after - equilibrium).abs() > (before - equilibrium).abs() + TOLERANCE
                    || ((before - equilibrium) * (after - equilibrium) < 0.0
                        && (after - equilibrium).abs() > TOLERANCE)
                    || pair[1].cumulative_infections < pair[0].cumulative_infections - TOLERANCE
            });
            let incidence_error = (run.metrics.cumulative_incidence_fraction
                - reference.metrics.cumulative_incidence_fraction)
                .abs();
            let final_error = (run.trajectory.last().unwrap().state.infected
                - reference.trajectory.last().unwrap().state.infected)
                .abs()
                / POPULATION;
            if wrong_direction
                || incidence_error > MAX_INCIDENCE_ERROR
                || final_error > MAX_INCIDENCE_ERROR
            {
                evidence.flagged.push(format!(
                    "SIS {parameters:?} {initial:?} duration {duration}: \
                     wrong direction {wrong_direction}, incidence error {incidence_error}, \
                     final error {final_error}"
                ));
            }
        }

        let parameters = SirsParameters {
            beta_per_day: beta,
            gamma_per_day: gamma,
            omega_per_day: omega,
        };
        let total = a + b + c;
        let initial = SirState {
            susceptible: POPULATION * a / total,
            infected: POPULATION * b / total,
            recovered: POPULATION * c / total,
        };
        if let (Ok(run), Ok(reference)) = (
            simulate_sirs(parameters, initial, coarse),
            simulate_sirs(parameters, initial, fine),
        ) {
            evidence.accepted_sirs += 1;
            let wrong_direction = run.trajectory.windows(2).any(|pair| {
                let step = pair[1].time_days - pair[0].time_days;
                let local = SimulationSettings::new(step, step / 64.0).unwrap();
                // A start with a tolerated negative value is not a valid state.
                let Ok(exact) = simulate_sirs(parameters, pair[0].state, local) else {
                    return false;
                };
                let exact: Vec<_> = exact.trajectory.iter().map(sir_entries).collect();
                against_exact_flow(sir_entries(&pair[0]), sir_entries(&pair[1]), &exact)
            });
            let (last, fine_last) = (
                sir_entries(run.trajectory.last().unwrap()),
                sir_entries(reference.trajectory.last().unwrap()),
            );
            // The last entry is the cumulative infections, so this covers the incidence.
            let final_error = (0..4)
                .map(|entry| (last[entry] - fine_last[entry]).abs() / POPULATION)
                .fold(0.0, f64::max);
            if wrong_direction || final_error > MAX_INCIDENCE_ERROR {
                evidence.flagged.push(format!(
                    "SIRS {parameters:?} {initial:?} duration {duration}: \
                     wrong direction {wrong_direction}, final error {final_error}"
                ));
            }
        }
    }
    evidence
}

#[test]
#[ignore = "about a minute in release mode; run with --release -- --ignored --nocapture"]
fn full_stability_evidence() {
    let evidence = gather(250_000);
    println!(
        "accepted SIR runs: {}, accepted SEIR runs: {}, accepted SIS runs: {}, \
         accepted SIRS runs: {}, flagged: {}",
        evidence.accepted_sir,
        evidence.accepted_seir,
        evidence.accepted_sis,
        evidence.accepted_sirs,
        evidence.flagged.len()
    );
    assert!(
        evidence.accepted_sir >= 200_000,
        "{}",
        evidence.accepted_sir
    );
    assert!(
        evidence.accepted_seir >= 80_000,
        "{}",
        evidence.accepted_seir
    );
    assert!(evidence.accepted_sis >= 80_000, "{}", evidence.accepted_sis);
    assert!(
        evidence.accepted_sirs >= 80_000,
        "{}",
        evidence.accepted_sirs
    );
    assert!(evidence.flagged.is_empty(), "{:#?}", evidence.flagged);
}
