//! Evidence for the timestep stability limit.
//!
//! `simulate_sir` and `simulate_seir` reject a run when the model's fastest
//! decay rate times the step exceeds 2. That limit was chosen from
//! measurement, not derived, so this file keeps the measurement runnable.
//! It draws random configurations, runs each one at a step of one day and
//! again at 1/64 day, and counts accepted runs whose coarse result is
//! qualitatively wrong.
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
    SeirParameters, SeirState, SimulationSettings, SirParameters, SirState, simulate_seir,
    simulate_sir,
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
    /// One line per accepted run that is qualitatively wrong.
    flagged: Vec<String>,
}

/// Runs `configurations` random cases of each model and collects the evidence.
fn gather(configurations: u32) -> Evidence {
    let mut generator = Generator(42);
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
    }
    evidence
}

#[test]
#[ignore = "about a minute in release mode; run with --release -- --ignored --nocapture"]
fn full_stability_evidence() {
    let evidence = gather(250_000);
    println!(
        "accepted SIR runs: {}, accepted SEIR runs: {}, flagged: {}",
        evidence.accepted_sir,
        evidence.accepted_seir,
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
    assert!(evidence.flagged.is_empty(), "{:#?}", evidence.flagged);
}
