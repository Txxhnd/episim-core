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
//! The SIRS and SEIRS rows also reuse the SIR and SEIR draws, and take
//! `omega` from a second generator. No compartment of either model is monotone everywhere, so a run is
//! flagged if, over some coarse step, a compartment or the cumulative
//! infections move against the exact flow from the same starting state,
//! where that flow is monotone over the step. The exact flow is a run of 64
//! substeps from that state. A run is also flagged if its incidence or any
//! final compartment, as a fraction of `N`, misses the fine run by more than
//! 0.1.
//!
//! The measurement is ignored by default because it needs a release build:
//! about four minutes for both tests together. `docs/numerical_methods.md`
//! quotes its result.
//! The workflow `.github/workflows/stability-evidence.yml` runs it every
//! Monday and on demand.
//!
//! ```sh
//! cargo test --release --test stability_evidence_tests -- --ignored --nocapture
//! ```
//!
//! A second ignored test, `stability_evidence_with_schedules`, repeats the
//! measurement for runs with a transmission schedule; see its comment.
//!
//! A smaller sample is not kept in the normal suite: with 500 configurations
//! it could not detect either a loosened limit or a missing coupled SEIR rate,
//! while this one detects both.

use episim_core::{
    Observation, SeirParameters, SeirState, SeirsParameters, SimulationError, SimulationSettings,
    SirParameters, SirState, SirsParameters, SisParameters, SisState, TransmissionSchedule,
    simulate_seir, simulate_seirs, simulate_sir, simulate_sirs, simulate_sis,
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
    accepted_seirs: u32,
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

/// SEIR compartments and cumulative infections of an observation.
fn seir_entries(observation: &Observation<SeirState>) -> [f64; 5] {
    let state = observation.state;
    [
        state.susceptible,
        state.exposed,
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
            simulate_sir(parameters, initial, coarse.clone()),
            simulate_sir(parameters, initial, fine.clone()),
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
            simulate_seir(parameters, initial, coarse.clone()),
            simulate_seir(parameters, initial, fine.clone()),
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
            simulate_sis(parameters, initial, coarse.clone()),
            simulate_sis(parameters, initial, fine.clone()),
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
            simulate_sirs(parameters, initial, coarse.clone()),
            simulate_sirs(parameters, initial, fine.clone()),
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

        let parameters = SeirsParameters {
            beta_per_day: beta,
            sigma_per_day: sigma,
            gamma_per_day: gamma,
            omega_per_day: omega,
        };
        let total = a + b + c + d;
        let initial = SeirState {
            susceptible: POPULATION * a / total,
            exposed: POPULATION * d / total,
            infected: POPULATION * b / total,
            recovered: POPULATION * c / total,
        };
        if let (Ok(run), Ok(reference)) = (
            simulate_seirs(parameters, initial, coarse.clone()),
            simulate_seirs(parameters, initial, fine.clone()),
        ) {
            evidence.accepted_seirs += 1;
            let wrong_direction = run.trajectory.windows(2).any(|pair| {
                let step = pair[1].time_days - pair[0].time_days;
                let local = SimulationSettings::new(step, step / 64.0).unwrap();
                // A start with a tolerated negative value is not a valid state.
                let Ok(exact) = simulate_seirs(parameters, pair[0].state, local) else {
                    return false;
                };
                let exact: Vec<_> = exact.trajectory.iter().map(seir_entries).collect();
                against_exact_flow(seir_entries(&pair[0]), seir_entries(&pair[1]), &exact)
            });
            let (last, fine_last) = (
                seir_entries(run.trajectory.last().unwrap()),
                seir_entries(reference.trajectory.last().unwrap()),
            );
            // The last entry is the cumulative infections, so this covers the incidence.
            let final_error = (0..5)
                .map(|entry| (last[entry] - fine_last[entry]).abs() / POPULATION)
                .fold(0.0, f64::max);
            if wrong_direction || final_error > MAX_INCIDENCE_ERROR {
                evidence.flagged.push(format!(
                    "SEIRS {parameters:?} {initial:?} duration {duration}: \
                     wrong direction {wrong_direction}, final error {final_error}"
                ));
            }
        }
    }
    evidence
}

#[test]
#[ignore = "minutes in release mode; run with --release -- --ignored --nocapture"]
fn full_stability_evidence() {
    let evidence = gather(250_000);
    println!(
        "accepted SIR runs: {}, accepted SEIR runs: {}, accepted SIS runs: {}, \
         accepted SIRS runs: {}, accepted SEIRS runs: {}, flagged: {}",
        evidence.accepted_sir,
        evidence.accepted_seir,
        evidence.accepted_sis,
        evidence.accepted_sirs,
        evidence.accepted_seirs,
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
    // The SIRS limit, beta + gamma + omega, accepts about one configuration
    // in seven at a one-day step.
    assert!(
        evidence.accepted_sirs >= 30_000,
        "{}",
        evidence.accepted_sirs
    );
    // The SEIRS limit is conservative; it accepts about one configuration in 29.
    assert!(
        evidence.accepted_seirs >= 6_000,
        "{}",
        evidence.accepted_seirs
    );
    assert!(evidence.flagged.is_empty(), "{:#?}", evidence.flagged);
}

/// Whether a scheduled run goes wrong: over some step, an entry moves against
/// the exact flow from the same state, or a final entry misses the fine run
/// by more than 0.1 N. Every step lies in one segment of the schedule, so the
/// exact flow over it is a constant-parameter run at that step's multiplier,
/// `local(state, multiplier, step)`, with 64 substeps.
fn scheduled_run_flagged<S: Copy, const K: usize>(
    run: &[Observation<S>],
    fine_last: &Observation<S>,
    schedule: &TransmissionSchedule,
    entries: fn(&Observation<S>) -> [f64; K],
    local: impl Fn(S, f64, f64) -> Result<Vec<Observation<S>>, SimulationError>,
) -> bool {
    let wrong_direction = run.windows(2).any(|pair| {
        let step = pair[1].time_days - pair[0].time_days;
        let multiplier = schedule.multiplier_at(pair[0].time_days);
        // A start with a tolerated negative value is not a valid state.
        let Ok(exact) = local(pair[0].state, multiplier, step) else {
            return false;
        };
        let exact: Vec<_> = exact.iter().map(entries).collect();
        against_exact_flow(entries(&pair[0]), entries(&pair[1]), &exact)
    });
    let (last, reference) = (entries(run.last().unwrap()), entries(fine_last));
    // The last entry is the cumulative infections, so this covers the incidence.
    let final_error = (0..K)
        .map(|entry| (last[entry] - reference[entry]).abs() / POPULATION)
        .fold(0.0, f64::max);
    wrong_direction || final_error > MAX_INCIDENCE_ERROR
}

/// SIS compartments and cumulative infections of an observation.
fn sis_entries(observation: &Observation<SisState>) -> [f64; 3] {
    let state = observation.state;
    [
        state.susceptible,
        state.infected,
        observation.cumulative_infections,
    ]
}

/// Runs a model coarse and fine with `$schedule`, counts the run as accepted
/// when both succeed, and records it as flagged per `scheduled_run_flagged`.
/// `$scaled` builds the model's parameters with beta multiplied by its argument.
macro_rules! scheduled_row {
    ($evidence:expr, $index:expr, $label:expr, $simulate:ident, $initial:expr, $scaled:expr,
     $schedule:expr, $duration:expr, $entries:expr) => {{
        let with = |step: f64| {
            SimulationSettings::new($duration, step)
                .unwrap()
                .with_transmission_schedule($schedule.clone())
                .unwrap()
        };
        if let (Ok(run), Ok(reference)) = (
            $simulate($scaled(1.0), $initial, with(1.0)),
            $simulate($scaled(1.0), $initial, with(1.0 / 64.0)),
        ) {
            $evidence.accepted[$index] += 1;
            let local = |state, multiplier: f64, step: f64| {
                $simulate(
                    $scaled(multiplier),
                    state,
                    SimulationSettings::new(step, step / 64.0).unwrap(),
                )
                .map(|result| result.trajectory)
            };
            if scheduled_run_flagged(
                &run.trajectory,
                reference.trajectory.last().unwrap(),
                &$schedule,
                $entries,
                local,
            ) {
                $evidence.flagged.push(format!(
                    "{} {:?} {:?} duration {}",
                    $label,
                    $scaled(1.0),
                    $schedule,
                    $duration
                ));
            }
        }
    }};
}

#[derive(Debug, Default)]
struct ScheduledEvidence {
    /// Accepted runs of SIR, SEIR, SIS, SIRS and SEIRS, in that order.
    accepted: [u32; 5],
    flagged: Vec<String>,
}

/// As `gather`, for runs with a transmission schedule. The configurations
/// repeat `gather`'s draws (generators 42 and 7), and a third generator
/// (seed 11) adds a three-segment schedule: changes on whole days, so on the
/// one-day and 1/64-day grids, and multipliers in `[0, 2)`. The stability
/// limit must then hold for the largest beta in force.
fn gather_scheduled(configurations: u32) -> ScheduledEvidence {
    let mut generator = Generator(42);
    let mut omega_generator = Generator(7);
    let mut schedule_generator = Generator(11);
    let mut evidence = ScheduledEvidence::default();
    for _ in 0..configurations {
        let beta = generator.unit() * 2.1;
        let sigma = (generator.unit() * 2.1).max(1e-3);
        let gamma = (generator.unit() * 2.1).max(1e-3);
        let mut duration = 2.0 + (generator.unit() * 40.0).floor();
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
        let first_change = 1.0 + (schedule_generator.unit() * (duration - 1.0)).floor();
        let second_change = first_change + 1.0 + (schedule_generator.unit() * 40.0).floor();
        let schedule = TransmissionSchedule::new(&[
            (first_change, schedule_generator.unit() * 2.0),
            (second_change, schedule_generator.unit() * 2.0),
        ])
        .unwrap();

        let total = a + b + c;
        let sir0 = SirState {
            susceptible: POPULATION * a / total,
            infected: POPULATION * b / total,
            recovered: POPULATION * c / total,
        };
        let total = a + b + c + d;
        let seir0 = SeirState {
            susceptible: POPULATION * a / total,
            exposed: POPULATION * d / total,
            infected: POPULATION * b / total,
            recovered: POPULATION * c / total,
        };
        let sis0 = SisState {
            susceptible: POPULATION * a / (a + b),
            infected: POPULATION * b / (a + b),
        };

        let sir = |m: f64| SirParameters {
            beta_per_day: beta * m,
            gamma_per_day: gamma,
        };
        scheduled_row!(
            evidence,
            0,
            "SIR",
            simulate_sir,
            sir0,
            sir,
            schedule,
            duration,
            sir_entries
        );
        let seir = |m: f64| SeirParameters {
            beta_per_day: beta * m,
            sigma_per_day: sigma,
            gamma_per_day: gamma,
        };
        scheduled_row!(
            evidence,
            1,
            "SEIR",
            simulate_seir,
            seir0,
            seir,
            schedule,
            duration,
            seir_entries
        );
        let sis = |m: f64| SisParameters {
            beta_per_day: beta * m,
            gamma_per_day: gamma,
        };
        scheduled_row!(
            evidence,
            2,
            "SIS",
            simulate_sis,
            sis0,
            sis,
            schedule,
            duration,
            sis_entries
        );
        let sirs = |m: f64| SirsParameters {
            beta_per_day: beta * m,
            gamma_per_day: gamma,
            omega_per_day: omega,
        };
        scheduled_row!(
            evidence,
            3,
            "SIRS",
            simulate_sirs,
            sir0,
            sirs,
            schedule,
            duration,
            sir_entries
        );
        let seirs = |m: f64| SeirsParameters {
            beta_per_day: beta * m,
            sigma_per_day: sigma,
            gamma_per_day: gamma,
            omega_per_day: omega,
        };
        scheduled_row!(
            evidence,
            4,
            "SEIRS",
            simulate_seirs,
            seir0,
            seirs,
            schedule,
            duration,
            seir_entries
        );
    }
    evidence
}

#[test]
#[ignore = "minutes in release mode; run with --release -- --ignored --nocapture"]
fn stability_evidence_with_schedules() {
    let evidence = gather_scheduled(250_000);
    let [sir, seir, sis, sirs, seirs] = evidence.accepted;
    println!(
        "with schedules: accepted SIR runs: {sir}, accepted SEIR runs: {seir}, \
         accepted SIS runs: {sis}, accepted SIRS runs: {sirs}, \
         accepted SEIRS runs: {seirs}, flagged: {}",
        evidence.flagged.len()
    );
    // Enough accepted runs per model that a pass means something.
    assert!(
        evidence.accepted.iter().all(|&count| count >= 5_000),
        "{evidence:?}"
    );
    assert!(evidence.flagged.is_empty(), "{:#?}", evidence.flagged);
}
