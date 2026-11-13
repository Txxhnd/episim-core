//! Shows how the size and timing of an SIR outbreak depend on R0.
//!
//! Run with `cargo run --example r0_sweep`. Each row simulates one year in a
//! population of 10,000 with one infected person at the start. People stay
//! infectious for 1 / gamma = 10 days on average; only the transmission rate
//! beta changes, so `R0 = beta / gamma`.
//!
//! Below R0 = 1 the single case never spreads. Just above it the epidemic is
//! so slow that it is still going at day 365, so the infected share by then
//! is below the final size. Further up, the peak comes sooner and higher, and
//! the final size approaches the whole population.

use episim_core::{SimulationError, SimulationSettings, SirParameters, SirState, simulate_sir};

const GAMMA_PER_DAY: f64 = 0.1;
const POPULATION: f64 = 10_000.0;

fn main() -> Result<(), SimulationError> {
    println!("  R0 | herd immunity | peak infected | peak day | infected by day 365 | final size");
    // R0 from 0.5 to 4.0 in steps of 0.25.
    for index in 0..15_u8 {
        let r0 = 0.5 + 0.25 * f64::from(index);
        let parameters = SirParameters {
            beta_per_day: r0 * GAMMA_PER_DAY,
            gamma_per_day: GAMMA_PER_DAY,
        };
        let initial_state = SirState {
            susceptible: POPULATION - 1.0,
            infected: 1.0,
            recovered: 0.0,
        };
        // A quarter-day step, recorded once per day.
        let settings = SimulationSettings::new(365.0, 0.25)?.with_record_every(4)?;
        let metrics = simulate_sir(parameters, initial_state, settings)?.metrics;
        // SIR reports both; only models with waning immunity leave them out.
        let (Some(threshold), Some(final_size)) =
            (metrics.herd_immunity_threshold, metrics.final_size_fraction)
        else {
            unreachable!("SIR reports a herd-immunity threshold and a final size");
        };

        println!(
            "{r0:4.2} | {:12.1}% | {:12.1}% | {:8} | {:18.1}% | {:9.1}%",
            100.0 * threshold,
            100.0 * metrics.peak_infectious / POPULATION,
            metrics.peak_infectious_time_days,
            100.0 * metrics.cumulative_incidence_fraction,
            100.0 * final_size
        );
    }
    Ok(())
}
