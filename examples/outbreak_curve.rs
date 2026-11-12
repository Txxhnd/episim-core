//! Prints one SIR outbreak as CSV, one row per simulated day, for a
//! spreadsheet or plotting tool.
//!
//! Run with `cargo run --example outbreak_curve > outbreak.csv`, open the
//! file in a spreadsheet, and chart the susceptible, infected and recovered
//! columns against the day. The last column is the model's effective
//! reproduction number, `R0 * S / N`: the infected peak on the day it falls
//! below 1.

use episim_core::{SimulationError, SimulationSettings, SirParameters, SirState, simulate_sir};

fn main() -> Result<(), SimulationError> {
    // R0 = 3, and people stay infectious for 1 / gamma = 10 days on average.
    let parameters = SirParameters {
        beta_per_day: 0.3,
        gamma_per_day: 0.1,
    };
    let initial_state = SirState {
        susceptible: 990.0,
        infected: 10.0,
        recovered: 0.0,
    };
    // A quarter-day step for accuracy, recorded once per day.
    let settings = SimulationSettings::new(160.0, 0.25)?.with_record_every(4)?;
    let result = simulate_sir(parameters, initial_state, settings)?;

    println!("day,susceptible,infected,recovered,effective_r");
    for observation in &result.trajectory {
        let state = observation.state;
        println!(
            "{},{:.3},{:.3},{:.3},{:.4}",
            observation.time_days,
            state.susceptible,
            state.infected,
            state.recovered,
            parameters.effective_reproduction_number(&state)?
        );
    }
    Ok(())
}
