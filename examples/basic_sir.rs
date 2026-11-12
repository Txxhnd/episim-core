//! Simulates a basic SIR outbreak and prints a summary.
//!
//! Run with `cargo run --example basic_sir`.

use episim_core::{SimulationError, SimulationSettings, SirParameters, SirState, simulate_sir};

fn main() -> Result<(), SimulationError> {
    let parameters = SirParameters {
        beta_per_day: 0.3,
        gamma_per_day: 0.1,
    };
    let initial_state = SirState {
        susceptible: 990.0,
        infected: 10.0,
        recovered: 0.0,
    };
    let settings = SimulationSettings::new(160.0, 0.25)?;

    let result = simulate_sir(parameters, initial_state, settings)?;
    let metrics = result.metrics;
    let diagnostics = result.diagnostics;

    println!("SIR model");
    println!(
        "  beta = {} /day, gamma = {} /day",
        parameters.beta_per_day, parameters.gamma_per_day
    );
    println!(
        "  initial S = {}, I = {}, R = {}",
        initial_state.susceptible, initial_state.infected, initial_state.recovered
    );
    println!(
        "  duration = {} days, step = {} days",
        settings.duration_days(),
        settings.step_days()
    );
    println!("Metrics");
    println!(
        "  basic reproduction number R0 = {:.2}",
        metrics.basic_reproduction_number
    );
    println!(
        "  herd-immunity threshold = {:.1}% of the population",
        100.0 * metrics.herd_immunity_threshold
    );
    println!(
        "  peak infectious = {:.1} on day {}",
        metrics.peak_infectious, metrics.peak_infectious_time_days
    );
    println!(
        "  cumulative infections since start = {:.1} ({:.1}% of the population)",
        metrics.cumulative_infections_since_start,
        100.0 * metrics.cumulative_incidence_fraction
    );
    println!(
        "  final size when the epidemic is over = {:.1}% of the population",
        100.0 * metrics.final_size_fraction
    );
    println!("Diagnostics");
    println!("  steps executed = {}", diagnostics.steps_executed);
    println!(
        "  maximum relative population drift = {:e}",
        diagnostics.maximum_relative_population_drift
    );
    println!(
        "  tolerated negative values = {}",
        diagnostics.tolerated_negative_value_count
    );
    Ok(())
}
