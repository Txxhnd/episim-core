//! Simulates a basic SEIR outbreak and prints a summary.
//!
//! Run with `cargo run --example basic_seir`.

use episim_core::{SeirParameters, SeirState, SimulationError, SimulationSettings, simulate_seir};

fn main() -> Result<(), SimulationError> {
    let parameters = SeirParameters {
        beta_per_day: 0.5,
        sigma_per_day: 0.2,
        gamma_per_day: 0.1,
    };
    let initial_state = SeirState {
        susceptible: 980.0,
        exposed: 10.0,
        infected: 10.0,
        recovered: 0.0,
    };
    let settings = SimulationSettings::new(160.0, 0.25)?;

    let result = simulate_seir(parameters, initial_state, settings)?;
    let metrics = result.metrics;
    let diagnostics = result.diagnostics;

    println!("SEIR model");
    println!(
        "  beta = {} /day, sigma = {} /day, gamma = {} /day",
        parameters.beta_per_day, parameters.sigma_per_day, parameters.gamma_per_day
    );
    println!(
        "  initial S = {}, E = {}, I = {}, R = {}",
        initial_state.susceptible,
        initial_state.exposed,
        initial_state.infected,
        initial_state.recovered
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
        "  peak infectious = {:.1} on day {}",
        metrics.peak_infectious, metrics.peak_infectious_time_days
    );
    println!(
        "  cumulative infections since start = {:.1} ({:.1}% of the population)",
        metrics.cumulative_infections_since_start,
        100.0 * metrics.cumulative_incidence_fraction
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
