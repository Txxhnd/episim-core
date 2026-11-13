use crate::metrics::EpidemicMetrics;

/// Model state recorded at one simulation time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation<S> {
    /// Simulation time in days since the start of the simulation.
    pub time_days: f64,
    /// New infections between day zero and this time, in individuals; 0 at
    /// day zero. It is integrated with the compartments as
    /// `dC/dt = beta * S * I / N`, so in a model where people can be infected
    /// again it counts every infection, not every person.
    pub cumulative_infections: f64,
    /// Compartment values at that time.
    pub state: S,
}

/// Numerical health indicators collected over a whole simulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SimulationDiagnostics {
    /// Number of integration steps executed.
    pub steps_executed: usize,
    /// Largest absolute difference, in individuals, between the total
    /// population after any accepted step and the initial population.
    pub maximum_absolute_population_drift: f64,
    /// Largest absolute population drift divided by the initial population.
    pub maximum_relative_population_drift: f64,
    /// Number of compartment values, over all accepted steps whether recorded
    /// or not, that were negative but within the roundoff tolerance
    /// `64 * f64::EPSILON * max(N0, 1)`. Such values are reported unchanged
    /// in the trajectory, never clipped to zero.
    pub tolerated_negative_value_count: usize,
}

/// Complete outcome of a simulation.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulationResult<S> {
    /// The initial observation at day zero followed by one observation per
    /// recorded integration step: every step by default, or every N-th step
    /// with [`SimulationSettings::with_record_every`](crate::SimulationSettings::with_record_every).
    /// The last one is always the final state, at the requested duration.
    pub trajectory: Vec<Observation<S>>,
    /// Epidemiological metrics derived from the trajectory.
    pub metrics: EpidemicMetrics,
    /// Numerical diagnostics collected during integration.
    pub diagnostics: SimulationDiagnostics,
}
