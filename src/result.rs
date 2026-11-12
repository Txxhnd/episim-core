/// Model state recorded at one simulation time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation<S> {
    /// Simulation time in days since the start of the simulation.
    pub time_days: f64,
    /// Compartment values at that time.
    pub state: S,
}

/// Numerical health indicators collected over a whole simulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SimulationDiagnostics {
    /// Number of integration steps executed.
    pub steps_executed: usize,
    /// Largest absolute difference, in individuals, between the total
    /// population of a recorded state and the initial population.
    pub maximum_absolute_population_drift: f64,
    /// Largest absolute population drift divided by the initial population.
    pub maximum_relative_population_drift: f64,
    /// Number of recorded compartment values that were negative but within
    /// the roundoff tolerance `64 * f64::EPSILON * max(N0, 1)`. Such values
    /// are reported unchanged in the trajectory, never clipped to zero.
    pub tolerated_negative_value_count: usize,
}
