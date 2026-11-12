use crate::error::SimulationError;

const DEFAULT_MAX_STEPS: usize = 1_000_000;

/// Validated simulation duration, fixed timestep, integration step limit,
/// and recording interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SimulationSettings {
    duration_days: f64,
    step_days: f64,
    max_steps: usize,
    record_every: usize,
}

impl SimulationSettings {
    /// Creates settings with a duration and fixed timestep, both in days.
    ///
    /// The maximum number of integration steps defaults to 1,000,000, and
    /// every step is recorded.
    ///
    /// # Errors
    ///
    /// Returns [`SimulationError::InvalidDuration`] or
    /// [`SimulationError::InvalidTimeStep`] if the corresponding value is not
    /// finite and strictly positive.
    pub fn new(duration_days: f64, step_days: f64) -> Result<Self, SimulationError> {
        if !(duration_days.is_finite() && duration_days > 0.0) {
            return Err(SimulationError::InvalidDuration);
        }
        if !(step_days.is_finite() && step_days > 0.0) {
            return Err(SimulationError::InvalidTimeStep);
        }
        Ok(Self {
            duration_days,
            step_days,
            max_steps: DEFAULT_MAX_STEPS,
            record_every: 1,
        })
    }

    /// Replaces the maximum number of integration steps.
    ///
    /// # Errors
    ///
    /// Returns [`SimulationError::InvalidMaximumSteps`] if `max_steps` is zero.
    pub fn with_max_steps(self, max_steps: usize) -> Result<Self, SimulationError> {
        if max_steps == 0 {
            return Err(SimulationError::InvalidMaximumSteps);
        }
        Ok(Self { max_steps, ..self })
    }

    /// Records only every `steps`-th integration step in the trajectory.
    ///
    /// The initial state and the final state are always recorded. Every step
    /// is still integrated and checked; this only reduces how many
    /// observations are stored. Peak metrics are read from the recorded
    /// observations, so a larger interval makes them coarser.
    ///
    /// # Errors
    ///
    /// Returns [`SimulationError::InvalidRecordInterval`] if `steps` is zero.
    pub fn with_record_every(self, steps: usize) -> Result<Self, SimulationError> {
        if steps == 0 {
            return Err(SimulationError::InvalidRecordInterval);
        }
        Ok(Self {
            record_every: steps,
            ..self
        })
    }

    /// The same run at half the timestep, recording the same times.
    ///
    /// Halving the step doubles the number of steps, so the step limit and
    /// the recording interval are doubled with it.
    pub(crate) fn halved(&self) -> Self {
        Self {
            duration_days: self.duration_days,
            step_days: self.step_days / 2.0,
            max_steps: self.max_steps.saturating_mul(2),
            record_every: self.record_every.saturating_mul(2),
        }
    }

    /// Total simulation duration in days.
    pub fn duration_days(&self) -> f64 {
        self.duration_days
    }

    /// Fixed integration timestep in days.
    pub fn step_days(&self) -> f64 {
        self.step_days
    }

    /// Maximum number of integration steps a simulation may execute.
    pub fn max_steps(&self) -> usize {
        self.max_steps
    }

    /// Number of integration steps between recorded observations.
    pub fn record_every(&self) -> usize {
        self.record_every
    }
}
