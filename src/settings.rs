use crate::error::SimulationError;
use crate::schedule::TransmissionSchedule;

const DEFAULT_MAX_STEPS: usize = 1_000_000;

/// Validated simulation duration, fixed timestep, integration step limit,
/// recording interval, and optional transmission schedule.
///
/// Settings are `Clone` but not `Copy`, because they can hold a schedule.
/// Clone them to reuse them after passing them to a simulation.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulationSettings {
    duration_days: f64,
    step_days: f64,
    max_steps: usize,
    record_every: usize,
    transmission_schedule: Option<TransmissionSchedule>,
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
            transmission_schedule: None,
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

    /// Multiplies the transmission rate by `schedule` over the run.
    ///
    /// Every change must fall on the step grid: its day must be a whole
    /// number of steps, within four units in the last place. Then every step
    /// lies inside one segment of the schedule and keeps the solver's
    /// fourth-order accuracy. Changes on or after the end of the run have no
    /// effect.
    ///
    /// # Errors
    ///
    /// Returns [`SimulationError::ScheduleChangeOffStepGrid`] with the first
    /// day that is not a whole number of steps.
    ///
    /// # Examples
    ///
    /// ```
    /// use episim_core::{SimulationError, SimulationSettings, TransmissionSchedule};
    ///
    /// let lockdown = TransmissionSchedule::new(&[(20.0, 0.4), (50.0, 1.0)])?;
    /// let settings = SimulationSettings::new(160.0, 0.25)?.with_transmission_schedule(lockdown)?;
    /// assert!(settings.transmission_schedule().is_some());
    ///
    /// // A change at day 20.3 does not fall on a grid of quarter days.
    /// let off_grid = TransmissionSchedule::new(&[(20.3, 0.4)])?;
    /// assert_eq!(
    ///     SimulationSettings::new(160.0, 0.25)?.with_transmission_schedule(off_grid),
    ///     Err(SimulationError::ScheduleChangeOffStepGrid { day: 20.3 })
    /// );
    /// # Ok::<(), SimulationError>(())
    /// ```
    pub fn with_transmission_schedule(
        self,
        schedule: TransmissionSchedule,
    ) -> Result<Self, SimulationError> {
        for &(day, _) in schedule.changes() {
            let on_grid = (day / self.step_days).round() * self.step_days;
            if (on_grid - day).abs() > 4.0 * f64::EPSILON * day {
                return Err(SimulationError::ScheduleChangeOffStepGrid { day });
            }
        }
        Ok(Self {
            transmission_schedule: Some(schedule),
            ..self
        })
    }

    /// The transmission schedule, if one was attached.
    pub fn transmission_schedule(&self) -> Option<&TransmissionSchedule> {
        self.transmission_schedule.as_ref()
    }

    /// The schedule's changes that take effect before the end of the run.
    fn changes_in_run(&self) -> impl Iterator<Item = &(f64, f64)> {
        self.transmission_schedule
            .iter()
            .flat_map(TransmissionSchedule::changes)
            .filter(|&&(day, _)| day < self.duration_days)
    }

    /// Multiplier for the step that starts at `start_index` steps. Changes
    /// are matched by step index, not by comparing times, so rounding in
    /// `index * step` cannot move a step into the wrong segment.
    pub(crate) fn multiplier_for_step(&self, start_index: usize) -> f64 {
        self.transmission_schedule
            .iter()
            .flat_map(TransmissionSchedule::changes)
            .rev()
            .find(|&&(day, _)| (day / self.step_days).round() <= start_index as f64)
            .map_or(1.0, |&(_, multiplier)| multiplier)
    }

    /// Largest multiplier in force at any time during the run, at least 1.
    pub(crate) fn largest_multiplier(&self) -> f64 {
        self.changes_in_run()
            .map(|&(_, multiplier)| multiplier)
            .fold(1.0, f64::max)
    }

    /// Multiplier in force at the end of the run.
    pub(crate) fn final_multiplier(&self) -> f64 {
        self.changes_in_run()
            .last()
            .map_or(1.0, |&(_, multiplier)| multiplier)
    }

    /// Whether the schedule sets any multiplier other than 1.
    pub(crate) fn varies_transmission(&self) -> bool {
        self.transmission_schedule
            .iter()
            .flat_map(TransmissionSchedule::changes)
            .any(|&(_, multiplier)| multiplier != 1.0)
    }

    /// The same run at half the timestep, recording the same times.
    ///
    /// Halving the step doubles the number of steps, so the step limit and
    /// the recording interval are doubled with it. A schedule on the original
    /// step grid is also on the halved one.
    pub(crate) fn halved(&self) -> Self {
        Self {
            duration_days: self.duration_days,
            step_days: self.step_days / 2.0,
            max_steps: self.max_steps.saturating_mul(2),
            record_every: self.record_every.saturating_mul(2),
            transmission_schedule: self.transmission_schedule.clone(),
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
