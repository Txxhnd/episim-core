use crate::error::SimulationError;

/// A piecewise-constant multiplier on the transmission rate `beta`, for
/// interventions such as lockdowns and school closures, and for their end.
///
/// The multiplier is 1 until the first change. Each change sets a new
/// multiplier from its day on, so the transmission rate in force on day `t`
/// is `beta * multiplier_at(t)`. A multiplier below 1 models fewer contacts,
/// above 1 more. Attach a schedule to a run with
/// [`SimulationSettings::with_transmission_schedule`](crate::SimulationSettings::with_transmission_schedule).
///
/// # Examples
///
/// ```
/// use episim_core::TransmissionSchedule;
///
/// // Contacts at 40% from day 20, back to normal from day 50.
/// let lockdown = TransmissionSchedule::new(&[(20.0, 0.4), (50.0, 1.0)])?;
/// assert_eq!(lockdown.multiplier_at(10.0), 1.0);
/// assert_eq!(lockdown.multiplier_at(20.0), 0.4);
/// assert_eq!(lockdown.multiplier_at(60.0), 1.0);
/// # Ok::<(), episim_core::SimulationError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct TransmissionSchedule {
    changes: Vec<(f64, f64)>,
}

impl TransmissionSchedule {
    /// Creates a schedule from `(day, multiplier)` changes.
    ///
    /// # Errors
    ///
    /// Returns [`SimulationError::InvalidSchedule`] unless every day is
    /// finite and greater than zero, the days strictly increase, and every
    /// multiplier is finite and non-negative.
    pub fn new(changes: &[(f64, f64)]) -> Result<Self, SimulationError> {
        let mut previous_day = 0.0;
        for &(day, multiplier) in changes {
            if !(day.is_finite()
                && day > previous_day
                && multiplier.is_finite()
                && multiplier >= 0.0)
            {
                return Err(SimulationError::InvalidSchedule);
            }
            previous_day = day;
        }
        Ok(Self {
            changes: changes.to_vec(),
        })
    }

    /// The multiplier in force on `day`: that of the latest change on or
    /// before `day`, or 1 before the first change. A change counts as on or
    /// before `day` within four units in the last place of the change day,
    /// the tolerance `with_transmission_schedule` allows for the step grid,
    /// so a recorded time such as `3 * 0.3 = 0.8999999999999999` reports the
    /// multiplier the run applied from day 0.9 on.
    ///
    /// Multiply the model-derived effective reproduction number of a state
    /// (for example
    /// [`SirParameters::effective_reproduction_number`](crate::SirParameters::effective_reproduction_number))
    /// by this value to get the effective reproduction number on that day.
    pub fn multiplier_at(&self, day: f64) -> f64 {
        let in_force = self.changes.partition_point(|&(change_day, _)| {
            change_day - day <= 4.0 * f64::EPSILON * change_day
        });
        in_force
            .checked_sub(1)
            .map_or(1.0, |last| self.changes[last].1)
    }

    /// The `(day, multiplier)` changes, in order.
    pub(crate) fn changes(&self) -> &[(f64, f64)] {
        &self.changes
    }
}
