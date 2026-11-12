use std::fmt;

/// Errors returned by simulation configuration and execution.
///
/// New variants may be added in later versions, so a `match` on this type
/// needs a wildcard arm.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SimulationError {
    /// An epidemiological parameter is non-finite or outside its valid range.
    InvalidParameter {
        /// Name of the rejected parameter.
        name: &'static str,
    },
    /// An initial compartment value is non-finite or negative.
    InvalidInitialState {
        /// Name of the rejected compartment.
        compartment: &'static str,
    },
    /// The total initial population is non-finite or not strictly positive.
    InvalidPopulation,
    /// The simulation duration is non-finite or not strictly positive.
    InvalidDuration,
    /// The integration timestep is non-finite or not strictly positive.
    InvalidTimeStep,
    /// The maximum number of integration steps is zero.
    InvalidMaximumSteps,
    /// The simulation requires more integration steps than the configured maximum.
    StepLimitExceeded,
    /// The timestep is too small to advance the floating-point simulation clock.
    TimeStepTooSmall,
    /// The timestep is too large for the solver to integrate stably with one
    /// of the rate parameters: their product exceeds the stability limit of
    /// the classical Runge-Kutta method, about 2.785.
    TimeStepTooLarge {
        /// Name of the rate parameter the timestep is too large for.
        parameter: &'static str,
    },
    /// A solver stage produced a non-finite intermediate value.
    NonFiniteIntermediate {
        /// Solver stage in which the non-finite value appeared.
        stage: &'static str,
    },
    /// An accepted integration step produced a non-finite compartment value.
    NonFiniteNumericalState {
        /// Name of the affected compartment.
        compartment: &'static str,
    },
    /// A compartment became negative beyond the roundoff tolerance.
    NegativeCompartment {
        /// Name of the affected compartment.
        compartment: &'static str,
    },
    /// The total population drifted beyond the conservation tolerance.
    PopulationConservationViolation,
    /// A derived metric is not representable as a finite value.
    MetricOverflow,
    /// Memory for the simulation trajectory could not be reserved.
    AllocationFailure,
}

impl fmt::Display for SimulationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter { name } => write!(f, "invalid parameter `{name}`"),
            Self::InvalidInitialState { compartment } => {
                write!(f, "invalid initial value for compartment `{compartment}`")
            }
            Self::InvalidPopulation => {
                f.write_str("total initial population must be finite and greater than zero")
            }
            Self::InvalidDuration => {
                f.write_str("simulation duration must be finite and greater than zero")
            }
            Self::InvalidTimeStep => f.write_str("timestep must be finite and greater than zero"),
            Self::InvalidMaximumSteps => f.write_str("maximum steps must be greater than zero"),
            Self::StepLimitExceeded => {
                f.write_str("required integration steps exceed the configured maximum")
            }
            Self::TimeStepTooSmall => {
                f.write_str("timestep is too small to advance the simulation clock")
            }
            Self::TimeStepTooLarge { parameter } => write!(
                f,
                "timestep is too large for stable integration with `{parameter}`"
            ),
            Self::NonFiniteIntermediate { stage } => {
                write!(f, "non-finite intermediate value in solver stage `{stage}`")
            }
            Self::NonFiniteNumericalState { compartment } => {
                write!(f, "non-finite value in compartment `{compartment}`")
            }
            Self::NegativeCompartment { compartment } => {
                write!(
                    f,
                    "compartment `{compartment}` is negative beyond roundoff tolerance"
                )
            }
            Self::PopulationConservationViolation => {
                f.write_str("population drift exceeds the conservation tolerance")
            }
            Self::MetricOverflow => f.write_str("derived metric is not finite"),
            Self::AllocationFailure => f.write_str("failed to reserve memory for the trajectory"),
        }
    }
}

impl std::error::Error for SimulationError {}
