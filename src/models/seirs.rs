use crate::accuracy::StepErrorEstimate;
use crate::error::SimulationError;
use crate::metrics::{
    EpidemicMetrics, basic_reproduction_number, effective_reproduction_number,
    herd_immunity_threshold,
};
use crate::models::seir::SeirState;
use crate::models::{Model, estimate, simulate};
use crate::result::SimulationResult;
use crate::settings::SimulationSettings;
use crate::validation;

/// Rate parameters of the SEIRS model, constant for the whole simulation.
///
/// The model uses [`SeirState`] for its compartments. With
/// `omega_per_day = 0` it is the SEIR model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeirsParameters {
    /// Transmission rate in inverse days; must be finite and non-negative.
    pub beta_per_day: f64,
    /// Rate of progression from exposed to infectious in inverse days; must
    /// be finite and strictly positive.
    pub sigma_per_day: f64,
    /// Recovery rate in inverse days; must be finite and strictly positive.
    pub gamma_per_day: f64,
    /// Rate in inverse days at which recovered individuals lose their
    /// immunity and become susceptible again; must be finite and
    /// non-negative. `1 / omega_per_day` is the mean duration of immunity.
    pub omega_per_day: f64,
}

impl SeirsParameters {
    pub(crate) fn validate(&self) -> Result<(), SimulationError> {
        validation::non_negative_rate("beta_per_day", self.beta_per_day)?;
        validation::positive_rate("sigma_per_day", self.sigma_per_day)?;
        validation::positive_rate("gamma_per_day", self.gamma_per_day)?;
        validation::non_negative_rate("omega_per_day", self.omega_per_day)
    }

    /// Model-derived effective reproduction number in `state`: `R0 * S / N`,
    /// where `N` is the total of `state`.
    ///
    /// Because `d(E + I)/dt = gamma * I * (effective R - 1)`, the exposed and
    /// infected together grow while it is above 1 and decline while it is
    /// below 1. As immunity wanes the susceptibles grow again, so it can rise
    /// back above 1; it tends to 1 as the model settles at its endemic
    /// equilibrium.
    ///
    /// This is a property of the model. It is not an estimate of the effective
    /// reproduction number `Rt` that public-health agencies infer from case
    /// data.
    ///
    /// # Errors
    ///
    /// Returns [`SimulationError::InvalidParameter`] for an invalid rate,
    /// [`SimulationError::InvalidPopulation`] if the total of `state` is not
    /// finite and positive, and [`SimulationError::MetricOverflow`] if the
    /// result is not finite.
    pub fn effective_reproduction_number(&self, state: &SeirState) -> Result<f64, SimulationError> {
        self.validate()?;
        let r0 = basic_reproduction_number(self.beta_per_day, self.gamma_per_day)?;
        effective_reproduction_number(
            r0,
            state.susceptible,
            state.susceptible + state.exposed + state.infected + state.recovered,
        )
    }
}

/// Simulates the deterministic SEIRS model with a fixed-step classical RK4 integrator.
///
/// The model is SEIR with waning immunity: recovered individuals become
/// susceptible again at rate `omega`. It assumes a closed, homogeneously
/// mixing population with constant parameters. Transmission is
/// frequency-dependent, with `N` fixed at the initial total population:
///
/// ```text
/// dS/dt = -beta * S * I / N + omega * R
/// dE/dt =  beta * S * I / N - sigma * E
/// dI/dt =  sigma * E - gamma * I
/// dR/dt =  gamma * I - omega * R
/// ```
///
/// With `R0 > 1` and `omega > 0` the epidemic does not end: after damped
/// waves it settles at an endemic equilibrium, whose infected fraction is
/// reported as `endemic_prevalence`. People can be infected more than once,
/// so `cumulative_infections_since_start` counts infections, not people.
/// There is no final size, so `final_size_fraction` is `None`.
///
/// # Errors
///
/// Returns a [`SimulationError`] if a parameter or the initial state is
/// invalid, if the timestep is too large for stable integration with these
/// rates, if the settings require more steps than their limit or a
/// timestep too small to advance the clock, if the integration produces a
/// non-finite value, a materially negative compartment, or excessive
/// population drift, if a metric is not finite, or if the trajectory cannot
/// be allocated.
///
/// # Examples
///
/// ```
/// use episim_core::{SeirState, SeirsParameters, SimulationSettings, simulate_seirs};
///
/// // Immunity lasts 100 days on average.
/// let parameters = SeirsParameters {
///     beta_per_day: 0.5,
///     sigma_per_day: 0.2,
///     gamma_per_day: 0.1,
///     omega_per_day: 0.01,
/// };
/// let initial_state = SeirState { susceptible: 980.0, exposed: 10.0, infected: 10.0, recovered: 0.0 };
/// let result = simulate_seirs(parameters, initial_state, SimulationSettings::new(5000.0, 0.25)?)?;
///
/// // (1 - 1/R0) / (1 + gamma/sigma + gamma/omega) = 0.8 / 11.5 of the population.
/// let endemic = result.metrics.endemic_prevalence.unwrap();
/// assert!((endemic - 0.8 / 11.5).abs() < 1e-12);
/// let last = result.trajectory.last().unwrap();
/// assert!((last.state.infected / 1000.0 - endemic).abs() < 1e-6);
/// # Ok::<(), episim_core::SimulationError>(())
/// ```
pub fn simulate_seirs(
    parameters: SeirsParameters,
    initial_state: SeirState,
    settings: SimulationSettings,
) -> Result<SimulationResult<SeirState>, SimulationError> {
    simulate(&parameters, initial_state, settings)
}

/// Estimates the error that the timestep introduces into [`simulate_seirs`]
/// for the same arguments.
///
/// The simulation is run twice, once as given and once at half the timestep,
/// and the two trajectories, cumulative infections included, are compared at
/// every recorded time. This costs about three times as much as one run. The
/// half-step run is allowed twice the configured step limit.
///
/// # Errors
///
/// Returns any [`SimulationError`] that [`simulate_seirs`] returns for either run.
pub fn estimate_seirs_step_error(
    parameters: SeirsParameters,
    initial_state: SeirState,
    settings: SimulationSettings,
) -> Result<StepErrorEstimate, SimulationError> {
    estimate(&parameters, initial_state, settings)
}

impl Model<5> for SeirsParameters {
    type State = SeirState;
    const COMPARTMENTS: &'static [&'static str] =
        &["susceptible", "exposed", "infected", "recovered"];

    fn validate(&self) -> Result<(), SimulationError> {
        SeirsParameters::validate(self)
    }

    fn population(state: &SeirState) -> Result<f64, SimulationError> {
        state.initial_population()
    }

    fn beta(&self) -> f64 {
        self.beta_per_day
    }

    fn gamma(&self) -> f64 {
        self.gamma_per_day
    }

    fn fastest_rate(&self) -> f64 {
        // The sum of all four rates. Tighter formulas built from SEIR's coupled
        // rate let through runs whose recovered compartment rises over a step
        // where it should fall. Measured; see `tests/stability_evidence_tests.rs`.
        self.beta_per_day + self.sigma_per_day + self.gamma_per_day + self.omega_per_day
    }

    fn to_array(state: &SeirState) -> [f64; 5] {
        [
            state.susceptible,
            state.exposed,
            state.infected,
            state.recovered,
            0.0,
        ]
    }

    fn from_array([susceptible, exposed, infected, recovered, _]: [f64; 5]) -> SeirState {
        SeirState {
            susceptible,
            exposed,
            infected,
            recovered,
        }
    }

    fn derivative(
        &self,
        population: f64,
        [susceptible, exposed, infected, recovered, _]: [f64; 5],
    ) -> [f64; 5] {
        // Dividing first avoids multiplying two large population values.
        let infection_flow = self.beta_per_day * (susceptible / population) * infected;
        let progression_flow = self.sigma_per_day * exposed;
        let recovery_flow = self.gamma_per_day * infected;
        let waning_flow = self.omega_per_day * recovered;
        [
            waning_flow - infection_flow,
            infection_flow - progression_flow,
            progression_flow - recovery_flow,
            recovery_flow - waning_flow,
            infection_flow,
        ]
    }

    fn infectious(state: &SeirState) -> f64 {
        state.infected
    }

    fn model_metrics(&self, metrics: &mut EpidemicMetrics, initial: &SeirState, _population: f64) {
        let (sigma_per_day, gamma_per_day, omega_per_day) =
            (self.sigma_per_day, self.gamma_per_day, self.omega_per_day);
        let r0 = metrics.basic_reproduction_number;
        metrics.herd_immunity_threshold = Some(herd_immunity_threshold(r0));
        // At the endemic equilibrium S* = N / R0, sigma * E* = gamma * I* and
        // gamma * I* = omega * R*, so the N (1 - 1/R0) not susceptible split as
        // E* : I* : R* = gamma/sigma : 1 : gamma/omega. Written with omega in the
        // numerator, so that omega = 0 gives 0 without dividing by zero. With
        // nobody exposed or infected at the start the infection can never take hold.
        metrics.endemic_prevalence = Some(if initial.exposed + initial.infected > 0.0 {
            herd_immunity_threshold(r0) * omega_per_day
                / (omega_per_day + gamma_per_day + gamma_per_day * omega_per_day / sigma_per_day)
        } else {
            0.0
        });
    }
}
