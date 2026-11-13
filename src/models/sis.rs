use crate::accuracy::{StepErrorEstimate, estimate_step_error};
use crate::error::SimulationError;
use crate::metrics::{
    basic_reproduction_number, calculate_metrics, effective_reproduction_number,
    herd_immunity_threshold,
};
use crate::result::SimulationResult;
use crate::settings::SimulationSettings;
use crate::simulation::driver::run_simulation;
use crate::validation;

/// Rate parameters of the SIS model, constant for the whole simulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SisParameters {
    /// Transmission rate in inverse days; must be finite and non-negative.
    pub beta_per_day: f64,
    /// Recovery rate in inverse days, at which the infected become
    /// susceptible again; must be finite and strictly positive.
    pub gamma_per_day: f64,
}

/// Population counts in each SIS compartment.
///
/// Counts are expected values and may be fractional. Every compartment must
/// be finite and non-negative, and the total must be finite and strictly
/// positive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SisState {
    /// Individuals who can become infected, including those who have
    /// recovered: the SIS model has no immunity.
    pub susceptible: f64,
    /// Individuals who are currently infectious.
    pub infected: f64,
}

impl SisParameters {
    pub(crate) fn validate(&self) -> Result<(), SimulationError> {
        validation::non_negative_rate("beta_per_day", self.beta_per_day)?;
        validation::positive_rate("gamma_per_day", self.gamma_per_day)
    }

    /// Model-derived effective reproduction number in `state`: `R0 * S / N`,
    /// where `N` is the total of `state`.
    ///
    /// Because `dI/dt = gamma * I * (effective R - 1)`, the infected grow while
    /// it is above 1 and decline while it is below 1. It tends to 1 as the
    /// model settles at its endemic equilibrium.
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
    pub fn effective_reproduction_number(&self, state: &SisState) -> Result<f64, SimulationError> {
        self.validate()?;
        let r0 = basic_reproduction_number(self.beta_per_day, self.gamma_per_day)?;
        effective_reproduction_number(r0, state.susceptible, state.susceptible + state.infected)
    }
}

impl SisState {
    /// Validates the state and returns the total population it defines.
    pub(crate) fn initial_population(&self) -> Result<f64, SimulationError> {
        validation::initial_population(&[
            ("susceptible", self.susceptible),
            ("infected", self.infected),
        ])
    }
}

/// Simulates the deterministic SIS model with a fixed-step classical RK4 integrator.
///
/// The model assumes a closed, homogeneously mixing population with constant
/// parameters, in which recovery gives no immunity. Transmission is
/// frequency-dependent, with `N` fixed at the initial total population:
///
/// ```text
/// dS/dt = -beta * S * I / N + gamma * I
/// dI/dt =  beta * S * I / N - gamma * I
/// ```
///
/// Instead of ending, an epidemic with `R0 > 1` settles at an endemic
/// equilibrium with a fraction `1 - 1/R0` of the population infected, reported
/// as `endemic_prevalence`. People can be infected many times, so
/// `cumulative_infections_since_start` counts infections, not people, and
/// `cumulative_incidence_fraction` can exceed 1. There is no herd-immunity
/// threshold and no final size: both metrics are `None`.
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
/// use episim_core::{SimulationSettings, SisParameters, SisState, simulate_sis};
///
/// let parameters = SisParameters { beta_per_day: 0.3, gamma_per_day: 0.1 };
/// let initial_state = SisState { susceptible: 990.0, infected: 10.0 };
/// let result = simulate_sis(parameters, initial_state, SimulationSettings::new(200.0, 0.25)?)?;
///
/// // R0 is 3, so two thirds of the population end up infected at any time.
/// let endemic = result.metrics.endemic_prevalence.unwrap();
/// assert!((endemic - 2.0 / 3.0).abs() < 1e-12);
/// let last = result.trajectory.last().unwrap();
/// assert!((last.state.infected / 1000.0 - endemic).abs() < 1e-6);
/// // Each person has been infected several times on average.
/// assert!(result.metrics.cumulative_incidence_fraction > 1.0);
/// # Ok::<(), episim_core::SimulationError>(())
/// ```
pub fn simulate_sis(
    parameters: SisParameters,
    initial_state: SisState,
    settings: SimulationSettings,
) -> Result<SimulationResult<SisState>, SimulationError> {
    parameters.validate()?;
    let population = initial_state.initial_population()?;
    let SisParameters {
        beta_per_day,
        gamma_per_day,
    } = parameters;
    // beta + gamma is the decay rate of a small departure from the state with
    // everyone infected. The 1.25 * gamma floor keeps gamma * step at most
    // 1.6: with little transmission and many infected, a larger step makes
    // the RK4 increment of the cumulative infections negative. Both were
    // measured; see `tests/stability_evidence_tests.rs`.
    validation::stable_step(
        settings.step_days().min(settings.duration_days()),
        (beta_per_day + gamma_per_day).max(1.25 * gamma_per_day),
    )?;
    // Checked now so that an unusable R0 is reported before any integration.
    basic_reproduction_number(beta_per_day, gamma_per_day)?;

    let output = run_simulation(
        [initial_state.susceptible, initial_state.infected, 0.0],
        population,
        settings,
        |[susceptible, infected, _]| {
            // Dividing first avoids multiplying two large population values.
            let infection_flow = beta_per_day * (susceptible / population) * infected;
            let recovery_flow = gamma_per_day * infected;
            [
                recovery_flow - infection_flow,
                infection_flow - recovery_flow,
                infection_flow,
            ]
        },
        |[susceptible, infected, _]| SisState {
            susceptible,
            infected,
        },
        &["susceptible", "infected"],
    )?;
    let mut metrics = calculate_metrics(
        &output.trajectory,
        population,
        beta_per_day,
        gamma_per_day,
        |state| state.infected,
    )?;
    // At the endemic equilibrium S* = N / R0, so I* / N = 1 - 1/R0.
    metrics.endemic_prevalence = Some(herd_immunity_threshold(metrics.basic_reproduction_number));
    Ok(SimulationResult {
        trajectory: output.trajectory,
        metrics,
        diagnostics: output.diagnostics,
    })
}

/// Estimates the error that the timestep introduces into [`simulate_sis`] for
/// the same arguments.
///
/// The simulation is run twice, once as given and once at half the timestep,
/// and the two trajectories, cumulative infections included, are compared at
/// every recorded time. This costs about three times as much as one run. The
/// half-step run is allowed twice the configured step limit.
///
/// # Errors
///
/// Returns any [`SimulationError`] that [`simulate_sis`] returns for either run.
pub fn estimate_sis_step_error(
    parameters: SisParameters,
    initial_state: SisState,
    settings: SimulationSettings,
) -> Result<StepErrorEstimate, SimulationError> {
    let run = simulate_sis(parameters, initial_state, settings)?;
    let half_step_run = simulate_sis(parameters, initial_state, settings.halved())?;
    Ok(estimate_step_error(
        &run.trajectory,
        &half_step_run.trajectory,
        initial_state.initial_population()?,
        |state| [state.susceptible, state.infected],
    ))
}
