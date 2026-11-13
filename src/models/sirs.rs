use crate::accuracy::{StepErrorEstimate, estimate_step_error};
use crate::error::SimulationError;
use crate::metrics::{
    basic_reproduction_number, calculate_metrics, effective_reproduction_number,
    herd_immunity_threshold,
};
use crate::models::sir::SirState;
use crate::result::SimulationResult;
use crate::settings::SimulationSettings;
use crate::simulation::driver::run_simulation;
use crate::validation;

/// Rate parameters of the SIRS model, constant for the whole simulation.
///
/// The model uses [`SirState`] for its compartments. With
/// `omega_per_day = 0` it is the SIR model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SirsParameters {
    /// Transmission rate in inverse days; must be finite and non-negative.
    pub beta_per_day: f64,
    /// Recovery rate in inverse days; must be finite and strictly positive.
    pub gamma_per_day: f64,
    /// Rate in inverse days at which recovered individuals lose their
    /// immunity and become susceptible again; must be finite and
    /// non-negative. `1 / omega_per_day` is the mean duration of immunity.
    pub omega_per_day: f64,
}

impl SirsParameters {
    pub(crate) fn validate(&self) -> Result<(), SimulationError> {
        validation::non_negative_rate("beta_per_day", self.beta_per_day)?;
        validation::positive_rate("gamma_per_day", self.gamma_per_day)?;
        validation::non_negative_rate("omega_per_day", self.omega_per_day)
    }

    /// Model-derived effective reproduction number in `state`: `R0 * S / N`,
    /// where `N` is the total of `state`.
    ///
    /// Because `dI/dt = gamma * I * (effective R - 1)`, the infected grow while
    /// it is above 1 and decline while it is below 1. As immunity wanes the
    /// susceptibles grow again, so it can rise back above 1; it tends to 1 as
    /// the model settles at its endemic equilibrium.
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
    pub fn effective_reproduction_number(&self, state: &SirState) -> Result<f64, SimulationError> {
        self.validate()?;
        let r0 = basic_reproduction_number(self.beta_per_day, self.gamma_per_day)?;
        effective_reproduction_number(
            r0,
            state.susceptible,
            state.susceptible + state.infected + state.recovered,
        )
    }
}

/// Simulates the deterministic SIRS model with a fixed-step classical RK4 integrator.
///
/// The model is SIR with waning immunity: recovered individuals become
/// susceptible again at rate `omega`. It assumes a closed, homogeneously
/// mixing population with constant parameters. Transmission is
/// frequency-dependent, with `N` fixed at the initial total population:
///
/// ```text
/// dS/dt = -beta * S * I / N + omega * R
/// dI/dt =  beta * S * I / N - gamma * I
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
/// use episim_core::{SimulationSettings, SirState, SirsParameters, simulate_sirs};
///
/// // Immunity lasts 100 days on average.
/// let parameters = SirsParameters { beta_per_day: 0.3, gamma_per_day: 0.1, omega_per_day: 0.01 };
/// let initial_state = SirState { susceptible: 990.0, infected: 10.0, recovered: 0.0 };
/// let result = simulate_sirs(parameters, initial_state, SimulationSettings::new(5000.0, 0.25)?)?;
///
/// // (1 - 1/R0) * omega / (omega + gamma) = (2/3) / 11 of the population stays infected.
/// let endemic = result.metrics.endemic_prevalence.unwrap();
/// assert!((endemic - 2.0 / 33.0).abs() < 1e-12);
/// let last = result.trajectory.last().unwrap();
/// assert!((last.state.infected / 1000.0 - endemic).abs() < 1e-6);
/// # Ok::<(), episim_core::SimulationError>(())
/// ```
pub fn simulate_sirs(
    parameters: SirsParameters,
    initial_state: SirState,
    settings: SimulationSettings,
) -> Result<SimulationResult<SirState>, SimulationError> {
    parameters.validate()?;
    let population = initial_state.initial_population()?;
    let SirsParameters {
        beta_per_day,
        gamma_per_day,
        omega_per_day,
    } = parameters;
    // With waning immunity the largest single rate is not enough: beta + gamma
    // is the decay rate of a small departure from the state with everyone
    // infected, and recovered individuals leave R at omega. Measured; see
    // `tests/stability_evidence_tests.rs`.
    validation::stable_step(
        settings.step_days().min(settings.duration_days()),
        (beta_per_day + gamma_per_day).max(omega_per_day),
    )?;
    // Checked now so that an unusable R0 is reported before any integration.
    basic_reproduction_number(beta_per_day, gamma_per_day)?;

    let output = run_simulation(
        [
            initial_state.susceptible,
            initial_state.infected,
            initial_state.recovered,
            0.0,
        ],
        population,
        settings,
        |[susceptible, infected, recovered, _]| {
            // Dividing first avoids multiplying two large population values.
            let infection_flow = beta_per_day * (susceptible / population) * infected;
            let recovery_flow = gamma_per_day * infected;
            let waning_flow = omega_per_day * recovered;
            [
                waning_flow - infection_flow,
                infection_flow - recovery_flow,
                recovery_flow - waning_flow,
                infection_flow,
            ]
        },
        |[susceptible, infected, recovered, _]| SirState {
            susceptible,
            infected,
            recovered,
        },
        &["susceptible", "infected", "recovered"],
    )?;
    let mut metrics = calculate_metrics(
        &output.trajectory,
        population,
        beta_per_day,
        gamma_per_day,
        |state| state.infected,
    )?;
    let r0 = metrics.basic_reproduction_number;
    metrics.herd_immunity_threshold = Some(herd_immunity_threshold(r0));
    // At the endemic equilibrium S* = N / R0, and gamma * I* = omega * R*
    // with R* = N - S* - I*, which gives
    // I* / N = (1 - 1/R0) * omega / (omega + gamma). Without waning
    // (omega = 0) the infection still dies out, so the prevalence is 0.
    metrics.endemic_prevalence =
        Some(herd_immunity_threshold(r0) * omega_per_day / (omega_per_day + gamma_per_day));
    Ok(SimulationResult {
        trajectory: output.trajectory,
        metrics,
        diagnostics: output.diagnostics,
    })
}

/// Estimates the error that the timestep introduces into [`simulate_sirs`]
/// for the same arguments.
///
/// The simulation is run twice, once as given and once at half the timestep,
/// and the two trajectories, cumulative infections included, are compared at
/// every recorded time. This costs about three times as much as one run. The
/// half-step run is allowed twice the configured step limit.
///
/// # Errors
///
/// Returns any [`SimulationError`] that [`simulate_sirs`] returns for either run.
pub fn estimate_sirs_step_error(
    parameters: SirsParameters,
    initial_state: SirState,
    settings: SimulationSettings,
) -> Result<StepErrorEstimate, SimulationError> {
    let run = simulate_sirs(parameters, initial_state, settings)?;
    let half_step_run = simulate_sirs(parameters, initial_state, settings.halved())?;
    Ok(estimate_step_error(
        &run.trajectory,
        &half_step_run.trajectory,
        initial_state.initial_population()?,
        |state| [state.susceptible, state.infected, state.recovered],
    ))
}
