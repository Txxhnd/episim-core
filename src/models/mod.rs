use crate::accuracy::{StepErrorEstimate, estimate_step_error};
use crate::error::SimulationError;
use crate::metrics::{EpidemicMetrics, basic_reproduction_number, calculate_metrics};
use crate::result::SimulationResult;
use crate::settings::SimulationSettings;
use crate::simulation::driver::run_simulation;
use crate::validation;

pub(crate) mod seir;
pub(crate) mod seirs;
pub(crate) mod sir;
pub(crate) mod sirs;
pub(crate) mod sis;

/// What a compartmental model contributes to the shared pipeline in
/// [`simulate`]. Implemented by each model's parameter type.
///
/// The state array holds the compartments followed by the cumulative
/// infection count, so `N` is one more than the number of compartments.
/// This trait is private: the public API stays one concrete function per model.
pub(crate) trait Model<const N: usize> {
    /// The public state type.
    type State: Copy;
    /// Compartment names in array order, used in errors.
    const COMPARTMENTS: &'static [&'static str];

    /// Rejects invalid parameters.
    fn validate(&self) -> Result<(), SimulationError>;
    /// Validates an initial state and returns the population it defines.
    fn population(state: &Self::State) -> Result<f64, SimulationError>;
    /// Transmission rate in inverse days.
    fn beta(&self) -> f64;
    /// Recovery rate in inverse days.
    fn gamma(&self) -> f64;
    /// Fastest rate at which any state of the model can decay; sets the
    /// stability limit.
    fn fastest_rate(&self) -> f64;
    /// The compartments of `state` followed by a cumulative infection count of 0.
    fn to_array(state: &Self::State) -> [f64; N];
    /// The state held in the compartments of `values`.
    fn from_array(values: [f64; N]) -> Self::State;
    /// Derivative of the compartments and of the cumulative infections.
    fn derivative(&self, population: f64, values: [f64; N]) -> [f64; N];
    /// The infectious compartment, for the peak metrics.
    fn infectious(state: &Self::State) -> f64;
    /// Fills in the metrics that only some models report.
    fn model_metrics(&self, metrics: &mut EpidemicMetrics, initial: &Self::State, population: f64);
}

/// Validates the inputs, applies the stability limit, integrates, and
/// derives the metrics for any model.
pub(crate) fn simulate<const N: usize, M: Model<N>>(
    model: &M,
    initial_state: M::State,
    settings: SimulationSettings,
) -> Result<SimulationResult<M::State>, SimulationError> {
    model.validate()?;
    let population = M::population(&initial_state)?;
    // A duration shorter than the step is integrated in one shorter step.
    validation::stable_step(
        settings.step_days().min(settings.duration_days()),
        model.fastest_rate(),
    )?;
    // Checked now so that an unusable R0 is reported before any integration.
    basic_reproduction_number(model.beta(), model.gamma())?;

    let output = run_simulation(
        M::to_array(&initial_state),
        population,
        settings,
        |values| model.derivative(population, values),
        M::from_array,
        M::COMPARTMENTS,
    )?;
    let mut metrics = calculate_metrics(
        &output.trajectory,
        population,
        model.beta(),
        model.gamma(),
        M::infectious,
    )?;
    model.model_metrics(&mut metrics, &initial_state, population);
    Ok(SimulationResult {
        trajectory: output.trajectory,
        metrics,
        diagnostics: output.diagnostics,
    })
}

/// Runs [`simulate`] at the configured step and at half of it, and compares
/// the two trajectories.
pub(crate) fn estimate<const N: usize, M: Model<N>>(
    model: &M,
    initial_state: M::State,
    settings: SimulationSettings,
) -> Result<StepErrorEstimate, SimulationError> {
    let run = simulate(model, initial_state, settings)?;
    let half_step_run = simulate(model, initial_state, settings.halved())?;
    // The last array entry is the cumulative infection count, 0 in every
    // converted state; the comparison reads the recorded count separately.
    Ok(estimate_step_error(
        &run.trajectory,
        &half_step_run.trajectory,
        M::population(&initial_state)?,
        M::to_array,
    ))
}
