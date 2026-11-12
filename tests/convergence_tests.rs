//! Timestep-refinement checks for the nonlinear SIR and SEIR systems.
//!
//! Neither system has a closed-form solution in time once transmission is
//! active, so the reference is a run with a much finer step. It is a
//! numerical reference, not an exact solution: its own error is about
//! 16^4 = 65536 times smaller than that of the finest step under test.

use episim_core::{
    SeirParameters, SeirState, SimulationSettings, SirParameters, SirState, simulate_seir,
    simulate_sir,
};

/// Every step is a power of two, so the time grids coincide exactly.
const STEPS: [f64; 3] = [1.0, 0.5, 0.25];
const REFERENCE_STEP: f64 = 1.0 / 64.0;

/// Error-reduction factors for halving the step twice, measured as the
/// largest compartment difference from the reference over the whole trajectory.
fn refinement_ratios(trajectory_for: impl Fn(SimulationSettings) -> Vec<[f64; 4]>) -> [f64; 2] {
    let run = |step| trajectory_for(SimulationSettings::new(40.0, step).unwrap());
    let reference = run(REFERENCE_STEP);
    let errors = STEPS.map(|step| {
        let stride = (step / REFERENCE_STEP) as usize;
        run(step)
            .iter()
            .enumerate()
            .flat_map(|(index, state)| {
                let expected = reference[index * stride];
                (0..4).map(move |compartment| (state[compartment] - expected[compartment]).abs())
            })
            .fold(0.0, f64::max)
    });
    assert!(
        errors[2] > 1e-9,
        "finest error {} is too close to the floating-point floor to measure a ratio",
        errors[2]
    );
    [errors[0] / errors[1], errors[1] / errors[2]]
}

#[test]
fn sir_timestep_refinement_reduces_reference_error() {
    let parameters = SirParameters {
        beta_per_day: 0.5,
        gamma_per_day: 0.1,
    };
    let initial = SirState {
        susceptible: 990.0,
        infected: 10.0,
        recovered: 0.0,
    };
    let ratios = refinement_ratios(|settings| {
        simulate_sir(parameters, initial, settings)
            .unwrap()
            .trajectory
            .iter()
            .map(|o| {
                [
                    o.state.susceptible,
                    o.state.infected,
                    o.state.recovered,
                    0.0,
                ]
            })
            .collect()
    });
    for ratio in ratios {
        assert!(ratio >= 8.0, "error-reduction ratios {ratios:?}");
    }
}

#[test]
fn seir_timestep_refinement_reduces_reference_error() {
    let parameters = SeirParameters {
        beta_per_day: 0.5,
        sigma_per_day: 0.2,
        gamma_per_day: 0.1,
    };
    let initial = SeirState {
        susceptible: 980.0,
        exposed: 10.0,
        infected: 10.0,
        recovered: 0.0,
    };
    let ratios = refinement_ratios(|settings| {
        simulate_seir(parameters, initial, settings)
            .unwrap()
            .trajectory
            .iter()
            .map(|o| {
                [
                    o.state.susceptible,
                    o.state.exposed,
                    o.state.infected,
                    o.state.recovered,
                ]
            })
            .collect()
    });
    for ratio in ratios {
        assert!(ratio >= 8.0, "error-reduction ratios {ratios:?}");
    }
}
