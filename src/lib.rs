//! A lightweight epidemiological simulation engine.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![deny(missing_debug_implementations)]

mod error;
mod metrics;
// The SEIR validators have no caller until `simulate_seir` exists; `expect`
// starts warning once they do, so this attribute cannot outlive its reason.
#[cfg_attr(not(test), expect(dead_code))]
mod models;
mod result;
mod settings;
mod simulation;
mod solver;
mod validation;

pub use error::SimulationError;
pub use metrics::EpidemicMetrics;
pub use models::seir::{SeirParameters, SeirState};
pub use models::sir::{SirParameters, SirState, simulate_sir};
pub use result::{Observation, SimulationDiagnostics, SimulationResult};
pub use settings::SimulationSettings;
