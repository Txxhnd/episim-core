//! A lightweight epidemiological simulation engine.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![deny(missing_debug_implementations)]

mod error;
// The validators have no caller until the simulation functions exist;
// `expect` starts warning once they do, so these attributes cannot outlive their reason.
#[cfg_attr(not(test), expect(dead_code))]
mod models;
mod settings;
#[cfg_attr(not(test), expect(dead_code))]
mod solver;
#[cfg_attr(not(test), expect(dead_code))]
mod validation;

pub use error::SimulationError;
pub use models::seir::{SeirParameters, SeirState};
pub use models::sir::{SirParameters, SirState};
pub use settings::SimulationSettings;
