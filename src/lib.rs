//! A lightweight epidemiological simulation engine.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![deny(missing_debug_implementations)]

mod error;
mod settings;

pub use error::SimulationError;
pub use settings::SimulationSettings;
