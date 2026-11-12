//! A lightweight epidemiological simulation engine.
//!
//! The crate integrates two deterministic compartmental models with a
//! fixed-step classical fourth-order Runge-Kutta solver:
//!
//! - [`simulate_sir`]: susceptible, infected, recovered.
//! - [`simulate_seir`]: susceptible, exposed, infected, recovered.
//!
//! Both take rate parameters in inverse days, an initial state in population
//! counts, and [`SimulationSettings`] in days. They return a
//! [`SimulationResult`] holding the full trajectory, [`EpidemicMetrics`], and
//! [`SimulationDiagnostics`], or a [`SimulationError`].
//!
//! The models assume a closed, homogeneously mixing population, constant
//! parameters, permanent immunity, and frequency-dependent transmission. The
//! results describe the model under those assumptions; they are not
//! forecasts or estimates of a real outbreak.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![deny(missing_debug_implementations)]

mod error;
mod metrics;
mod models;
mod result;
mod settings;
mod simulation;
mod solver;
mod validation;

pub use error::SimulationError;
pub use metrics::EpidemicMetrics;
pub use models::seir::{SeirParameters, SeirState, simulate_seir};
pub use models::sir::{SirParameters, SirState, simulate_sir};
pub use result::{Observation, SimulationDiagnostics, SimulationResult};
pub use settings::SimulationSettings;

// Runs the Rust examples in the Markdown documentation as doctests.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
#[doc = include_str!("../docs/api_design.md")]
struct MarkdownExamples;
