//! A lightweight epidemiological simulation engine.
//!
//! The crate integrates deterministic compartmental models with a fixed-step
//! classical fourth-order Runge-Kutta solver:
//!
//! - [`simulate_sir`]: susceptible, infected, recovered.
//! - [`simulate_seir`]: susceptible, exposed, infected, recovered.
//! - [`simulate_sis`]: susceptible, infected, and susceptible again, with no
//!   immunity.
//! - [`simulate_sirs`]: SIR with immunity that wanes.
//! - [`simulate_seirs`]: SEIR with immunity that wanes.
//!
//! Each takes rate parameters in inverse days, an initial state in population
//! counts, and [`SimulationSettings`] in days. They return a
//! [`SimulationResult`] holding the full trajectory, [`EpidemicMetrics`], and
//! [`SimulationDiagnostics`], or a [`SimulationError`].
//!
//! A fixed-step solver does not report how accurate a run is. To estimate
//! the error the timestep introduces, use [`estimate_sir_step_error`],
//! [`estimate_seir_step_error`], [`estimate_sis_step_error`],
//! [`estimate_sirs_step_error`] or [`estimate_seirs_step_error`].
//!
//! The models assume a closed, homogeneously mixing population, constant
//! parameters, and frequency-dependent transmission. Recovery gives permanent
//! immunity in SIR and SEIR, immunity that wanes in SIRS and SEIRS, and none in SIS. The
//! results describe the model under those assumptions; they are not
//! forecasts or estimates of a real outbreak.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![deny(missing_debug_implementations)]

mod accuracy;
mod error;
mod metrics;
mod models;
mod result;
mod schedule;
mod settings;
mod simulation;
mod solver;
mod validation;

pub use accuracy::StepErrorEstimate;
pub use error::SimulationError;
pub use metrics::EpidemicMetrics;
pub use models::seir::{SeirParameters, SeirState, estimate_seir_step_error, simulate_seir};
pub use models::seirs::{SeirsParameters, estimate_seirs_step_error, simulate_seirs};
pub use models::sir::{SirParameters, SirState, estimate_sir_step_error, simulate_sir};
pub use models::sirs::{SirsParameters, estimate_sirs_step_error, simulate_sirs};
pub use models::sis::{SisParameters, SisState, estimate_sis_step_error, simulate_sis};
pub use result::{Observation, SimulationDiagnostics, SimulationResult};
pub use schedule::TransmissionSchedule;
pub use settings::SimulationSettings;

// Runs the Rust examples in the Markdown documentation as doctests.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
#[doc = include_str!("../docs/api_design.md")]
struct MarkdownExamples;
