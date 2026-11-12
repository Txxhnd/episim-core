# EPISIM-CORE
## Product Vision & MVP Product Requirements Document

**Document status:** Draft for Review  
**Version:** 0.1  
**Date:** November 11, 2026  
**Product type:** Rust Computational Library  
**Development target:** 5 weeks  
**Release target:** v0.1.0  
**Working name:** EPISIM-CORE

---

# PART I. PRODUCT VISION

## 1. Product Overview

EPISIM-CORE is a lightweight, modular epidemiological simulation engine written entirely in Rust.

Its purpose is to provide a reusable computational foundation for simulating infectious disease dynamics through mathematical compartmental models.

The library allows developers and computational researchers to define population states, configure epidemiological parameters, execute simulations, and retrieve structured results without implementing numerical integration and epidemiological calculations from scratch.

EPISIM-CORE is designed as an independent computational library. It does not require a graphical interface, web server, database, or command-line application.

The initial product focuses on deterministic compartmental models, with SIR and SEIR as its foundational implementations.

## 2. Vision Statement

**To become a reliable, extensible, and scientifically transparent Rust foundation for epidemiological simulation and computational infectious disease modeling.**

EPISIM-CORE aims to make epidemiological simulation accessible to developers while maintaining clear mathematical assumptions, reproducible numerical computation, and well-defined software interfaces.

The long-term vision is to support a growing collection of epidemiological models that can be integrated into research tools, educational applications, analytical software, and other simulation systems.

## 3. Mission

The product mission is to:

1. Provide accurate and reproducible implementations of established epidemiological models.
2. Simplify the implementation of compartmental simulations for Rust developers.
3. Make mathematical assumptions and numerical behavior transparent.
4. Establish a modular foundation for introducing additional epidemiological models.
5. Provide a lightweight alternative to building a complete epidemiological simulation framework from scratch.

## 4. Target Audience

### Primary audience

**Rust developers and software engineers**

Developers who need epidemiological simulation capabilities within their own applications and do not want to implement the underlying mathematical models independently.

### Secondary audience

**Computational epidemiology researchers and students**

People working with mathematical disease models who want a reusable Rust-based computational component.

### Additional audience

**Educational software developers**

Developers building tools for teaching epidemiology, mathematical modeling, or infectious disease dynamics.

## 5. Core Product Principles

### Scientific transparency

Every implemented model must have documented equations, parameter definitions, units, and assumptions.

### Numerical reliability

Simulation results must be validated through mathematical invariants, numerical convergence tests, and reference scenarios.

### Minimalism

The core library should avoid unnecessary dependencies, infrastructure, and application-level functionality.

### Reusability

The public API should allow other applications to consume simulation results without depending on EPISIM's internal implementation details.

### Extensibility

The architecture should permit additional models and numerical methods to be introduced without requiring a complete redesign.

### Reproducibility

Identical inputs and execution environments must produce consistent deterministic simulation results.

## 6. Long-Term Product Direction

The following capabilities represent the broader product vision. They are not commitments for the first release.

| Development area | Potential capabilities |
|---|---|
| Compartmental models | SIR, SEIR, SIRS, SIS, SEIRS |
| Advanced population models | Age-stratified and heterogeneous populations |
| Transmission dynamics | Time-varying transmission parameters |
| Intervention modeling | Vaccination, isolation, and social interventions |
| Numerical methods | Adaptive integration and alternative solvers |
| Stochastic modeling | Stochastic compartmental simulation |
| Individual simulation | Agent-based epidemiological models |
| Analytical capabilities | Sensitivity analysis and parameter estimation |
| Ecosystem integration | WebAssembly and downstream application bindings |

The long-term direction is to evolve from a small deterministic simulation library into a modular epidemiological computational engine.

Expansion will depend on demonstrated requirements and validation needs rather than simply increasing the number of supported features.

---

# PART II. MVP PRODUCT REQUIREMENTS DOCUMENT

## 1. Executive Summary

EPISIM-CORE v0.1.0 will establish the foundational functionality of the library through two deterministic compartmental epidemiological models:

- SIR: Susceptible, Infected, Recovered.
- SEIR: Susceptible, Exposed, Infected, Recovered.

The MVP will provide a reusable Rust API for configuring, executing, and retrieving epidemiological simulations.

The first release prioritizes mathematical correctness, numerical reliability, predictable API behavior, and minimal dependencies.

The MVP will not include a graphical interface, CLI, web application, stochastic simulation, or agent-based modeling.

**The primary objective is to deliver a working and well-documented Rust library that developers can integrate into their own applications.**

## 2. Problem Statement

Developers implementing epidemiological simulation software frequently need to handle several responsibilities:

- Translating epidemiological equations into executable algorithms.
- Implementing numerical integration methods.
- Validating initial population conditions and model parameters.
- Maintaining population conservation.
- Detecting numerical instability.
- Calculating epidemiological indicators.
- Structuring and exposing simulation results.

These responsibilities can lead to duplicated implementations, inconsistent assumptions, and insufficient numerical validation.

Existing Rust projects already demonstrate that SIR and SEIR simulation are feasible within the ecosystem. However, there remains an opportunity to build a focused library with an explicit emphasis on model assumptions, numerical validation, consistent APIs, and reproducible computational behavior.

EPISIM-CORE addresses this problem by providing a small, independently usable simulation engine.

## 3. MVP Goal

The MVP must allow a Rust developer to:

1. Configure a basic SIR or SEIR simulation.
2. Define the initial population distribution.
3. Specify epidemiological parameters.
4. Configure the simulation duration and time step.
5. Execute the simulation through a public Rust API.
6. Retrieve the resulting epidemiological time series.
7. Access basic simulation metrics.
8. Receive meaningful errors when inputs or numerical results are invalid.

The MVP is successful when these capabilities work reliably without requiring external application infrastructure.

## 4. Target Users and User Stories

### User Story 1: Execute a simulation

As a Rust developer, I want to configure an epidemiological model and execute it through a simple API so that I can incorporate disease simulations into my application.

**Expected outcome:** A valid configuration produces a complete simulation result.

### User Story 2: Inspect simulation trajectories

As a computational researcher, I want to retrieve the number of individuals in each compartment at every simulation time point so that I can analyze disease dynamics.

**Expected outcome:** The library returns a structured time-series representation.

### User Story 3: Evaluate epidemic characteristics

As a researcher, I want to retrieve basic epidemiological metrics so that I can evaluate the simulated epidemic trajectory.

**Expected outcome:** The library provides clearly defined metrics derived from the simulation.

### User Story 4: Validate model configuration

As a developer, I want invalid parameters and initial conditions to generate descriptive errors so that I can identify configuration problems.

**Expected outcome:** Invalid input is rejected without unexpected panics.

### User Story 5: Verify numerical reliability

As a computational researcher, I want predictable numerical behavior and documented model assumptions so that I can assess whether simulation results are trustworthy within the model's limitations.

**Expected outcome:** The library passes mathematical invariant and numerical convergence tests.

## 5. MVP Scope

### 5.1 Included Features

| Feature | Priority | Description |
|---|---|---|
| SIR model | P0 | Deterministic SIR simulation |
| SEIR model | P0 | Deterministic SEIR simulation |
| RK4 solver | P0 | Fixed-step fourth-order Runge-Kutta |
| Population initialization | P0 | User-defined initial compartment values |
| Parameter validation | P0 | Validate parameters and initial conditions |
| Simulation configuration | P0 | Configure duration and time step |
| Time-series output | P0 | Return compartment values over time |
| Basic epidemiological metrics | P0 | Calculate selected metrics |
| Numerical diagnostics | P0 | Detect invalid numerical states |
| Typed error handling | P0 | Structured and descriptive errors |
| Unit tests | P0 | Test mathematical behavior |
| Integration tests | P0 | Test complete simulations |
| API documentation | P0 | Document public library interfaces |
| Usage examples | P0 | Provide executable SIR and SEIR examples |
| Performance benchmarks | P1 | Record baseline computation performance |

P0 indicates a mandatory release requirement. P1 indicates a supporting feature that should be completed where feasible.

### 5.2 Explicitly Excluded Features

The following are outside the MVP scope:

- Command-line interface.
- Graphical user interface.
- Desktop application.
- Web application or backend server.
- WebAssembly bindings.
- Stochastic simulation.
- Agent-based simulation.
- Spatial epidemiology.
- Age-stratified populations.
- Birth and death dynamics.
- Vaccination and intervention modules.
- Time-varying transmission parameters.
- Parameter estimation and calibration.
- Forecasting from real-world epidemiological data.
- Database integration.
- Built-in JSON or CSV file export.

These exclusions are intentional. The first release must establish a dependable computational foundation before expanding its functionality.

## 6. Mathematical Model Requirements

### 6.1 SIR Model

The SIR model divides the population into three compartments:

- Susceptible (S): individuals who can become infected.
- Infected (I): individuals who are currently infectious.
- Recovered (R): individuals who have recovered and are assumed immune.

The governing equations are:

    dS/dt = -βSI/N

    dI/dt = βSI/N - γI

    dR/dt = γI

Where:

- β = transmission rate per day.
- γ = recovery rate per day.
- N = total population.

Population conservation:

    S + I + R = N

### 6.2 SEIR Model

The SEIR model introduces an additional exposed compartment.

- Susceptible (S): individuals susceptible to infection.
- Exposed (E): infected individuals who are not yet infectious.
- Infected (I): infectious individuals.
- Recovered (R): individuals who have recovered and are assumed immune.

The governing equations are:

    dS/dt = -βSI/N

    dE/dt = βSI/N - σE

    dI/dt = σE - γI

    dR/dt = γI

Where:

- β = transmission rate per day.
- σ = progression rate from exposed to infectious per day.
- γ = recovery rate per day.
- N = total population.

Population conservation:

    S + E + I + R = N

The exposed compartment represents the latent, non-infectious period in this simplified model. It should not automatically be interpreted as the full clinical incubation period.

### 6.3 Shared Modeling Assumptions

Both MVP models must explicitly document the following assumptions:

1. The population is closed.
2. No births, deaths, or migration occur.
3. All individuals experience homogeneous mixing.
4. Epidemiological parameters remain constant throughout a simulation.
5. Recovered individuals do not return to the susceptible compartment.
6. Only the infected compartment transmits infection.
7. The model is deterministic.

The transmission formulation uses frequency-dependent transmission, with infection pressure proportional to I/N.

The total population is derived from the initial state rather than independently configured.

## 7. Parameter and Data Requirements

### 7.1 Initial Population State

The initial state must contain finite, non-negative floating-point values.

SIR:

    initial_state = {
        susceptible,
        infected,
        recovered
    }

SEIR:

    initial_state = {
        susceptible,
        exposed,
        infected,
        recovered
    }

Population values are represented as f64.

Fractional values are permitted because deterministic compartmental simulations may produce expected population counts rather than integer individual counts.

The total initial population must be greater than zero and finite.

### 7.2 Epidemiological Parameters

| Parameter | Requirement |
|---|---|
| β | Finite and greater than or equal to zero |
| γ | Finite and strictly greater than zero |
| σ | Finite and strictly greater than zero for SEIR |

The time convention for the MVP is days. Epidemiological rate parameters must therefore use inverse days.

### 7.3 Simulation Settings

The user must provide:

- Simulation duration in days.
- Fixed integration step in days.

Requirements:

- Duration must be finite and greater than zero.
- Time step must be finite and greater than zero.
- Time step must not exceed the solver's stability limit: β, γ, or (for SEIR) σ multiplied by the time step must not exceed 2.785.
- The final integration step must be shortened when necessary to reach the requested duration.
- The initial state and final state must both appear in the resulting trajectory.

## 8. Numerical Solver Requirements

The MVP will use the classical fourth-order Runge-Kutta method (RK4).

The solver must:

1. Integrate the SIR and SEIR differential equations.
2. Support fixed time steps.
3. Support fractional time steps.
4. Handle a final partial time interval.
5. Preserve the simulation time associated with each state.
6. Detect non-finite numerical results.
7. Detect materially negative compartment values.
8. Avoid silently clipping invalid values to zero.
9. Reject, before the simulation starts, a time step too large for stable integration with the configured rates.

The numerical implementation must not assume that RK4 guarantees non-negative compartment values for every possible time-step size.

The solver must document its numerical limitations and return a typed error when the simulation becomes numerically invalid beyond its documented tolerance.

Beyond its stability limit RK4 can turn decay into growth without producing a negative value or any population drift, so this case cannot be left to checks made during the simulation. It must be rejected from the configuration. The limit guards stability only; it does not guarantee accuracy.

## 9. Simulation Output Requirements

The library must return a structured simulation result.

### 9.1 Time Series

Each observation must contain:

- Simulation time.
- Susceptible population.
- Exposed population, when applicable.
- Infected population.
- Recovered population.

The trajectory must contain the initial state and the final state.

If the requested duration is not evenly divisible by the selected time step, the final interval must be shorter so the output reaches the requested simulation duration.

### 9.2 Epidemiological Metrics

The MVP will calculate:

**Basic reproduction number (R0)**

For the specified simple SIR and SEIR models:

    R0 = β / γ

This is the theoretical model-derived reproduction number under the stated assumptions. It is not an estimate derived from observed case data.

**Peak infectious population**

The maximum recorded number of infectious individuals in the trajectory.

**Peak infectious time**

The simulation time at which the maximum recorded infectious population occurs.

Peak metrics are based on the returned time-series observations. Their temporal precision depends on the selected time step.

**New infections since simulation start**

    cumulative_infections(t) = S(0) - S(t)

**Cumulative incidence fraction**

    cumulative_incidence_fraction(t) =
        [S(0) - S(t)] / N(0)

This definition measures the proportion of the initial population newly infected during the simulation. Initial prevalent infections are not counted as newly occurring infections.

The MVP must not present this value as a data-derived estimate of incidence in a real population.

## 10. Functional Requirements

| ID | Requirement | Priority | Acceptance Condition |
|---|---|---|---|
| FR-01 | SIR simulation | P0 | Valid SIR configuration produces a trajectory |
| FR-02 | SEIR simulation | P0 | Valid SEIR configuration produces a trajectory |
| FR-03 | Initial state validation | P0 | Invalid population values are rejected |
| FR-04 | Parameter validation | P0 | Invalid epidemiological rates are rejected |
| FR-05 | Time validation | P0 | Invalid duration and timestep, including a timestep beyond the solver's stability limit, are rejected |
| FR-06 | Numerical integration | P0 | RK4 produces valid results for supported scenarios |
| FR-07 | Time-series generation | P0 | Output contains initial and final states |
| FR-08 | Epidemiological metrics | P0 | Defined metrics are calculated correctly |
| FR-09 | Numerical diagnostics | P0 | Materially invalid numerical states produce errors |
| FR-10 | Deterministic execution | P0 | Identical inputs produce equivalent results |
| FR-11 | Documentation | P0 | Public interfaces have explanatory documentation |
| FR-12 | Runnable examples | P0 | SIR and SEIR examples compile and execute |

## 11. Non-Functional Requirements

### Performance

- The library must run efficiently on ordinary x86_64 CPU hardware.
- No GPU dependency is permitted.
- Simulation performance must be benchmarked using optimized release builds.
- The benchmark environment and workload must be documented.
- The benchmark must cover representative trajectory lengths, including 10,000, 100,000, and 1,000,000 integration steps where feasible.

No absolute runtime guarantee will be claimed before measurement.

### Reliability

- Invalid user input must return structured errors.
- Core simulation code must not panic under normal user input.
- Numerical instability must be distinguishable from invalid configuration.
- No silent correction of materially invalid model states is permitted.

### Reproducibility

- No random-number generation is required in the MVP.
- Repeated runs using identical configurations must produce equivalent results.
- Cross-platform bitwise identity is not a release requirement.

### Maintainability

- The project must use idiomatic Rust conventions.
- Public interfaces must have clear documentation.
- Model equations and solver implementation must remain logically separated.
- Core functionality must not depend on application-specific infrastructure.

### Dependencies

- No mandatory runtime dependency is preferred.
- Serialization dependencies are deferred.
- Benchmark dependencies must remain development-only.
- Unsafe Rust is prohibited in the initial release.

### Compatibility

- The project must compile on the target Intel macOS development environment.
- The declared minimum supported Rust version must be documented.
- Standard Rust compilation should make additional platform integrations unnecessary.

## 12. Success Metrics

The MVP will be evaluated through the following measurable release criteria.

| Metric | Target |
|---|---|
| SIR model functionality | All required SIR tests pass |
| SEIR model functionality | All required SEIR tests pass |
| Population conservation | Relative drift no greater than 1e-9 in designated reference fixtures |
| Numerical positivity | No materially negative state values in supported reference scenarios |
| Numerical convergence | Halving the timestep reduces trajectory error by at least a factor of 8 in selected smooth reference scenarios, before the floating-point error floor |
| Input validation | All defined invalid-input test cases return typed errors |
| Determinism | Repeated same-environment runs produce equivalent results |
| Documentation | All public API items are documented |
| Documentation tests | All included runnable documentation examples pass |
| Benchmarking | Baseline benchmark results are recorded and reproducible |
| Runtime dependencies | Zero mandatory runtime dependencies |

The numerical thresholds above are proposed release targets. They must be verified against independent reference solutions and documented test fixtures before release.

### Developer Experience Metric

A developer familiar with basic Rust should be able to:

1. Add the crate to a project.
2. Construct a valid SIR configuration.
3. Run a simulation.
4. Retrieve the trajectory.

The quick-start example should demonstrate this workflow without requiring additional infrastructure or custom application code.

## 13. MVP Release Acceptance Criteria

The v0.1.0 release is considered ready when:

- [ ] SIR implementation is complete.
- [ ] SEIR implementation is complete.
- [ ] RK4 integration is validated.
- [ ] Initial state and parameter validation are implemented.
- [ ] Time-series output is correctly generated.
- [ ] Epidemiological metrics pass their reference tests.
- [ ] Population conservation tests pass.
- [ ] Numerical convergence tests pass.
- [ ] Invalid numerical states are handled explicitly.
- [ ] Unit and integration tests pass.
- [ ] Documentation examples compile.
- [ ] Public API documentation is complete.
- [ ] Release benchmarks are recorded.
- [ ] Release packaging is verified.

## 14. Risks and Mitigation

| Risk | Potential impact | Mitigation |
|---|---|---|
| Incorrect mathematical formulation | Scientifically invalid output | Validate equations against established references |
| Numerical instability | Invalid compartment values | Step-size tests and numerical diagnostics |
| Ambiguous epidemiological metrics | Misinterpretation of results | Document precise definitions and units |
| Overengineering | Delayed MVP delivery | Keep public API concrete and limited |
| Excessive dependencies | Reduced portability | Keep core dependency-free |
| Scope expansion | Missed delivery target | Explicitly defer advanced models |
| API redesign | Downstream compatibility problems | Review public API before release |

## 15. Development Roadmap

**Target duration: 5 weeks**

| Week | Focus | Deliverable |
|---|---|---|
| 1 | Foundation and API | Crate structure, public types, validation |
| 2 | SIR implementation | Working SIR engine and RK4 solver |
| 3 | SEIR implementation | Working SEIR model and integration tests |
| 4 | Metrics and numerical validation | Complete result structures and analytical tests |
| 5 | Quality assurance and documentation | Benchmark results, examples, and v0.1.0 release |

The timeline may extend to six weeks if numerical validation identifies issues that require additional investigation.

## 16. Research and Competitive Context

EPISIM-CORE enters an existing ecosystem of epidemiological modeling software.

Relevant references include:

- **EpiModel:** An established R framework for deterministic and stochastic epidemic modeling.
- **sir-epidemic-model:** A Rust implementation of deterministic SIR with fixed-step RK4 integration.
- **jivanu:** A broader Rust microbiology library that includes SIR, SEIR, and SIRS functionality.

EPISIM-CORE will not compete by supporting the greatest number of models in its first release.

Its intended differentiation is a narrowly scoped library emphasizing:

1. Explicit mathematical assumptions.
2. Tested numerical behavior.
3. Clear and consistent Rust interfaces.
4. Predictable input and output contracts.
5. Minimal dependencies.
6. Reproducible computation.

The competitive landscape must be reviewed again before public release.

## 17. Reference Materials

### Epidemiological Foundations

1. Kermack, W. O., & McKendrick, A. G. (1927). *Contributions to the Mathematical Theory of Epidemics*. [PubMed](https://pubmed.ncbi.nlm.nih.gov/2059741/).

2. *Structural identifiability and observability of compartmental models of the COVID-19 pandemic*. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7752088/).

3. CDC. (2026). *Rt: Estimating the Direction of Disease Transmission*. [CDC Modeling Handbook](https://www.cdc.gov/cfa-modeling-and-forecasting/modeling-handbook/mh-rt.html).

4. Efendiyev (2026). *Numerical Investigation of a Diffusive SIR Model: Focus on Positivity Preservation*. [DOI](https://doi.org/10.1002/mma.70610).

### Existing Software

5. [sir-epidemic-model, Docs.rs](https://docs.rs/sir-epidemic-model/latest/sir_epidemic_model/).

6. [jivanu epidemiology module, Docs.rs](https://docs.rs/jivanu/latest/jivanu/epidemiology/).

7. [EpiModel official website](https://epimodel.org/).

### Rust Engineering

8. [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/).

9. [The Rustdoc Book](https://doc.rust-lang.org/rustdoc/how-to-write-documentation.html).

10. [The Cargo Book: Features](https://doc.rust-lang.org/cargo/reference/features.html).

---

# PART III. PRODUCT DECISIONS AND APPROVAL

## 18. Proposed MVP Decisions

The following decisions define the proposed baseline for implementation.

| Decision | Proposed selection |
|---|---|
| Product architecture | Pure Rust library |
| Initial models | SIR and SEIR |
| Numerical integration | Classical RK4 |
| Integration strategy | Fixed timestep |
| Numerical precision | f64 |
| Population representation | Floating-point expected counts |
| Model parameters | Constant over simulation duration |
| Population dynamics | Closed population |
| Public model API | Concrete SIR and SEIR interfaces |
| Runtime dependencies | None |
| Serialization | Deferred |
| CLI and GUI | Excluded |
| Initial release | v0.1.0 |

The public API should initially expose concrete model-specific functions. Internal abstractions may be shared between models, but a generic public model-extension trait is deferred until a real extension requirement has been established.

## 19. Definition of Success

EPISIM-CORE v0.1.0 will be successful when it demonstrates that a small Rust library can provide reliable and reusable epidemiological simulation functionality without introducing unnecessary application infrastructure.

The first release is not expected to become a comprehensive epidemiological modeling platform.

It must instead establish a trustworthy computational foundation that can support future development.

**Product Vision:** A reusable Rust foundation for epidemiological simulation.

**MVP Mission:** Deliver a validated, deterministic SIR and SEIR simulation library with a clear public API and reproducible numerical behavior.

**Release Objective:** Establish EPISIM-CORE v0.1.0 as a technically sound foundation for subsequent model expansion.