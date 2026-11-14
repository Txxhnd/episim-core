# Competitive review, November 2026

Redone on 2026-11-14 before the first crates.io release, as section 16 of
the (frozen) product requirements document requires. It covers `episim-core`
0.2.0.

## Method

- **crates.io:** the search API for "epidemic", "epidemiology", "sir model",
  "seir", "compartmental", "disease model" and "infection", then the
  metadata of every relevant crate. Downloads and dates are as of 2026-11-14.
- **Close reading:** the published source of the three closest Rust crates
  (`sir-epidemic-model` 0.6.0, `jivanu` 1.0.0, `yatrosci-med-epidemiology`
  0.1.0), downloaded from crates.io and read, not run.
- **Other ecosystems:** web searches for R, Python, Julia and C++ tools.
  Versions, licences and required dependencies were then read from CRAN,
  r-universe and PyPI metadata, and `{epidemics}`' release notes were read
  for its current solver. Capabilities come from documentation, not
  hands-on use.

## Name

`episim-core`, `episim` and `epi-sim` are all free on crates.io.

Similar names exist elsewhere: EpiSimR, an R Shiny app for epidemic
simulation on CRAN, and the long-standing EpiSims agent-based simulator. The
crates.io name is unambiguous, but searches for "episim" will also find
those.

## Rust

| Crate | Version (date) | Downloads | License | What it is |
|---|---|---|---|---|
| `jivanu` | 1.0.0 (2026-03-31) | 674 | GPL-3.0-only | A microbiology engine whose `epidemiology` module has single-step SIR, SEIR and SIRS functions, R0, herd-immunity threshold, vaccination coverage and case-fatality helpers. It uses RK4 only with the optional `hisab` feature and forward Euler otherwise, and clips negative compartments to 0 after each step. |
| `sir-epidemic-model` | 0.6.0 (2026-07-01) | 24 | Unlicense | SIR only, on population fractions, with fixed-step RK4, peak and final size. Its final size is a fixed-point iteration that returns its last value after 200 rounds without reporting non-convergence. No repository link. |
| `yatrosci-med-epidemiology` | 0.1.0 (2026-07-03) | 20 | AGPL-3.0-or-later | 96 lines: SIR, SEIR and SIS derivative helpers, with no integrator. |
| `henad`, `henad-models` | 0.3.0 (2026-10-10) | 13, 44 | MIT OR Apache-2.0 | A parallel agent-based modelling engine with example models. |
| `epirust` | 0.0.1 (2020-01-27) | 1,787 | Apache-2.0/MIT | Agent-based epidemiological simulation; one release, six years ago. |
| `sir_ddft` | 2.0.0 (2024-05-14) | 4,099 | AGPL-3.0 | SIR dynamical density functional theory: a spatial PDE model, a research code. |
| `lau-network-science` | 0.1.0 (2026-06-08) | 26 | MIT | General network science, including epidemic spreading on networks. |
| `rebop` | 0.10.0 (2026-08-24) | 25,024 | MIT | Fast stochastic (Gillespie) simulation of reaction networks; a stochastic SIR is one use. |
| `pflow`, `tpt-sci-reaction-network` | 0.3.0, 0.1.1 | 75, 106 | MIT, MIT OR Apache-2.0 | Generic modelling languages (Petri nets, reaction networks) that can express compartmental models. |

Beyond these, any general ODE solver crate can integrate an SIR system.

**Reading.** No Rust crate offers what `episim-core` does in one place:

- five deterministic compartmental models behind one checked API;
- a stability limit backed by kept evidence;
- typed errors instead of clipping;
- a step-error estimate;
- transmission schedules;
- metrics with exact references.

The two closest, `jivanu` and `sir-epidemic-model`, are narrower. `jivanu`
has the broadest overlap (SIR, SEIR, SIRS), but its GPL licence, Euler
default and silent clipping are exactly what this project chose against.
The niche is small: the closest crates have tens to hundreds of downloads.

## Other ecosystems

The mature tools are in R and Python. Versions and licences below come from
CRAN, r-universe and PyPI metadata, read on 2026-11-14. Dependency counts are
the packages each declares as required (R: Imports, Depends and LinkingTo;
Python: `requires_dist`, which can include optional extras). Capabilities
come from each project's documentation and release notes, not from running
it. Where a fact was not confirmed, the table says so.

### Comparison

| | `episim-core` 0.2.0 | EpiModel 2.7.0 | `{epidemics}` 0.5.0 | odin 1.5.12 | epimodels 1.5.0 | epydemix 1.5.0 | epiworldR 0.14 | Starsim 3.7.4 |
|---|---|---|---|---|---|---|---|---|
| Language | Rust | R | R | R DSL, compiled to C | Python | Python | C++ core, R and Python wrappers | Python |
| Licence | MIT | GPL-3 | MIT | MIT | GPL-3.0-or-later | GPL-3 | MIT | MIT |
| Required dependencies | none at run time | 25 | 11 | 8 | 24 | 19 | 4 | 65 |
| Kind of model | deterministic ODE | deterministic (`dcm`), stochastic individual-contact, stochastic network | deterministic ODE with demographic groups | any ODE or discrete-time system the user writes | ODE, discrete-time and continuous-time Markov chain | stochastic compartmental | agent-based | agent-based, with metapopulation and compartmental options |
| Built-in compartmental models | SIR, SEIR, SIS, SIRS, SEIRS | SI, SIR, SIS (`dcm`), others as user functions | a default SEIR-type model, Vacamole, diphtheria and other published models | none: a language for writing them | 27 classes, including SIR, SIS, SIRS, SEIR, vector-borne and multi-strain | defined by the user | SIS, SIR, SEIR and others | SIR and others |
| Integrator | fixed-step RK4 | deSolve | odin, so deSolve (since 0.5.0; previously Boost odeint) | deSolve | SciPy (RK45, LSODA and others) or JAX | stochastic simulation | agent-based steps | agent-based steps |
| Accuracy control | none automatic; the user chooses the step. A stability limit rejects unstable steps, and a step-error estimate is on request | adaptive, from deSolve's tolerances | adaptive, from deSolve | adaptive, from deSolve | adaptive, from SciPy or JAX | not applicable | not applicable | not applicable |
| Interventions over time | piecewise multiplier on `beta` | not confirmed for `dcm` | non-pharmaceutical interventions and vaccination regimes | whatever the user writes | not confirmed | interventions and policy scenarios | not confirmed | yes, in the agent-based framework |
| Waning or no immunity | SIS, SIRS, SEIRS | SIS | not confirmed | whatever the user writes | SIS, SIRS and others | defined by the user | SIS and others | not confirmed |
| Age or contact structure | none | two groups | age groups with contact matrices | whatever the user writes | not confirmed | contact matrices and layers | networks | networks |
| Fitting to data | none | network-model estimation | none (scenarios) | through the newer odin2, dust2 and monty stack | yes | approximate Bayesian computation | not confirmed | not confirmed |
| Failure behaviour | typed errors for invalid input, instability, negative states and population drift; never clips | R errors and warnings from the package and deSolve | R errors | deSolve diagnostics | Python exceptions | Python exceptions | not confirmed | not confirmed |

Smaller tools are left out of the table: `finalsize` (R, MIT, final size in
heterogeneous populations), Compyrtment (Python, MIT, last released 2021),
epidemik (Python, MIT, 0.2.0), the Julia sir-julia recipes and Epidemics.jl,
and Covasim, now part of Starsim.

### What the comparison shows

- **Breadth:** every mature tool covers more than `episim-core`. They offer
  stochastic or agent-based simulation, demographic or network structure, or
  fitting to data. `{epidemics}` is the closest in intent: published
  compartmental models with interventions and vaccination for scenario work.
- **Accuracy control is the main technical gap.** The deterministic tools in
  R and Python use adaptive solvers that hold the error to a tolerance by
  themselves. `episim-core` uses a fixed step. It rejects unstable steps and
  can estimate the error, but it does not control it, so the user must pick
  a step. This was a deliberate choice for a predictable time grid and no
  dependencies, but it is a real cost for users who want accuracy without
  thinking about steps.
- **Where `episim-core` is distinctive:**
  - It runs in-process from Rust with no runtime dependencies. The
    alternatives need an R or Python runtime and 4 to 65 required packages.
  - Its numerical behaviour is pinned by tests against exact references, and
    its step limits by kept evidence.
  - Invalid states are typed errors, never clipped.
  - It is MIT, as are `{epidemics}`, odin, epiworld and Starsim. EpiModel,
    epimodels and epydemix are GPL.
- **Silent failure is not only a Rust-crate problem.** `{epidemics}` 0.5.0's
  release notes fix a bug that silently dropped or merged demographic groups
  in the output when a model had 10 or more groups. That is the class of
  failure this project's checks are built to turn into errors. It is not a
  claim that `episim-core` is free of bugs: its own issues #1 to #3 were
  silent wrong values.

## Where `episim-core` 0.2.0 stands

What it has that the reviewed Rust crates do not, and, for the dependency
footprint and failure behaviour, the tools in other languages too:

- **Numerical guarantees that are tested and kept in the repository:**
  - a stability limit per model, backed by a 250,000-run evidence harness
    that CI runs weekly, including with schedules;
  - typed errors instead of clipping;
  - population conservation checked every step;
  - a step-error estimate.
- **Exact references in the tests:** the logistic SIS solution, endemic
  equilibria, the final-size equation, and the S–R relation.
- **Five models and transmission schedules** behind one API, with no
  runtime dependencies and an MIT licence. That licence is easier to adopt
  than the GPL and AGPL of three of the Rust neighbours.

Gaps against the wider field:

- **Adaptive step control.** This is not excluded by any decision, so it is
  the one gap worth reconsidering. The product document's long-term
  direction already lists adaptive integration.
- **Outside the current scope by decision:**
  - stochastic and agent-based models;
  - network and age structure;
  - vaccination as a process;
  - fitting to data.

## Recommendations before publishing

1. **Publish under `episim-core`.** The name is free, and the description
   and keywords are adequate.
2. **Consider the keywords.** They are `epidemiology`, `simulation`, `sir`,
   `seir` and `ode`. crates.io allows five. Nothing found by this review
   competes for `sirs` or `seirs`, but `ode` brings in readers looking for
   numerics, and that matches the project's strength. No change recommended.
3. **Consider adaptive step control for a later release.** An embedded
   pair such as Dormand-Prince 5(4) would remove the main technical gap
   against the R and Python tools. It needs its own design study: the
   stability limit, the evidence harness, the recording interval and the
   schedule's on-grid rule all assume a fixed step.
4. **State the positioning in the README in one sentence:** a small,
   dependency-free Rust library for deterministic compartmental models with
   checked numerics, not a replacement for EpiModel or `{epidemics}`. The
   README's scope section already says what is missing; it does not yet say
   what sets the library apart.
5. **Expect few users.** The closest crates have tens to hundreds of
   downloads. The library's value is as a correct, well-documented reference
   and teaching tool, as the product document intends, not as a widely used
   dependency.

## Sources

- crates.io search and crate API, queried 2026-11-14; crate sources from
  static.crates.io.
- CRAN, r-universe and PyPI package metadata (versions, licences,
  dependencies) for EpiModel, epidemics, odin, epiworldR, finalsize,
  epimodels, epydemix, Starsim, Covasim, Compyrtment and epidemik, read
  2026-11-14.
- [EpiModel `dcm` reference](https://epimodel.github.io/EpiModel/reference/dcm.html),
  [New DCMs with EpiModel](https://rdrr.io/cran/EpiModel/f/inst/web/NewDCMs.Rmd),
  [EpiModel on CRAN](https://cran.rstudio.com/web/packages/EpiModel/index.html),
  [EpiModel on r-universe](https://epimodel.r-universe.dev/EpiModel/DESCRIPTION),
  [MIDAS catalog entry](https://catalog.midasnetwork.us/collection/1923).
- [epidemics on r-universe](https://epiverse-trace.r-universe.dev/epidemics),
  [epidemics releases](https://www.github.com/epiverse-trace/epidemics/releases),
  [epidemics v0.1.0 announcement](https://epiverse-trace.github.io/posts/epidemics_v0.1.0/index.html),
  [epidemics NEWS](https://epiverse-trace.r-universe.dev/epidemics/NEWS),
  [epidemics design principles](https://epiverse-trace.r-universe.dev/epidemics/doc/design-principles.Rmd).
- [odin on GitHub](https://www.github.com/mrc-ide/odin),
  [odin and dust tutorial](https://mrc-ide.github.io/odin-dust-tutorial/odin.html),
  [odin2](https://mrc-ide.r-universe.dev/odin2/doc/odin2.Rmd).
- [epiworldR on GitHub](https://github.com/UofUEpi/epiworldR),
  [epiworldR manual](https://cran.r-universe.dev/epiworldR/doc/manual.pdf),
  [epiworldpy](https://uofuepibio.github.io/epiworldpy).
- [epimodels on PyPI](https://pypi.org/project/epimodels),
  [epimodels wiki](https://github.com/fccoelho/epimodels/wiki),
  [Epydemix paper](https://www.networkscienceinstitute.org/publications/epydemix-an-open-source-python-package-for-epidemic-modeling-with-integrated-approximate-bayesian-calibration),
  [Starsim on PyPI](https://pypi.org/project/starsim/3.5.2/),
  [Covasim](https://github.com/starsimhub/covasim),
  [epydemix documentation](https://epydemix.readthedocs.io/),
  [Compyrtment on PyPI](https://pypi.org/project/compyrtment/),
  [popdynamics](https://github.com/popdynamics/popdynamics),
  [seir](https://seir.readthedocs.io/en/latest/readme.html),
  [epidemik on PyPI](https://pypi.org/project/epidemik/).
- [sir-julia tutorials](https://github.com/epirecipes/sir-julia/wiki),
  [Epidemics.jl](https://pratikunterwegs.github.io/Epidemics.jl/dev/).
- [EpiSimR on CRAN](https://cran.r-project.org/web/packages/EpiSimR/index.html).
