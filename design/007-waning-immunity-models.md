# Design 007: SIS, SIRS and SEIRS (waning immunity)

Spike for plan 007, on branch `spike/007-waning-immunity` from `9468df0`.
The prototype on this branch is not for merging. Only this document, or a
build plan written from it, moves on. Plan 005 (threshold metrics) had not
landed at the starting commit. Q3 therefore says how its metrics would apply
once it lands.

All commands run from the repository root. The prototype is in these files:

| File | What it is |
|---|---|
| `src/simulation/driver.rs` | `run_simulation_with_accumulators(initial, compartments, ...)`; `run_simulation` now calls it with `compartments = N`. `check_state` takes slices. |
| `src/models/sis.rs` | `SisParameters`, `SisState` (with `cumulative_infections`), `simulate_sis`, guard, metrics |
| `src/spike.rs` | `#[doc(hidden)]` SIS, SIRS and SEIRS derivatives with an accumulator, and an unguarded `integrate` |
| `tests/sis_prototype_tests.rs` | logistic closed form, conservation, accumulator checks |
| `tests/waning_stability_evidence_tests.rs` | stability harness for every candidate formula (ignored) |
| `tests/waning_equilibrium_tests.rs` | endemic equilibria and the shape of the infectious curve |

## Summary

1. Count infections with an accumulator `dC/dt = beta S I / N` for every model. The driver takes a count of leading compartments, and the solver is unchanged. It equals `S(0) - S(T)` for SIR to 1.3e-12 and the exact SIS integral to 4.1e-5 people.
2. Guards, measured over 250,000 configurations each with zero flagged runs: SIS `max(beta + gamma, 1.25 gamma)`, SIRS `max(beta + gamma, omega)`, SEIRS `beta + sigma + gamma + omega`. Every SIR-style "largest rate" extension fails.
3. Metrics: replace final size with endemic prevalence `I*/N` (all three formulas confirmed to six decimals). Keep `R0 S / N` unchanged (exact). Report `1 - 1/R0` as a herd-immunity threshold only for SIRS and SEIRS. Give the SIS peak time a tolerance, because roundoff decides it.
4. API: `SisParameters`/`SisState`/`simulate_sis`. `SirsParameters` and `SeirsParameters` reuse `SirState` and `SeirState`. One `#[non_exhaustive]` `EpidemicMetrics` with `Option` fields, and `cumulative_infections` on `Observation`.
5. Build in five plans: the accumulator for SIR and SEIR first, then metrics shape (after 005), then SIS, SIRS and SEIRS.

## Q1 Counting infections

**Recommendation: candidate A, for every model including SIR and SEIR.**
Integrate an accumulator `C` with `dC/dt = beta S I / N` as the last entry
of the state array. Give the driver a count of how many leading entries form
the population: `run_simulation_with_accumulators(initial, compartments, ...)`
in the prototype. Compute `cumulative_infections_since_start` as
`C(T) - C(0)` for every model.

### How the driver tells compartments from accumulators

| Option | Result in the prototype |
|---|---|
| A runtime count of leading compartments | Works. The driver passes `&state[..compartments]` to `check_state`, which now takes slices, so its unit tests compile unchanged. `rk4_step` is untouched. It still checks every entry, `C` included, for finiteness. |
| Two arrays, compartments and accumulators | Rejected. The driver is generic over `const N`, and stable Rust cannot write `[f64; P + A]` (that needs the unstable `generic_const_exprs`). Two arrays would need either the solver to integrate a pair of arrays, which changes `src/solver/rk4.rs` (the plan's STOP line), or concatenating and splitting on every derivative call. |

The count could be a `const` generic in a build. A runtime `usize` with an
`assert!` was enough for the spike, and it costs nothing measurable. The
accumulator is not checked for sign, because it only grows. It is not
summed into `N`. Its initial value is validated as finite and non-negative.

### Measurements

| Check | Result | Command |
|---|---|---|
| SIR with an accumulator: largest `|C(t) - (S(0) - S(t))|` over 160 days at 1/4 day, for `(beta, gamma)` = (0.3, 0.1), (0.5, 0.1), (1.9, 0.4) | 6.8e-13, 1.3e-12, 6.8e-13 people: roundoff | `cargo test --test sis_prototype_tests -- --nocapture` (`sir_accumulator_equals_the_susceptible_drop_within_roundoff`) |
| SIS: largest `|C(t) - exact|` over 60 days at 1/4 day, where exact is `I(t) - I(0) + gamma * integral of the logistic I` in closed form; `(beta, gamma, I0)` = (0.3, 0.1, 10), (1.5, 1.4, 10), (0.1, 0.3, 500) | 4.1e-5, 1.0e-6, 2.6e-5 people | same command (`sis_cumulative_infections_match_the_exact_integral_of_incidence`) |
| SIS: `C(T)/N` | 3.257 and 3.898 at 60 days; 12.59 at 200 days (`beta = 0.3`, `gamma = 0.1`), within 3.7e-5 people of exact | same command (`sis_metrics_count_reinfections`) |

The exact integral used is
`integral of I over [0, T] = K T + (K / r) ln((1 + A e^{-rT}) / (1 + A))`,
with `r = beta - gamma`, `K = N (1 - gamma/beta)` and `A = K / I(0) - 1`.
It holds for `R0 < 1` as well, where `K < 0`.

### What it means for the two metrics

- `cumulative_infections_since_start` keeps its name and its meaning: new
  infections between 0 and `T`. With waning immunity it counts
  reinfections, which is what incidence means. For SIR and SEIR, `C(T) - C(0)`
  equals `S(0) - S(T)` to about 1e-12 people, so computing it the new way
  changes nothing a user can see. One computation for all models is
  simpler than two.
- `cumulative_incidence_fraction` keeps its formula, `(C(T) - C(0)) / N(0)`,
  but its doc must change. It is infections per head of the initial
  population, not the fraction of people infected. It exceeds 1 for SIS
  (measured 12.59 after 200 days) and can for SIRS and SEIRS. For SIR and
  SEIR the two readings coincide. Renaming it, for example to
  `cumulative_infections_per_capita`, is an open question for the owner.
- `calculate_metrics` would take an `infections_of` closure returning
  `C`, in place of `susceptible_of`. The prototype passes `-C` as
  `susceptible_of` (marked `shortcut:`) to avoid touching the SIR and SEIR
  callers.
- The step-error estimate should compare `C` as well. It is integrated
  state with its own discretisation error.

## Q2 Stability guard

**Recommendation** (fastest rate times step at most 2, the library's
`MAX_RATE_TIMES_STEP`, unchanged):

| Model | Fastest decay rate |
|---|---|
| SIS | `max(beta + gamma, 1.25 * gamma)` |
| SIRS | `max(beta + gamma, omega)` |
| SEIRS | `beta + sigma + gamma + omega` |

Each is the simplest formula tried that gave zero flagged runs. None of the
SIR-style "largest single rate" extensions is enough. The advisor's
inconclusive probe is now settled: SIS needs more than `max(beta, gamma)`.

### Harness

`tests/waning_stability_evidence_tests.rs` copies the generator, seed (42),
sample size (250,000 configurations), steps (1 day against 1/64 day),
tolerance (`1e-9 N`) and threshold (0.1) of
`tests/stability_evidence_tests.rs`. `omega` comes from a second generator
(seed 7) so every other value repeats that harness exactly. Every
configuration runs once without a guard. Each candidate then accepts the
runs whose rate times step is at most 2. Checks, per run:

- **compartment / accumulator direction**: over a coarse step, an entry
  moves the opposite way, by more than the tolerance, to the exact flow
  from the same coarse state, where that flow is monotone over the step.
  "Exact" is 64 RK4 substeps from that state. This generalises the SIR
  checks (S never rises, R never falls, I falls below threshold) to models
  where no compartment is monotone everywhere.
- **overshoot** (SIS): `I` crosses `I*`, which the logistic never does.
- **incidence error**: `|C(T) - C_fine(T)| / N > 0.1`, the original check
  with `C` in place of `S(0) - S(T)`.
- **final-state error**: any compartment at `T` differs from the fine run
  by more than 0.1 N. This stands in for "the equilibrium from a fine run",
  since many durations (2 to 42 days) end before the run settles.

**Control rows.** SIR and SEIR, run through the same checks with
`omega = 0`, give **226,590 and 107,633 accepted runs with 0 flagged**. Those
are exactly the counts in `docs/numerical_methods.md`, so the generator copy
is faithful. The stricter direction checks and the accumulator checks also
do not flag anything the library already accepts.

### Results

Command, for every row:
`cargo test --release --test waning_stability_evidence_tests -- --ignored --nocapture`
(348 s on this machine). "Smallest flagged rate x step" is over all
flagged runs, accepted or not. A value above 2 means the candidate has
margin. "Errors" are accepted runs where the coarse run returned a
`SimulationError` (for example `NegativeCompartment`). Those are loud
failures, not counted as accepted, the same as in the existing harness.
The check columns count accepted flagged runs per check. One run can fail
several checks.

| Model | Candidate | Accepted | Flagged | Errors | Smallest flagged rate x step | Compartment dir. / accumulator dir. / overshoot / incidence / final state |
|---|---|---|---|---|---|---|
| SIS | none | 249221 | 18074 | 779 | 0.000 | 0 / 5705 / 318 / 16111 / 25 |
| SIS | `max(beta, gamma)` | 226310 | 12943 | 295 | 1.152 | 0 / 3529 / 120 / 11831 / 14 |
| SIS | `beta + gamma` | 113098 | 374 | 0 | 1.688 (2.106 ignoring accumulator direction) | 0 / 374 / 0 / 0 / 0 |
| SIS | **`max(beta + gamma, 1.25 gamma)`** | 108560 | **0** | 0 | 2.049 | 0 / 0 / 0 / 0 / 0 |
| SIS | `1.25 (beta + gamma)` | 72329 | 0 | 0 | 2.110 | 0 / 0 / 0 / 0 / 0 |
| SIRS | none | 243224 | 4516 | 6776 | 0.000 | 0 / 41 / 0 / 4493 / 35 |
| SIRS | `max(beta, gamma, omega)` | 212468 | 2410 | 3304 | 1.563 | 0 / 15 / 0 / 2402 / 16 |
| SIRS | **`max(beta + gamma, omega)`** | 106450 | **0** | 1215 | 2.042 | 0 / 0 / 0 / 0 / 0 |
| SIRS | spectral `max(-Re lambda)` | 182444 | 306 | 1919 | 1.781 | 0 / 0 / 0 / 305 / 1 |
| SIRS | spectral `max abs(lambda)` | 158215 | 62 | 1395 | 1.874 | 0 / 0 / 0 / 61 / 1 |
| SIRS | `beta + gamma + omega` | 35922 | 0 | 0 | 2.869 | 0 / 0 / 0 / 0 / 0 |
| SEIRS | none | 201123 | 20870 | 48877 | 0.000 | 17708 / 0 / 0 / 7417 / 2509 |
| SEIRS | `max(beta, seir_coupled, omega)` | 100952 | 1185 | 1674 | 1.212 | 845 / 0 / 0 / 343 / 17 |
| SEIRS | spectral `max(-Re lambda)` | 81243 | 392 | 497 | 1.295 | 344 / 0 / 0 / 48 / 7 |
| SEIRS | spectral `max abs(lambda)` | 81243 | 392 | 497 | 1.295 | 344 / 0 / 0 / 48 / 7 |
| SEIRS | `max(beta + gamma, seir_coupled, omega)` | 70532 | 563 | 527 | 1.212 | 562 / 0 / 0 / 0 / 2 |
| SEIRS | `max(beta, seir_coupled) + omega` | 25638 | 2 | 7 | 1.894 | 2 / 0 / 0 / 0 / 0 |
| SEIRS | `max(beta + gamma, seir_coupled + omega)` | 23620 | 2 | 7 | 1.894 | 2 / 0 / 0 / 0 / 0 |
| SEIRS | `max(beta + gamma, seir_coupled, gamma + omega)` | 47369 | 179 | 39 | 1.400 | 179 / 0 / 0 / 0 / 0 |
| SEIRS | `max(beta, seir_coupled, sigma + gamma, omega)` | 84626 | 936 | 930 | 1.221 | 638 / 0 / 0 / 299 / 14 |
| SEIRS | `max(beta + gamma, seir_coupled, sigma + gamma, omega)` | 58046 | 392 | 130 | 1.221 | 392 / 0 / 0 / 0 / 0 |
| SEIRS | `max(beta, seir_coupled, sigma + gamma) + omega` | 21275 | 0 | 0 | 2.049 | 0 / 0 / 0 / 0 / 0 |
| SEIRS | **`beta + sigma + gamma + omega`** | 8501 | **0** | 0 | 2.401 | 0 / 0 / 0 / 0 / 0 |
| SIR control (`omega = 0`) | `max(beta, gamma)` (library) | 226590 | 0 | 15 | none flagged | 0 / 0 / 0 / 0 / 0 |
| SEIR control (`omega = 0`) | `max(beta, seir_coupled)` (library) | 107633 | 0 | 234 | 2.035 | 0 / 0 / 0 / 0 / 0 |

`seir_coupled` is the SEIR term in `src/models/seir.rs`. The spectral rows
are the largest `-Re lambda`, or `abs(lambda)`, of the Jacobian over a grid
of the state space. For SIRS that is the 2x2 Jacobian in `(I, R)` with
`S = N - I - R`, on steps of N/20:

```text
J = [ beta (1 - 2I/N - R/N) - gamma    -beta I/N ]
    [ gamma                            -omega    ]
```

For SEIRS it is the 3x3 Jacobian in `(E, I, R)` on steps of N/8. The test
`spectrum_helpers_reproduce_known_rates` (not ignored) checks the helpers:
with `omega = 0` the SIRS bound is exactly `max(beta, gamma)`, the library's
SIR formula, and the SEIRS bound is at least `seir_coupled`.

### What the numbers say

- **Linear stability is not what binds.** For SIR the spectral bound is
  exactly the library formula. For SIRS and SEIRS the spectral bounds still
  let through hundreds of flagged runs, almost all incidence errors (SIRS)
  or wrong directions (SEIRS). The guard has to be set by measurement, as the
  existing limit of 2 was.
- **SIS, two separate effects.** `beta + gamma`, the linearised rate at
  `I = N`, removes every compartment-level failure (smallest flagged
  product 2.106 without the accumulator check). The remaining 374 flags are
  all the accumulator stepping backwards. With `I` near `N` and `gamma * h`
  above about 1.6, RK4's fourth stage evaluates `I` below zero and `S`
  above `N`, so that stage's incidence `beta S I / N` is negative and
  outweighs the others. An example from the run: `C` changed by -7.0 where
  the exact change was +17.6. Evaluating the stages by hand for `beta = 0`
  and `I(0) = N`, the stage-weighted incidence first turns negative between
  `gamma * h = 1.62` and `1.64`. Hence the `1.25 gamma` term
  (`2 / 1.25 = 1.6`). This affects only models where `S` has an inflow.
  In SIR and SEIR, `dC/dt = -dS/dt`, so the RK4 increment of `C` is exactly
  the negated increment of `S`, which the existing checks already cover.
- **SIRS.** `max(beta + gamma, omega)` is the advisor's lower bound and is
  enough, with a thin margin (2.042, like SEIR's 2.035). It still leaves
  1,215 accepted runs that end in an error, against 15 for SIR. Those
  failures are loud, but users will see them. `beta + gamma + omega` has
  none, but accepts only a third as many runs.
- **SEIRS.** Every flag left for the tighter formulas is `R` rising over a
  step where the exact `R` falls (for example +3.4 against -6.1). That
  happens when `sigma + gamma` is large: the error in `I` along the fast
  `E -> I` chain lands in `R`. SEIR has the same error, but there `R` can
  only rise, so the check cannot see it. Two formulas reach zero:
  `beta + sigma + gamma + omega`, the simplest, with the larger margin (2.401)
  and 8,501 accepted, and `max(beta, seir_coupled, sigma + gamma) + omega`,
  with 21,275 accepted and margin 2.049. I recommend the plain sum, by the
  plan's "simplest" rule and its margin. The cost is that it accepts about one
  configuration in 29 at a one-day step, against one in 2.3 for SEIR. At typical rates
  (`beta = 0.5, sigma = 0.2, gamma = 0.1, omega = 0.01`) it still allows a
  2.47-day step. The choice is question 6 below.
- **Evidence, not proof.** As for SIR and SEIR, zero flags in 250,000
  configurations does not prove every accepted run correct. The thin margins
  (2.04 to 2.05) mean a larger sample could find flags near the limit.

## Q3 Metrics

### Endemic equilibria, measured

Each run starts from `S = 990, I = 10` (SIS, SIRS) or
`S = 980, E = 10, I = 10` (SEIRS) in `N = 1000`. It runs 5000 days at
1/64 day. Produced by
`cargo test --release --test waning_equilibrium_tests -- --nocapture`;
the test asserts agreement within 1e-6 N.

| Model and parameters | predicted I* | measured I(5000) | predicted S* | measured S(5000) | \|I(5000) - I(4900)\| | local maxima of I | first maximum is global |
|---|---|---|---|---|---|---|---|
| SIS beta 0.3 gamma 0.1 | 666.666667 | 666.666667 | 333.333333 | 333.333333 | 0.0e0 | 1 | true |
| SIS beta 1.5 gamma 1.4 | 66.666667 | 66.666667 | 933.333333 | 933.333333 | 3.6e-10 | 0 | false |
| SIS beta 0.1 gamma 0.3 | 0.000000 | 0.000000 | 1000.000000 | 1000.000000 | 0.0e0 | 0 | false |
| SIRS beta 0.3 gamma 0.1 omega 0.01 | 60.606061 | 60.606061 | 333.333333 | 333.333333 | 0.0e0 | 15 | true |
| SIRS beta 0.5 gamma 0.2 omega 0.05 | 120.000000 | 120.000000 | 400.000000 | 400.000000 | 0.0e0 | 10 | true |
| SIRS beta 1.8 gamma 0.3 omega 1.2 | 666.666667 | 666.666667 | 166.666667 | 166.666667 | 0.0e0 | 3 | true |
| SIRS beta 0.1 gamma 0.2 omega 0.05 | 0.000000 | 0.000000 | 1000.000000 | 1000.000000 | 1.5e-212 | 0 | false |
| SEIRS beta 0.5 sigma 0.2 gamma 0.1 omega 0.01 | 69.565217 | 69.565217 | 200.000000 | 200.000000 | 0.0e0 | 13 | true |
| SEIRS beta 0.8 sigma 0.5 gamma 0.25 omega 0.05 | 105.769231 | 105.769231 | 312.500000 | 312.500000 | 0.0e0 | 11 | true |
| SEIRS beta 1.8 sigma 1 gamma 0.3 omega 1.2 | 537.634409 | 537.634409 | 166.666667 | 166.666667 | 0.0e0 | 4 | true |
| SEIRS beta 0.1 sigma 0.2 gamma 0.2 omega 0.05 | 0.000000 | 0.000000 | 1000.000000 | 1000.000000 | 2.4e-124 | 0 | false |

All of the advisor's equilibrium formulas hold, including `I* = 0` and
`S* = N` when `R0 <= 1`: three or four parameter sets per model, agreement
to the six printed decimals. "First maximum is global: false" means there
is no interior maximum: `I` falls from the start, or rises monotonically
(SIS, `beta = 1.5`, `gamma = 1.4`).

### Metric by metric

| Metric | SIS | SIRS | SEIRS |
|---|---|---|---|
| `basic_reproduction_number` | keeps meaning, `beta / gamma` | keeps meaning | keeps meaning |
| `peak_infectious` | keeps meaning, but is usually `I*` approached from below | keeps meaning: first epidemic wave | keeps meaning |
| `peak_infectious_time_days` | **needs a new definition** (below) | keeps "earliest global maximum" | keeps it |
| `cumulative_infections_since_start` | keeps name and meaning, computed from `C` (Q1) | same | same |
| `cumulative_incidence_fraction` | keeps formula; doc becomes "per head", can exceed 1 | same | same |
| final size (plan 005 `final_size_fraction`) | **replace** with endemic prevalence `I*/N` | replace | replace |
| herd-immunity threshold `1 - 1/R0` (plan 005) | **do not report** under that name | keep, with new wording | keep, with new wording |
| effective R `R0 S / N` (plan 005) | keeps meaning exactly | keeps meaning exactly | keeps meaning exactly, for `E + I` |

**Peak.** In SIRS and SEIRS the infectious curve oscillates towards `I*`.
The fine runs above had 3 to 15 local maxima. In all six endemic runs the
first maximum was the global one. "Earliest global maximum" is therefore
still the right definition: it picks the first wave, which is what people
mean by the peak. Damping is not guaranteed for every start, though. A run
that starts with many recovered people can have a later wave above an
earlier one. The definition still says exactly what is returned, so it
stays.

For SIS the definition breaks in practice. `I` rises monotonically to `I*`
and then sits on a plateau where only roundoff moves it. The global maximum
of the recorded values is then decided by roundoff. In the 200-day SIS run
(`sis_metrics_count_reinfections`), `peak_infectious_time_days` came out as
**191.0**, although `I` had been within roundoff of `I*` for over 100 days.
The 5000-day run also counted one spurious "local maximum" on the plateau.
Recommendation: take the earliest recorded time at which `I` is within a
small tolerance of the maximum, for example `64 * EPSILON * N` as the driver
uses. That changes nothing measurable for SIR and SEIR. An alternative is to
document that for SIS the peak time is the time the plateau is reached, to
roundoff. That is an open question for the owner.

**Final size.** It does not exist: there is no final state, only an endemic
equilibrium. Replace it with `endemic_prevalence = I*/N`, where
`I*/N = 1 - 1/R0` for SIS, `(1 - 1/R0) omega / (omega + gamma)` for SIRS,
and `(1 - 1/R0) / (1 + gamma/sigma + gamma/omega)` for SEIRS, or 0 when
`R0 <= 1`. All are checked above. A companion `endemic_susceptible_fraction
= min(1, 1/R0)` is the same for all three. These are theoretical, like
`basic_reproduction_number`, and need no trajectory.

**Herd-immunity threshold.** `1 - 1/R0` is the non-susceptible fraction
above which `R_eff < 1` (derived below). In SIRS and SEIRS that is still a
fraction immune or infected, so the quantity keeps its meaning at any one
moment. But immunity wanes, so it is not a state the population stays in.
At the endemic equilibrium the non-susceptible fraction is exactly
`1 - 1/R0`. Suggested wording for SIRS and SEIRS: "the fraction of the
population that must be non-susceptible at a given time for infections to
decline at that time; with waning immunity it is not maintained without
continued infection or vaccination". For SIS nobody is ever immune, and
`1 - 1/R0` equals the endemic prevalence `I*/N`. Report it as
`endemic_prevalence`, not as a herd-immunity threshold. Recommendation:
plan 005's `herd_immunity_threshold` becomes `Option<f64>`, `None` for SIS.

**Effective R.** Waning only adds terms to the `S` and `R` equations, so
the infection equations are those of SIR and SEIR:

```text
SIS, SIRS:  dI/dt       = beta S I / N - gamma I = gamma I (R0 S/N - 1)
SEIRS:      d(E + I)/dt = beta S I / N - gamma I = gamma I (R0 S/N - 1)
```

So `R_eff = R0 S / N` sets the direction of `I` (SIS, SIRS) and of `E + I`
(SEIRS) exactly, whenever `I > 0`, as it does in SIR and SEIR. That is
derivation, not measurement. What changes is that `S` can grow again, so
`R_eff` can rise back above 1: that is what drives repeated waves. It
tends to 1 at the endemic equilibrium (`S* = N / R0`). Plan 005's
`effective_reproduction_number` method applies to all three models
unchanged.

## Q4 Public API

Following `SirParameters`, `SirState`, `simulate_sir` and
`estimate_sir_step_error`:

| Model | Parameters | State | Functions |
|---|---|---|---|
| SIS | `SisParameters { beta_per_day, gamma_per_day }` | `SisState { susceptible, infected }` | `simulate_sis`, `estimate_sis_step_error` |
| SIRS | `SirsParameters { beta_per_day, gamma_per_day, omega_per_day }` | **reuse `SirState`** | `simulate_sirs`, `estimate_sirs_step_error` |
| SEIRS | `SeirsParameters { beta_per_day, sigma_per_day, gamma_per_day, omega_per_day }` | **reuse `SeirState`** | `simulate_seirs`, `estimate_seirs_step_error` |

`omega_per_day` must be finite and non-negative. With `omega = 0` SIRS is
SIR exactly. The harness's SIR and SEIR control rows run that way.

- **States.** SIRS has the same compartments as SIR, and SEIRS as SEIR.
  Reusing the state types means one result type per compartment layout,
  and users can feed a SIR state into SIRS directly. Only the doc on
  `recovered` changes, from "assumed immune" to "immune (permanently in SIR
  and SEIR, temporarily in SIRS and SEIRS)". Parameters get their own types,
  because adding a field to `SirParameters` would break every SIR caller and
  hide the model choice in a number.
- **Where the accumulator lives.** The prototype put
  `cumulative_infections` into `SisState`. That works, but it makes the
  user supply a 0 that is not part of the population, and it would stop
  SIRS from reusing `SirState`. Recommendation: add
  `cumulative_infections: f64` to `Observation<S>`, next to `time_days`, for
  every model. `S` then stays pure compartments. That is a breaking change
  to a public struct, which is cheap now. The driver would fill it from the
  accumulator entry. That is an open question for the owner. Not prototyped,
  because changing `Observation` edits existing tests that construct it.
- **Metrics.** Recommendation: one `EpidemicMetrics` for all models, marked
  `#[non_exhaustive]` now, before the first release. Model-specific fields
  are `Option<f64>`: plan 005's `final_size_fraction` is `None` for
  waning models, and a new `endemic_prevalence` is `None` for SIR and SEIR.
  `herd_immunity_threshold` is `None` for SIS. Per-family metrics types
  would make `SimulationResult` generic over a second parameter, or need two
  result types, for a difference of two or three fields.
  `#[non_exhaustive]` means later additions, such as births, are not
  breaking. Its cost is that users cannot build `EpidemicMetrics` with a
  struct literal outside the crate. No test or example in the crate does
  today.

## Q5 Validation

Each new model should ship with these. Marks: **G** = an existing test
generalises with new parameters, **N** = needs a new reference.

| Test | SIS | SIRS | SEIRS | Reference |
|---|---|---|---|---|
| Parameter and state validation (negative, NaN, infinite, zero `gamma`, `omega >= 0`) | G | G | G | copy of `sir_rejects_*` |
| Closed form | N: logistic `I(t)` at 1/4 day (prototype measured at most 2.1e-5 people, error falling 15.7 to 16.4x per halving; bound 1e-4) | N: none exists; use `omega = 0` against SIR's existing references, plus the equilibrium | N: as SIRS, against SEIR | `sis_matches_the_logistic_closed_form_at_a_quarter_day_step` |
| Endemic equilibrium after a long fine run, within 1e-6 N, `R0 > 1` and `R0 <= 1` | N | N | N | `waning_equilibrium_tests.rs` |
| Conservation of `S + I (+ E) + R` with the accumulator excluded | G | G | G | `sis_conserves_the_population_and_ignores_the_accumulator` |
| Accumulator: `C(T)` equals `S(0) - S(T)` within roundoff for SIR and SEIR (and for SIRS and SEIRS at `omega = 0`) | n/a | N | N | `sir_accumulator_equals_the_susceptible_drop_within_roundoff` |
| Accumulator: SIS `C(T)` against the closed-form integral | N | n/a | n/a | `sis_cumulative_infections_match_the_exact_integral_of_incidence` |
| Guard boundary: accepted at `rate * h = 2`, `TimeStepTooLarge` just above, `max_step_days` reported | G | G | G | existing `*_surfaces_*` tests |
| Stability evidence harness, zero flagged runs | N: extend `stability_evidence_tests.rs` | N | N | `waning_stability_evidence_tests.rs` |
| Step-error estimate falls about sixteenfold per halving and matches a fine reference | G | G | G | `*_step_error_estimate_*` |
| `R_eff` direction (`I`, or `E + I`, falls when `R0 S/N < 1`) | N | N | N | derivation in Q3 |
| Metrics: endemic prevalence, peak-time tolerance, per-head incidence above 1 | N | N | N | Q3 |
| Convergence order (fourth) | G | G | G | `convergence_tests.rs` |

The SIS closed form is the strongest test available. SIRS and SEIRS have no
closed form, so they rely on `omega = 0` reducing to SIR and SEIR, on the
equilibrium, and on convergence order.

## Open questions for the owner

1. **Where does `cumulative_infections` live?** On `Observation`, in every
   state type, or in a separate field of the result? *Recommended: on
   `Observation<S>`, for all models.*
2. **Rename `cumulative_incidence_fraction`?** It can exceed 1 once
   reinfection exists. *Recommended: keep the name, document it as "per
   head of the initial population". Rename only if users read it as a
   proportion of people.*
3. **One metrics type with `Option` fields, plus `#[non_exhaustive]`?**
   *Recommended: yes, now, before publishing.*
4. **SIS peak time on a plateau.** Use a tolerance on the maximum, or
   document that the time is decided by roundoff? *Recommended: tolerance
   `64 * EPSILON * N`, applied to all models.*
5. **Should SIR and SEIR switch to the accumulator for
   `cumulative_infections_since_start`?** Their values would change by
   about 1e-12 people. *Recommended: yes, one code path.*
6. **SEIRS guard: simple or permissive?** `beta + sigma + gamma + omega`
   (8,501 of 250,000 accepted, margin 2.401) or
   `max(beta, seir_coupled, sigma + gamma) + omega` (21,275 accepted, margin
   2.049). Both give zero flags. *Recommended: the plain sum.*
7. **Must cumulative infections never step backwards?** If a small
   backward step of `C` is acceptable when the end result is within 0.1 N,
   the SIS guard can be `beta + gamma` (113,098 accepted instead of
   108,560). *Recommended: no backward steps, so keep `1.25 gamma`. A
   cumulative count that falls is visibly wrong.*
8. **SIRS accepted-run errors.** 1,215 runs that the guard accepts end in
   a `SimulationError`, against 15 for SIR. Is that acceptable, or should
   the guard be `beta + gamma + omega` (no errors, a third of the accepted
   runs)? *Recommended: accept them. They are loud, and the step is
   inaccurate there anyway. The error type was not recorded in this spike.*

## Proposed build plans

In this order. Each one leaves `cargo test` green and ships its docs
(`docs/mathematical_models.md`, `docs/numerical_methods.md`,
`docs/api_design.md`, `CHANGELOG.md`).

1. **Incidence accumulator for SIR and SEIR** (effort S). Scope: the driver
   learns a compartment count. SIR and SEIR integrate `C`.
   `cumulative_infections_since_start` comes from `C`, and
   `cumulative_infections` is added to `Observation` if the owner agrees
   (question 1). No new model. Files: `src/simulation/driver.rs`,
   `src/metrics.rs`, `src/models/sir.rs`, `src/models/seir.rs`,
   `src/result.rs`, driver and metrics unit tests (constructor changes only),
   and a new test that `C(T)` equals `S(0) - S(T)` within roundoff. Check:
   `stability_evidence_tests` still reports 226,590 and 107,633 accepted runs
   with none flagged.
2. **Metrics shape** (effort S, after plan 005). Scope: `#[non_exhaustive]`
   on `EpidemicMetrics`, `Option` for model-specific fields,
   `endemic_prevalence`, the peak-time tolerance (question 4), and the doc of
   `cumulative_incidence_fraction` (question 2). Files: `src/metrics.rs`,
   `src/models/*.rs`, docs.
3. **SIS** (effort M). Scope: `SisParameters`, `SisState`, `simulate_sis`,
   `estimate_sis_step_error`, guard from the Q2 table, and the Q5 tests,
   including SIS rows in `tests/stability_evidence_tests.rs`. Files:
   `src/models/sis.rs`, `src/lib.rs`, tests, docs.
4. **SIRS** (effort M). Scope: `SirsParameters`, `simulate_sirs` and
   `estimate_sirs_step_error` returning `SirState`, guard
   `max(beta + gamma, omega)`, and the Q5 tests. Files: `src/models/sirs.rs`,
   `src/models/sir.rs` (state doc), `src/lib.rs`, tests, docs.
5. **SEIRS** (effort M). As SIRS, with `SeirState`, and the SEIRS guard from
   Q2. If the owner wants a less conservative SEIRS guard (question 6), run a
   short follow-up spike first. The harness here makes that a matter of
   adding candidate rows.

The prototype's `src/spike.rs`, `tests/waning_*` and the SIS module do not
move on. Plans 3 to 5 should rebuild from the harness design in
`tests/waning_stability_evidence_tests.rs`: local exact steps for the
direction check, all checks reported per run, and SIR and SEIR control
rows.
