# Mathematical models

`episim-core` implements two deterministic compartmental models. This page
gives their equations, the units and valid ranges of every quantity, the
assumptions behind them, and the exact definition of each reported metric.

## Units

| Quantity | Unit |
|---|---|
| Time, duration, timestep | days |
| `beta_per_day`, `sigma_per_day`, `gamma_per_day` | per day |
| Compartments, peak infectious, cumulative infections | individuals |
| `R0`, cumulative incidence fraction, relative drift | dimensionless |

Compartment values are `f64` and may be fractional. A deterministic model
tracks expected counts, not whole individuals.

## Shared assumptions

1. The population is closed: no births, deaths, or migration.
2. Individuals mix homogeneously.
3. Parameters are constant over the simulation.
4. Recovered individuals are permanently immune.
5. Only the infected compartment transmits infection.
6. The dynamics are deterministic.

### Frequency-dependent transmission

New infections occur at the rate

```text
infection flow = beta * (S / N) * I
```

`N` is the sum of the initial compartments. It is computed once and held
fixed for the whole run; it is not recomputed from the evolving state. The
code evaluates `S / N` first, which avoids multiplying two large population
values together.

## SIR

Compartments: susceptible `S`, infected (infectious) `I`, recovered `R`.

```text
dS/dt = -beta * S * I / N
dI/dt =  beta * S * I / N - gamma * I
dR/dt =  gamma * I
```

| Parameter | Meaning | Valid range |
|---|---|---|
| `beta_per_day` | transmission rate | finite, `>= 0` |
| `gamma_per_day` | recovery rate | finite, `> 0` |

## SEIR

Compartments: susceptible `S`, exposed `E`, infected (infectious) `I`,
recovered `R`. Exposed individuals are infected but not yet infectious. `E`
models a latent period; it should not be read as the clinical incubation
period.

```text
dS/dt = -beta * S * I / N
dE/dt =  beta * S * I / N - sigma * E
dI/dt =  sigma * E - gamma * I
dR/dt =  gamma * I
```

| Parameter | Meaning | Valid range |
|---|---|---|
| `beta_per_day` | transmission rate | finite, `>= 0` |
| `sigma_per_day` | rate of progression from exposed to infectious | finite, `> 0` |
| `gamma_per_day` | recovery rate | finite, `> 0` |

## Initial state

Every compartment must be finite and non-negative. Their sum `N` must be
finite and strictly positive. The sum is checked separately because finite
compartments can still add up to infinity.

## Population conservation

In both models the derivatives sum to zero, so `N(t) = N(0)` for the exact
solution. The library checks this after every integration step; see
[numerical methods](numerical_methods.md).

## Metrics

All metrics are model-derived quantities. None is an estimate from observed
case data.

| Metric | Definition |
|---|---|
| `basic_reproduction_number` | `R0 = beta / gamma`, for both models. The theoretical number of secondary infections from one case in a fully susceptible population. |
| `peak_infectious` | The largest `I` among the recorded observations. |
| `peak_infectious_time_days` | The time of the earliest recorded observation at that peak. If `I` never rises above its initial value, this is day zero. |
| `cumulative_infections_since_start` | `S(0) - S(T)`, where `T` is the final simulation time. Infections already present at day zero are not counted. |
| `cumulative_incidence_fraction` | `(S(0) - S(T)) / N(0)`. |

The cumulative metrics refer to the final time `T` only. For any earlier
time `t`, compute `S(0) - S(t)` from the trajectory.

There is no effective reproduction number (`Rt`) in this version.

## Closed-form reference cases

With `beta = 0` there is no transmission and both systems are linear. These
solutions are used as independent references in the test suite.

SIR:

```text
S(t) = S0
I(t) = I0 * exp(-gamma * t)
R(t) = R0 + I0 * (1 - exp(-gamma * t))
```

SEIR, for `gamma != sigma`:

```text
S(t) = S0
E(t) = E0 * exp(-sigma * t)
I(t) = I0 * exp(-gamma * t)
       + sigma * E0 / (gamma - sigma) * (exp(-sigma * t) - exp(-gamma * t))
```

Here `R0` in the SIR solution is the initial recovered count, not the
reproduction number.
