# Mathematical models

`episim-core` implements five deterministic compartmental models. This page
gives their equations, the units and valid ranges of every quantity, the
assumptions behind them, and the exact definition of each reported metric.

## Units

| Quantity | Unit |
|---|---|
| Time, duration, timestep | days |
| `beta_per_day`, `sigma_per_day`, `gamma_per_day`, `omega_per_day` | per day |
| Compartments, peak infectious, cumulative infections | individuals |
| `R0`, cumulative incidence fraction, relative drift | dimensionless |

Compartment values are `f64` and may be fractional. A deterministic model
tracks expected counts, not whole individuals.

## Shared assumptions

1. The population is closed: no births, deaths, or migration.
2. Individuals mix homogeneously.
3. Parameters are constant over the simulation.
4. In SIR and SEIR, recovered individuals are permanently immune. In SIRS
   and SEIRS their immunity wanes, and in SIS recovery gives no immunity.
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

### Cumulative infections

Every model also integrates the number of new infections since day zero,

```text
dC/dt = beta * S * I / N,    C(0) = 0
```

`C` is not a compartment: nobody moves into or out of it, and it is not part
of `N`. Each observation reports it as `cumulative_infections`. In SIR and
SEIR everyone who leaves `S` is infected, so `C(t) = S(0) - S(t)` exactly. In
SIS, SIRS and SEIRS people can be infected again; `C` counts each infection, so it
can exceed `N`.

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

## SIS

Compartments: susceptible `S` and infected (infectious) `I`. Recovered
individuals return to `S` at once, so nobody is ever immune and people can
be infected again and again.

```text
dS/dt = -beta * S * I / N + gamma * I
dI/dt =  beta * S * I / N - gamma * I
```

| Parameter | Meaning | Valid range |
|---|---|---|
| `beta_per_day` | transmission rate | finite, `>= 0` |
| `gamma_per_day` | recovery rate, back to susceptible | finite, `> 0` |

With `S = N - I` this is the logistic equation. With `r = beta - gamma`,
`I* = N (1 - 1/R0)` and `A = I* / I(0) - 1`, its solution is

```text
I(t) = I* / (1 + A * exp(-r t))
```

`I` moves monotonically towards `I*` and never crosses it. When `R0 > 1` the
infection becomes endemic at `I*`; when `R0 <= 1` it dies out. The formula
holds in both cases, with `I* < 0` when `R0 < 1`.

## SIRS

SIR with waning immunity: recovered individuals become susceptible again at
rate `omega`, after an average of `1 / omega` days. The compartments are
those of SIR, and the model uses `SirState`.

```text
dS/dt = -beta * S * I / N + omega * R
dI/dt =  beta * S * I / N - gamma * I
dR/dt =  gamma * I - omega * R
```

| Parameter | Meaning | Valid range |
|---|---|---|
| `beta_per_day` | transmission rate | finite, `>= 0` |
| `gamma_per_day` | recovery rate | finite, `> 0` |
| `omega_per_day` | rate at which immunity wanes | finite, `>= 0` |

With `omega = 0` it is SIR exactly. With `omega > 0` and `R0 > 1` the
infected rise and fall in damped waves and settle at the endemic
equilibrium

```text
S* = N / R0
I* = N * (1 - 1/R0) * omega / (omega + gamma)
R* = N - S* - I*
```

There is no closed form for the path to it. When `R0 <= 1` the infection
dies out.

## SEIRS

SEIR with waning immunity, in the same way as SIRS. The model uses
`SeirState`.

```text
dS/dt = -beta * S * I / N + omega * R
dE/dt =  beta * S * I / N - sigma * E
dI/dt =  sigma * E - gamma * I
dR/dt =  gamma * I - omega * R
```

The parameters are those of SEIR plus `omega_per_day`, finite and `>= 0`.
With `omega = 0` it is SEIR exactly. With `omega > 0` and `R0 > 1` it settles
at

```text
S* = N / R0
I* = N * (1 - 1/R0) / (1 + gamma/sigma + gamma/omega)
E* = (gamma / sigma) * I*
R* = (gamma / omega) * I*
```

## Initial state

Every compartment must be finite and non-negative. Their sum `N` must be
finite and strictly positive. The sum is checked separately because finite
compartments can still add up to infinity.

## Population conservation

In every model the derivatives of the compartments sum to zero, so `N(t) = N(0)` for the exact
solution. The library checks this after every integration step; see
[numerical methods](numerical_methods.md).

## Metrics

All metrics are model-derived quantities. None is an estimate from observed
case data.

| Metric | Definition |
|---|---|
| `basic_reproduction_number` | `R0 = beta / gamma`, for every model. The theoretical number of secondary infections from one case in a fully susceptible population. |
| `herd_immunity_threshold` | `1 - 1/R0`, or 0 when `R0 <= 1`: the fraction of the population that must be immune for the infected to decline. An `Option`, `Some` for SIR, SEIR, SIRS and SEIRS. |
| `peak_infectious` | The largest `I` among the recorded observations. |
| `peak_infectious_time_days` | The time of the earliest recorded observation within `64 * f64::EPSILON * max(N0, 1)` individuals of that peak. If `I` never rises above its initial value by more than that, this is day zero. The tolerance matters only where `I` levels off at its peak; there roundoff alone would decide which observation is largest. |
| `cumulative_infections_since_start` | `C(T)`, where `T` is the final simulation time; in SIR and SEIR this is `S(0) - S(T)`. Infections already present at day zero are not counted. |
| `cumulative_incidence_fraction` | `C(T) / N(0)`: new infections per head of the initial population. In SIR and SEIR it is also the fraction of people newly infected. In SIS, SIRS and SEIRS it counts reinfections and can exceed 1. |
| `final_size_fraction` | `(S(0) - S(inf)) / N(0)`: the fraction newly infected by the time the epidemic is over, from the final size equation below. An `Option`, `Some` for SIR and SEIR. |
| `endemic_prevalence` | `I* / N`, or 0 when `R0 <= 1`: the fraction infected once the model settles at its endemic equilibrium. It is `1 - 1/R0` for SIS and `(1 - 1/R0) * omega / (omega + gamma)` for SIRS, and `(1 - 1/R0) / (1 + gamma/sigma + gamma/omega)` for SEIRS. It follows from the parameters alone. An `Option`, `Some` for SIS, SIRS and SEIRS. |

The cumulative metrics refer to the final time `T` only. For any earlier
time `t`, read `cumulative_infections` from the observation at `t`. As `T` grows,
`cumulative_incidence_fraction` approaches `final_size_fraction`; comparing
the two shows whether the simulated period covers the whole epidemic.

SIS has no final size and no herd-immunity threshold, because nobody becomes
immune; both are `None`. Its long-run quantity is `endemic_prevalence`. In
SIS, `1 - 1/R0` is the endemic prevalence: at the equilibrium a fraction
`1/R0` is susceptible and the rest are infected. The infected rise towards
`I*` and level off there, so `peak_infectious_time_days` is the time they
come within roundoff of it.

SIRS and SEIRS have no final size either, so `final_size_fraction` is
`None`, also when `omega = 0`; use `simulate_sir` or `simulate_seir` for that
case. Its herd-immunity threshold
keeps the formula `1 - 1/R0`: it is the fraction of the population that must
be non-susceptible at a given time for infections to decline at that time.
With waning immunity that state is not maintained without continued
infection or vaccination. At the endemic equilibrium the non-susceptible
fraction is exactly `1 - 1/R0`. The peak time is that of the first wave,
which in every measured run was also the highest.

## Effective reproduction number

```text
effective R = R0 * S / N
```

`SirParameters::effective_reproduction_number`,
`SeirParameters::effective_reproduction_number`,
`SisParameters::effective_reproduction_number`,
`SirsParameters::effective_reproduction_number` and
`SeirsParameters::effective_reproduction_number` compute it for any state,
with `N` taken as the total of that state. It is the average number of
people one infected person infects while a fraction `S / N` of the
population is still susceptible. It sets the direction of the epidemic
exactly:

```text
SIR, SIS, SIRS:  dI/dt       = gamma * I * (effective R - 1)
SEIR, SEIRS:     d(E + I)/dt = gamma * I * (effective R - 1)
```

In SIS, SIRS and SEIRS the susceptibles grow again as people recover or lose
their immunity, so effective R does not keep falling. In SIRS and SEIRS it
can rise back above 1, which starts the next wave; in both it tends to 1 as the infection
settles at `I*`.

So in SIR the infected peak exactly when effective R falls to 1, which is
when the susceptible fraction falls to `1 / R0`. The herd-immunity threshold
`1 - 1/R0` is the same condition seen from the immune side.

This is a model-derived quantity. It is not the effective reproduction
number `Rt` that public-health agencies estimate from case data, which
depends on reporting, testing, and statistical assumptions this library does
not model.

## Final size

Dividing `dS/dt = -beta * S * I / N` by `dR/dt = gamma * I` gives
`dS/dR = -R0 * S / N`, in SIR and SEIR. Every exact solution therefore
satisfies

```text
S(t) = S(0) * exp(-R0 * (R(t) - R(0)) / N)
```

at all times. When the epidemic is over, the exposed and infected are gone,
so `R(inf) = N - S(inf)`, and the susceptibles left over solve the final size
equation

```text
S(inf) = S(0) * exp(-R0 * (S(0) + seed - S(inf)) / N)
```

where `seed` is the number infected (plus, for SEIR, exposed) at the start.
With a positive seed the equation has exactly one solution between 0 and
`S(0)`; the library finds it by bisection. With no seed, no transmission, or
no susceptibles, the final size is 0.

For a vanishing seed in a fully susceptible population this is the textbook
relation `z = 1 - exp(-R0 * z)` for the fraction infected `z`: about 80%
when `R0 = 2`, and about 94% when `R0 = 3`.

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
