//! Baseline workloads: SIR and SEIR at 10,000, 100,000, and 1,000,000
//! integration steps. Run with `cargo bench`.

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use episim_core::{
    SeirParameters, SeirState, SimulationSettings, SirParameters, SirState, simulate_seir,
    simulate_sir,
};

/// (integration steps, timestep in days). Each timestep is a power of two,
/// so the duration is exact and the step count unambiguous, and each run
/// covers 120 to 200 simulated days, so every size exercises an active
/// outbreak. Stretching the horizon instead would mostly measure arithmetic on
/// subnormal values long after the outbreak has ended.
const WORKLOADS: [(usize, f64); 3] = [
    (10_000, 1.0 / 64.0),
    (100_000, 1.0 / 512.0),
    (1_000_000, 1.0 / 8192.0),
];

fn settings(steps: usize, step_days: f64) -> SimulationSettings {
    SimulationSettings::new(steps as f64 * step_days, step_days).unwrap()
}

fn bench_sir(criterion: &mut Criterion) {
    let parameters = SirParameters {
        beta_per_day: 0.3,
        gamma_per_day: 0.1,
    };
    let initial = SirState {
        susceptible: 990.0,
        infected: 10.0,
        recovered: 0.0,
    };
    let mut group = criterion.benchmark_group("sir");
    group.sample_size(10);
    for (steps, step_days) in WORKLOADS {
        let settings = settings(steps, step_days);
        let result = simulate_sir(parameters, initial, settings).unwrap();
        assert_eq!(result.diagnostics.steps_executed, steps);
        group.bench_function(BenchmarkId::from_parameter(steps), |bencher| {
            bencher.iter(|| simulate_sir(black_box(parameters), black_box(initial), settings))
        });
    }
    group.finish();
}

fn bench_seir(criterion: &mut Criterion) {
    let parameters = SeirParameters {
        beta_per_day: 0.5,
        sigma_per_day: 0.2,
        gamma_per_day: 0.1,
    };
    let initial = SeirState {
        susceptible: 980.0,
        exposed: 10.0,
        infected: 10.0,
        recovered: 0.0,
    };
    let mut group = criterion.benchmark_group("seir");
    group.sample_size(10);
    for (steps, step_days) in WORKLOADS {
        let settings = settings(steps, step_days);
        let result = simulate_seir(parameters, initial, settings).unwrap();
        assert_eq!(result.diagnostics.steps_executed, steps);
        group.bench_function(BenchmarkId::from_parameter(steps), |bencher| {
            bencher.iter(|| simulate_seir(black_box(parameters), black_box(initial), settings))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_sir, bench_seir);
criterion_main!(benches);
