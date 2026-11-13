use episim_core::{
    Observation, SeirParameters, SeirState, SeirsParameters, SimulationError, SimulationSettings,
    SirParameters, SirState, SirsParameters, SisParameters, SisState, TransmissionSchedule,
    estimate_sir_step_error, simulate_seir, simulate_seirs, simulate_sir, simulate_sirs,
    simulate_sis,
};

const SIR: SirParameters = SirParameters {
    beta_per_day: 0.3,
    gamma_per_day: 0.1,
};

const SIR0: SirState = SirState {
    susceptible: 990.0,
    infected: 10.0,
    recovered: 0.0,
};

/// Contacts at 40% from day 20, back to normal from day 50.
fn lockdown() -> TransmissionSchedule {
    TransmissionSchedule::new(&[(20.0, 0.4), (50.0, 1.0)]).unwrap()
}

fn settings(duration_days: f64, step_days: f64) -> SimulationSettings {
    SimulationSettings::new(duration_days, step_days).unwrap()
}

fn scheduled(
    duration_days: f64,
    step_days: f64,
    schedule: TransmissionSchedule,
) -> SimulationSettings {
    settings(duration_days, step_days)
        .with_transmission_schedule(schedule)
        .unwrap()
}

#[test]
fn schedule_rejects_invalid_changes() {
    for changes in [
        &[(0.0, 0.5)][..],
        &[(-1.0, 0.5)],
        &[(f64::NAN, 0.5)],
        &[(f64::INFINITY, 0.5)],
        &[(5.0, 0.5), (5.0, 0.8)],
        &[(5.0, 0.5), (3.0, 0.8)],
        &[(5.0, -0.1)],
        &[(5.0, f64::NAN)],
        &[(5.0, f64::INFINITY)],
    ] {
        assert_eq!(
            TransmissionSchedule::new(changes),
            Err(SimulationError::InvalidSchedule),
            "{changes:?}"
        );
    }
    // No changes at all, and a multiplier of 0 (no contacts), are valid.
    assert!(TransmissionSchedule::new(&[]).is_ok());
    assert!(TransmissionSchedule::new(&[(1.0, 0.0)]).is_ok());
}

#[test]
fn schedule_multiplier_is_one_before_the_first_change_and_holds_until_the_next() {
    let schedule = lockdown();
    for (day, multiplier) in [
        (0.0, 1.0),
        (19.99, 1.0),
        (20.0, 0.4),
        (49.99, 0.4),
        (50.0, 1.0),
        (1e6, 1.0),
    ] {
        assert_eq!(schedule.multiplier_at(day), multiplier, "day {day}");
    }
}

#[test]
fn settings_reject_changes_off_the_step_grid() {
    let off_grid = TransmissionSchedule::new(&[(20.0, 0.4), (20.3, 1.0)]).unwrap();
    assert_eq!(
        settings(100.0, 0.25).with_transmission_schedule(off_grid),
        Err(SimulationError::ScheduleChangeOffStepGrid { day: 20.3 })
    );
    // Decimal steps: 20 / 0.1 and 20.3 / 0.1 are not exact in binary, but
    // within the tolerance of whole numbers of steps.
    let decimal = TransmissionSchedule::new(&[(20.0, 0.4), (20.3, 1.0)]).unwrap();
    assert!(
        settings(100.0, 0.1)
            .with_transmission_schedule(decimal)
            .is_ok()
    );
    // A change after the end of the run is allowed and has no effect.
    let late = TransmissionSchedule::new(&[(500.0, 0.4)]).unwrap();
    let result = simulate_sir(SIR, SIR0, scheduled(100.0, 0.5, late)).unwrap();
    let plain = simulate_sir(SIR, SIR0, settings(100.0, 0.5)).unwrap();
    assert_eq!(result.trajectory, plain.trajectory);
}

#[test]
fn schedule_without_a_change_in_transmission_leaves_the_run_unchanged() {
    let plain = simulate_sir(SIR, SIR0, settings(160.0, 0.25)).unwrap();
    for schedule in [
        TransmissionSchedule::new(&[]).unwrap(),
        TransmissionSchedule::new(&[(30.0, 1.0)]).unwrap(),
    ] {
        let result = simulate_sir(SIR, SIR0, scheduled(160.0, 0.25, schedule)).unwrap();
        assert_eq!(result, plain);
    }
}

/// Largest difference in any compartment or the cumulative infections, at
/// the coarse run's times, from a run whose step is `stride` times smaller.
fn error_against(
    coarse: &[Observation<SirState>],
    fine: &[Observation<SirState>],
    stride: usize,
) -> f64 {
    coarse
        .iter()
        .enumerate()
        .flat_map(|(index, observation)| {
            let exact = fine[index * stride];
            let (state, reference) = (observation.state, exact.state);
            [
                (state.susceptible - reference.susceptible).abs(),
                (state.infected - reference.infected).abs(),
                (state.recovered - reference.recovered).abs(),
                (observation.cumulative_infections - exact.cumulative_infections).abs(),
            ]
        })
        .fold(0.0, f64::max)
}

#[test]
fn lockdown_keeps_fourth_order_accuracy() {
    // A step 1024 times finer than a day stands in for the exact solution;
    // the changes on days 20 and 50 lie on every grid used here.
    let run = |step| simulate_sir(SIR, SIR0, scheduled(160.0, step, lockdown())).unwrap();
    let fine = run(1.0 / 1024.0);
    let errors = [1usize, 2, 4].map(|per_day| {
        error_against(
            &run(1.0 / per_day as f64).trajectory,
            &fine.trajectory,
            1024 / per_day,
        )
    });
    for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
        assert!((12.0..20.0).contains(&ratio), "{errors:?}");
    }
    assert!(errors[2] < 1e-4, "{errors:?}");
}

#[test]
fn lockdown_lowers_the_epidemic_and_has_no_final_size() {
    let plain = simulate_sir(SIR, SIR0, settings(400.0, 1.0 / 16.0)).unwrap();
    let locked = simulate_sir(SIR, SIR0, scheduled(400.0, 1.0 / 16.0, lockdown())).unwrap();
    let (plain, locked) = (plain.metrics, locked.metrics);
    // Without the lockdown the run reaches the final size of R0 = 3.
    assert!(
        (plain.cumulative_incidence_fraction - 0.9312).abs() < 1e-3,
        "{plain:?}"
    );
    // With it, fewer are infected, and the infected peak when it starts.
    assert!(
        (locked.cumulative_incidence_fraction - 0.8171).abs() < 1e-3,
        "{locked:?}"
    );
    assert!(locked.peak_infectious < plain.peak_infectious - 50.0);
    assert_eq!(locked.peak_infectious_time_days, 20.0);
    // The final size equation assumes a constant beta; R0 and the threshold
    // are the values without the schedule.
    assert_eq!(locked.final_size_fraction, None);
    assert_eq!(
        locked.basic_reproduction_number,
        plain.basic_reproduction_number
    );
    assert_eq!(
        locked.herd_immunity_threshold,
        plain.herd_immunity_threshold
    );
}

#[test]
fn stability_limit_uses_the_largest_beta_in_force() {
    let fast = SirParameters {
        beta_per_day: 1.0,
        gamma_per_day: 0.1,
    };
    let run = |changes: &[(f64, f64)]| {
        let schedule = TransmissionSchedule::new(changes).unwrap();
        simulate_sir(fast, SIR0, scheduled(10.0, 1.0, schedule))
    };
    // Four times the contacts from day 5: beta * 4 = 4 per day, so the
    // largest accepted step is 2 / 4 = 0.5 days.
    assert_eq!(
        run(&[(5.0, 4.0)]).unwrap_err(),
        SimulationError::TimeStepTooLarge { max_step_days: 0.5 }
    );
    // A change at the end of the run never takes effect, and fewer contacts
    // never tighten the limit.
    assert!(run(&[(10.0, 4.0)]).is_ok());
    assert!(run(&[(5.0, 0.5)]).is_ok());

    // The shared pipeline applies the same rule to every model: SEIRS sums
    // its rates, 4 * 0.4 + 0.2 + 0.1 + 0.05 = 1.95, which a one-day step
    // accepts, while 4 * 0.5 + 0.35 = 2.35 it does not.
    let seirs = |beta_per_day| SeirsParameters {
        beta_per_day,
        sigma_per_day: 0.2,
        gamma_per_day: 0.1,
        omega_per_day: 0.05,
    };
    let seir0 = SeirState {
        susceptible: 980.0,
        exposed: 10.0,
        infected: 10.0,
        recovered: 0.0,
    };
    let four_times = || TransmissionSchedule::new(&[(5.0, 4.0)]).unwrap();
    assert!(simulate_seirs(seirs(0.4), seir0, scheduled(10.0, 1.0, four_times())).is_ok());
    assert!(matches!(
        simulate_seirs(seirs(0.5), seir0, scheduled(10.0, 1.0, four_times())),
        Err(SimulationError::TimeStepTooLarge { .. })
    ));
}

#[test]
fn endemic_prevalence_uses_the_multiplier_the_run_ends_under() {
    // R0 = 3, halved from day 10: the run settles at the endemic level of
    // R0 = 1.5, I* / N = 1 - 1/1.5 = 1/3.
    let parameters = SisParameters {
        beta_per_day: 0.3,
        gamma_per_day: 0.1,
    };
    let initial = SisState {
        susceptible: 990.0,
        infected: 10.0,
    };
    let halved = TransmissionSchedule::new(&[(10.0, 0.5)]).unwrap();
    let result = simulate_sis(parameters, initial, scheduled(800.0, 0.25, halved)).unwrap();
    let endemic = result.metrics.endemic_prevalence.unwrap();
    assert!((endemic - 1.0 / 3.0).abs() < 1e-12, "{endemic}");
    let last = result.trajectory.last().unwrap().state;
    assert!((last.infected / 1000.0 - endemic).abs() < 1e-6, "{last:?}");
}

#[test]
fn step_error_estimate_works_with_a_schedule() {
    // The half-step run keeps the schedule, which is on its grid too.
    let estimate = |step| {
        estimate_sir_step_error(SIR, SIR0, scheduled(160.0, step, lockdown()))
            .unwrap()
            .max_absolute_error
    };
    let halving = estimate(1.0) / estimate(0.5);
    assert!((12.0..20.0).contains(&halving), "{halving}");
}

#[test]
fn every_model_runs_a_lockdown_with_fewer_infections() {
    let run = |schedule: Option<TransmissionSchedule>| {
        let settings = match schedule {
            Some(schedule) => scheduled(160.0, 0.25, schedule),
            None => settings(160.0, 0.25),
        };
        let seir0 = SeirState {
            susceptible: 980.0,
            exposed: 10.0,
            infected: 10.0,
            recovered: 0.0,
        };
        [
            simulate_sir(SIR, SIR0, settings.clone()).unwrap().metrics,
            simulate_seir(
                SeirParameters {
                    beta_per_day: 0.5,
                    sigma_per_day: 0.2,
                    gamma_per_day: 0.1,
                },
                seir0,
                settings.clone(),
            )
            .unwrap()
            .metrics,
            simulate_sis(
                SisParameters {
                    beta_per_day: 0.3,
                    gamma_per_day: 0.1,
                },
                SisState {
                    susceptible: 990.0,
                    infected: 10.0,
                },
                settings.clone(),
            )
            .unwrap()
            .metrics,
            simulate_sirs(
                SirsParameters {
                    beta_per_day: 0.3,
                    gamma_per_day: 0.1,
                    omega_per_day: 0.01,
                },
                SIR0,
                settings.clone(),
            )
            .unwrap()
            .metrics,
            simulate_seirs(
                SeirsParameters {
                    beta_per_day: 0.5,
                    sigma_per_day: 0.2,
                    gamma_per_day: 0.1,
                    omega_per_day: 0.01,
                },
                seir0,
                settings,
            )
            .unwrap()
            .metrics,
        ]
    };
    for (model, (plain, locked)) in ["SIR", "SEIR", "SIS", "SIRS", "SEIRS"]
        .iter()
        .zip(run(None).iter().zip(run(Some(lockdown()))))
    {
        assert!(
            locked.cumulative_infections_since_start < plain.cumulative_infections_since_start,
            "{model}: {locked:?} vs {plain:?}"
        );
    }
}

#[test]
fn multiplier_times_effective_r_sets_the_direction_of_the_infected() {
    // dI/dt = gamma * I * (m(t) * R0 * S / N - 1): the infected rise exactly
    // while the scheduled effective R is above 1.
    let schedule = lockdown();
    let result = simulate_sir(SIR, SIR0, scheduled(120.0, 0.05, schedule.clone())).unwrap();
    for pair in result.trajectory.windows(2) {
        let r = schedule.multiplier_at(pair[0].time_days)
            * SIR.effective_reproduction_number(&pair[0].state).unwrap();
        if (r - 1.0).abs() > 1e-2 {
            let change = pair[1].state.infected - pair[0].state.infected;
            assert_eq!(change > 0.0, r > 1.0, "day {}: R {r}", pair[0].time_days);
        }
    }
}
