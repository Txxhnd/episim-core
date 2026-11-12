use episim_core::{SimulationError, SimulationSettings};

#[test]
fn settings_accepts_positive_finite_duration_and_step() {
    let settings = SimulationSettings::new(10.0, 0.25).unwrap();
    assert_eq!(settings.duration_days(), 10.0);
    assert_eq!(settings.step_days(), 0.25);
}

#[test]
fn settings_rejects_zero_negative_nan_and_infinite_duration() {
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            SimulationSettings::new(duration, 1.0),
            Err(SimulationError::InvalidDuration),
            "duration {duration}"
        );
    }
}

#[test]
fn settings_rejects_zero_negative_nan_and_infinite_step() {
    for step in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            SimulationSettings::new(10.0, step),
            Err(SimulationError::InvalidTimeStep),
            "step {step}"
        );
    }
}

#[test]
fn settings_defaults_to_one_million_max_steps() {
    let settings = SimulationSettings::new(10.0, 1.0).unwrap();
    assert_eq!(settings.max_steps(), 1_000_000);
}

#[test]
fn settings_rejects_zero_max_steps() {
    let settings = SimulationSettings::new(10.0, 1.0).unwrap();
    assert_eq!(
        settings.with_max_steps(0),
        Err(SimulationError::InvalidMaximumSteps)
    );
}

#[test]
fn settings_accepts_explicit_positive_max_steps() {
    let settings = SimulationSettings::new(10.0, 0.25)
        .unwrap()
        .with_max_steps(5_000_000)
        .unwrap();
    assert_eq!(settings.max_steps(), 5_000_000);
    assert_eq!(settings.duration_days(), 10.0);
    assert_eq!(settings.step_days(), 0.25);
}
