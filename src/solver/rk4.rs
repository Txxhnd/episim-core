use crate::error::SimulationError;

/// Advances `state` by one classical fourth-order Runge-Kutta step of size `step`.
///
/// Every stage input, stage slope, and the combined result is checked for
/// finiteness. Negative values are returned unchanged; judging them is the
/// caller's responsibility.
pub(crate) fn rk4_step<const N: usize, F>(
    state: [f64; N],
    step: f64,
    derivative: F,
) -> Result<[f64; N], SimulationError>
where
    F: Fn([f64; N]) -> [f64; N],
{
    if !step.is_finite() {
        return Err(SimulationError::InvalidTimeStep);
    }
    let slope_at = |stage, input| finite(stage, derivative(finite(stage, input)?));
    let offset =
        |slope: [f64; N], factor: f64| std::array::from_fn(|i| state[i] + factor * slope[i]);

    let k1 = slope_at("k1", state)?;
    let k2 = slope_at("k2", offset(k1, step / 2.0))?;
    let k3 = slope_at("k3", offset(k2, step / 2.0))?;
    let k4 = slope_at("k4", offset(k3, step))?;
    finite(
        "final",
        std::array::from_fn(|i| {
            state[i] + step / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        }),
    )
}

fn finite<const N: usize>(
    stage: &'static str,
    values: [f64; N],
) -> Result<[f64; N], SimulationError> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(values)
    } else {
        Err(SimulationError::NonFiniteIntermediate { stage })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Integrates dy/dt = -y from y(0) = 1 to t = 1 and returns the absolute error.
    fn decay_error(steps: usize) -> f64 {
        let step = 1.0 / steps as f64;
        let mut state = [1.0];
        for _ in 0..steps {
            state = rk4_step(state, step, |[y]| [-y]).unwrap();
        }
        (state[0] - (-1.0_f64).exp()).abs()
    }

    #[test]
    fn rk4_integrates_scalar_exponential_decay_with_fourth_order_accuracy() {
        let coarse = decay_error(10);
        let medium = decay_error(20);
        let fine = decay_error(40);
        assert!(coarse < 1e-6, "coarse error {coarse}");
        // Halving the step divides a fourth-order error by about 2^4 = 16.
        for ratio in [coarse / medium, medium / fine] {
            assert!((14.0..18.0).contains(&ratio), "error ratio {ratio}");
        }
    }

    #[test]
    fn rk4_integrates_constant_derivative_exactly_within_roundoff() {
        let [next] = rk4_step([1.0], 0.5, |_| [2.0]).unwrap();
        assert!((next - 2.0).abs() <= 4.0 * f64::EPSILON, "next {next}");
    }

    #[test]
    fn rk4_supports_three_and_four_component_arrays() {
        let three = rk4_step([1.0, 2.0, 3.0], 2.0, |_| [1.0, 0.0, -1.0]).unwrap();
        assert_eq!(three, [3.0, 2.0, 1.0]);

        let four = rk4_step([1.0, 2.0, 3.0, 4.0], 2.0, |_| [1.0, 0.0, -1.0, 0.5]).unwrap();
        assert_eq!(four, [3.0, 2.0, 1.0, 5.0]);
    }

    #[test]
    fn rk4_rejects_non_finite_step() {
        for step in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                rk4_step([1.0], step, |[y]| [-y]),
                Err(SimulationError::InvalidTimeStep),
                "step {step}"
            );
        }
    }

    #[test]
    fn rk4_rejects_non_finite_intermediate_stage() {
        let stage = |stage| Err(SimulationError::NonFiniteIntermediate { stage });
        // With y(0) = 0, a unit slope and a unit step, the stages are
        // evaluated at y = 0, 0.5, 0.5, and 1.
        let nan_above = |limit: f64| move |[y]: [f64; 1]| [if y > limit { f64::NAN } else { 1.0 }];

        assert_eq!(rk4_step([0.0], 1.0, |_| [f64::NAN]), stage("k1"));
        assert_eq!(rk4_step([0.0], 1.0, nan_above(0.25)), stage("k2"));
        assert_eq!(rk4_step([0.0], 1.0, nan_above(0.75)), stage("k4"));
        // A slope of 2 at y = 0.5 moves the third stage, and only it, to y = 1.
        let nan_at_one = |[y]: [f64; 1]| match y {
            0.0 => [1.0],
            0.5 => [2.0],
            _ => [f64::NAN],
        };
        assert_eq!(rk4_step([0.0], 1.0, nan_at_one), stage("k3"));
        // A finite slope can still overflow the state the next stage is evaluated at.
        assert_eq!(rk4_step([f64::MAX], 2.0, |_| [f64::MAX]), stage("k2"));
    }

    #[test]
    fn rk4_rejects_non_finite_final_state() {
        // Every stage is finite, but the weighted sum of the slopes overflows.
        assert_eq!(
            rk4_step([0.0], 1.0, |_| [f64::MAX / 4.0]),
            Err(SimulationError::NonFiniteIntermediate { stage: "final" })
        );
    }
}
