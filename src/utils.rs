// Margin of error when comparing floats.
pub const EPSILON: f64 = 0.00001;

pub fn is_almost_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < EPSILON
}
