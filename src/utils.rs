// Margin of error when comparing floats.
pub const EPSILON: f32 = 0.00001;

pub fn is_almost_equal(a: f32, b: f32) -> bool {
    (a - b).abs() < EPSILON
}
