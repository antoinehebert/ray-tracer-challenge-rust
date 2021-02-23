pub fn is_almost_equal(a: f32, b: f32) -> bool {
    (a - b).abs() < f32::EPSILON
}
