#[macro_export]
macro_rules! assert_almost_eq {
    ($lhs: expr, $rhs: expr) => {
        assert!(is_almost_equal($lhs, $rhs));
    };
}
