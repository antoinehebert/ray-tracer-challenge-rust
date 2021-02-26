use crate::utils::*;

use std::ops;
#[derive(Debug)]
pub struct Tuple {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32,
}

impl Tuple {
    pub fn new(x: f32, y: f32, z: f32, w: f32) -> Self {
        Tuple { x, y, z, w }
    }

    pub fn zero() -> Self {
        Tuple {
            x: 0.,
            y: 0.,
            z: 0.,
            w: 0.,
        }
    }

    fn is_point(&self) -> bool {
        self.w == 1.0
    }

    fn is_vector(&self) -> bool {
        self.w == 0.0
    }

    fn point(x: f32, y: f32, z: f32) -> Tuple {
        Tuple::new(x, y, z, 1.0)
    }

    fn vector(x: f32, y: f32, z: f32) -> Tuple {
        Tuple::new(x, y, z, 0.0)
    }

    fn magnitude(&self) -> f32 {
        assert!(self.is_vector());

        // .w really?!
        (self.x.powf(2.) + self.y.powf(2.) + self.z.powf(2.) + self.w.powf(2.)).sqrt()
    }

    fn normalize(&self) -> Self {
        assert!(self.is_vector());

        let magnitude = self.magnitude();
        // .w really?!
        Self::new(
            self.x / magnitude,
            self.y / magnitude,
            self.z / magnitude,
            self.w / magnitude,
        )
    }

    fn dot(&self, other: &Self) -> f32 {
        assert!(self.is_vector());

        // .w really?
        self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w
    }

    fn cross(&self, other: &Self) -> Self {
        assert!(self.is_vector() && other.is_vector());

        Tuple::vector(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
}

impl PartialEq for Tuple {
    fn eq(&self, other: &Self) -> bool {
        is_almost_equal(self.x, other.x)
            && is_almost_equal(self.y, other.y)
            && is_almost_equal(self.z, other.z)
            && is_almost_equal(self.w, other.w)
    }
}

impl ops::Add<Tuple> for Tuple {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(
            self.x + rhs.x,
            self.y + rhs.y,
            self.z + rhs.z,
            self.w + rhs.w,
        )
    }
}

impl ops::Sub<Tuple> for Tuple {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(
            self.x - rhs.x,
            self.y - rhs.y,
            self.z - rhs.z,
            self.w - rhs.w,
        )
    }
}

impl ops::Neg for Tuple {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z, -self.w)
    }
}

impl ops::Mul<f32> for Tuple {
    // The multiplication of rational numbers is a closed operation.
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs, self.w * rhs)
    }
}

impl ops::Div<f32> for Tuple {
    // The multiplication of rational numbers is a closed operation.
    type Output = Self;

    fn div(self, rhs: f32) -> Self {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs, self.w / rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_almost_eq;

    #[test]
    fn a_tuple_with_w_eq_1_is_a_point() {
        let t = Tuple::new(4.3, -4.2, 3.1, 1.0);

        assert_eq!(t.x, 4.3);
        assert_eq!(t.y, -4.2);
        assert_eq!(t.z, 3.1);
        assert_eq!(t.w, 1.0);

        assert!(t.is_point());
        assert!(!t.is_vector());
    }

    #[test]
    fn tuple_zero() {
        let t = Tuple::zero();

        assert_eq!(t.x, 0.);
        assert_eq!(t.y, 0.);
        assert_eq!(t.z, 0.);
        assert_eq!(t.w, 0.);
    }

    #[test]
    fn a_tuple_with_w_eq_0_is_a_vector() {
        let t = Tuple::new(4.3, -4.2, 3.1, 0.0);

        assert_eq!(t.x, 4.3);
        assert_eq!(t.y, -4.2);
        assert_eq!(t.z, 3.1);
        assert_eq!(t.w, 0.0);

        assert!(!t.is_point());
        assert!(t.is_vector());
    }

    #[test]
    fn point_creates_tuples_with_w_eq_1() {
        let point = Tuple::point(4.0, -4.0, 3.0);
        assert_eq!(point, Tuple::new(4.0, -4.0, 3.0, 1.0));
    }

    #[test]
    fn vector_creates_tuples_with_w_0() {
        let vector = Tuple::vector(4.0, -4.0, 3.0);
        assert_eq!(vector, Tuple::new(4.0, -4.0, 3.0, 0.0));
    }

    #[test]
    fn adding_two_tuples() {
        let a = Tuple::new(3.0, -2.0, 5.0, 1.0);
        let b = Tuple::new(-2.0, 3.0, 1.0, 0.0);

        assert_eq!(a + b, Tuple::new(1.0, 1.0, 6.0, 1.0));
    }

    #[test]
    fn subtracting_two_points() {
        let a = Tuple::point(3., 2., 1.);
        let b = Tuple::point(5., 6., 7.);

        assert_eq!(a - b, Tuple::vector(-2., -4., -6.));
    }

    #[test]
    fn subtracting_vector_from_point() {
        let a = Tuple::point(3., 2., 1.);
        let b = Tuple::vector(5., 6., 7.);

        assert_eq!(a - b, Tuple::point(-2., -4., -6.));
    }
    #[test]
    fn subtracting_vector_from_vector() {
        let a = Tuple::vector(3., 2., 1.);
        let b = Tuple::vector(5., 6., 7.);

        assert_eq!(a - b, Tuple::vector(-2., -4., -6.));
    }

    #[test]
    fn subtracting_a_vector_from_the_zero_vector() {
        let zero = Tuple::vector(0., 0., 0.);
        let v = Tuple::vector(1., -2., 3.);
        assert_eq!(zero - v, Tuple::vector(-1., 2., -3.));
    }

    #[test]
    fn negating_a_tuple() {
        let a = Tuple::new(1., -2., 3., -4.);

        assert_eq!(-a, Tuple::new(-1., 2., -3., 4.));
    }
    #[test]
    fn multiplying_a_tuple_by_a_scalar() {
        let a = Tuple::new(1., -2., 3., -4.);

        assert_eq!(a * 3.5, Tuple::new(3.5, -7., 10.5, -14.));
    }

    #[test]
    fn multiplying_a_tuple_by_a_fraction() {
        let a = Tuple::new(1., -2., 3., -4.);

        assert_eq!(a * 0.5, Tuple::new(0.5, -1., 1.5, -2.));
    }

    #[test]
    fn dividing_a_tuple_by_a_scalar() {
        let a = Tuple::new(1., -2., 3., -4.);

        assert_eq!(a / 2., Tuple::new(0.5, -1., 1.5, -2.));
    }
    #[test]
    fn computing_the_magnitude_of_vector_1_0_0() {
        let v = Tuple::vector(1., 0., 0.);
        assert_eq!(v.magnitude(), 1.);
    }
    #[test]
    fn computing_the_magnitude_of_vector_0_1_0() {
        let v = Tuple::vector(0., 1., 0.);
        assert_eq!(v.magnitude(), 1.);
    }
    #[test]
    fn computing_the_magnitude_of_vector_0_0_1() {
        let v = Tuple::vector(0., 0., 1.);
        assert_eq!(v.magnitude(), 1.);
    }
    #[test]
    fn computing_the_magnitude_of_vector_1_2_3() {
        let v = Tuple::vector(1., 2., 3.);
        assert_eq!(v.magnitude(), (14. as f32).sqrt());
    }
    #[test]
    fn computing_the_magnitude_of_neg_vector_1_2_3() {
        let v = Tuple::vector(-1., -2., -3.);
        assert_eq!(v.magnitude(), (14. as f32).sqrt());
    }

    #[test]
    fn normalizing_vector_4_0_0_gives_1_0_0() {
        let v = Tuple::vector(4., 0., 0.);
        assert_eq!(v.normalize(), Tuple::vector(1., 0., 0.));
    }

    #[test]
    fn normalizing_vector_1_2_3() {
        let v = Tuple::vector(1., 2., 3.);
        let norm = v.normalize();

        assert_eq!(norm, Tuple::vector(0.26726124, 0.5345225, 0.8017837));
        assert_almost_eq!(norm.magnitude(), 1.);
    }

    #[test]
    fn the_magnitude_of_a_normalized_vector() {
        let v = Tuple::vector(1., 2., 3.);
        let norm = v.normalize();

        assert_almost_eq!(norm.magnitude(), 1.);
    }

    #[test]
    fn the_dot_product_of_two_tuples() {
        let a = Tuple::vector(1., 2., 3.);
        let b = Tuple::vector(2., 3., 4.);

        assert_almost_eq!(a.dot(&b), 20.);
    }

    #[test]
    fn the_cross_product_of_two_vectors() {
        let a = Tuple::vector(1., 2., 3.);
        let b = Tuple::vector(2., 3., 4.);

        assert_eq!(a.cross(&b), Tuple::vector(-1., 2., -1.));
        assert_eq!(b.cross(&a), Tuple::vector(1., -2., 1.));
    }
}
