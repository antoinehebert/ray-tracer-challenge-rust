use crate::matrix::Matrix;
use crate::tuple::Tuple;

fn translation(x: f32, y: f32, z: f32) -> Matrix<4> {
    let mut result = Matrix::<4>::identity();
    result[0][3] = x;
    result[1][3] = y;
    result[2][3] = z;

    result
}

fn scaling(x: f32, y: f32, z: f32) -> Matrix<4> {
    let mut result = Matrix::<4>::identity();

    result[0][0] = x;
    result[1][1] = y;
    result[2][2] = z;

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multiplying_by_a_translation_matrix() {
        let transform = translation(5., -3., 2.);
        let p = Tuple::point(-3., 4., 5.);

        assert_eq!(transform * p, Tuple::point(2., 1., 7.));
    }

    #[test]
    fn multiplying_by_the_inverse_of_a_translation_matrix() {
        let transform = translation(5., -3., 2.);
        let inv = transform.inverse().expect("invertible");
        let p = Tuple::point(-3., 4., 5.);

        assert_eq!(inv * p, Tuple::point(-8., 7., 3.));
    }

    #[test]
    fn translation_does_not_affect_vectors() {
        let transform = translation(5., -3., 2.);
        let v = Tuple::vector(-3., 4., 5.);
        assert_eq!(transform * v, v);
    }

    #[test]
    fn a_scaling_matrix_applied_to_a_point() {
        let transform = scaling(2., 3., 4.);
        let p = Tuple::point(-4., 6., 8.);
        assert_eq!(transform * p, Tuple::point(-8., 18., 32.));
    }

    #[test]
    fn a_scaling_matrix_applied_to_a_vector() {
        let transform = scaling(2., 3., 4.);
        let v = Tuple::vector(-4., 6., 8.);
        assert_eq!(transform * v, Tuple::vector(-8., 18., 32.));
    }

    #[test]
    fn multiplying_by_the_inverse_of_a_scaling_matrix() {
        let transform = scaling(2., 3., 4.);
        let inv = transform.inverse().expect("invertible");
        let v = Tuple::vector(-4., 6., 8.);
        assert_eq!(inv * v, Tuple::vector(-2., 2., 2.));
    }

    #[test]
    fn reflection_is_scaling_by_a_negative_value() {
        let transform = scaling(-1., 1., 1.);
        let p = Tuple::point(2., 3., 4.);
        assert_eq!(transform * p, Tuple::point(-2., 3., 4.));
    }

    // # Scenario: Rotating a point around the x axis
    // #   Given p â† point(0, 1, 0)
    // #     And half_quarter â† rotation_x(Ï€ / 4)
    // #     And full_quarter â† rotation_x(Ï€ / 2)
    // #   Then half_quarter * p = point(0, âˆš2/2, âˆš2/2)
    // #     And full_quarter * p = point(0, 0, 1)

    // # Scenario: The inverse of an x-rotation rotates in the opposite direction
    // #   Given p â† point(0, 1, 0)
    // #     And half_quarter â† rotation_x(Ï€ / 4)
    // #     And inv â† inverse(half_quarter)
    // #   Then inv * p = point(0, âˆš2/2, -âˆš2/2)

    // # Scenario: Rotating a point around the y axis
    // #   Given p â† point(0, 0, 1)
    // #     And half_quarter â† rotation_y(Ï€ / 4)
    // #     And full_quarter â† rotation_y(Ï€ / 2)
    // #   Then half_quarter * p = point(âˆš2/2, 0, âˆš2/2)
    // #     And full_quarter * p = point(1, 0, 0)

    // # Scenario: Rotating a point around the z axis
    // #   Given p â† point(0, 1, 0)
    // #     And half_quarter â† rotation_z(Ï€ / 4)
    // #     And full_quarter â† rotation_z(Ï€ / 2)
    // #   Then half_quarter * p = point(-âˆš2/2, âˆš2/2, 0)
    // #     And full_quarter * p = point(-1, 0, 0)

    // # Scenario: A shearing transformation moves x in proportion to y
    // #   Given transform â† shearing(1, 0, 0, 0, 0, 0)
    // #     And p â† point(2, 3, 4)
    // #   Then transform * p = point(5, 3, 4)

    // # Scenario: A shearing transformation moves x in proportion to z
    // #   Given transform â† shearing(0, 1, 0, 0, 0, 0)
    // #     And p â† point(2, 3, 4)
    // #   Then transform * p = point(6, 3, 4)

    // # Scenario: A shearing transformation moves y in proportion to x
    // #   Given transform â† shearing(0, 0, 1, 0, 0, 0)
    // #     And p â† point(2, 3, 4)
    // #   Then transform * p = point(2, 5, 4)

    // # Scenario: A shearing transformation moves y in proportion to z
    // #   Given transform â† shearing(0, 0, 0, 1, 0, 0)
    // #     And p â† point(2, 3, 4)
    // #   Then transform * p = point(2, 7, 4)

    // # Scenario: A shearing transformation moves z in proportion to x
    // #   Given transform â† shearing(0, 0, 0, 0, 1, 0)
    // #     And p â† point(2, 3, 4)
    // #   Then transform * p = point(2, 3, 6)

    // # Scenario: A shearing transformation moves z in proportion to y
    // #   Given transform â† shearing(0, 0, 0, 0, 0, 1)
    // #     And p â† point(2, 3, 4)
    // #   Then transform * p = point(2, 3, 7)

    // # Scenario: Individual transformations are applied in sequence
    // #   Given p â† point(1, 0, 1)
    // #     And A â† rotation_x(Ï€ / 2)
    // #     And B â† scaling(5, 5, 5)
    // #     And C â† translation(10, 5, 7)
    // #   # apply rotation first
    // #   When p2 â† A * p
    // #   Then p2 = point(1, -1, 0)
    // #   # then apply scaling
    // #   When p3 â† B * p2
    // #   Then p3 = point(5, -5, 0)
    // #   # then apply translation
    // #   When p4 â† C * p3
    // #   Then p4 = point(15, 0, 7)

    // # Scenario: Chained transformations must be applied in reverse order
    // #   Given p â† point(1, 0, 1)
    // #     And A â† rotation_x(Ï€ / 2)
    // #     And B â† scaling(5, 5, 5)
    // #     And C â† translation(10, 5, 7)
    // #   When T â† C * B * A
    // #   Then T * p = point(15, 0, 7)

    // # Scenario: The transformation matrix for the default orientation
    // #   Given from â† point(0, 0, 0)
    // #     And to â† point(0, 0, -1)
    // #     And up â† vector(0, 1, 0)
    // #   When t â† view_transform(from, to, up)
    // #   Then t = identity_matrix

    // # Scenario: A view transformation matrix looking in positive z direction
    // #   Given from â† point(0, 0, 0)
    // #     And to â† point(0, 0, 1)
    // #     And up â† vector(0, 1, 0)
    // #   When t â† view_transform(from, to, up)
    // #   Then t = scaling(-1, 1, -1)

    // # Scenario: The view transformation moves the world
    // #   Given from â† point(0, 0, 8)
    // #     And to â† point(0, 0, 0)
    // #     And up â† vector(0, 1, 0)
    // #   When t â† view_transform(from, to, up)
    // #   Then t = translation(0, 0, -8)

    // # Scenario: An arbitrary view transformation
    // #   Given from â† point(1, 3, 2)
    // #     And to â† point(4, -2, 8)
    // #     And up â† vector(1, 1, 0)
    // #   When t â† view_transform(from, to, up)
    // #   Then t is the following 4x4 matrix:
    // #       | -0.50709 | 0.50709 |  0.67612 | -2.36643 |
    // #       |  0.76772 | 0.60609 |  0.12122 | -2.82843 |
    // #       | -0.35857 | 0.59761 | -0.71714 |  0.00000 |
    // #       |  0.00000 | 0.00000 |  0.00000 |  1.00000 |
}
