use crate::matrix::Matrix;
use crate::tuple::Tuple;

pub fn translation(x: f64, y: f64, z: f64) -> Matrix<4> {
    let mut result = Matrix::<4>::identity();
    result[0][3] = x;
    result[1][3] = y;
    result[2][3] = z;

    result
}

pub fn scaling(x: f64, y: f64, z: f64) -> Matrix<4> {
    let mut result = Matrix::<4>::identity();

    result[0][0] = x;
    result[1][1] = y;
    result[2][2] = z;

    result
}

pub fn rotation_x(rad: f64) -> Matrix<4> {
    let mut result = Matrix::<4>::identity();

    let cos = rad.cos();
    result[1][1] = cos;
    result[2][2] = cos;

    let sin = rad.sin();
    result[1][2] = -sin;
    result[2][1] = sin;

    result
}

pub fn rotation_y(rad: f64) -> Matrix<4> {
    let mut result = Matrix::<4>::identity();

    let cos = rad.cos();
    result[0][0] = cos;
    result[2][2] = cos;

    let sin = rad.sin();
    result[0][2] = sin;
    result[2][0] = -sin;

    result
}

pub fn rotation_z(rad: f64) -> Matrix<4> {
    let mut result = Matrix::<4>::identity();

    let cos = rad.cos();
    result[0][0] = cos;
    result[1][1] = cos;

    let sin = rad.sin();
    result[0][1] = -sin;
    result[1][0] = sin;

    result
}

pub fn shearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Matrix<4> {
    let mut result = Matrix::<4>::identity();

    result[0][1] = xy;
    result[0][2] = xz;

    result[1][0] = yx;
    result[1][2] = yz;

    result[2][0] = zx;
    result[2][1] = zy;

    result
}

pub fn view_transform(from: Tuple, to: Tuple, up: Tuple) -> Matrix<4> {
    let forward = (to - from).normalize();
    let upn = up.normalize();
    let left = forward.cross(&upn);
    let true_up = left.cross(&forward);

    let orientation = Matrix::<4>::new([
        [left.x, left.y, left.z, 0.],
        [true_up.x, true_up.y, true_up.z, 0.],
        [-forward.x, -forward.y, -forward.z, 0.],
        [0., 0., 0., 1.],
    ]);
    orientation * translation(-from.x, -from.y, -from.z)
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

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

    #[test]
    fn rotating_a_point_around_the_x_axis() {
        let p = Tuple::point(0., 1., 0.);
        let half_quarter = rotation_x(PI / 4.);
        let full_quarter = rotation_x(PI / 2.);

        assert_eq!(
            half_quarter * p,
            Tuple::point(0., (2. as f64).sqrt() / 2., (2. as f64).sqrt() / 2.)
        );
        assert_eq!(full_quarter * p, Tuple::point(0., 0., 1.));
    }

    #[test]
    fn the_inverse_of_an_x_rotation_rotates_in_the_opposite_direction() {
        let p = Tuple::point(0., 1., 0.);
        let half_quarter = rotation_x(PI / 4.);
        let inv = half_quarter.inverse().expect("invertible");
        assert_eq!(
            inv * p,
            Tuple::point(0., (2. as f64).sqrt() / 2., -(2. as f64).sqrt() / 2.)
        );
    }

    #[test]
    fn rotating_a_point_around_the_y_axis() {
        let p = Tuple::point(0., 0., 1.);
        let half_quarter = rotation_y(PI / 4.);
        let full_quarter = rotation_y(PI / 2.);

        assert_eq!(
            half_quarter * p,
            Tuple::point((2. as f64).sqrt() / 2., 0., (2. as f64).sqrt() / 2.)
        );
        assert_eq!(full_quarter * p, Tuple::point(1., 0., 0.));
    }

    #[test]
    fn rotating_a_point_around_the_z_axis() {
        let p = Tuple::point(0., 1., 0.);
        let half_quarter = rotation_z(PI / 4.);
        let full_quarter = rotation_z(PI / 2.);

        assert_eq!(
            half_quarter * p,
            Tuple::point(-(2. as f64).sqrt() / 2., (2. as f64).sqrt() / 2., 0.)
        );
        assert_eq!(full_quarter * p, Tuple::point(-1., 0., 0.));
    }

    #[test]
    fn a_shearing_transformation_moves_x_in_proportion_to_y() {
        let transform = shearing(1., 0., 0., 0., 0., 0.);
        let p = Tuple::point(2., 3., 4.);
        assert_eq!(transform * p, Tuple::point(5., 3., 4.));
    }

    #[test]
    fn a_shearing_transformation_moves_x_in_proportion_to_z() {
        let transform = shearing(0., 1., 0., 0., 0., 0.);
        let p = Tuple::point(2., 3., 4.);
        assert_eq!(transform * p, Tuple::point(6., 3., 4.));
    }

    #[test]
    fn a_shearing_transformation_moves_y_in_proportion_to_x() {
        let transform = shearing(0., 0., 1., 0., 0., 0.);
        let p = Tuple::point(2., 3., 4.);
        assert_eq!(transform * p, Tuple::point(2., 5., 4.));
    }

    #[test]
    fn a_shearing_transformation_moves_y_in_proportion_to_z() {
        let transform = shearing(0., 0., 0., 1., 0., 0.);
        let p = Tuple::point(2., 3., 4.);
        assert_eq!(transform * p, Tuple::point(2., 7., 4.));
    }

    #[test]
    fn a_shearing_transformation_moves_z_in_proportion_to_x() {
        let transform = shearing(0., 0., 0., 0., 1., 0.);
        let p = Tuple::point(2., 3., 4.);
        assert_eq!(transform * p, Tuple::point(2., 3., 6.));
    }

    #[test]
    fn a_shearing_transformation_moves_z_in_proportion_to_y() {
        let transform = shearing(0., 0., 0., 0., 0., 1.);
        let p = Tuple::point(2., 3., 4.);
        assert_eq!(transform * p, Tuple::point(2., 3., 7.));
    }

    #[test]
    fn individual_transformations_are_applied_in_sequence() {
        let p = Tuple::point(1., 0., 1.);
        let a = rotation_x(PI / 2.);
        let b = scaling(5., 5., 5.);
        let c = translation(10., 5., 7.);

        // apply rotation first
        let p2 = a * p;
        assert_eq!(p2, Tuple::point(1., -1., 0.));

        // then apply scaling
        let p3 = b * p2;
        assert_eq!(p3, Tuple::point(5., -5., 0.));

        // then apply translation
        let p4 = c * p3;
        assert_eq!(p4, Tuple::point(15., 0., 7.));
    }

    #[test]
    fn chained_transformations_must_be_applied_in_reverse_order() {
        let p = Tuple::point(1., 0., 1.);
        let a = rotation_x(PI / 2.);
        let b = scaling(5., 5., 5.);
        let c = translation(10., 5., 7.);

        let t = c * b * a;
        assert_eq!(t * p, Tuple::point(15., 0., 7.));
    }

    #[test]
    fn the_transformation_matrix_for_the_default_orientation() {
        let from = Tuple::point(0., 0., 0.);
        let to = Tuple::point(0., 0., -1.);
        let up = Tuple::vector(0., 1., 0.);
        let t = view_transform(from, to, up);
        assert_eq!(t, Matrix::<4>::identity());
    }

    #[test]
    fn a_view_transformation_matrix_looking_in_positive_z_direction() {
        let from = Tuple::point(0., 0., 0.);
        let to = Tuple::point(0., 0., 1.);
        let up = Tuple::vector(0., 1., 0.);
        let t = view_transform(from, to, up);
        assert_eq!(t, scaling(-1., 1., -1.));
    }

    #[test]
    fn the_view_transformation_moves_the_world() {
        let from = Tuple::point(0., 0., 8.);
        let to = Tuple::point(0., 0., 0.);
        let up = Tuple::vector(0., 1., 0.);
        let t = view_transform(from, to, up);
        assert_eq!(t, translation(0., 0., -8.));
    }

    #[test]
    fn an_arbitrary_view_transformation() {
        let from = Tuple::point(1., 3., 2.);
        let to = Tuple::point(4., -2., 8.);
        let up = Tuple::vector(1., 1., 0.);
        let t = view_transform(from, to, up);
        assert_eq!(
            t,
            Matrix::<4>::new([
                [-0.50709, 0.50709, 0.67612, -2.36643],
                [0.76772, 0.60609, 0.12122, -2.82843],
                [-0.35857, 0.59761, -0.71714, 0.00000],
                [0.00000, 0.00000, 0.00000, 1.00000]
            ])
        );
    }
}
