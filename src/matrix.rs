// use crate::assert_almost_eq;
use crate::tuple::*;
use crate::utils::*;
use std::ops;

////////////////////////////////////////////////////////////////////////////////
// Matrix4
////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Copy, Clone)]
struct Matrix4 {
    values: [[f32; 4]; 4],
}

impl Matrix4 {
    fn new(v0: [f32; 4], v1: [f32; 4], v2: [f32; 4], v3: [f32; 4]) -> Self {
        Self {
            values: [
                [v0[0], v0[1], v0[2], v0[3]],
                [v1[0], v1[1], v1[2], v1[3]],
                [v2[0], v2[1], v2[2], v2[3]],
                [v3[0], v3[1], v3[2], v3[3]],
            ],
        }
    }

    fn zero() -> Self {
        Matrix4::new(
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
        )
    }

    fn identity() -> Self {
        Matrix4::new(
            [1., 0., 0., 0.],
            [0., 1., 0., 0.],
            [0., 0., 1., 0.],
            [0., 0., 0., 1.],
        )
    }

    fn transpose(&self) -> Self {
        let mut result = Self::zero();

        for row in 0..4 {
            for col in 0..4 {
                result[col][row] = self[row][col]
            }
        }

        result
    }
}

impl ops::Index<usize> for Matrix4 {
    type Output = [f32; 4];

    fn index(&self, idx: usize) -> &Self::Output {
        &self.values[idx]
    }
}

impl ops::IndexMut<usize> for Matrix4 {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.values[idx]
    }
}

impl PartialEq for Matrix4 {
    fn eq(&self, other: &Self) -> bool {
        for y in 0..4 {
            for x in 0..4 {
                if !is_almost_equal(self[x][y], other[x][y]) {
                    return false;
                }
            }
        }
        true
    }
}

impl ops::Mul for Matrix4 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = Matrix4::zero();

        for row in 0..4 {
            for col in 0..4 {
                result[row][col] = self[row][0] * rhs[0][col]
                    + self[row][1] * rhs[1][col]
                    + self[row][2] * rhs[2][col]
                    + self[row][3] * rhs[3][col];
            }
        }

        result
    }
}

impl ops::Mul<Tuple> for Matrix4 {
    type Output = Tuple;

    fn mul(self, rhs: Tuple) -> Self::Output {
        let mut result = Tuple::zero();

        result.x =
            self[0][0] * rhs.x + self[0][1] * rhs.y + self[0][2] * rhs.z + self[0][3] * rhs.w;

        result.y =
            self[1][0] * rhs.x + self[1][1] * rhs.y + self[1][2] * rhs.z + self[1][3] * rhs.w;

        result.z =
            self[2][0] * rhs.x + self[2][1] * rhs.y + self[2][2] * rhs.z + self[2][3] * rhs.w;

        result.w =
            self[3][0] * rhs.x + self[3][1] * rhs.y + self[3][2] * rhs.z + self[3][3] * rhs.w;

        result
    }
}

////////////////////////////////////////////////////////////////////////////////
// Matrix2
////////////////////////////////////////////////////////////////////////////////

struct Matrix2 {
    values: [[f32; 2]; 2],
}

impl Matrix2 {
    fn new(v0: [f32; 2], v1: [f32; 2]) -> Self {
        Self {
            values: [[v0[0], v0[1]], [v1[0], v1[1]]],
        }
    }
}

impl ops::Index<usize> for Matrix2 {
    type Output = [f32; 2];

    fn index(&self, idx: usize) -> &Self::Output {
        &self.values[idx]
    }
}

////////////////////////////////////////////////////////////////////////////////
// Matrix3
////////////////////////////////////////////////////////////////////////////////

struct Matrix3 {
    values: [[f32; 3]; 3],
}

impl Matrix3 {
    fn new(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> Self {
        Self {
            values: [
                [v0[0], v0[1], v0[2]],
                [v1[0], v1[1], v1[2]],
                [v2[0], v2[1], v2[2]],
            ],
        }
    }
}

impl ops::Index<usize> for Matrix3 {
    type Output = [f32; 3];

    fn index(&self, idx: usize) -> &Self::Output {
        &self.values[idx]
    }
}

////////////////////////////////////////////////////////////////////////////////
// Tests
////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructing_and_inspecting_a_4x4_matrix() {
        let m = Matrix4::new(
            [1.0, 2.0, 3.0, 4.0],
            [5.5, 6.5, 7.5, 8.5],
            [9.0, 10.0, 11.0, 12.0],
            [13.5, 14.5, 15.5, 16.5],
        );

        assert_eq!(m[0][0], 1.);
        assert_eq!(m[0][3], 4.);
        assert_eq!(m[1][0], 5.5);
        assert_eq!(m[1][2], 7.5);
        assert_eq!(m[2][2], 11.);
        assert_eq!(m[3][0], 13.5);
        assert_eq!(m[3][2], 15.5);
    }

    #[test]
    fn a_2x2_matrix_ought_to_be_representable() {
        let m = Matrix2::new([-3., 5.], [1., -2.]);
        assert_eq!(m[0][0], -3.);
        assert_eq!(m[0][1], 5.);
        assert_eq!(m[1][0], 1.);
        assert_eq!(m[1][1], -2.);
    }

    #[test]
    fn a_3x3_matrix_ought_to_be_representable() {
        let m = Matrix3::new([-3., 5., 0.], [1., -2., -7.], [0., 1., 1.]);

        assert_eq!(m[0][0], -3.);
        assert_eq!(m[1][1], -2.);
        assert_eq!(m[2][2], 1.);
    }

    #[test]
    fn matrix_equality_with_identical_matrices() {
        let a = Matrix4::new(
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        );
        let b = Matrix4::new(
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        );
        assert!(a == b);
    }

    #[test]
    fn matrix_equality_with_different_matrices() {
        let a = Matrix4::new(
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        );
        let b = Matrix4::new(
            [2., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 1.],
        );
        assert!(a != b);
    }

    #[test]
    fn multiplying_two_matrices() {
        let a = Matrix4::new(
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        );
        let b = Matrix4::new(
            [-2., 1., 2., 3.],
            [3., 2., 1., -1.],
            [4., 3., 6., 5.],
            [1., 2., 7., 8.],
        );

        let expectation = Matrix4::new(
            [20., 22., 50., 48.],
            [44., 54., 114., 108.],
            [40., 58., 110., 102.],
            [16., 26., 46., 42.],
        );

        assert_eq!(a * b, expectation);
    }

    #[test]
    fn a_matrix_multiplied_by_a_tuple() {
        let a = Matrix4::new(
            [1., 2., 3., 4.],
            [2., 4., 4., 2.],
            [8., 6., 4., 1.],
            [0., 0., 0., 1.],
        );
        let b = Tuple::new(1., 2., 3., 1.);
        assert_eq!(a * b, Tuple::new(18., 24., 33., 1.));
    }

    #[test]
    fn multiplying_a_matrix_by_the_identity_matrix() {
        let a = Matrix4::new(
            [0., 1., 2., 4.],
            [1., 2., 4., 8.],
            [2., 4., 8., 16.],
            [4., 8., 16., 32.],
        );

        assert_eq!(a * Matrix4::identity(), a);
    }

    #[test]
    fn multiplying_the_identity_matrix_by_a_tuple() {
        let a = Tuple::new(1., 2., 3., 4.);

        assert_eq!(Matrix4::identity() * a, a);
    }

    #[test]
    fn transposing_a_matrix() {
        let a = Matrix4::new(
            [0., 9., 3., 0.],
            [9., 8., 0., 8.],
            [1., 8., 5., 3.],
            [0., 0., 5., 8.],
        );

        let expected = Matrix4::new(
            [0., 9., 1., 0.],
            [9., 8., 8., 0.],
            [3., 0., 5., 5.],
            [0., 8., 3., 8.],
        );

        assert_eq!(a.transpose(), expected);
    }

    #[test]
    fn transposing_the_identity_matrix() {
        assert_eq!(Matrix4::identity().transpose(), Matrix4::identity());
    }
    //
    // Scenario: Calculating the determinant of a 2x2 matrix
    //   Given the following 2x2 matrix A:
    //     |  1 | 5 |
    //     | -3 | 2 |
    //   Then determinant(A) = 17
    //
    // Scenario: A submatrix of a 3x3 matrix is a 2x2 matrix
    //   Given the following 3x3 matrix A:
    //     |  1 | 5 |  0 |
    //     | -3 | 2 |  7 |
    //     |  0 | 6 | -3 |
    //   Then submatrix(A, 0, 2) is the following 2x2 matrix:
    //     | -3 | 2 |
    //     |  0 | 6 |
    //
    // Scenario: A submatrix of a 4x4 matrix is a 3x3 matrix
    //   Given the following 4x4 matrix A:
    //     | -6 |  1 |  1 |  6 |
    //     | -8 |  5 |  8 |  6 |
    //     | -1 |  0 |  8 |  2 |
    //     | -7 |  1 | -1 |  1 |
    //   Then submatrix(A, 2, 1) is the following 3x3 matrix:
    //     | -6 |  1 | 6 |
    //     | -8 |  8 | 6 |
    //     | -7 | -1 | 1 |
    //
    // Scenario: Calculating a minor of a 3x3 matrix
    //   Given the following 3x3 matrix A:
    //       |  3 |  5 |  0 |
    //       |  2 | -1 | -7 |
    //       |  6 | -1 |  5 |
    //     And B â† submatrix(A, 1, 0)
    //   Then determinant(B) = 25
    //     And minor(A, 1, 0) = 25
    //
    // Scenario: Calculating a cofactor of a 3x3 matrix
    //   Given the following 3x3 matrix A:
    //       |  3 |  5 |  0 |
    //       |  2 | -1 | -7 |
    //       |  6 | -1 |  5 |
    //   Then minor(A, 0, 0) = -12
    //     And cofactor(A, 0, 0) = -12
    //     And minor(A, 1, 0) = 25
    //     And cofactor(A, 1, 0) = -25
    //
    // Scenario: Calculating the determinant of a 3x3 matrix
    //   Given the following 3x3 matrix A:
    //     |  1 |  2 |  6 |
    //     | -5 |  8 | -4 |
    //     |  2 |  6 |  4 |
    //   Then cofactor(A, 0, 0) = 56
    //     And cofactor(A, 0, 1) = 12
    //     And cofactor(A, 0, 2) = -46
    //     And determinant(A) = -196
    //
    // Scenario: Calculating the determinant of a 4x4 matrix
    //   Given the following 4x4 matrix A:
    //     | -2 | -8 |  3 |  5 |
    //     | -3 |  1 |  7 |  3 |
    //     |  1 |  2 | -9 |  6 |
    //     | -6 |  7 |  7 | -9 |
    //   Then cofactor(A, 0, 0) = 690
    //     And cofactor(A, 0, 1) = 447
    //     And cofactor(A, 0, 2) = 210
    //     And cofactor(A, 0, 3) = 51
    //     And determinant(A) = -4071
    //
    // Scenario: Testing an invertible matrix for invertibility
    //   Given the following 4x4 matrix A:
    //     |  6 |  4 |  4 |  4 |
    //     |  5 |  5 |  7 |  6 |
    //     |  4 | -9 |  3 | -7 |
    //     |  9 |  1 |  7 | -6 |
    //   Then determinant(A) = -2120
    //     And A is invertible
    //
    // Scenario: Testing a noninvertible matrix for invertibility
    //   Given the following 4x4 matrix A:
    //     | -4 |  2 | -2 | -3 |
    //     |  9 |  6 |  2 |  6 |
    //     |  0 | -5 |  1 | -5 |
    //     |  0 |  0 |  0 |  0 |
    //   Then determinant(A) = 0
    //     And A is not invertible
    //
    // Scenario: Calculating the inverse of a matrix
    //   Given the following 4x4 matrix A:
    //       | -5 |  2 |  6 | -8 |
    //       |  1 | -5 |  1 |  8 |
    //       |  7 |  7 | -6 | -7 |
    //       |  1 | -3 |  7 |  4 |
    //     And B â† inverse(A)
    //   Then determinant(A) = 532
    //     And cofactor(A, 2, 3) = -160
    //     And B[3,2] = -160/532
    //     And cofactor(A, 3, 2) = 105
    //     And B[2,3] = 105/532
    //     And B is the following 4x4 matrix:
    //       |  0.21805 |  0.45113 |  0.24060 | -0.04511 |
    //       | -0.80827 | -1.45677 | -0.44361 |  0.52068 |
    //       | -0.07895 | -0.22368 | -0.05263 |  0.19737 |
    //       | -0.52256 | -0.81391 | -0.30075 |  0.30639 |
    //
    // Scenario: Calculating the inverse of another matrix
    //   Given the following 4x4 matrix A:
    //     |  8 | -5 |  9 |  2 |
    //     |  7 |  5 |  6 |  1 |
    //     | -6 |  0 |  9 |  6 |
    //     | -3 |  0 | -9 | -4 |
    //   Then inverse(A) is the following 4x4 matrix:
    //     | -0.15385 | -0.15385 | -0.28205 | -0.53846 |
    //     | -0.07692 |  0.12308 |  0.02564 |  0.03077 |
    //     |  0.35897 |  0.35897 |  0.43590 |  0.92308 |
    //     | -0.69231 | -0.69231 | -0.76923 | -1.92308 |
    //
    // Scenario: Calculating the inverse of a third matrix
    //   Given the following 4x4 matrix A:
    //     |  9 |  3 |  0 |  9 |
    //     | -5 | -2 | -6 | -3 |
    //     | -4 |  9 |  6 |  4 |
    //     | -7 |  6 |  6 |  2 |
    //   Then inverse(A) is the following 4x4 matrix:
    //     | -0.04074 | -0.07778 |  0.14444 | -0.22222 |
    //     | -0.07778 |  0.03333 |  0.36667 | -0.33333 |
    //     | -0.02901 | -0.14630 | -0.10926 |  0.12963 |
    //     |  0.17778 |  0.06667 | -0.26667 |  0.33333 |
    //
    // Scenario: Multiplying a product by its inverse
    //   Given the following 4x4 matrix A:
    //       |  3 | -9 |  7 |  3 |
    //       |  3 | -8 |  2 | -9 |
    //       | -4 |  4 |  4 |  1 |
    //       | -6 |  5 | -1 |  1 |
    //     And the following 4x4 matrix B:
    //       |  8 |  2 |  2 |  2 |
    //       |  3 | -1 |  7 |  0 |
    //       |  7 |  0 |  5 |  4 |
    //       |  6 | -2 |  0 |  5 |
    //     And C â† A * B
    //   Then C * inverse(B) = A
}
