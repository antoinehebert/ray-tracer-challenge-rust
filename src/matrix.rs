use crate::assert_almost_eq;
use crate::tuple::*;
use crate::utils::*;
use std::ops;

// TODO: Use const generics when 1.51 is released, to dry the matrix implementation and only keep one!

////////////////////////////////////////////////////////////////////////////////
// Matrix4
////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Copy, Clone)]
struct Matrix<const SIZE: usize> {
    values: [[f32; SIZE]; SIZE],
}

impl<const SIZE: usize> Matrix<SIZE> {
    fn new(values: [[f32; SIZE]; SIZE]) -> Self {
        Self { values }
    }

    fn zero() -> Self {
        Self::new([[0.; SIZE]; SIZE])
    }

    fn identity() -> Self {
        let mut result = Self::zero();

        for v in 0..SIZE {
            result[v][v] = 1.;
        }

        result
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

    fn determinant(&self) -> f32 {
        let mut result = 0.;
        if SIZE == 2 {
            result = self[0][0] * self[1][1] - self[0][1] * self[1][0];
        } else {
            for column in 0..SIZE {
                result += self.values[0][column] * self.cofactor(0, column);
            }
        }

        result
    }

    // FIXME: This is a hack because SIZE-1 doesn't work with const generics yet...
    fn submatrix3(&self, row: usize, col: usize) -> Matrix<3> {
        assert!(SIZE == 4, "Called submatrix3 on a different size than 4!");
        let mut result = Matrix::<3>::zero();

        let mut res_x = 0;
        let mut res_y = 0;

        for x in 0..SIZE {
            if x == row {
                continue;
            }

            for y in 0..SIZE {
                if y == col {
                    continue;
                }

                result[res_x][res_y] = self[x][y];

                res_y += 1;
            }

            res_x += 1;
            res_y = 0;
        }

        result
    }

    // FIXME: This is a hack because SIZE-1 doesn't work with const generics yet...
    fn submatrix2(&self, row: usize, col: usize) -> Matrix<2> {
        assert!(SIZE == 3, "Called submatrix2 on a different size than 3!");

        let mut result = Matrix::<2>::zero();

        let mut res_x = 0;
        let mut res_y = 0;

        for x in 0..SIZE {
            if x == row {
                continue;
            }

            for y in 0..SIZE {
                if y == col {
                    continue;
                }

                result[res_x][res_y] = self[x][y];

                res_y += 1;
            }

            res_x += 1;
            res_y = 0;
        }

        result
    }

    fn minor(&self, row: usize, col: usize) -> f32 {
        if SIZE == 4 {
            self.submatrix3(row, col).determinant()
        } else if SIZE == 3 {
            self.submatrix2(row, col).determinant()
        } else {
            panic!(
                "Unsupported SIZE={} used in Matrix::minor, supported values are: 3, 4.",
                SIZE
            )
        }
    }

    // Using RES_SIZE here instead of SIZE-1 because it doesn't work yet...
    fn cofactor(&self, row: usize, col: usize) -> f32 {
        let result = self.minor(row, col);
        if (row + col) % 2 == 0 {
            result
        } else {
            -result
        }
    }
}

impl<const SIZE: usize> ops::Index<usize> for Matrix<SIZE> {
    type Output = [f32; SIZE];

    fn index(&self, idx: usize) -> &Self::Output {
        &self.values[idx]
    }
}

impl<const SIZE: usize> ops::IndexMut<usize> for Matrix<SIZE> {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.values[idx]
    }
}

impl<const SIZE: usize> PartialEq for Matrix<SIZE> {
    fn eq(&self, other: &Self) -> bool {
        for y in 0..SIZE {
            for x in 0..SIZE {
                if !is_almost_equal(self[x][y], other[x][y]) {
                    return false;
                }
            }
        }
        true
    }
}

impl<const SIZE: usize> ops::Mul for Matrix<SIZE> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = Self::zero();

        for row in 0..SIZE {
            for col in 0..SIZE {
                let mut val = 0.;
                for n in 0..SIZE {
                    val += self[row][n] * rhs[n][col];
                }
                result[row][col] = val;
            }
        }

        result
    }
}

impl ops::Mul<Tuple> for Matrix<4> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructing_and_inspecting_a_4x4_matrix() {
        let m = Matrix::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.5, 6.5, 7.5, 8.5],
            [9.0, 10.0, 11.0, 12.0],
            [13.5, 14.5, 15.5, 16.5],
        ]);

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
        let m = Matrix::new([[-3., 5.], [1., -2.]]);
        assert_eq!(m[0][0], -3.);
        assert_eq!(m[0][1], 5.);
        assert_eq!(m[1][0], 1.);
        assert_eq!(m[1][1], -2.);
    }

    #[test]
    fn a_3x3_matrix_ought_to_be_representable() {
        let m = Matrix::new([[-3., 5., 0.], [1., -2., -7.], [0., 1., 1.]]);

        assert_eq!(m[0][0], -3.);
        assert_eq!(m[1][1], -2.);
        assert_eq!(m[2][2], 1.);
    }

    #[test]
    fn matrix_equality_with_identical_matrices() {
        let a = Matrix::new([
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        ]);
        let b = Matrix::new([
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        ]);
        assert!(a == b);
    }

    #[test]
    fn matrix_equality_with_different_matrices() {
        let a = Matrix::new([
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        ]);
        let b = Matrix::new([
            [2., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 1.],
        ]);
        assert!(a != b);
    }

    #[test]
    fn multiplying_two_matrices() {
        let a = Matrix::new([
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        ]);
        let b = Matrix::new([
            [-2., 1., 2., 3.],
            [3., 2., 1., -1.],
            [4., 3., 6., 5.],
            [1., 2., 7., 8.],
        ]);

        let expectation = Matrix::new([
            [20., 22., 50., 48.],
            [44., 54., 114., 108.],
            [40., 58., 110., 102.],
            [16., 26., 46., 42.],
        ]);

        assert_eq!(a * b, expectation);
    }

    #[test]
    fn a_matrix_multiplied_by_a_tuple() {
        let a = Matrix::new([
            [1., 2., 3., 4.],
            [2., 4., 4., 2.],
            [8., 6., 4., 1.],
            [0., 0., 0., 1.],
        ]);
        let b = Tuple::new(1., 2., 3., 1.);
        assert_eq!(a * b, Tuple::new(18., 24., 33., 1.));
    }

    #[test]
    fn multiplying_a_matrix_by_the_identity_matrix() {
        let a = Matrix::new([
            [0., 1., 2., 4.],
            [1., 2., 4., 8.],
            [2., 4., 8., 16.],
            [4., 8., 16., 32.],
        ]);

        assert_eq!(a * Matrix::identity(), a);
    }

    #[test]
    fn multiplying_the_identity_matrix_by_a_tuple() {
        let a = Tuple::new(1., 2., 3., 4.);

        assert_eq!(Matrix::identity() * a, a);
    }

    #[test]
    fn transposing_a_matrix() {
        let a = Matrix::new([
            [0., 9., 3., 0.],
            [9., 8., 0., 8.],
            [1., 8., 5., 3.],
            [0., 0., 5., 8.],
        ]);

        let expected = Matrix::new([
            [0., 9., 1., 0.],
            [9., 8., 8., 0.],
            [3., 0., 5., 5.],
            [0., 8., 3., 8.],
        ]);

        assert_eq!(a.transpose(), expected);
    }

    #[test]
    fn transposing_the_identity_matrix() {
        assert_eq!(Matrix::identity().transpose(), Matrix::<4>::identity());
    }

    #[test]
    fn calculating_the_determinant_of_a_2x2_matrix() {
        let a = Matrix::new([[1., 5.], [-3., 2.]]);
        assert_almost_eq!(a.determinant(), 17.);
    }

    #[test]
    fn a_submatrix_of_a_3x3_matrix_is_a_2x2_matrix() {
        let a = Matrix::new([[1., 5., 0.], [-3., 2., 7.], [0., 6., -3.]]);

        let expectation = Matrix::new([[-3., 2.], [0., 6.]]);

        assert_eq!(a.submatrix2(0, 2), expectation);
    }

    #[test]
    fn a_submatrix_of_a_4x4_matrix_is_a_3x3_matrix() {
        let a = Matrix::new([
            [-6., 1., 1., 6.],
            [-8., 5., 8., 6.],
            [-1., 0., 8., 2.],
            [-7., 1., -1., 1.],
        ]);

        let expectation = Matrix::new([[-6., 1., 6.], [-8., 8., 6.], [-7., -1., 1.]]);
        assert_eq!(a.submatrix3(2, 1), expectation);
    }

    #[test]
    fn calculating_a_minor_of_a_3x3_matrix() {
        let a = Matrix::new([[3., 5., 0.], [2., -1., -7.], [6., -1., 5.]]);
        let b = a.submatrix2(1, 0);

        assert_almost_eq!(b.determinant(), 25.);
        assert_almost_eq!(a.minor(1, 0), 25.);
    }

    #[test]
    fn calculating_a_cofactor_of_a_3x3_matrix() {
        let a = Matrix::new([[3., 5., 0.], [2., -1., -7.], [6., -1., 5.]]);

        assert_almost_eq!(a.minor(0, 0), -12.);
        assert_almost_eq!(a.cofactor(0, 0), -12.);
        assert_almost_eq!(a.minor(1, 0), 25.);
        assert_almost_eq!(a.cofactor(1, 0), -25.);
    }

    #[test]
    fn calculating_the_determinant_of_a_3x3_matrix() {
        let a = Matrix::<3>::new([[1., 2., 6.], [-5., 8., -4.], [2., 6., 4.]]);
        assert_almost_eq!(a.cofactor(0, 0), 56.);
        assert_almost_eq!(a.cofactor(0, 1), 12.);
        assert_almost_eq!(a.cofactor(0, 2), -46.);
        assert_almost_eq!(a.determinant(), -196.);
    }

    #[test]
    fn calculating_the_determinant_of_a_4x4_matrix() {
        let a = Matrix::<4>::new([
            [-2., -8., 3., 5.],
            [-3., 1., 7., 3.],
            [1., 2., -9., 6.],
            [-6., 7., 7., -9.],
        ]);

        assert_almost_eq!(a.cofactor(0, 0), 690.);
        assert_almost_eq!(a.cofactor(0, 1), 447.);
        assert_almost_eq!(a.cofactor(0, 2), 210.);
        assert_almost_eq!(a.cofactor(0, 3), 51.);
        assert_almost_eq!(a.determinant(), -4071.);
    }

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
