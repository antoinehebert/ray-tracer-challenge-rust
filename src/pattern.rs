use crate::{color::*, matrix::Matrix, shape::Shape, tuple::Tuple};

#[derive(Debug, PartialEq, Copy, Clone)]
enum PatternKind {
    Stripe(Color, Color),
    Gradient(Color, Color),
    Ring(Color, Color), // Like a target...
    Checkers(Color, Color),

    // TODO: Yikes, test induced damage :(.
    Test,
}

#[derive(Debug, PartialEq, Clone)]
pub struct Pattern {
    pub transform: Matrix<4>,
    kind: PatternKind,
}

impl Pattern {
    pub fn stripe(a: Color, b: Color) -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            kind: PatternKind::Stripe(a, b),
        }
    }

    pub fn gradient(from: Color, to: Color) -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            kind: PatternKind::Gradient(from, to),
        }
    }

    pub fn ring(a: Color, b: Color) -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            kind: PatternKind::Ring(a, b),
        }
    }

    pub fn checkers(a: Color, b: Color) -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            kind: PatternKind::Checkers(a, b),
        }
    }

    // TODO: Yikes! Test induced damage.
    pub fn test_pattern() -> Pattern {
        Pattern {
            transform: Matrix::<4>::identity(),
            kind: PatternKind::Test,
        }
    }

    fn color_at(&self, point: &Tuple) -> Color {
        match self.kind {
            PatternKind::Stripe(a, b) => {
                if point.x.floor() % 2.0 == 0.0 {
                    a
                } else {
                    b
                }
            }
            PatternKind::Gradient(from, to) => from + (to - from) * (point.x - point.x.floor()),
            PatternKind::Ring(a, b) => {
                if (point.x.powf(2.0) + point.z.powf(2.0)).sqrt().floor() % 2.0 == 0.0 {
                    a
                } else {
                    b
                }
            }
            PatternKind::Checkers(a, b) => {
                if (point.x.floor() + point.y.floor() + point.z.floor()) % 2.0 == 0.0 {
                    a
                } else {
                    b
                }
            }
            // Used to test that we're properly transforming world coordinates into local ones.
            PatternKind::Test => Color::new(point.x, point.y, point.z),
        }
    }

    // TODO: move on Shape?
    pub fn color_at_shape(&self, object: &Shape, world_point: &Tuple) -> Color {
        let object_point =
            object.transform().inverse().expect("should be invertible") * *world_point;
        let pattern_point = self.transform.inverse().expect("should be invertible") * object_point;

        self.color_at(&pattern_point)
    }
}

// Note: Some of these tests are assuming we're using some form of an abstract class... We're using a mix of struc
// and enum here so it doesn't really fit this idiom really well...
#[cfg(test)]
mod tests {
    use super::*;
    use crate::transformations::*;

    #[test]
    fn creating_a_stripe_pattern() {
        let pattern = Pattern::stripe(WHITE, BLACK);

        if let PatternKind::Stripe(a, b) = pattern.kind {
            assert_eq!(a, WHITE);
            assert_eq!(b, BLACK);
        } else {
            assert!(false); // You should not be here!
        }
    }

    #[test]
    fn a_stripe_pattern_is_constant_in_y() {
        let pattern = Pattern::stripe(WHITE, BLACK);

        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 1.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 2.0, 0.0)), WHITE);
    }

    #[test]
    fn a_stripe_pattern_is_constant_in_z() {
        let pattern = Pattern::stripe(WHITE, BLACK);

        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 1.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 2.0)), WHITE);
    }

    #[test]
    fn a_stripe_pattern_alternates_in_x() {
        let pattern = Pattern::stripe(WHITE, BLACK);

        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.9, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(1.0, 0.0, 0.0)), BLACK);
        assert_eq!(pattern.color_at(&Tuple::point(-0.1, 0.0, 0.0)), BLACK);
        assert_eq!(pattern.color_at(&Tuple::point(-1.0, 0.0, 0.0)), BLACK);
        assert_eq!(pattern.color_at(&Tuple::point(-1.1, 0.0, 0.0)), WHITE);
    }

    #[test]
    fn stripes_with_an_object_transformation() {
        let mut object = Shape::sphere();
        object.set_transform(scaling(2.0, 2.0, 2.0));
        let pattern = Pattern::stripe(WHITE, BLACK);
        let c = pattern.color_at_shape(&object, &Tuple::point(1.5, 0.0, 0.0));
        assert_eq!(c, WHITE);
    }

    #[test]
    fn stripes_with_a_pattern_transformation() {
        let object = Shape::sphere();
        let mut pattern = Pattern::stripe(WHITE, BLACK);
        pattern.transform = scaling(2.0, 2.0, 2.0);
        let c = pattern.color_at_shape(&object, &Tuple::point(1.5, 0.0, 0.0));
        assert_eq!(c, WHITE);
    }

    #[test]
    fn stripes_with_both_an_object_and_a_pattern_transformation() {
        let mut object = Shape::sphere();
        object.set_transform(scaling(2.0, 2.0, 2.0));
        let mut pattern = Pattern::stripe(WHITE, BLACK);
        pattern.transform = translation(0.5, 0.0, 0.0);
        let c = pattern.color_at_shape(&object, &Tuple::point(2.5, 0.0, 0.0));
        assert_eq!(c, WHITE);
    }

    #[test]
    fn the_default_pattern_transformation() {
        let pattern = Pattern::stripe(WHITE, BLACK);
        assert_eq!(pattern.transform, Matrix::<4>::identity());
    }

    #[test]
    fn assigning_a_transformation() {
        let mut pattern = Pattern::test_pattern();
        pattern.transform = translation(1.0, 2.0, 3.0);

        assert_eq!(pattern.transform, translation(1.0, 2.0, 3.0));
    }

    #[test]
    fn a_pattern_with_an_object_transformation() {
        let mut shape = Shape::sphere();
        shape.set_transform(scaling(2.0, 2.0, 2.0));
        let pattern = Pattern::test_pattern();

        let c = pattern.color_at_shape(&shape, &Tuple::point(2.0, 3.0, 4.0));

        assert_eq!(c, Color::new(1.0, 1.5, 2.0));
    }

    #[test]
    fn a_pattern_with_a_pattern_transformation() {
        let shape = Shape::sphere();
        let mut pattern = Pattern::test_pattern();
        pattern.transform = scaling(2.0, 2.0, 2.0);

        let c = pattern.color_at_shape(&shape, &Tuple::point(2.0, 3.0, 4.0));

        assert_eq!(c, Color::new(1.0, 1.5, 2.0));
    }

    #[test]
    fn a_pattern_with_both_an_object_and_a_pattern_transformation() {
        let mut shape = Shape::sphere();
        shape.set_transform(scaling(2.0, 2.0, 2.0));
        let mut pattern = Pattern::test_pattern();
        pattern.transform = translation(0.5, 1.0, 1.5);

        let c = pattern.color_at_shape(&shape, &Tuple::point(2.5, 3.0, 3.5));

        assert_eq!(c, Color::new(0.75, 0.5, 0.25));
    }

    #[test]
    fn a_gradient_linearly_interpolates_between_colors() {
        let pattern = Pattern::gradient(WHITE, BLACK);

        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(
            pattern.color_at(&Tuple::point(0.25, 0.0, 0.0)),
            Color::new(0.75, 0.75, 0.75)
        );
        assert_eq!(
            pattern.color_at(&Tuple::point(0.5, 0.0, 0.0)),
            Color::new(0.5, 0.5, 0.5)
        );
        assert_eq!(
            pattern.color_at(&Tuple::point(0.75, 0.0, 0.0)),
            Color::new(0.25, 0.25, 0.25)
        );
    }

    #[test]
    fn a_ring_should_extend_in_both_x_and_z() {
        let pattern = Pattern::ring(WHITE, BLACK);

        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(1.0, 0.0, 0.0)), BLACK);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 1.0)), BLACK);
        // 0.708 = just slightly more than √2/2.0
        assert_eq!(pattern.color_at(&Tuple::point(0.708, 0.0, 0.708)), BLACK);
    }

    #[test]
    fn checkers_should_repeat_in_x() {
        let pattern = Pattern::checkers(WHITE, BLACK);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.99, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(1.01, 0.0, 0.0)), BLACK);
    }

    #[test]
    fn checkers_should_repeat_in_y() {
        let pattern = Pattern::checkers(WHITE, BLACK);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.99, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 1.01, 0.0)), BLACK);
    }
    #[test]
    fn checkers_should_repeat_in_z() {
        let pattern = Pattern::checkers(WHITE, BLACK);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 0.99)), WHITE);
        assert_eq!(pattern.color_at(&Tuple::point(0.0, 0.0, 1.01)), BLACK);
    }
}
