use crate::{color::*, matrix::Matrix, shape::Shape, tuple::Tuple};

#[derive(Debug, PartialEq, Clone)]
pub struct Stripe {
    a: Color,
    b: Color,
    transform: Matrix<4>,
}

impl Stripe {
    pub fn new(a: Color, b: Color) -> Self {
        Self {
            a,
            b,
            transform: Matrix::<4>::identity(),
        }
    }

    pub fn stripe_at(&self, point: &Tuple) -> Color {
        if point.x.floor() % 2.0 == 0.0 {
            self.a
        } else {
            self.b
        }
    }

    pub fn stripe_at_object(&self, object: &Shape, world_point: &Tuple) -> Color {
        let object_point = object.transform.inverse().expect("should be invertible") * *world_point;
        let pattern_point = self.transform.inverse().expect("should be invertible") * object_point;

        self.stripe_at(&pattern_point)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transformations::*;

    #[test]
    fn creating_a_stripe_pattern() {
        let pattern = Stripe::new(WHITE, BLACK);
        assert_eq!(pattern.a, WHITE);
        assert_eq!(pattern.b, BLACK);
    }

    #[test]
    fn a_stripe_pattern_is_constant_in_y() {
        let pattern = Stripe::new(WHITE, BLACK);

        assert_eq!(pattern.stripe_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.stripe_at(&Tuple::point(0.0, 1.0, 0.0)), WHITE);
        assert_eq!(pattern.stripe_at(&Tuple::point(0.0, 2.0, 0.0)), WHITE);
    }

    #[test]
    fn a_stripe_pattern_is_constant_in_z() {
        let pattern = Stripe::new(WHITE, BLACK);

        assert_eq!(pattern.stripe_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.stripe_at(&Tuple::point(0.0, 0.0, 1.0)), WHITE);
        assert_eq!(pattern.stripe_at(&Tuple::point(0.0, 0.0, 2.0)), WHITE);
    }

    #[test]
    fn a_stripe_pattern_alternates_in_x() {
        let pattern = Stripe::new(WHITE, BLACK);

        assert_eq!(pattern.stripe_at(&Tuple::point(0.0, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.stripe_at(&Tuple::point(0.9, 0.0, 0.0)), WHITE);
        assert_eq!(pattern.stripe_at(&Tuple::point(1.0, 0.0, 0.0)), BLACK);
        assert_eq!(pattern.stripe_at(&Tuple::point(-0.1, 0.0, 0.0)), BLACK);
        assert_eq!(pattern.stripe_at(&Tuple::point(-1.0, 0.0, 0.0)), BLACK);
        assert_eq!(pattern.stripe_at(&Tuple::point(-1.1, 0.0, 0.0)), WHITE);
    }

    #[test]
    fn stripes_with_an_object_transformation() {
        let mut object = Shape::sphere();
        object.transform = scaling(2.0, 2.0, 2.0);
        let pattern = Stripe::new(WHITE, BLACK);
        let c = pattern.stripe_at_object(&object, &Tuple::point(1.5, 0.0, 0.0));
        assert_eq!(c, WHITE);
    }

    #[test]
    fn stripes_with_a_pattern_transformation() {
        let object = Shape::sphere();
        let mut pattern = Stripe::new(WHITE, BLACK);
        pattern.transform = scaling(2.0, 2.0, 2.0);
        let c = pattern.stripe_at_object(&object, &Tuple::point(1.5, 0.0, 0.0));
        assert_eq!(c, WHITE);
    }

    #[test]
    fn stripes_with_both_an_object_and_a_pattern_transformation() {
        let mut object = Shape::sphere();
        object.transform = scaling(2.0, 2.0, 2.0);
        let mut pattern = Stripe::new(WHITE, BLACK);
        pattern.transform = translation(0.5, 0.0, 0.0);
        let c = pattern.stripe_at_object(&object, &Tuple::point(2.5, 0.0, 0.0));
        assert_eq!(c, WHITE);
    }

    // Scenario: The default pattern transformation
    //   Given pattern â† test_pattern()
    //   Then pattern.transform = identity_matrix

    // Scenario: Assigning a transformation
    //   Given pattern â† test_pattern()
    //   When set_pattern_transform(pattern, translation(1.0, 2.0, 3.0))
    //   Then pattern.transform = translation(1.0, 2.0, 3.0)

    // Scenario: A pattern with an object transformation
    //   Given shape â† sphere()
    //     And set_transform(shape, scaling(2.0, 2.0, 2.0))
    //     And pattern â† test_pattern()
    //   When c â† pattern_at_shape(pattern, shape, point(2.0, 3.0, 4.0))
    //   Then c = color(1.0, 1.5, 2.0)

    // Scenario: A pattern with a pattern transformation
    //   Given shape â† sphere()
    //     And pattern â† test_pattern()
    //     And set_pattern_transform(pattern, scaling(2.0, 2.0, 2.0))
    //   When c â† pattern_at_shape(pattern, shape, point(2.0, 3.0, 4.0))
    //   Then c = color(1.0, 1.5, 2.0)

    // Scenario: A pattern with both an object and a pattern transformation
    //   Given shape â† sphere()
    //     And set_transform(shape, scaling(2.0, 2.0, 2.0))
    //     And pattern â† test_pattern()
    //     And set_pattern_transform(pattern, translation(0.5, 1.0, 1.5))
    //   When c â† pattern_at_shape(pattern, shape, point(2.5, 3.0, 3.5))
    //   Then c = color(0.75, 0.5, 0.25)

    // Scenario: A gradient linearly interpolates between colors
    //   Given pattern â† gradient_pattern(white, black)
    //   Then pattern_at(pattern, point(0.0, 0.0, 0.0)) = white
    //     And pattern_at(pattern, point(0.25, 0.0, 0.0)) = color(0.75, 0.75, 0.75)
    //     And pattern_at(pattern, point(0.5, 0.0, 0.0)) = color(0.5, 0.5, 0.5)
    //     And pattern_at(pattern, point(0.75, 0.0, 0.0)) = color(0.25, 0.25, 0.25)

    // Scenario: A ring should extend in both x and z
    //   Given pattern â† ring_pattern(white, black)
    //   Then pattern_at(pattern, point(0.0, 0.0, 0.0)) = white
    //     And pattern_at(pattern, point(1.0, 0.0, 0.0)) = black
    //     And pattern_at(pattern, point(0.0, 0.0, 1.0)) = black
    //     # 0.708 = just slightly more than âˆš2/2.0
    //     And pattern_at(pattern, point(0.708, 0.0, 0.708)) = black

    // Scenario: Checkers should repeat in x
    //   Given pattern â† checkers_pattern(white, black)
    //   Then pattern_at(pattern, point(0.0, 0.0, 0.0)) = white
    //     And pattern_at(pattern, point(0.99, 0.0, 0.0)) = white
    //     And pattern_at(pattern, point(1.01, 0.0, 0.0)) = black

    // Scenario: Checkers should repeat in y
    //   Given pattern â† checkers_pattern(white, black)
    //   Then pattern_at(pattern, point(0.0, 0.0, 0.0)) = white
    //     And pattern_at(pattern, point(0.0, 0.99, 0.0)) = white
    //     And pattern_at(pattern, point(0.0, 1.01, 0.0)) = black

    // Scenario: Checkers should repeat in z
    //   Given pattern â† checkers_pattern(white, black)
    //   Then pattern_at(pattern, point(0.0, 0.0, 0.0)) = white
    //     And pattern_at(pattern, point(0.0, 0.0, 0.99)) = white
    // And pattern_at(pattern, point(0.0, 0.0, 1.01)) = black
}
