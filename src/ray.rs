use crate::tuple::Tuple;
use crate::{matrix::Matrix, transformations::*};

pub struct Ray {
    pub origin: Tuple,
    pub direction: Tuple,
}

impl Ray {
    pub fn new(origin: Tuple, direction: Tuple) -> Self {
        Self { origin, direction }
    }

    fn position(&self, time: f32) -> Tuple {
        self.origin + self.direction * time
    }

    pub fn transform(&self, matrix: Matrix<4>) -> Self {
        Ray {
            origin: matrix * self.origin,
            direction: matrix * self.direction,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creating_and_querying_a_ray() {
        let origin = Tuple::point(1., 2., 3.);
        let direction = Tuple::vector(4., 5., 6.);

        let r = Ray::new(origin, direction);

        assert_eq!(r.origin, origin);
        assert_eq!(r.direction, direction);
    }

    #[test]
    fn computing_a_point_from_a_distance() {
        let r = Ray::new(Tuple::point(2., 3., 4.), Tuple::vector(1., 0., 0.));
        assert_eq!(r.position(0.), Tuple::point(2., 3., 4.));
        assert_eq!(r.position(1.), Tuple::point(3., 3., 4.));
        assert_eq!(r.position(-1.), Tuple::point(1., 3., 4.));
        assert_eq!(r.position(2.5), Tuple::point(4.5, 3., 4.));
    }

    #[test]
    fn translating_a_ray() {
        let r = Ray::new(Tuple::point(1., 2., 3.), Tuple::vector(0., 1., 0.));
        let m = translation(3., 4., 5.);
        let r2 = r.transform(m);
        assert_eq!(r2.origin, Tuple::point(4., 6., 8.));
        assert_eq!(r2.direction, Tuple::vector(0., 1., 0.));
    }

    #[test]
    fn scaling_a_ray() {
        let r = Ray::new(Tuple::point(1., 2., 3.), Tuple::vector(0., 1., 0.));
        let m = scaling(2., 3., 4.);
        let r2 = r.transform(m);
        assert_eq!(r2.origin, Tuple::point(2., 6., 12.));
        assert_eq!(r2.direction, Tuple::vector(0., 3., 0.));
    }
}
