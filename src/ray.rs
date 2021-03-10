use crate::tuple::Tuple;

struct Ray {
    origin: Tuple,
    direction: Tuple,
}

impl Ray {
    fn new(origin: Tuple, direction: Tuple) -> Self {
        Self { origin, direction }
    }

    fn position(&self, time: f32) -> Tuple {
        self.origin + self.direction * time
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

    // Scenario: Translating a ray
    //   Given r â† ray(point(1., 2., 3.), vector(0., 1., 0.))
    //     And m â† translation(3., 4., 5.)
    //   When r2. â† transform(r, m)
    //   Then r2..origin = point(4., 6., 8.)
    //     And r2..direction = vector(0., 1., 0.)

    // Scenario: Scaling a ray
    //   Given r â† ray(point(1., 2., 3.), vector(0., 1., 0.))
    //     And m â† scaling(2., 3., 4.)
    //   When r2. â† transform(r, m)
    //   Then r2..origin = point(2., 6., 12.)
    //     And r2.direction = vector(0., 3., 0.)
}
