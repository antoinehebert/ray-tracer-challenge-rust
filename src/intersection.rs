use crate::ray::Ray;
use crate::sphere::*;
use crate::tuple::Tuple;
use crate::utils::*;
#[derive(Debug, Clone, PartialEq)]
pub struct Intersection<'a> {
    pub t: f64, // intersection "time"
    pub object: &'a Sphere,
}

impl<'a> Intersection<'a> {
    pub fn new(t: f64, object: &'a Sphere) -> Self {
        Self { t, object }
    }

    pub fn prepare_computations(&self, ray: &Ray) -> Computations {
        let point = ray.position(self.t);
        let eyev = -ray.direction;
        let mut normalv = self.object.normal_at(point);

        let inside = normalv.dot(&eyev) < 0.0;
        if inside {
            normalv = -normalv;
        }

        Computations {
            t: self.t,
            object: self.object,
            point,
            over_point: point + normalv * EPSILON,
            eyev,
            inside,
            normalv,
        }
    }
}

pub type Intersections<'a> = Vec<Intersection<'a>>;

pub struct Computations<'a> {
    pub t: f64,
    pub object: &'a Sphere,
    pub point: Tuple,
    pub over_point: Tuple, // Prevent acne effect.
    pub eyev: Tuple,
    pub inside: bool,
    pub normalv: Tuple,
}

// TODO: make this a method on Intersections.
pub fn hit<'a>(xs: &'a Intersections) -> Option<&'a Intersection<'a>> {
    xs.iter()
        .filter(|x| x.t >= 0.)
        .min_by(|x, y| x.t.partial_cmp(&y.t).unwrap_or(std::cmp::Ordering::Equal))
}

#[cfg(test)]
mod tests {
    use crate::{assert_almost_eq, transformations::translation};

    use super::*;

    #[test]
    fn an_intersection_encapsulates_t_and_object() {
        let s = Sphere::new();
        let i = Intersection::new(3.5, &s);

        assert_eq!(i.t, 3.5);
        assert_eq!(i.object, &s);
    }
    #[test]
    fn aggregating_intersections() {
        let s = Sphere::new();
        let i1 = Intersection::new(1., &s);
        let i2 = Intersection::new(2., &s);
        let xs = vec![i1, i2];

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn the_hit_when_all_intersections_have_positive_t() {
        let s = Sphere::new();
        let i1 = Intersection::new(1., &s);
        let i2 = Intersection::new(2., &s);
        let xs = vec![i1.clone(), i2];

        let i = hit(&xs).unwrap();
        assert_eq!(i, &i1);
    }

    #[test]
    fn the_hit_when_some_intersections_have_negative_t() {
        let s = Sphere::new();
        let i1 = Intersection::new(-1., &s);
        let i2 = Intersection::new(1., &s);
        let xs = vec![i1, i2.clone()];

        let i = hit(&xs).unwrap();
        assert_eq!(i, &i2);
    }

    #[test]
    fn the_hit_when_all_intersections_have_negative_t() {
        let s = Sphere::new();
        let i1 = Intersection::new(-2., &s);
        let i2 = Intersection::new(-1., &s);
        let xs = vec![i1, i2];

        assert!(hit(&xs).is_none());
    }

    #[test]
    fn the_hit_is_always_the_lowest_nonnegative_intersection() {
        let s = Sphere::new();
        let i1 = Intersection::new(5., &s);
        let i2 = Intersection::new(7., &s);
        let i3 = Intersection::new(-3., &s);
        let i4 = Intersection::new(2., &s);
        let xs = vec![i1, i2, i3, i4.clone()];

        let i = hit(&xs).unwrap();
        assert_eq!(i, &i4);
    }

    #[test]
    fn precomputing_the_state_of_an_intersection() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let i = Intersection::new(4.0, &shape);
        let comps = i.prepare_computations(&r);

        assert_eq!(comps.t, i.t);
        assert_eq!(comps.object, i.object);
        assert_eq!(comps.point, Tuple::point(0.0, 0.0, -1.0));
        assert_eq!(comps.eyev, Tuple::vector(0.0, 0.0, -1.0));
        assert_eq!(comps.normalv, Tuple::vector(0.0, 0.0, -1.0));
    }

    // Scenario: Precomputing the reflection vector
    //   Given shape â† plane()
    //     And r â† ray(point(0, 1, -1), vector(0, -âˆš2/2, âˆš2/2))
    //     And i â† intersection(âˆš2, shape)
    //   When comps â† prepare_computations(i, r)
    //   Then comps.reflectv = vector(0, âˆš2/2, âˆš2/2)

    #[test]
    fn the_hit_when_an_intersection_occurs_on_the_outside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let i = Intersection::new(4.0, &shape);
        let comps = i.prepare_computations(&r);
        assert!(!comps.inside);
    }

    #[test]
    fn the_hit_when_an_intersection_occurs_on_the_inside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let i = Intersection::new(1.0, &shape);
        let comps = i.prepare_computations(&r);

        assert_eq!(comps.t, i.t);
        assert_eq!(comps.object, i.object);
        assert_eq!(comps.point, Tuple::point(0.0, 0.0, 1.0));
        assert_eq!(comps.eyev, Tuple::vector(0.0, 0.0, -1.0));
        assert!(comps.inside);
        // Normal would have been (0, 0, 1), but is inverted!
        assert_eq!(comps.normalv, Tuple::vector(0.0, 0.0, -1.0));
    }

    #[test]
    fn the_hit_should_offset_the_point() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Sphere::new();
        shape.transform(translation(0.0, 0.0, 1.0));
        let i = Intersection::new(5.0, &shape);
        let comps = i.prepare_computations(&r);
        assert!(comps.over_point.z < -EPSILON / 2.0);
        assert!(comps.point.z > comps.over_point.z);
    }
    // Scenario: The under point is offset below the surface
    //   Given r â† ray(point(0, 0, -5), vector(0, 0, 1))
    //     And shape â† glass_sphere() with:
    //       | transform | translation(0, 0, 1) |
    //     And i â† intersection(5, shape)
    //     And xs â† intersections(i)
    //   When comps â† prepare_computations(i, r, xs)
    //   Then comps.under_point.z > EPSILON/2
    //     And comps.point.z < comps.under_point.z

    // Scenario Outline: Finding n1 and n2 at various intersections
    //   Given A â† glass_sphere() with:
    //       | transform                 | scaling(2, 2, 2) |
    //       | material.refractive_index | 1.5              |
    //     And B â† glass_sphere() with:
    //       | transform                 | translation(0, 0, -0.25) |
    //       | material.refractive_index | 2.0                      |
    //     And C â† glass_sphere() with:
    //       | transform                 | translation(0, 0, 0.25) |
    //       | material.refractive_index | 2.5                     |
    //     And r â† ray(point(0, 0, -4), vector(0, 0, 1))
    //     And xs â† intersections(2:A, 2.75:B, 3.25:C, 4.75:B, 5.25:C, 6:A)
    //   When comps â† prepare_computations(xs[<index>], r, xs)
    //   Then comps.n1 = <n1>
    //     And comps.n2 = <n2>

    //   Examples:
    //     | index | n1  | n2  |
    //     | 0     | 1.0 | 1.5 |
    //     | 1     | 1.5 | 2.0 |
    //     | 2     | 2.0 | 2.5 |
    //     | 3     | 2.5 | 2.5 |
    //     | 4     | 2.5 | 1.5 |
    //     | 5     | 1.5 | 1.0 |

    // Scenario: The Schlick approximation under total internal reflection
    //   Given shape â† glass_sphere()
    //     And r â† ray(point(0, 0, âˆš2/2), vector(0, 1, 0))
    //     And xs â† intersections(-âˆš2/2:shape, âˆš2/2:shape)
    //   When comps â† prepare_computations(xs[1], r, xs)
    //     And reflectance â† schlick(comps)
    //   Then reflectance = 1.0

    // Scenario: The Schlick approximation with a perpendicular viewing angle
    //   Given shape â† glass_sphere()
    //     And r â† ray(point(0, 0, 0), vector(0, 1, 0))
    //     And xs â† intersections(-1:shape, 1:shape)
    //   When comps â† prepare_computations(xs[1], r, xs)
    //     And reflectance â† schlick(comps)
    //   Then reflectance = 0.04

    // Scenario: The Schlick approximation with small angle and n2 > n1
    //   Given shape â† glass_sphere()
    //     And r â† ray(point(0, 0.99, -2), vector(0, 0, 1))
    //     And xs â† intersections(1.8589:shape)
    //   When comps â† prepare_computations(xs[0], r, xs)
    //     And reflectance â† schlick(comps)
    //   Then reflectance = 0.48873

    // Scenario: An intersection can encapsulate `u` and `v`
    //   Given s â† triangle(point(0, 1, 0), point(-1, 0, 0), point(1, 0, 0))
    //   When i â† intersection_with_uv(3.5, s, 0.2, 0.4)
    //   Then i.u = 0.2
    // And i.v = 0.4
}
