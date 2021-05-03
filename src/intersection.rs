use crate::ray::Ray;
use crate::shape::Shape;
use crate::tuple::Tuple;
use crate::utils::*;

#[derive(Debug, Clone, PartialEq)]
pub struct Intersection<'a> {
    pub t: f64, // intersection "time"
    pub object: &'a Shape,
}

impl<'a> Intersection<'a> {
    pub fn new(t: f64, object: &'a Shape) -> Self {
        Self { t, object }
    }

    pub fn prepare_computations(&self, ray: &Ray, xs: &Intersections) -> Computations {
        let point = ray.position(self.t);
        let eyev = -ray.direction;
        let mut normalv = self.object.normal_at(&point);

        let inside = normalv.dot(&eyev) < 0.0;
        if inside {
            normalv = -normalv;
        }

        let reflectv = ray.direction.reflect(&normalv);

        let mut containers: Vec<&Shape> = vec![];
        let mut n1 = 1.0;
        let mut n2 = 1.0;
        'refraction_index_loop: for i in xs {
            if i == self {
                if !containers.is_empty() {
                    n1 = containers
                        .last()
                        .expect("there should be a last element in the list")
                        .material
                        .refractive_index
                };
            }

            match containers.iter().position(|&o| o == i.object) {
                Some(index) => {
                    containers.remove(index);
                }
                None => {
                    containers.push(i.object);
                }
            }

            if i == self {
                if !containers.is_empty() {
                    n2 = containers
                        .last()
                        .expect("there should be a last element in the list")
                        .material
                        .refractive_index
                };
                break 'refraction_index_loop;
            }
        }

        Computations {
            t: self.t,
            object: self.object,
            point,
            over_point: point + normalv * EPSILON,
            under_point: point - normalv * EPSILON,
            eyev,
            inside,
            normalv,
            reflectv,
            n1,
            n2,
        }
    }
}

pub type Intersections<'a> = Vec<Intersection<'a>>;

pub struct Computations<'a> {
    pub t: f64,
    pub object: &'a Shape,
    pub point: Tuple,
    pub over_point: Tuple,  // Prevent objects from shadowing themselves
    pub under_point: Tuple, // Where refracted rays will originate
    pub eyev: Tuple,
    pub inside: bool,
    pub normalv: Tuple,
    pub reflectv: Tuple,
    pub n1: f64, // refractive_index from
    pub n2: f64, // refractive_index to
}

// TODO: make this a method on Intersections.
pub fn hit<'a>(xs: &'a Intersections) -> Option<&'a Intersection<'a>> {
    xs.iter()
        .filter(|x| x.t >= 0.)
        .min_by(|x, y| x.t.partial_cmp(&y.t).unwrap_or(std::cmp::Ordering::Equal))
}

#[cfg(test)]
mod tests {
    use crate::transformations::*;

    use super::*;

    #[test]
    fn an_intersection_encapsulates_t_and_object() {
        let s = Shape::sphere();
        let i = Intersection::new(3.5, &s);

        assert_eq!(i.t, 3.5);
        assert_eq!(i.object, &s);
    }
    #[test]
    fn aggregating_intersections() {
        let s = Shape::sphere();
        let i1 = Intersection::new(1., &s);
        let i2 = Intersection::new(2., &s);
        let xs = vec![i1, i2];

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn the_hit_when_all_intersections_have_positive_t() {
        let s = Shape::sphere();
        let i1 = Intersection::new(1., &s);
        let i2 = Intersection::new(2., &s);
        let xs = vec![i1.clone(), i2];

        let i = hit(&xs).unwrap();
        assert_eq!(i, &i1);
    }

    #[test]
    fn the_hit_when_some_intersections_have_negative_t() {
        let s = Shape::sphere();
        let i1 = Intersection::new(-1., &s);
        let i2 = Intersection::new(1., &s);
        let xs = vec![i1, i2.clone()];

        let i = hit(&xs).unwrap();
        assert_eq!(i, &i2);
    }

    #[test]
    fn the_hit_when_all_intersections_have_negative_t() {
        let s = Shape::sphere();
        let i1 = Intersection::new(-2., &s);
        let i2 = Intersection::new(-1., &s);
        let xs = vec![i1, i2];

        assert!(hit(&xs).is_none());
    }

    #[test]
    fn the_hit_is_always_the_lowest_nonnegative_intersection() {
        let s = Shape::sphere();
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
        let shape = Shape::sphere();
        let i = Intersection::new(4.0, &shape);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);

        assert_eq!(comps.t, i.t);
        assert_eq!(comps.object, i.object);
        assert_eq!(comps.point, Tuple::point(0.0, 0.0, -1.0));
        assert_eq!(comps.eyev, Tuple::vector(0.0, 0.0, -1.0));
        assert_eq!(comps.normalv, Tuple::vector(0.0, 0.0, -1.0));
    }

    #[test]
    fn precomputing_the_reflection_vector() {
        let shape = Shape::plane();
        let r = Ray::new(
            Tuple::point(0.0, 1.0, -1.0),
            Tuple::vector(0.0, -(2.0 as f64).sqrt() / 2.0, (2.0 as f64).sqrt() / 2.0),
        );
        let i = Intersection::new((2 as f64).sqrt(), &shape);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        assert_eq!(
            comps.reflectv,
            Tuple::vector(0.0, (2.0 as f64).sqrt() / 2.0, (2.0 as f64).sqrt() / 2.0)
        );
    }
    #[test]
    fn the_hit_when_an_intersection_occurs_on_the_outside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::sphere();
        let i = Intersection::new(4.0, &shape);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        assert!(!comps.inside);
    }

    #[test]
    fn the_hit_when_an_intersection_occurs_on_the_inside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::sphere();
        let i = Intersection::new(1.0, &shape);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);

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
        let mut shape = Shape::sphere();
        shape.transform = translation(0.0, 0.0, 1.0);
        let i = Intersection::new(5.0, &shape);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        assert!(comps.over_point.z < -EPSILON / 2.0);
        assert!(comps.point.z > comps.over_point.z);
    }

    //                    -------------
    //                ---/             \---
    //             --/          A          \--
    //           -/                           \-
    //          /                               \
    //         /                                 \
    //        /          ------- -------          \
    //       /         -/      -X-      \-         \
    //      /         /       /   \       \         \
    //     0|       1/      2/    3\      4\       5 |
    //------+--------+-------+-----+-------+--------+------------>
    //      |        \       \     /       /        |
    //      \         \       \   /       /         /
    //       \         -\   B  -X-   C   /-         /
    //        \          ------- -------          /
    //         \                                 /
    //          \                               /
    //           -\                           /-
    //             --\                     /--
    //                ---\             /---
    //                    -------------
    #[test]
    fn finding_n1_and_n2_at_various_intersections() {
        let mut a = Shape::glass_sphere();
        a.transform = scaling(2.0, 2.0, 2.0);
        a.material.refractive_index = 1.5;

        let mut b = Shape::glass_sphere();
        b.transform = translation(0.0, 0.0, -0.25);
        b.material.refractive_index = 2.0;

        let mut c = Shape::glass_sphere();
        c.transform = translation(0.0, 0.0, 0.25);
        c.material.refractive_index = 2.5;

        let r = Ray::new(Tuple::point(0.0, 0.0, -4.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs: Intersections = vec![
            Intersection::new(2.0, &a),
            Intersection::new(2.75, &b),
            Intersection::new(3.25, &c),
            Intersection::new(4.75, &b),
            Intersection::new(5.25, &c),
            Intersection::new(6.0, &a),
        ];

        let examples = [
            [1.0, 1.5],
            [1.5, 2.0],
            [2.0, 2.5],
            [2.5, 2.5],
            [2.5, 1.5],
            [1.5, 1.0],
        ];

        for (index, ns) in examples.iter().enumerate() {
            let comps = (xs[index]).prepare_computations(&r, &xs);
            assert_eq!(comps.n1, ns[0]);
            assert_eq!(comps.n2, ns[1]);
        }
    }

    #[test]
    fn the_under_point_is_offset_below_the_surface() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut shape = Shape::glass_sphere();
        shape.transform = translation(0., 0., 1.);
        let i = Intersection::new(5., &shape);
        let xs: Intersections = vec![i];
        let comps = xs[0].prepare_computations(&r, &xs);
        assert!(comps.under_point.z > (EPSILON / 2.));
        assert!(comps.point.z < comps.under_point.z);
    }

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
