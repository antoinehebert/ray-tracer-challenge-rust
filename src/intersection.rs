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
        Self { t, object: &object }
    }

    pub fn prepare_computations(&self, ray: &Ray, xs: &Intersections) -> Computations {
        let point = ray.position(self.t);
        let eyev = -ray.direction;
        let mut normalv = Shape::normal_at(&self.object, &point);

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
                        .get_material()
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
                        .get_material()
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

    pub fn hit(xs: &'a Intersections) -> Option<&'a Intersection<'a>> {
        xs.iter()
            .filter(|x| x.t >= 0.)
            .min_by(|x, y| x.t.partial_cmp(&y.t).unwrap_or(std::cmp::Ordering::Equal))
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

impl<'a> Computations<'a> {
    // Approximation of Fresnel Effect
    //
    // Returns a number between 0 and 1, inclusive. This number is called the reflectance and represents what fraction
    // of the light is reflected, given the surface information at the hit.
    pub fn schlick(&self) -> f64 {
        // find the cosine of the angle between the eye and normal vectors
        let mut cos = self.eyev.dot(&self.normalv);

        // Total internal reflection can only occur if n1 > n2
        if self.n1 > self.n2 {
            // DRY: we have the same code in World#refracted_color
            let n = self.n1 / self.n2;
            let sin2_t = n.powi(2) * (1.0 - cos.powi(2));
            if sin2_t > 1.0 {
                return 1.0;
            }

            // compute cosine of theta_t using trig identity
            let cos_t = (1.0 - sin2_t).sqrt();
            // when n1 > n2, use cos(theta_t) instead
            cos = cos_t;
        }

        let r0 = ((self.n1 - self.n2) / (self.n1 + self.n2)).powi(2);
        r0 + (1.0 - r0) * (1.0 - cos).powi(5)
    }
}

#[cfg(test)]
mod tests {
    use crate::{assert_almost_eq, transformations::*};

    use super::*;

    #[test]
    fn an_intersection_encapsulates_t_and_object() {
        let s = Shape::sphere();
        let i = Intersection::new(3.5, &s);

        assert_eq!(i.t, 3.5);
        assert_eq!(*i.object, s);
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

        let i = Intersection::hit(&xs).unwrap();
        assert_eq!(i, &i1);
    }

    #[test]
    fn the_hit_when_some_intersections_have_negative_t() {
        let s = Shape::sphere();
        let i1 = Intersection::new(-1., &s);
        let i2 = Intersection::new(1., &s);
        let xs = vec![i1, i2.clone()];

        let i = Intersection::hit(&xs).unwrap();
        assert_eq!(i, &i2);
    }

    #[test]
    fn the_hit_when_all_intersections_have_negative_t() {
        let s = Shape::sphere();
        let i1 = Intersection::new(-2., &s);
        let i2 = Intersection::new(-1., &s);
        let xs = vec![i1, i2];

        assert!(Intersection::hit(&xs).is_none());
    }

    #[test]
    fn the_hit_is_always_the_lowest_nonnegative_intersection() {
        let s = Shape::sphere();
        let i1 = Intersection::new(5., &s);
        let i2 = Intersection::new(7., &s);
        let i3 = Intersection::new(-3., &s);
        let i4 = Intersection::new(2., &s);
        let xs = vec![i1, i2, i3, i4.clone()];

        let i = Intersection::hit(&xs).unwrap();
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
        shape.set_transform(translation(0.0, 0.0, 1.0));
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
        a.set_transform(scaling(2.0, 2.0, 2.0));
        a.get_material_mut().refractive_index = 1.5;

        let mut b = Shape::glass_sphere();
        b.set_transform(translation(0.0, 0.0, -0.25));
        b.get_material_mut().refractive_index = 2.0;

        let mut c = Shape::glass_sphere();
        c.set_transform(translation(0.0, 0.0, 0.25));
        c.get_material_mut().refractive_index = 2.5;

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
        shape.set_transform(translation(0., 0., 1.));
        let i = Intersection::new(5., &shape);
        let xs: Intersections = vec![i];
        let comps = xs[0].prepare_computations(&r, &xs);
        assert!(comps.under_point.z > (EPSILON / 2.));
        assert!(comps.point.z < comps.under_point.z);
    }

    #[test]
    fn the_schlick_approximation_under_total_internal_reflection() {
        let shape = Shape::glass_sphere();
        let r = Ray::new(
            Tuple::point(0.0, 0.0, 2.0_f64.sqrt() / 2.0),
            Tuple::vector(0.0, 1.0, 0.0),
        );
        let xs = vec![
            Intersection::new(-2.0_f64.sqrt() / 2.0, &shape),
            Intersection::new(2.0_f64.sqrt() / 2.0, &shape),
        ];
        let comps = xs[1].prepare_computations(&r, &xs);
        let reflectance = comps.schlick();

        assert_eq!(reflectance, 1.0);
    }

    #[test]
    fn the_schlick_approximation_with_a_perpendicular_viewing_angle() {
        let shape = Shape::glass_sphere();
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let xs = vec![
            Intersection::new(-1.0, &shape),
            Intersection::new(1.0, &shape),
        ];
        let comps = xs[1].prepare_computations(&r, &xs);
        let reflectance = comps.schlick();

        assert_almost_eq!(reflectance, 0.04);
    }

    #[test]
    fn the_schlick_approximation_with_small_angle_and_n2_gt_n1() {
        let shape = Shape::glass_sphere();
        let r = Ray::new(Tuple::point(0.0, 0.99, -2.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = vec![Intersection::new(1.8589, &shape)];
        let comps = xs[0].prepare_computations(&r, &xs);
        let reflectance = comps.schlick();

        assert_almost_eq!(reflectance, 0.48873);
    }

    // Scenario: An intersection can encapsulate `u` and `v`
    //   Given s â† triangle(point(0, 1, 0), point(-1, 0, 0), point(1, 0, 0))
    //   When i â† intersection_with_uv(3.5, s, 0.2, 0.4)
    //   Then i.u = 0.2
    // And i.v = 0.4
}
