use crate::intersection::*;
use crate::matrix::Matrix;
use crate::ray::*;
use crate::transformations::*;
use crate::tuple::Tuple;
use crate::utils::*;

#[derive(Debug)]
pub struct Sphere {
    transform: Matrix<4>,
}

impl Sphere {
    pub fn new() -> Self {
        Sphere {
            transform: Matrix::<4>::identity(),
        }
    }

    // Returns intersection points (time) along `ray`.
    pub fn intersect(&self, ray: &Ray) -> Intersections {
        // We want the the sphere to always be centered at the world origin, so we move the ray by the inverse of the
        // sphere transformation.
        let ray = ray.transform(
            self.transform
                .inverse()
                .expect("sphere transform should be invertible"),
        );

        // The sphere is always centered at the world origin...
        let sphere_to_ray = ray.origin - Tuple::point(0., 0., 0.);
        let a = ray.direction.dot(&ray.direction);
        let b = 2. * ray.direction.dot(&sphere_to_ray);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.;

        // Classic quadratic formula!
        let discriminant = b.powf(2.) - 4. * a * c;

        let mut result = Vec::new();

        if discriminant >= 0. {
            let sqrt = discriminant.sqrt();
            result.push(Intersection::new((-b - sqrt) / (2. * a), &self));
            result.push(Intersection::new((-b + sqrt) / (2. * a), &self));
        }

        result
    }

    fn transform(&mut self, matrix: Matrix<4>) {
        self.transform = self.transform * matrix;
    }
}

impl PartialEq for Sphere {
    fn eq(&self, _other: &Self) -> bool {
        // all spheres are equal for now... since they are all centered at the origin and they all have the same radius
        // of 1.
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_almost_eq;

    #[test]
    fn a_ray_intersects_a_sphere_at_two_points() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, 4.0);
        assert_almost_eq!(xs[1].t, 6.0);
    }

    #[test]
    fn a_ray_intersects_a_sphere_at_a_tangen() {
        let r = Ray::new(Tuple::point(0., 1., -5.), Tuple::vector(0., 0., 1.));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, 5.0);
        assert_almost_eq!(xs[1].t, 5.0);
    }

    #[test]
    fn a_ray_misses_a_sphere() {
        let r = Ray::new(Tuple::point(0., 2., -5.), Tuple::vector(0., 0., 1.));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn a_ray_originates_inside_a_spher() {
        let r = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, -1.0);
        assert_almost_eq!(xs[1].t, 1.0);
    }

    #[test]
    fn a_sphere_is_behind_a_ray() {
        let r = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, -6.0);
        assert_almost_eq!(xs[1].t, -4.0);
    }

    #[test]
    fn intersect_sets_the_object_on_the_intersection() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].object, &s);
        assert_eq!(xs[1].object, &s);
    }

    #[test]
    fn a_sphere_s_default_transformations() {
        let s = Sphere::new();
        assert_eq!(s.transform, Matrix::<4>::identity())
    }

    #[test]
    fn changing_a_sphere_s_transformations() {
        let mut s = Sphere::new();
        let t = translation(2., 3., 4.);
        s.transform(t);

        assert_eq!(s.transform, t)
    }

    #[test]
    fn intersecting_a_scaled_sphere_with_a_ray() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut s = Sphere::new();

        s.transform(scaling(2., 2., 2.));
        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, 3.);
        assert_almost_eq!(xs[1].t, 7.);
    }

    #[test]
    fn intersecting_a_translated_sphere_with_a_ray() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut s = Sphere::new();

        s.transform(translation(5., 0., 0.));
        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 0);
    }
    // Scenario: The normal on a sphere at a point on the x axis
    //   Given s â† sphere()
    //   When n â† normal_at(s, point(1., 0., 0.))
    //   Then n = vector(1., 0., 0.)

    // Scenario: The normal on a sphere at a point on the y axis
    //   Given s â† sphere()
    //   When n â† normal_at(s, point(0., 1., 0.))
    //   Then n = vector(0., 1., 0.)

    // Scenario: The normal on a sphere at a point on the z axis
    //   Given s â† sphere()
    //   When n â† normal_at(s, point(0., 0., 1.))
    //   Then n = vector(0., 0., 1.)

    // Scenario: The normal on a sphere at a nonaxial point
    //   Given s â† sphere()
    //   When n â† normal_at(s, point(âˆš3./3., âˆš3./3., âˆš3./3.))
    //   Then n = vector(âˆš3./3., âˆš3./3., âˆš3./3.)

    // Scenario: The normal is a normalized vector
    //   Given s â† sphere()
    //   When n â† normal_at(s, point(âˆš3./3., âˆš3./3., âˆš3./3.))
    //   Then n = normalize(n)

    // Scenario: Computing the normal on a translated sphere
    //   Given s â† sphere()
    //     And set_transform(s, translation(0., 1., 0.))
    //   When n â† normal_at(s, point(0., 1.70711, -0.70711))
    //   Then n = vector(0., 0.70711, -0.70711)

    // Scenario: Computing the normal on a transformed sphere
    //   Given s â† sphere()
    //     And m â† scaling(1., 0.5, 1.) * rotation_z(Ï€/5.)
    //     And set_transform(s, m)
    //   When n â† normal_at(s, point(0., âˆš2./2., -âˆš2./2.))
    //   Then n = vector(0., 0.97014, -0.24254)

    // Scenario: A sphere has a default material
    //   Given s â† sphere()
    //   When m â† s.material
    //   Then m = material()

    // Scenario: A sphere may be assigned a material
    //   Given s â† sphere()
    //     And m â† material()
    //     And m.ambient â† 1
    //   When s.material â† m
    //   Then s.material = m

    // Scenario: A helper for producing a sphere with a glassy material
    //   Given s â† glass_sphere()
    //   Then s.transform = identity_matrix
    //     And s.material.transparency = 1.0
    //     And s.material.refractive_index = 1.5
}
