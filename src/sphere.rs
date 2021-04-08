use crate::intersection::*;
use crate::material::*;
use crate::matrix::Matrix;
use crate::ray::*;
use crate::shape::Shape;
use crate::tuple::Tuple;

#[derive(Debug, PartialEq)]
pub struct Sphere {
    transform: Matrix<4>,
    material: Material,
}

impl Sphere {
    pub fn new() -> Self {
        Sphere {
            transform: Matrix::<4>::identity(),
            material: Material::new(),
        }
    }
}

impl Shape for Sphere {
    fn transform(&self) -> &Matrix<4> {
        &self.transform
    }

    fn set_transform(&mut self, transform: Matrix<4>) {
        self.transform = transform;
    }

    fn material(&self) -> &Material {
        &self.material
    }

    fn set_material(&mut self, material: Material) {
        self.material = material
    }

    fn local_intersect(&self, local_ray: &Ray) -> Intersections {
        // The sphere is always centered at the world origin...
        let sphere_to_ray = local_ray.origin - Tuple::point(0., 0., 0.);
        let a = local_ray.direction.dot(&local_ray.direction);
        let b = 2. * local_ray.direction.dot(&sphere_to_ray);
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

    fn local_normal_at(&self, local_point: &Tuple) -> Tuple {
        *local_point - Tuple::point(0., 0., 0.)
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::assert_almost_eq;
    use crate::transformations::*;
    use crate::utils::*;

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
        s.set_transform(t);

        assert_eq!(s.transform, t)
    }

    #[test]
    fn intersecting_a_scaled_sphere_with_a_ray() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut s = Sphere::new();

        s.set_transform(scaling(2., 2., 2.));
        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, 3.);
        assert_almost_eq!(xs[1].t, 7.);
    }

    #[test]
    fn intersecting_a_translated_sphere_with_a_ray() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut s = Sphere::new();

        s.set_transform(translation(5., 0., 0.));
        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 0);
    }
    #[test]
    fn the_normal_on_a_sphere_at_a_point_on_the_x_axis() {
        let s = Sphere::new();
        let n = s.normal_at(&Tuple::point(1., 0., 0.));
        assert_eq!(n, Tuple::vector(1., 0., 0.));
    }

    #[test]
    fn the_normal_on_a_sphere_at_a_point_on_the_y_axis() {
        let s = Sphere::new();
        let n = s.normal_at(&Tuple::point(0., 1., 0.));
        assert_eq!(n, Tuple::vector(0., 1., 0.));
    }

    #[test]
    fn the_normal_on_a_sphere_at_a_point_on_the_z_axis() {
        let s = Sphere::new();
        let n = s.normal_at(&Tuple::point(0., 0., 1.));
        assert_eq!(n, Tuple::vector(0., 0., 1.));
    }

    #[test]
    fn the_normal_on_a_sphere_at_a_nonaxial_point() {
        let s = Sphere::new();
        let n = s.normal_at(&Tuple::point(
            (3. as f64).sqrt() / 3.,
            (3. as f64).sqrt() / 3.,
            (3. as f64).sqrt() / 3.,
        ));
        assert_eq!(
            n,
            Tuple::vector(
                (3. as f64).sqrt() / 3.,
                (3. as f64).sqrt() / 3.,
                (3. as f64).sqrt() / 3.
            )
        );
    }

    #[test]
    fn the_normal_is_a_normalized_vector() {
        let s = Sphere::new();
        let n = s.normal_at(&Tuple::point(
            (3. as f64).sqrt() / 3.,
            (3. as f64).sqrt() / 3.,
            (3. as f64).sqrt() / 3.,
        ));
        assert_eq!(n, n.normalize());
    }

    #[test]
    fn computing_the_normal_on_a_translated_sphere() {
        let mut s = Sphere::new();
        s.set_transform(translation(0., 1., 0.));

        let n = s.normal_at(&Tuple::point(0., 1.70711, -0.70711));
        assert_eq!(n, Tuple::vector(0., 0.70711, -0.70711));
    }

    #[test]
    fn computing_the_normal_on_a_transformed_sphere() {
        let mut s = Sphere::new();
        let m = scaling(1., 0.5, 1.) * rotation_z(PI / 5.);
        s.set_transform(m);
        let n = s.normal_at(&Tuple::point(
            0.,
            (2 as f64).sqrt() / 2.,
            -(2 as f64).sqrt() / 2.,
        ));
        assert_eq!(n, Tuple::vector(0., 0.97014, -0.24254));
    }

    #[test]
    fn a_sphere_has_a_default_material() {
        let s = Sphere::new();
        let m = s.material;

        assert_eq!(m, Material::new());
    }

    #[test]
    fn a_sphere_may_be_assigned_a_material() {
        let mut s = Sphere::new();
        let mut m = Material::new();
        m.ambient = 1.234;
        s.material = m;
        assert_eq!(s.material, m);
    }

    // Scenario: A helper for producing a sphere with a glassy material
    //   Given s â† glass_sphere()
    //   Then s.transform = identity_matrix
    //     And s.material.transparency = 1.0
    //     And s.material.refractive_index = 1.5
}
