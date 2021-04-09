use crate::{
    intersection::{Intersection, Intersections},
    material::Material,
    matrix::Matrix,
    ray::Ray,
    tuple::Tuple,
    utils::EPSILON,
};

#[derive(Debug, PartialEq)]
enum ShapeKind {
    Sphere, // The sphere is always centered at the world origin...
    Plane,  // Plane is in xy, with the normal pointing in the positive y direction.
}

#[derive(Debug, PartialEq)]
pub struct Shape {
    pub transform: Matrix<4>,
    pub material: Material,
    kind: ShapeKind,
}

impl Shape {
    pub fn sphere() -> Self {
        Shape {
            transform: Matrix::<4>::identity(),
            material: Material::new(),
            kind: ShapeKind::Sphere,
        }
    }

    pub fn plane() -> Self {
        Shape {
            transform: Matrix::<4>::identity(),
            material: Material::new(),
            kind: ShapeKind::Plane,
        }
    }

    // Returns intersection points (time) along `ray`.
    pub fn intersect(&self, world_ray: &Ray) -> Intersections {
        let local_ray = world_ray.transform(
            self.transform
                .inverse()
                .expect("shape transfor should be invertible"),
        );

        let mut result = Vec::new();
        match self.kind {
            ShapeKind::Sphere { .. } => {
                // The sphere is always centered at the world origin...
                let sphere_to_ray = local_ray.origin - Tuple::point(0., 0., 0.);
                let a = local_ray.direction.dot(&local_ray.direction);
                let b = 2. * local_ray.direction.dot(&sphere_to_ray);
                let c = sphere_to_ray.dot(&sphere_to_ray) - 1.;

                // Classic quadratic formula!
                let discriminant = b.powf(2.) - 4. * a * c;

                if discriminant >= 0. {
                    let sqrt = discriminant.sqrt();
                    result.push(Intersection::new((-b - sqrt) / (2. * a), &self));
                    result.push(Intersection::new((-b + sqrt) / (2. * a), &self));
                }
            }
            ShapeKind::Plane { .. } => {
                // Plane is in xy, with the normal pointing in the positive y direction.
                if local_ray.direction.y.abs() >= EPSILON {
                    result.push(Intersection::new(
                        -local_ray.origin.y / local_ray.direction.y,
                        &self,
                    ));
                }
            }
        };

        result
    }

    pub fn normal_at(&self, world_point: &Tuple) -> Tuple {
        let sphere_inverted_transform = self
            .transform
            .inverse()
            .expect("Transform should be invertible");
        let local_point = sphere_inverted_transform * *world_point;

        let local_normal = match self.kind {
            ShapeKind::Sphere => local_point - Tuple::point(0., 0., 0.),
            ShapeKind::Plane => Tuple::point(0.0, 1.0, 0.0),
        };

        let mut world_normal = sphere_inverted_transform.transpose() * local_normal;
        // Hack: Instead of removing any translation by taking a 3x3 submatrix of the transform, we just set w to 0.
        world_normal.w = 0.;

        world_normal.normalize()
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
    fn the_default_transformation() {
        let s = Shape::sphere();
        assert_eq!(s.transform, Matrix::identity());
    }

    #[test]
    fn assigning_a_transformation() {
        let mut s = Shape::sphere();
        s.transform = translation(2.0, 3.0, 4.0);

        assert_eq!(s.transform, translation(2.0, 3.0, 4.0));
    }

    #[test]
    fn the_default_material() {
        let s = Shape::sphere();
        let m = s.material;

        assert_eq!(m, Material::new());
    }

    #[test]
    fn assigning_a_material() {
        let mut s = Shape::sphere();
        let mut m = Material::new();
        m.ambient = 1.0;
        s.material = m;
        assert_eq!(s.material, m);
    }

    #[test]
    fn a_ray_intersects_a_sphere_at_two_points() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let s = Shape::sphere();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, 4.0);
        assert_almost_eq!(xs[1].t, 6.0);
    }

    #[test]
    fn a_ray_intersects_a_sphere_at_a_tangent() {
        let r = Ray::new(Tuple::point(0., 1., -5.), Tuple::vector(0., 0., 1.));
        let s = Shape::sphere();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, 5.0);
        assert_almost_eq!(xs[1].t, 5.0);
    }

    #[test]
    fn a_ray_misses_a_sphere() {
        let r = Ray::new(Tuple::point(0., 2., -5.), Tuple::vector(0., 0., 1.));
        let s = Shape::sphere();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn a_ray_originates_inside_a_sphere() {
        let r = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
        let s = Shape::sphere();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, -1.0);
        assert_almost_eq!(xs[1].t, 1.0);
    }

    #[test]
    fn a_sphere_is_behind_a_ray() {
        let r = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
        let s = Shape::sphere();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, -6.0);
        assert_almost_eq!(xs[1].t, -4.0);
    }

    #[test]
    fn intersect_sets_the_object_on_the_intersection() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let s = Shape::sphere();
        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].object, &s);
        assert_eq!(xs[1].object, &s);
    }

    #[test]
    fn a_sphere_s_default_transformations() {
        let s = Shape::sphere();
        assert_eq!(s.transform, Matrix::<4>::identity())
    }

    #[test]
    fn changing_a_sphere_s_transformations() {
        let mut s = Shape::sphere();
        let t = translation(2., 3., 4.);
        s.transform = t;

        assert_eq!(s.transform, t)
    }

    #[test]
    fn intersecting_a_scaled_sphere_with_a_ray() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut s = Shape::sphere();

        s.transform = scaling(2., 2., 2.);
        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);

        assert_almost_eq!(xs[0].t, 3.);
        assert_almost_eq!(xs[1].t, 7.);
    }

    #[test]
    fn intersecting_a_translated_sphere_with_a_ray() {
        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut s = Shape::sphere();

        s.transform = translation(5., 0., 0.);
        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 0);
    }
    #[test]
    fn the_normal_on_a_sphere_at_a_point_on_the_x_axis() {
        let s = Shape::sphere();
        let n = s.normal_at(&Tuple::point(1., 0., 0.));
        assert_eq!(n, Tuple::vector(1., 0., 0.));
    }

    #[test]
    fn the_normal_on_a_sphere_at_a_point_on_the_y_axis() {
        let s = Shape::sphere();
        let n = s.normal_at(&Tuple::point(0., 1., 0.));
        assert_eq!(n, Tuple::vector(0., 1., 0.));
    }

    #[test]
    fn the_normal_on_a_sphere_at_a_point_on_the_z_axis() {
        let s = Shape::sphere();
        let n = s.normal_at(&Tuple::point(0., 0., 1.));
        assert_eq!(n, Tuple::vector(0., 0., 1.));
    }

    #[test]
    fn the_normal_on_a_sphere_at_a_nonaxial_point() {
        let s = Shape::sphere();
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
        let s = Shape::sphere();
        let n = s.normal_at(&Tuple::point(
            (3. as f64).sqrt() / 3.,
            (3. as f64).sqrt() / 3.,
            (3. as f64).sqrt() / 3.,
        ));
        assert_eq!(n, n.normalize());
    }

    #[test]
    fn computing_the_normal_on_a_translated_sphere() {
        let mut s = Shape::sphere();
        s.transform = translation(0., 1., 0.);

        let n = s.normal_at(&Tuple::point(0., 1.70711, -0.70711));
        assert_eq!(n, Tuple::vector(0., 0.70711, -0.70711));
    }

    #[test]
    fn computing_the_normal_on_a_transformed_sphere() {
        let mut s = Shape::sphere();
        let m = scaling(1., 0.5, 1.) * rotation_z(PI / 5.);
        s.transform = m;
        let n = s.normal_at(&Tuple::point(
            0.,
            (2 as f64).sqrt() / 2.,
            -(2 as f64).sqrt() / 2.,
        ));
        assert_eq!(n, Tuple::vector(0., 0.97014, -0.24254));
    }

    #[test]
    fn a_sphere_has_a_default_material() {
        let s = Shape::sphere();
        let m = s.material;

        assert_eq!(m, Material::new());
    }

    #[test]
    fn a_sphere_may_be_assigned_a_material() {
        let mut s = Shape::sphere();
        let mut m = Material::new();
        m.ambient = 1.234;
        s.material = m;
        assert_eq!(s.material, m);
    }

    // --------------
    // Feature: Shape
    // --------------

    // Scenario: A shape has a parent attribute
    //   Given s â† test_shape()
    //   Then s.parent is nothing

    // Scenario: Converting a point from world to object space
    //   Given g1 â† group()
    //     And set_transform(g1, rotation_y(Ï€/2))
    //     And g2 â† group()
    //     And set_transform(g2, scaling(2, 2, 2))
    //     And add_child(g1, g2)
    //     And s â† sphere()
    //     And set_transform(s, translation(5, 0, 0))
    //     And add_child(g2, s)
    //   When p â† world_to_object(s, point(-2, 0, -10))
    //   Then p = point(0, 0, -1)

    // Scenario: Converting a normal from object to world space
    //   Given g1 â† group()
    //     And set_transform(g1, rotation_y(Ï€/2))
    //     And g2 â† group()
    //     And set_transform(g2, scaling(1, 2, 3))
    //     And add_child(g1, g2)
    //     And s â† sphere()
    //     And set_transform(s, translation(5, 0, 0))
    //     And add_child(g2, s)
    //   When n â† normal_to_world(s, vector(âˆš3/3, âˆš3/3, âˆš3/3))
    //   Then n = vector(0.2857, 0.4286, -0.8571)

    // Scenario: Finding the normal on a child object
    //   Given g1 â† group()
    //     And set_transform(g1, rotation_y(Ï€/2))
    //     And g2 â† group()
    //     And set_transform(g2, scaling(1, 2, 3))
    //     And add_child(g1, g2)
    //     And s â† sphere()
    //     And set_transform(s, translation(5, 0, 0))
    //     And add_child(g2, s)
    //   When n â† normal_at(s, point(1.7321, 1.1547, -5.5774))
    //   Then n = vector(0.2857, 0.4286, -0.8571)

    // ---------------
    // Feature: Sphere
    // ---------------

    // Scenario: A helper for producing a sphere with a glassy material
    //   Given s â† glass_sphere()
    //   Then s.transform = identity_matrix
    //     And s.material.transparency = 1.0
    //     And s.material.refractive_index = 1.5

    // ---------------
    // Feature: Planes
    // ---------------

    #[test]
    fn the_normal_of_a_plane_is_constant_everywhere() {
        let p = Shape::plane();

        let n1 = p.normal_at(&Tuple::point(0.0, 0.0, 0.0));
        let n2 = p.normal_at(&Tuple::point(10.0, 0.0, -10.0));
        let n3 = p.normal_at(&Tuple::point(-5.0, 0.0, 150.0));

        assert_eq!(n1, Tuple::vector(0.0, 1.0, 0.0));
        assert_eq!(n2, Tuple::vector(0.0, 1.0, 0.0));
        assert_eq!(n3, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn intersect_with_a_ray_parallel_to_the_plane() {
        let p = Shape::plane();
        let r = Ray::new(Tuple::point(0.0, 10.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = p.intersect(&r);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn intersect_with_a_coplanar_ray() {
        let p = Shape::plane();
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = p.intersect(&r);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn a_ray_intersecting_a_plane_from_above() {
        let p = Shape::plane();
        let r = Ray::new(Tuple::point(0.0, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        let xs = p.intersect(&r);
        assert_eq!(xs.len(), 1);
        assert_eq!(xs[0].t, 1.0);
        assert_eq!(xs[0].object, &p);
    }

    #[test]
    fn a_ray_intersecting_a_plane_from_below() {
        let p = Shape::plane();
        let r = Ray::new(Tuple::point(0.0, -1.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let xs = p.intersect(&r);
        assert_eq!(xs.len(), 1);
        assert_eq!(xs[0].t, 1.0);
        assert_eq!(xs[0].object, &p);
    }
}
