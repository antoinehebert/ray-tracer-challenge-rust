use std::mem::swap;

use crate::{
    intersection::{Intersection, Intersections},
    material::Material,
    matrix::Matrix,
    ray::Ray,
    tuple::Tuple,
    utils::{is_almost_equal, EPSILON},
};

#[derive(Debug, PartialEq, Clone)]
enum ShapeKind {
    Sphere, // The sphere is always centered at the world origin...
    Plane,  // Plane is in xy, with the normal pointing in the positive y direction.
    Cube,   // Centered at the world origin and going from -1 to 1.
    Cylinder {
        minimum: f64,
        maximum: f64,
        capped: bool,
    }, // Radius of 1, extending to infinity in both +y and -y unless it's capped.
    Cone {
        minimum: f64,
        maximum: f64,
        capped: bool,
    }, // Radius of 1, extending to infinity in both +y and -y unless it's capped.
}

#[derive(Debug, PartialEq, Clone)]
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

    pub fn glass_sphere() -> Self {
        let mut material = Material::new();
        material.transparency = 1.0;
        material.refractive_index = 1.5;

        Self {
            transform: Matrix::<4>::identity(),
            material: material,
            kind: ShapeKind::Sphere,
        }
    }

    pub fn plane() -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            material: Material::new(),
            kind: ShapeKind::Plane,
        }
    }

    pub fn cube() -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            material: Material::new(),
            kind: ShapeKind::Cube,
        }
    }

    pub fn infinite_cylinder() -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            material: Material::new(),
            kind: ShapeKind::Cylinder {
                minimum: -f64::INFINITY,
                maximum: f64::INFINITY,
                capped: false,
            },
        }
    }

    pub fn cylinder(min: f64, max: f64, capped: bool) -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            material: Material::new(),
            kind: ShapeKind::Cylinder {
                minimum: min,
                maximum: max,
                capped,
            },
        }
    }

    pub fn infinite_cone() -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            material: Material::new(),
            kind: ShapeKind::Cone {
                minimum: -f64::INFINITY,
                maximum: f64::INFINITY,
                capped: false,
            },
        }
    }

    pub fn cone(min: f64, max: f64, capped: bool) -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            material: Material::new(),
            kind: ShapeKind::Cone {
                minimum: min,
                maximum: max,
                capped,
            },
        }
    }

    // Returns intersection points (time) along `ray`.
    pub fn intersect(&self, world_ray: &Ray) -> Intersections {
        let local_ray = world_ray.transform(
            self.transform
                .inverse()
                .expect("shape transform should be invertible"),
        );

        let mut result = Vec::new();
        match self.kind {
            ShapeKind::Sphere => {
                // The sphere is always centered at the world origin...
                let sphere_to_ray = local_ray.origin - Tuple::point(0., 0., 0.);
                let a = local_ray.direction.dot(&local_ray.direction);
                let b = 2. * local_ray.direction.dot(&sphere_to_ray);
                let c = sphere_to_ray.dot(&sphere_to_ray) - 1.;

                // Classic quadratic formula!
                let discriminant = b.powi(2) - 4. * a * c;

                if discriminant >= 0. {
                    let sqrt = discriminant.sqrt();
                    result.push(Intersection::new((-b - sqrt) / (2. * a), &self));
                    result.push(Intersection::new((-b + sqrt) / (2. * a), &self));
                }
            }
            ShapeKind::Plane => {
                // Plane is in xy, with the normal pointing in the positive y direction.
                if local_ray.direction.y.abs() >= EPSILON {
                    result.push(Intersection::new(
                        -local_ray.origin.y / local_ray.direction.y,
                        &self,
                    ));
                }
            }
            ShapeKind::Cube => {
                // Treat the cube as if it were composed of 6 panes. The intersection of the ray with that square will
                // always be those two points: the largest minimum t value and the smallest maximum t value.
                //
                // 2D example: Here B is the larges mimimum t and C the smallest maximum t.
                //
                //                |                |  .>
                //                |                 .
                //                |               .D
                //        --------+------------.-C-+--------
                //                |          .     |
                //                |       .        |
                //                |    .           |
                //                | .              |
                //                .                |
                //              . B                |
                //         --A.---+----------------+--------
                //          .     |                |
                //        .       |                |
                //                |                |
                //
                // TODO: Optimize: You don't need to check all 6 planes if it's clear that the ray misses.
                let (xtmin, xtmax) = check_axis(local_ray.origin.x, local_ray.direction.x);
                let (ytmin, ytmax) = check_axis(local_ray.origin.y, local_ray.direction.y);
                let (ztmin, ztmax) = check_axis(local_ray.origin.z, local_ray.direction.z);

                let tmin = xtmin.max(ytmin).max(ztmin);
                let tmax = xtmax.min(ytmax).min(ztmax);

                if tmax >= tmin {
                    result.push(Intersection::new(tmin, &self));
                    result.push(Intersection::new(tmax, &self));
                }
            }
            ShapeKind::Cylinder {
                minimum, maximum, ..
            } => {
                let a = local_ray.direction.x.powi(2) + local_ray.direction.z.powi(2);

                // Intersect walls
                if !is_almost_equal(a, 0.0) {
                    let b = 2.0 * local_ray.origin.x * local_ray.direction.x
                        + 2.0 * local_ray.origin.z * local_ray.direction.z;
                    let c = local_ray.origin.x.powi(2) + local_ray.origin.z.powi(2) - 1.0;

                    let discriminant = b.powi(2) - 4.0 * a * c;

                    if discriminant >= 0.0 {
                        let sqrt = discriminant.sqrt();
                        let mut t0 = (-b - sqrt) / (2. * a);
                        let mut t1 = (-b + sqrt) / (2. * a);

                        if t0 > t1 {
                            swap(&mut t0, &mut t1);
                        }

                        let y0 = local_ray.origin.y + t0 * local_ray.direction.y;
                        if minimum < y0 && y0 < maximum {
                            result.push(Intersection::new(t0, &self));
                        }

                        let y1 = local_ray.origin.y + t1 * local_ray.direction.y;
                        if minimum < y1 && y1 < maximum {
                            result.push(Intersection::new(t1, &self));
                        }
                    }
                }

                self.intersect_caps(&mut result, &local_ray);
            }
            ShapeKind::Cone {
                minimum, maximum, ..
            } => {
                let a = local_ray.direction.x.powi(2) - local_ray.direction.y.powi(2)
                    + local_ray.direction.z.powi(2);
                let b = 2.0 * local_ray.origin.x * local_ray.direction.x
                    - 2.0 * local_ray.origin.y * local_ray.direction.y
                    + 2.0 * local_ray.origin.z * local_ray.direction.z;
                let c = local_ray.origin.x.powi(2) - local_ray.origin.y.powi(2)
                    + local_ray.origin.z.powi(2);

                // Intersect walls
                if is_almost_equal(a, 0.0) {
                    if !is_almost_equal(b, 0.0) {
                        let t = -c / (2.0 * b);
                        result.push(Intersection::new(t, &self));
                    }
                } else {
                    let discriminant = b.powi(2) - 4.0 * a * c;

                    if discriminant >= 0.0 {
                        let sqrt = discriminant.sqrt();
                        let mut t0 = (-b - sqrt) / (2. * a);
                        let mut t1 = (-b + sqrt) / (2. * a);

                        if t0 > t1 {
                            swap(&mut t0, &mut t1);
                        }

                        let y0 = local_ray.origin.y + t0 * local_ray.direction.y;
                        if minimum < y0 && y0 < maximum {
                            result.push(Intersection::new(t0, &self));
                        }

                        let y1 = local_ray.origin.y + t1 * local_ray.direction.y;
                        if minimum < y1 && y1 < maximum {
                            result.push(Intersection::new(t1, &self));
                        }
                    }
                }

                self.intersect_caps(&mut result, &local_ray);
            }
        };

        result
    }

    fn intersect_caps<'a>(&'a self, intersections: &mut Vec<Intersection<'a>>, local_ray: &Ray) {
        let (maximum, minimum, capped) = match self.kind {
            ShapeKind::Cylinder {
                maximum,
                minimum,
                capped,
            } => (maximum, minimum, capped),
            ShapeKind::Cone {
                maximum,
                minimum,
                capped,
            } => (maximum, minimum, capped),
            _ => panic!("Expected a cylinder or a cone."),
        };

        if !capped {
            return;
        }

        if is_almost_equal(local_ray.direction.y, 0.0) {
            return;
        }

        let t = (minimum - local_ray.origin.y) / local_ray.direction.y;
        // Check for an intersection with the lower end cap by intersecting​ the ray with the plane at
        // y=cyl.minimum​.
        if self.check_cap(&local_ray, t) {
            intersections.push(Intersection::new(t, &self));
        }

        let t = (maximum - local_ray.origin.y) / local_ray.direction.y;
        // Check for an intersection with the upper end cap by intersecting​ the ray with the plane at
        // y=cyl.maximum.
        if self.check_cap(&local_ray, t) {
            intersections.push(Intersection::new(t, &self));
        }
    }

    // A helper function to reduce duplication.​ Checks to see if the intersection at `t` is within a radius​ of y (the
    // radius of your cylinders/cone) from the y axis.​
    //
    // Note: y is always 1 for cylinders.
    fn check_cap(&self, ray: &Ray, t: f64) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let y = ray.origin.y + t * ray.direction.y;
        let z = ray.origin.z + t * ray.direction.z;

        x.powi(2) + z.powi(2) <= y.abs()
    }

    pub fn normal_at(&self, world_point: &Tuple) -> Tuple {
        let shape_inverted_transform = self
            .transform
            .inverse()
            .expect("Transform should be invertible");
        let local_point = shape_inverted_transform * *world_point;

        let local_normal = match self.kind {
            ShapeKind::Sphere => local_point - Tuple::point(0.0, 0.0, 0.0),
            ShapeKind::Plane => Tuple::vector(0.0, 1.0, 0.0),
            ShapeKind::Cube => {
                let x_abs = local_point.x.abs();
                let y_abs = local_point.y.abs();
                let z_abs = local_point.z.abs();

                let maxc = x_abs.max(y_abs).max(z_abs);

                if maxc == x_abs {
                    Tuple::vector(local_point.x, 0.0, 0.0)
                } else if maxc == y_abs {
                    Tuple::vector(0.0, local_point.y, 0.0)
                } else {
                    Tuple::vector(0.0, 0.0, local_point.z)
                }
            }
            ShapeKind::Cylinder {
                minimum, maximum, ..
            } => {
                // Square of the distance from the y axis.
                let dist = local_point.x.powi(2) + local_point.z.powi(2);

                if dist < 1.0 && local_point.y >= maximum - EPSILON {
                    Tuple::vector(0.0, 1.0, 0.0)
                } else if dist < 1.0 && local_point.y <= minimum + EPSILON {
                    Tuple::vector(0.0, -1.0, 0.0)
                } else {
                    Tuple::vector(local_point.x, 0.0, local_point.z)
                }
            }
            ShapeKind::Cone { .. } => {
                let mut y = (local_point.x.powi(2) + local_point.z.powi(2)).sqrt();
                if local_point.y > 0.0 {
                    y = -y;
                }
                Tuple::vector(local_point.x, y, local_point.z)
            }
        };
        assert!(local_normal.is_vector());

        let mut world_normal = shape_inverted_transform.transpose() * local_normal;
        // Hack: Instead of removing any translation by taking a 3x3 submatrix of the transform, we just set w to 0.
        world_normal.w = 0.;

        world_normal.normalize()
    }
}

fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
    let tmin_numerator = (-1.0) - origin;
    let tmax_numerator = 1.0 - origin;

    let mut tmin: f64;
    let mut tmax: f64;
    if direction.abs() >= EPSILON {
        tmin = tmin_numerator / direction;
        tmax = tmax_numerator / direction;
    } else {
        tmin = tmin_numerator * f64::INFINITY;
        tmax = tmax_numerator * f64::INFINITY;
    }

    if tmin > tmax {
        swap(&mut tmin, &mut tmax);
    }

    (tmin, tmax)
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::assert_almost_eq;
    use crate::transformations::*;

    //
    // Sphere
    //

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
        s.material = m.clone();
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
        s.material = m.clone();
        assert_eq!(s.material, m);
    }

    #[test]
    fn a_helper_for_producing_a_sphere_with_a_glassy_material() {
        let s = Shape::glass_sphere();
        assert_eq!(s.transform, Matrix::<4>::identity());
        assert_eq!(s.material.transparency, 1.0);
        assert_eq!(s.material.refractive_index, 1.5);
    }

    //
    // Shape
    //

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

    //
    // Planes
    //

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

    //
    // Cubes
    //

    #[test]
    fn a_ray_intersects_a_cube() {
        // +x
        internal_a_ray_intersects_a_cube(
            Tuple::point(5.0, 0.5, 0.0),
            Tuple::vector(-1.0, 0.0, 0.0),
            4.0,
            6.0,
        );
        // -x
        internal_a_ray_intersects_a_cube(
            Tuple::point(-5.0, 0.5, 0.0),
            Tuple::vector(1.0, 0.0, 0.0),
            4.0,
            6.0,
        );
        // +y
        internal_a_ray_intersects_a_cube(
            Tuple::point(0.5, 5.0, 0.0),
            Tuple::vector(0.0, -1.0, 0.0),
            4.0,
            6.0,
        );
        // -y
        internal_a_ray_intersects_a_cube(
            Tuple::point(0.5, -5.0, 0.0),
            Tuple::vector(0.0, 1.0, 0.0),
            4.0,
            6.0,
        );
        // +z
        internal_a_ray_intersects_a_cube(
            Tuple::point(0.5, 0.0, 5.0),
            Tuple::vector(0.0, 0.0, -1.0),
            4.0,
            6.0,
        );
        // -z
        internal_a_ray_intersects_a_cube(
            Tuple::point(0.5, 0.0, -5.0),
            Tuple::vector(0.0, 0.0, 1.0),
            4.0,
            6.0,
        );
        // inside
        internal_a_ray_intersects_a_cube(
            Tuple::point(0.0, 0.5, 0.0),
            Tuple::vector(0.0, 0.0, 1.0),
            -1.0,
            1.0,
        );
    }

    fn internal_a_ray_intersects_a_cube(origin: Tuple, direction: Tuple, t1: f64, t2: f64) {
        let c = Shape::cube();
        let r = Ray::new(origin, direction);
        let xs = c.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, t1);
        assert_almost_eq!(xs[1].t, t2);
    }

    fn internal_a_ray_misses_a_cube(origin: Tuple, direction: Tuple) {
        let c = Shape::cube();
        let r = Ray::new(origin, direction);
        let xs = c.intersect(&r);
        assert!(xs.is_empty());
    }

    #[test]
    fn a_ray_misses_a_cube() {
        internal_a_ray_misses_a_cube(
            Tuple::point(-2.0, 0.0, 0.0),
            Tuple::vector(0.2673, 0.5345, 0.8018),
        );
        internal_a_ray_misses_a_cube(
            Tuple::point(0.0, -2.0, 0.0),
            Tuple::vector(0.8018, 0.2673, 0.5345),
        );
        internal_a_ray_misses_a_cube(
            Tuple::point(0.0, 0.0, -2.0),
            Tuple::vector(0.5345, 0.8018, 0.2673),
        );
        internal_a_ray_misses_a_cube(Tuple::point(2.0, 0.0, 2.0), Tuple::vector(0.0, 0.0, -1.0));
        internal_a_ray_misses_a_cube(Tuple::point(0.0, 2.0, 2.0), Tuple::vector(0.0, -1.0, 0.0));
        internal_a_ray_misses_a_cube(Tuple::point(2.0, 2.0, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
    }

    fn internal_the_normal_on_the_surface_of_a_cube(point: Tuple, expected_normal: Tuple) {
        let c = Shape::cube();
        let normal = c.normal_at(&point);
        assert_eq!(normal, expected_normal);
    }

    #[test]
    fn the_normal_on_the_surface_of_a_cube() {
        internal_the_normal_on_the_surface_of_a_cube(
            Tuple::point(1.0, 0.5, -0.8),
            Tuple::vector(1.0, 0.0, 0.0),
        );
        internal_the_normal_on_the_surface_of_a_cube(
            Tuple::point(-1.0, -0.2, 0.9),
            Tuple::vector(-1.0, 0.0, 0.0),
        );
        internal_the_normal_on_the_surface_of_a_cube(
            Tuple::point(-0.4, 1.0, -0.1),
            Tuple::vector(0.0, 1.0, 0.0),
        );
        internal_the_normal_on_the_surface_of_a_cube(
            Tuple::point(0.3, -1.0, -0.7),
            Tuple::vector(0.0, -1.0, 0.0),
        );
        internal_the_normal_on_the_surface_of_a_cube(
            Tuple::point(-0.6, 0.3, 1.0),
            Tuple::vector(0.0, 0.0, 1.0),
        );
        internal_the_normal_on_the_surface_of_a_cube(
            Tuple::point(0.4, 0.4, -1.0),
            Tuple::vector(0.0, 0.0, -1.0),
        );
        internal_the_normal_on_the_surface_of_a_cube(
            Tuple::point(1.0, 1.0, 1.0),
            Tuple::vector(1.0, 0.0, 0.0),
        );
        internal_the_normal_on_the_surface_of_a_cube(
            Tuple::point(-1.0, -1.0, -1.0),
            Tuple::vector(-1.0, 0.0, 0.0),
        );
    }

    //
    // Feature: Cylinders
    //

    fn internal_a_ray_misses_a_cylinder(origin: Tuple, direction: Tuple) {
        assert!(origin.is_point());
        assert!(direction.is_vector());

        let r = Ray::new(origin, direction.normalize());
        let cyl = Shape::infinite_cylinder();

        let xs = cyl.intersect(&r);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn a_ray_misses_a_cylinder() {
        internal_a_ray_misses_a_cylinder(Tuple::point(1., 0., 0.), Tuple::vector(0., 1., 0.));
        internal_a_ray_misses_a_cylinder(Tuple::point(0., 0., 0.), Tuple::vector(0., 1., 0.));
        internal_a_ray_misses_a_cylinder(Tuple::point(0., 0., -5.), Tuple::vector(1., 1., 1.));
    }

    fn internal_a_ray_strikes_a_cylinder(origin: Tuple, direction: Tuple, t0: f64, t1: f64) {
        let cyl = Shape::infinite_cylinder();
        let ray = Ray::new(origin, direction.normalize());
        let xs = cyl.intersect(&ray);

        assert_eq!(xs.len(), 2);
        assert_almost_eq!(xs[0].t, t0);
        assert_almost_eq!(xs[1].t, t1);
    }

    #[test]
    fn a_ray_strikes_a_cylinder() {
        internal_a_ray_strikes_a_cylinder(
            Tuple::point(1.0, 0.0, -5.0),
            Tuple::vector(0.0, 0.0, 1.0),
            5.0,
            5.0,
        );
        internal_a_ray_strikes_a_cylinder(
            Tuple::point(0.0, 0.0, -5.0),
            Tuple::vector(0.0, 0.0, 1.0),
            4.0,
            6.0,
        );
        internal_a_ray_strikes_a_cylinder(
            Tuple::point(0.5, 0.0, -5.0),
            Tuple::vector(0.1, 1.0, 1.0),
            6.80798,
            7.08872,
        );
    }

    #[test]
    fn normal_vector_on_a_cylinder() {
        let cyl = Shape::infinite_cylinder();
        assert_eq!(
            cyl.normal_at(&Tuple::point(1.0, 0.0, 0.0)),
            Tuple::vector(1.0, 0.0, 0.0)
        );
        assert_eq!(
            cyl.normal_at(&Tuple::point(0.0, 5.0, -1.0)),
            Tuple::vector(0.0, 0.0, -1.0)
        );
        assert_eq!(
            cyl.normal_at(&Tuple::point(0.0, -2.0, 1.0)),
            Tuple::vector(0.0, 0.0, 1.0)
        );
        assert_eq!(
            cyl.normal_at(&Tuple::point(-1.0, 1.0, 0.0)),
            Tuple::vector(-1.0, 0.0, 0.0)
        );
    }

    #[test]
    fn the_default_minimum_and_maximum_for_a_cylinder() {
        let cyl = Shape::infinite_cylinder();

        if let ShapeKind::Cylinder {
            minimum, maximum, ..
        } = cyl.kind
        {
            assert_eq!(minimum, -f64::INFINITY);
            assert_eq!(maximum, f64::INFINITY);
        } else {
            assert!(false);
        }
    }

    fn internal_intersecting_a_constrained_cylinder(ray: Ray) -> usize {
        let cyl = Shape::cylinder(1.0, 2.0, false);

        let xs = cyl.intersect(&ray);
        xs.len()
    }

    #[test]
    fn intersecting_a_constrained_cylinder() {
        // Ray shooting from within the cylinder, but going out the unclosed top.
        assert_eq!(
            internal_intersecting_a_constrained_cylinder(Ray::new(
                Tuple::point(0.0, 1.5, 0.0),
                Tuple::vector(0.1, 1.0, 0.0)
            )),
            0
        );
        // Ray passing over
        assert_eq!(
            internal_intersecting_a_constrained_cylinder(Ray::new(
                Tuple::point(0.0, 3.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0)
            )),
            0
        );
        // Ray passing under
        assert_eq!(
            internal_intersecting_a_constrained_cylinder(Ray::new(
                Tuple::point(0.0, 0.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0)
            )),
            0
        );
        // Ray on top
        assert_eq!(
            internal_intersecting_a_constrained_cylinder(Ray::new(
                Tuple::point(0.0, 2.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0)
            )),
            0
        );
        // Ray on bottom
        assert_eq!(
            internal_intersecting_a_constrained_cylinder(Ray::new(
                Tuple::point(0.0, 1.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0)
            )),
            0
        );
        // Ray in the middle
        assert_eq!(
            internal_intersecting_a_constrained_cylinder(Ray::new(
                Tuple::point(0.0, 1.5, -2.0),
                Tuple::vector(0.0, 0.0, 1.0)
            )),
            2
        );
    }

    #[test]
    fn the_default_closed_value_for_a_cylinder() {
        let cyl = Shape::infinite_cylinder();

        if let ShapeKind::Cylinder { capped, .. } = cyl.kind {
            assert!(!capped);
        } else {
            assert!(false);
        }
    }

    fn internal_intersecting_the_caps_of_a_closed_cylinder(
        point: Tuple,
        direction: Tuple,
        count: usize,
    ) {
        let cyl = Shape::cylinder(1.0, 2.0, true);

        let r = Ray::new(point, direction.normalize());
        assert_eq!(cyl.intersect(&r).len(), count);
    }

    #[test]
    fn intersecting_the_caps_of_a_closed_cylinder() {
        internal_intersecting_the_caps_of_a_closed_cylinder(
            Tuple::point(0.0, 3.0, 0.0),
            Tuple::vector(0.0, -1.0, 0.0),
            2,
        );
        internal_intersecting_the_caps_of_a_closed_cylinder(
            Tuple::point(0.0, 3.0, -2.0),
            Tuple::vector(0.0, -1.0, 2.0),
            2,
        );
        // corner case
        internal_intersecting_the_caps_of_a_closed_cylinder(
            Tuple::point(0.0, 4.0, -2.0),
            Tuple::vector(0.0, -1.0, 1.0),
            2,
        );
        internal_intersecting_the_caps_of_a_closed_cylinder(
            Tuple::point(0.0, 0.0, -2.0),
            Tuple::vector(0.0, 1.0, 2.0),
            2,
        );
        // corner case
        internal_intersecting_the_caps_of_a_closed_cylinder(
            Tuple::point(0.0, -1.0, -2.0),
            Tuple::vector(0.0, 1.0, 1.0),
            2,
        );
    }

    #[test]
    fn the_normal_vector_on_a_cylinder_s_end_caps() {
        let test_it = |point: Tuple, normal: Tuple| {
            let cyl = Shape::cylinder(1.0, 2.0, true);
            let n = cyl.normal_at(&point);

            assert_eq!(n, normal);
        };

        test_it(Tuple::point(0.0, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        test_it(Tuple::point(0.5, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        test_it(Tuple::point(0.0, 1.0, 0.5), Tuple::vector(0.0, -1.0, 0.0));
        test_it(Tuple::point(0.0, 2.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        test_it(Tuple::point(0.5, 2.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        test_it(Tuple::point(0.0, 2.0, 0.5), Tuple::vector(0.0, 1.0, 0.0));
    }

    //
    // Feature: Cones
    //

    #[test]
    fn intersecting_a_cone_with_a_ray() {
        let test_it = |origin: Tuple, direction: Tuple, t0: f64, t1: f64| {
            let shape = Shape::infinite_cone();
            let direction = direction.normalize();
            let ray = Ray::new(origin, direction);
            let xs = shape.intersect(&ray);

            assert_eq!(xs.len(), 2);
            assert_almost_eq!(xs[0].t, t0);
            assert_almost_eq!(xs[1].t, t1);
        };

        test_it(
            Tuple::point(0.0, 0.0, -5.0),
            Tuple::vector(0.0, 0.0, 1.0),
            5.0,
            5.0,
        );
        test_it(
            Tuple::point(0.0, 0.0, -5.0),
            Tuple::vector(1.0, 1.0, 1.0),
            8.66025,
            8.66025,
        );
        test_it(
            Tuple::point(1.0, 1.0, -5.0),
            Tuple::vector(-0.5, -1.0, 1.0),
            4.55006,
            49.44994,
        );
    }

    #[test]
    fn intersecting_a_cone_with_a_ray_parallel_to_one_of_its_halves() {
        let shape = Shape::infinite_cone();
        let direction = Tuple::vector(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(0.0, 0.0, -1.0), direction);
        let xs = shape.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_almost_eq!(xs[0].t, 0.35355);
    }

    #[test]
    fn intersecting_a_cone_s_end_caps() {
        let test_it = |origin: Tuple, direction: Tuple, count: usize| {
            let shape = Shape::cone(-0.5, 0.5, true);
            let direction = direction.normalize();
            let r = Ray::new(origin, direction);
            let xs = shape.intersect(&r);
            assert_eq!(xs.len(), count);
        };

        test_it(
            Tuple::point(0.0, 0.0, -5.0),
            Tuple::vector(0.0, 1.0, 0.0),
            0,
        );
        test_it(
            Tuple::point(0.0, 0.0, -0.25),
            Tuple::vector(0.0, 1.0, 1.0),
            2,
        );
        test_it(
            Tuple::point(0.0, 0.0, -0.25),
            Tuple::vector(0.0, 1.0, 0.0),
            4,
        );
    }

    #[test]
    fn computing_the_normal_vector_on_a_cone() {
        let test_it = |point: Tuple, normal: Tuple| {
            let shape = Shape::infinite_cone();
            let n = shape.normal_at(&point);
            assert_eq!(n, normal.normalize()); // Normalize, because in the book it's testing local_normal_at.
        };

        test_it(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 0.0));
        test_it(
            Tuple::point(1.0, 1.0, 1.0),
            Tuple::vector(1.0, -(2.0 as f64).sqrt(), 1.0),
        );
        test_it(Tuple::point(-1.0, -1.0, 0.0), Tuple::vector(-1.0, 1.0, 0.0));
    }
}
