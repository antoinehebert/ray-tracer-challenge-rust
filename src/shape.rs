use std::mem::swap;

use uuid::Uuid;

use crate::{
    intersection::{Intersection, Intersections},
    material::Material,
    matrix::Matrix,
    ray::Ray,
    tuple::Tuple,
    utils::{is_almost_equal, EPSILON},
    world::World,
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
    Group {
        uuid: Uuid,
        shapes: Vec<Shape>,
    },
}
#[derive(Debug, PartialEq, Clone)]
pub struct Shape {
    kind: ShapeKind,

    pub parent_id: Option<Uuid>,
    pub transform: Matrix<4>,
    pub material: Material,
}

impl Shape {
    // Returns intersection points (time) along `ray`.
    pub fn intersect(&self, world_ray: &Ray) -> Intersections {
        let local_ray = world_ray.transform(
            self.transform
                .inverse()
                .expect("shape transform should be invertible"),
        );

        let mut result = Vec::new();

        match &self.kind {
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
                        if *minimum < y0 && y0 < *maximum {
                            result.push(Intersection::new(t0, &self));
                        }

                        let y1 = local_ray.origin.y + t1 * local_ray.direction.y;
                        if *minimum < y1 && y1 < *maximum {
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
                        if *minimum < y0 && y0 < *maximum {
                            result.push(Intersection::new(t0, &self));
                        }

                        let y1 = local_ray.origin.y + t1 * local_ray.direction.y;
                        if *minimum < y1 && y1 < *maximum {
                            result.push(Intersection::new(t1, &self));
                        }
                    }
                }

                self.intersect_caps(&mut result, &local_ray);
            }
            ShapeKind::Group { shapes, .. } => {
                let mut shape_results: Intersections = shapes
                    .into_iter()
                    .flat_map(|shape| shape.intersect(&local_ray))
                    .collect();
                shape_results
                    .sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));

                result.append(&mut shape_results);
            }
        };

        result
    }

    pub fn normal_at(&self, world_point: &Tuple, world: Option<&World>) -> Tuple {
        let local_point = world_to_object(self, world_point, world);

        let local_normal = match &self.kind {
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
            ShapeKind::Group { .. } => {
                unreachable!();
            }
        };
        assert!(local_normal.is_vector());

        let mut world_normal = normal_to_world(self, &local_normal, world);
        // Hack: Instead of removing any translation by taking a 3x3 submatrix of the transform, we just set w to 0.
        world_normal.w = 0.;

        world_normal.normalize()
    }

    pub fn sphere() -> Self {
        Self {
            kind: ShapeKind::Sphere,
            parent_id: None,
            transform: Matrix::<4>::identity(),
            material: Material::new(),
        }
    }

    pub fn glass_sphere() -> Self {
        let mut material = Material::new();
        material.transparency = 1.0;
        material.refractive_index = 1.5;

        Self {
            kind: ShapeKind::Sphere,
            parent_id: None,
            transform: Matrix::<4>::identity(),
            material: material,
        }
    }

    pub fn plane() -> Self {
        Self {
            kind: ShapeKind::Plane,
            parent_id: None,
            transform: Matrix::<4>::identity(),
            material: Material::new(),
        }
    }

    pub fn cube() -> Self {
        Self {
            kind: ShapeKind::Cube,
            parent_id: None,
            transform: Matrix::<4>::identity(),
            material: Material::new(),
        }
    }

    pub fn infinite_cylinder() -> Self {
        Self {
            kind: ShapeKind::Cylinder {
                minimum: -f64::INFINITY,
                maximum: f64::INFINITY,
                capped: false,
            },
            parent_id: None,
            transform: Matrix::<4>::identity(),
            material: Material::new(),
        }
    }

    pub fn cylinder(min: f64, max: f64, capped: bool) -> Self {
        Self {
            kind: ShapeKind::Cylinder {
                minimum: min,
                maximum: max,
                capped,
            },
            parent_id: None,
            transform: Matrix::<4>::identity(),
            material: Material::new(),
        }
    }

    pub fn infinite_cone() -> Self {
        Self {
            kind: ShapeKind::Cone {
                minimum: -f64::INFINITY,
                maximum: f64::INFINITY,
                capped: false,
            },
            parent_id: None,
            transform: Matrix::<4>::identity(),
            material: Material::new(),
        }
    }

    pub fn cone(min: f64, max: f64, capped: bool) -> Self {
        Self {
            kind: ShapeKind::Cone {
                minimum: min,
                maximum: max,
                capped,
            },
            parent_id: None,
            transform: Matrix::<4>::identity(),
            material: Material::new(),
        }
    }

    pub fn group() -> Self {
        Self {
            kind: ShapeKind::Group {
                shapes: vec![],
                uuid: Uuid::new_v4(),
            },
            parent_id: None,
            transform: Matrix::<4>::identity(),
            material: Material::new(),
        }
    }

    // Maybe passing shapes in the constructor would work too?
    pub fn add_child(&mut self, mut shape: Shape) {
        match &mut self.kind {
            ShapeKind::Group { shapes, uuid } => {
                shape.parent_id = Some(uuid.clone());
                shapes.push(shape);
            }
            _ => panic!(),
        }
    }

    pub fn get_shape_by_id(&self, id: Uuid) -> Option<&Self> {
        match &self.kind {
            ShapeKind::Group { uuid, shapes } => {
                if id == *uuid {
                    return Some(self);
                }

                for shape in shapes {
                    match shape.get_shape_by_id(id) {
                        Some(shape) => return Some(shape),
                        None => (),
                    }
                }
            }
            _ => (),
        }

        None
    }

    pub fn shapes(&self) -> Option<&Vec<Shape>> {
        match &self.kind {
            ShapeKind::Group { shapes, .. } => Some(&shapes),
            _ => None, // It would be nice to return an empty vector here instead, so callers wouldn't have to unwrap.
        }
    }

    fn intersect_caps<'b>(&'b self, intersections: &mut Vec<Intersection<'b>>, local_ray: &Ray) {
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

// TODO: HACK: Passing in the world here because this his how I can find parents. Maybe using Rc or unsafe pointers
// would make my life easier...
fn world_to_object(s: &Shape, point: &Tuple, world: Option<&World>) -> Tuple {
    let mut point = *point;

    match s.parent_id {
        Some(parent_id) => {
            let world = world.expect("You can't have a parent without passing in a world.");
            if let Some(parent) = world.get_shape_by_id(parent_id) {
                point = world_to_object(parent, &point, Some(world)).clone();
            }
        }
        None => (),
    };

    s.transform.inverse().expect("Should be invertible") * point
}

// TODO: HACK: Passing in the world here because this his how I can find parents. Maybe using Rc or unsafe pointers
// would make my life easier...
fn normal_to_world(shape: &Shape, normal: &Tuple, world: Option<&World>) -> Tuple {
    let mut normal = shape.transform.inverse().unwrap().transpose() * *normal;
    normal.w = 0.0;
    normal = normal.normalize();

    if let Some(uuid) = shape.parent_id {
        let world = world.expect("You can't have a parent without passing in a world.");

        let parent = world.get_shape_by_id(uuid);
        if let Some(parent) = parent {
            normal = normal_to_world(parent, &normal, Some(world));
        }
    }
    normal
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::assert_almost_eq;
    use crate::color::WHITE;
    use crate::light::Light;
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
        let n = s.normal_at(&Tuple::point(1., 0., 0.), None);
        assert_eq!(n, Tuple::vector(1., 0., 0.));
    }

    #[test]
    fn the_normal_on_a_sphere_at_a_point_on_the_y_axis() {
        let s = Shape::sphere();
        let n = s.normal_at(&Tuple::point(0., 1., 0.), None);
        assert_eq!(n, Tuple::vector(0., 1., 0.));
    }

    #[test]
    fn the_normal_on_a_sphere_at_a_point_on_the_z_axis() {
        let s = Shape::sphere();
        let n = s.normal_at(&Tuple::point(0., 0., 1.), None);
        assert_eq!(n, Tuple::vector(0., 0., 1.));
    }

    #[test]
    fn the_normal_on_a_sphere_at_a_nonaxial_point() {
        let s = Shape::sphere();
        let n = s.normal_at(
            &Tuple::point(
                (3. as f64).sqrt() / 3.,
                (3. as f64).sqrt() / 3.,
                (3. as f64).sqrt() / 3.,
            ),
            None,
        );
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
        let n = s.normal_at(
            &Tuple::point(
                (3. as f64).sqrt() / 3.,
                (3. as f64).sqrt() / 3.,
                (3. as f64).sqrt() / 3.,
            ),
            None,
        );
        assert_eq!(n, n.normalize());
    }

    #[test]
    fn computing_the_normal_on_a_translated_sphere() {
        let mut s = Shape::sphere();
        s.transform = translation(0., 1., 0.);

        let n = s.normal_at(&Tuple::point(0., 1.70711, -0.70711), None);
        assert_eq!(n, Tuple::vector(0., 0.70711, -0.70711));
    }

    #[test]
    fn computing_the_normal_on_a_transformed_sphere() {
        let mut s = Shape::sphere();
        let m = scaling(1., 0.5, 1.) * rotation_z(PI / 5.);
        s.transform = m;
        let n = s.normal_at(
            &Tuple::point(0., (2 as f64).sqrt() / 2., -(2 as f64).sqrt() / 2.),
            None,
        );
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

    #[test]
    fn a_shape_has_a_parent_attribute() {
        let s = Shape::sphere();
        assert!(s.parent_id.is_none());
    }

    #[test]
    fn converting_a_point_from_world_to_object_space() {
        let mut g1 = Shape::group();
        g1.transform = rotation_y(PI / 2.);
        let mut g2 = Shape::group();
        g2.transform = scaling(2., 2., 2.);
        let mut s = Shape::sphere();
        s.transform = translation(5., 0., 0.);
        g2.add_child(s);
        g1.add_child(g2);

        let light = Light::new(Tuple::point(0.0, 0.0, 0.0), WHITE);
        let mut w = World::new(light);
        w.objects = vec![g1];

        // FIXME: Yuk!
        let s = &w.objects[0].shapes().unwrap()[0].shapes().unwrap()[0];
        let p = world_to_object(&s, &Tuple::point(-2., 0., -10.), Some(&w));
        assert_eq!(p, Tuple::point(0., 0., -1.));
    }

    #[test]
    fn converting_a_normal_from_object_to_world_space() {
        let mut g2 = Shape::group();
        g2.transform = scaling(1., 2., 3.);
        let mut s = Shape::sphere();
        s.transform = translation(5., 0., 0.);
        g2.add_child(s);
        let mut g1 = Shape::group();
        g1.transform = rotation_y(PI / 2.);
        g1.add_child(g2);

        let light = Light::new(Tuple::point(0.0, 0.0, 0.0), WHITE);
        let mut w = World::new(light);
        w.objects = vec![g1];

        // FIXME: Yuk!
        let s = &w.objects[0].shapes().unwrap()[0].shapes().unwrap()[0];

        let n = normal_to_world(
            &s,
            &Tuple::vector(
                (3 as f64).sqrt() / 3.,
                (3 as f64).sqrt() / 3.,
                (3 as f64).sqrt() / 3.,
            ),
            Some(&w),
        );
        assert_eq!(n, Tuple::vector(0.28571, 0.42857, -0.85714));
    }

    #[test]
    fn finding_the_normal_on_a_child_object() {
        let mut g2 = Shape::group();
        g2.transform = scaling(1., 2., 3.);
        let mut s = Shape::sphere();
        s.transform = translation(5., 0., 0.);
        g2.add_child(s);
        let mut g1 = Shape::group();
        g1.transform = rotation_y(PI / 2.);
        g1.add_child(g2);

        // FIXME: Yuk!
        let light = Light::new(Tuple::point(0.0, 0.0, 0.0), WHITE);
        let mut w = World::new(light);
        w.objects = vec![g1];
        let s = &w.objects[0].shapes().unwrap()[0].shapes().unwrap()[0];

        let n = s.normal_at(&Tuple::point(1.7321, 1.1547, -5.5774), Some(&w));
        assert_eq!(n, Tuple::vector(0.28570, 0.42854, -0.85716));
    }
    //
    // Planes
    //

    #[test]
    fn the_normal_of_a_plane_is_constant_everywhere() {
        let p = Shape::plane();

        let n1 = p.normal_at(&Tuple::point(0.0, 0.0, 0.0), None);
        let n2 = p.normal_at(&Tuple::point(10.0, 0.0, -10.0), None);
        let n3 = p.normal_at(&Tuple::point(-5.0, 0.0, 150.0), None);

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
        let normal = c.normal_at(&point, None);
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
            cyl.normal_at(&Tuple::point(1.0, 0.0, 0.0), None),
            Tuple::vector(1.0, 0.0, 0.0)
        );
        assert_eq!(
            cyl.normal_at(&Tuple::point(0.0, 5.0, -1.0), None),
            Tuple::vector(0.0, 0.0, -1.0)
        );
        assert_eq!(
            cyl.normal_at(&Tuple::point(0.0, -2.0, 1.0), None),
            Tuple::vector(0.0, 0.0, 1.0)
        );
        assert_eq!(
            cyl.normal_at(&Tuple::point(-1.0, 1.0, 0.0), None),
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
            let n = cyl.normal_at(&point, None);

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
            let n = shape.normal_at(&point, None);
            assert_eq!(n, normal.normalize()); // Normalize, because in the book it's testing local_normal_at.
        };

        test_it(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 0.0));
        test_it(
            Tuple::point(1.0, 1.0, 1.0),
            Tuple::vector(1.0, -(2.0 as f64).sqrt(), 1.0),
        );
        test_it(Tuple::point(-1.0, -1.0, 0.0), Tuple::vector(-1.0, 1.0, 0.0));
    }

    //
    //Feature: Groups
    //

    #[test]
    fn creating_a_new_group() {
        let g = Shape::group();

        assert_eq!(g.transform, Matrix::<4>::identity());
        assert_eq!(g.shapes().unwrap().len(), 0);
    }

    #[test]
    fn adding_a_child_to_a_group() {
        let mut g = Shape::group();
        let s = Shape::sphere();
        g.add_child(s);

        let shapes = g.shapes().unwrap();
        assert_eq!(shapes.len(), 1);

        match g.kind {
            ShapeKind::Group { uuid, .. } => assert_eq!(shapes[0].parent_id.unwrap(), uuid),
            _ => panic!(),
        }
    }

    #[test]
    fn intersecting_a_ray_with_an_empty_group() {
        let g = Shape::group();
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 0.0));

        let xs = g.intersect(&r);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn intersecting_a_ray_with_a_nonempty_group() {
        let mut g = Shape::group();

        let s1 = Shape::sphere();
        let mut s2 = Shape::sphere();
        s2.transform = translation(0., 0., -3.);
        let mut s3 = Shape::sphere();
        s3.transform = translation(5., 0., 0.);

        g.add_child(s1);
        g.add_child(s2);
        g.add_child(s3);

        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));

        let xs = g.intersect(&r);
        assert_eq!(xs.len(), 4);
        assert_eq!(*xs[0].object, g.shapes().unwrap()[1]); // s2
        assert_eq!(*xs[1].object, g.shapes().unwrap()[1]); // s2
        assert_eq!(*xs[1].object, g.shapes().unwrap()[1]); // s1
        assert_eq!(*xs[3].object, g.shapes().unwrap()[0]); // s1
    }

    #[test]
    fn intersecting_a_transformed_group() {
        let mut g = Shape::group();
        g.transform = scaling(2., 2., 2.);
        let mut s = Shape::sphere();
        s.transform = translation(5., 0., 0.);
        g.add_child(s);
        let r = Ray::new(Tuple::point(10., 0., -10.), Tuple::vector(0., 0., 1.));

        let xs = g.intersect(&r);
        assert_eq!(xs.len(), 2);
    }
}
