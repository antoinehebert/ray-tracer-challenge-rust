use crate::{
    intersection::Intersections, material::Material, matrix::Matrix, ray::Ray, tuple::Tuple,
};

pub trait Shape {
    // TODO: Is there a way to have fields in traits, or to DRY because all shapes will have a transform and a material.
    fn transform(&self) -> &Matrix<4>;
    fn set_transform(&mut self, transform: Matrix<4>);

    fn material(&self) -> &Material;
    fn set_material(&mut self, material: Material);

    // Returns intersection points (time) along `ray`.
    fn intersect(&self, world_ray: &Ray) -> Intersections {
        let local_ray = world_ray.transform(
            self.transform()
                .inverse()
                .expect("shape transfor should be invertible"),
        );
        self.local_intersect(&local_ray)
    }

    fn normal_at(&self, world_point: &Tuple) -> Tuple {
        let sphere_inverted_transform = self
            .transform()
            .inverse()
            .expect("Transform should be invertible");
        let local_point = sphere_inverted_transform * *world_point;

        let local_normal = self.local_normal_at(&local_point);

        let mut world_normal = sphere_inverted_transform.transpose() * local_normal;
        // Hack: Instead of removing any translation by taking a 3x3 submatrix of the transform, we just set w to 0.
        world_normal.w = 0.;

        world_normal.normalize()
    }

    //
    // Private:
    // TODO: Can we have private methods on traits?
    //
    fn local_intersect(&self, local_ray: &Ray) -> Intersections;
    fn local_normal_at(&self, local_point: &Tuple) -> Tuple;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::Matrix;
    use crate::transformations::*;
    use crate::tuple::Tuple;
    use std::f64::consts::PI;

    // Can't set this in TestShape directly, because local_intersect can't mutate self...
    // TODO: Figure out if there's a better way of doing this in Rust.
    static mut saved_ray: Option<Ray> = None;

    struct TestShape {
        transform: Matrix<4>,
        material: Material,
    }

    impl Shape for TestShape {
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

        fn local_intersect(&self, ray: &Ray) -> Intersections {
            unsafe {
                saved_ray = Some(ray.clone());
            }
            Vec::new()
        }

        fn local_normal_at(&self, local_point: &Tuple) -> Tuple {
            local_point.clone()
        }
    }

    fn test_shape() -> TestShape {
        TestShape {
            transform: Matrix::identity(),
            material: Material::new(),
        }
    }

    // Is this really testing something meaningful? As far as I can tell it's impossible have an abstract parent class
    // like in C++, that would hold a common field and default implementation...
    #[test]
    fn the_default_transformation() {
        let s = test_shape();
        assert_eq!(*s.transform(), Matrix::identity());
    }

    #[test]
    fn assigning_a_transformation() {
        let mut s = test_shape();
        s.set_transform(translation(2.0, 3.0, 4.0));

        assert_eq!(*s.transform(), translation(2.0, 3.0, 4.0));
    }

    #[test]
    fn the_default_material() {
        let s = test_shape();
        let m = s.material();

        assert_eq!(*m, Material::new());
    }

    #[test]
    fn assigning_a_material() {
        let mut s = test_shape();
        let mut m = Material::new();
        m.ambient = 1.0;
        s.material = m;
        assert_eq!(*s.material(), m);
    }

    #[test]
    fn intersecting_a_scaled_shape_with_a_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut s = test_shape();
        s.set_transform(scaling(2.0, 2.0, 2.0));
        let _xs = s.intersect(&r);

        unsafe {
            assert_eq!(saved_ray.unwrap().origin, Tuple::point(0.0, 0.0, -2.5));
            assert_eq!(saved_ray.unwrap().direction, Tuple::vector(0.0, 0.0, 0.5));
        }
    }

    #[test]
    fn intersecting_a_translated_shape_with_a_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut s = test_shape();
        s.set_transform(translation(5.0, 0.0, 0.0));
        let _xs = s.intersect(&r);

        unsafe {
            assert_eq!(saved_ray.unwrap().origin, Tuple::point(-5.0, 0.0, -5.0));
            assert_eq!(saved_ray.unwrap().direction, Tuple::vector(0.0, 0.0, 1.0));
        }
    }

    #[test]
    fn computing_the_normal_on_a_translated_shape() {
        let mut s = test_shape();
        s.set_transform(translation(0.0, 1.0, 0.0));
        let n = s.normal_at(&Tuple::point(0.0, 1.70711, -0.70711));
        assert_eq!(n, Tuple::vector(0.0, 0.70711, -0.70711));
    }

    #[test]
    fn computing_the_normal_on_a_transformed_shape() {
        let mut s = test_shape();
        let m = scaling(1.0, 0.5, 1.0) * rotation_z(PI / 5.0);
        s.set_transform(m);
        let n = s.normal_at(&Tuple::point(
            0.0,
            (2.0 as f64).sqrt() / 2.0,
            -(2.0 as f64).sqrt() / 2.0,
        ));
        assert_eq!(n, Tuple::vector(0.0, 0.97014, -0.24254));
    }

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
}
