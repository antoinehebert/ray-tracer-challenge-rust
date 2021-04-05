use crate::{
    color::{Color, BLACK, WHITE},
    intersection::*,
    light::Light,
    material::Material,
    ray::Ray,
    sphere::Sphere,
    transformations::scaling,
    tuple::Tuple,
};

pub struct World {
    pub objects: Vec<Sphere>,
    pub light: Light,
}

impl World {
    pub fn new(light: Light) -> Self {
        Self {
            objects: Vec::new(),
            light: light,
        }
    }

    pub fn default_world() -> Self {
        let light = Light::new(Tuple::point(-10.0, 10.0, -10.0), WHITE);

        let mut material = Material::new();
        material.color = Color::new(0.8, 1.0, 0.6);
        material.diffuse = 0.7;
        material.specular = 0.2;
        let mut s1 = Sphere::new();
        s1.material = material;

        let mut s2 = Sphere::new();
        s2.transform(scaling(0.5, 0.5, 0.5));

        Self {
            objects: vec![s1, s2],
            light: light,
        }
    }

    fn intersect<'a, 'b>(&'a self, ray: &'b Ray) -> Intersections<'a> {
        let mut result = Intersections::new();

        for object in &self.objects {
            let mut xs = object.intersect(&ray);
            result.append(&mut xs);
        }

        result.sort_by(|x, y| x.t.partial_cmp(&y.t).unwrap_or(std::cmp::Ordering::Equal));

        result
    }

    fn shade_hit(&self, comps: &Computations) -> Color {
        comps
            .object
            .material
            .lighting(&self.light, &comps.point, &comps.eyev, &comps.normalv)
    }

    pub fn color_at(&self, ray: &Ray) -> Color {
        let intersections = self.intersect(&ray);
        let xs = hit(&intersections);

        match xs {
            None => BLACK,
            Some(intersection) => self.shade_hit(&intersection.prepare_computations(&ray)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_almost_eq;
    use crate::utils::*;

    #[test]
    fn creating_a_world() {
        let light = Light::new(Tuple::point(-0.0, 0.0, 0.0), WHITE);
        let w = World::new(light.clone());
        assert_eq!(w.objects.len(), 0);
        assert_eq!(w.light, light);
    }

    #[test]
    fn the_default_world() {
        let light = Light::new(Tuple::point(-10.0, 10.0, -10.0), WHITE);

        let mut material = Material::new();
        material.color = Color::new(0.8, 1.0, 0.6);
        material.diffuse = 0.7;
        material.specular = 0.2;
        let mut s1 = Sphere::new();
        s1.material = material;

        let mut s2 = Sphere::new();
        s2.transform(scaling(0.5, 0.5, 0.5));

        let w = World::default_world();

        assert_eq!(w.light, light);
        assert!(w.objects.contains(&s1));
        assert!(w.objects.contains(&s2));
    }

    #[test]
    fn intersect_a_world_with_a_ray() {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = w.intersect(&r);
        assert_eq!(xs.len(), 4);
        assert_almost_eq!(xs[0].t, 4.0);
        assert_almost_eq!(xs[1].t, 4.5);
        assert_almost_eq!(xs[2].t, 5.5);
        assert_almost_eq!(xs[3].t, 6.0);
    }

    #[test]
    fn shading_an_intersection() {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = &w.objects[0];
        let i = Intersection::new(4.0, shape);
        let comps = i.prepare_computations(&r);
        let c = w.shade_hit(&comps);

        assert_eq!(c, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shading_an_intersection_from_the_inside() {
        let mut w = World::default_world();
        w.light = Light::new(Tuple::point(0.0, 0.25, 0.0), WHITE);
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = &w.objects[1];
        let i = Intersection::new(0.5, shape);
        let comps = i.prepare_computations(&r);
        let c = w.shade_hit(&comps);

        assert_eq!(c, Color::new(0.90498, 0.90498, 0.90498));
    }

    #[test]
    fn the_color_when_a_ray_misses() {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 1.0, 0.0));
        let c = w.color_at(&r);
        assert_eq!(c, Color::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn the_color_when_a_ray_hits() {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let c = w.color_at(&r);
        assert_eq!(c, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn the_color_with_an_intersection_behind_the_ray() {
        let mut w = World::default_world();
        w.objects[0].material.ambient = 1.0;
        w.objects[1].material.ambient = 1.0;
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.75), Tuple::vector(0.0, 0.0, -1.0));
        let c = w.color_at(&r);
        assert_eq!(c, w.objects[1].material.color);
    }
    // Scenario: There is no shadow when nothing is collinear with point and light
    //   Given w â† default_world()
    //     And p â† point(0, 10, 0)
    //    Then is_shadowed(w, p) is false

    // Scenario: The shadow when an object is between the point and the light
    //   Given w â† default_world()
    //     And p â† point(10, -10, 10)
    //    Then is_shadowed(w, p) is true

    // Scenario: There is no shadow when an object is behind the light
    //   Given w â† default_world()
    //     And p â† point(-20, 20, -20)
    //    Then is_shadowed(w, p) is false

    // Scenario: There is no shadow when an object is behind the point
    //   Given w â† default_world()
    //     And p â† point(-2, 2, -2)
    //    Then is_shadowed(w, p) is false

    // Scenario: shade_hit() is given an intersection in shadow
    //   Given w â† world()
    //     And w.light â† point_light(point(0, 0, -10), color(1, 1, 1))
    //     And s1 â† sphere()
    //     And s1 is added to w
    //     And s2 â† sphere() with:
    //       | transform | translation(0, 0, 10) |
    //     And s2 is added to w
    //     And r â† ray(point(0, 0, 5), vector(0, 0, 1))
    //     And i â† intersection(4, s2)
    //   When comps â† prepare_computations(i, r)
    //     And c â† shade_hit(w, comps)
    //   Then c = color(0.1, 0.1, 0.1)

    // Scenario: The reflected color for a nonreflective material
    //   Given w â† default_world()
    //     And r â† ray(point(0, 0, 0), vector(0, 0, 1))
    //     And shape â† the second object in w
    //     And shape.material.ambient â† 1
    //     And i â† intersection(1, shape)
    //   When comps â† prepare_computations(i, r)
    //     And color â† reflected_color(w, comps)
    //   Then color = color(0, 0, 0)

    // Scenario: The reflected color for a reflective material
    //   Given w â† default_world()
    //     And shape â† plane() with:
    //       | material.reflective | 0.5                   |
    //       | transform           | translation(0, -1, 0) |
    //     And shape is added to w
    //     And r â† ray(point(0, 0, -3), vector(0, -âˆš2/2, âˆš2/2))
    //     And i â† intersection(âˆš2, shape)
    //   When comps â† prepare_computations(i, r)
    //     And color â† reflected_color(w, comps)
    //   Then color = color(0.19032, 0.2379, 0.14274)

    // Scenario: shade_hit() with a reflective material
    //   Given w â† default_world()
    //     And shape â† plane() with:
    //       | material.reflective | 0.5                   |
    //       | transform           | translation(0, -1, 0) |
    //     And shape is added to w
    //     And r â† ray(point(0, 0, -3), vector(0, -âˆš2/2, âˆš2/2))
    //     And i â† intersection(âˆš2, shape)
    //   When comps â† prepare_computations(i, r)
    //     And color â† shade_hit(w, comps)
    //   Then color = color(0.87677, 0.92436, 0.82918)

    // Scenario: color_at() with mutually reflective surfaces
    //   Given w â† world()
    //     And w.light â† point_light(point(0, 0, 0), color(1, 1, 1))
    //     And lower â† plane() with:
    //       | material.reflective | 1                     |
    //       | transform           | translation(0, -1, 0) |
    //     And lower is added to w
    //     And upper â† plane() with:
    //       | material.reflective | 1                    |
    //       | transform           | translation(0, 1, 0) |
    //     And upper is added to w
    //     And r â† ray(point(0, 0, 0), vector(0, 1, 0))
    //   Then color_at(w, r) should terminate successfully

    // Scenario: The reflected color at the maximum recursive depth
    //   Given w â† default_world()
    //     And shape â† plane() with:
    //       | material.reflective | 0.5                   |
    //       | transform           | translation(0, -1, 0) |
    //     And shape is added to w
    //     And r â† ray(point(0, 0, -3), vector(0, -âˆš2/2, âˆš2/2))
    //     And i â† intersection(âˆš2, shape)
    //   When comps â† prepare_computations(i, r)
    //     And color â† reflected_color(w, comps, 0)
    //   Then color = color(0, 0, 0)

    // Scenario: The refracted color with an opaque surface
    //   Given w â† default_world()
    //     And shape â† the first object in w
    //     And r â† ray(point(0, 0, -5), vector(0, 0, 1))
    //     And xs â† intersections(4:shape, 6:shape)
    //   When comps â† prepare_computations(xs[0], r, xs)
    //     And c â† refracted_color(w, comps, 5)
    //   Then c = color(0, 0, 0)

    // Scenario: The refracted color at the maximum recursive depth
    //   Given w â† default_world()
    //     And shape â† the first object in w
    //     And shape has:
    //       | material.transparency     | 1.0 |
    //       | material.refractive_index | 1.5 |
    //     And r â† ray(point(0, 0, -5), vector(0, 0, 1))
    //     And xs â† intersections(4:shape, 6:shape)
    //   When comps â† prepare_computations(xs[0], r, xs)
    //     And c â† refracted_color(w, comps, 0)
    //   Then c = color(0, 0, 0)

    // Scenario: The refracted color under total internal reflection
    //   Given w â† default_world()
    //     And shape â† the first object in w
    //     And shape has:
    //       | material.transparency     | 1.0 |
    //       | material.refractive_index | 1.5 |
    //     And r â† ray(point(0, 0, âˆš2/2), vector(0, 1, 0))
    //     And xs â† intersections(-âˆš2/2:shape, âˆš2/2:shape)
    //   # NOTE: this time you're inside the sphere, so you need
    //   # to look at the second intersection, xs[1], not xs[0]
    //   When comps â† prepare_computations(xs[1], r, xs)
    //     And c â† refracted_color(w, comps, 5)
    //   Then c = color(0, 0, 0)

    // Scenario: The refracted color with a refracted ray
    //   Given w â† default_world()
    //     And A â† the first object in w
    //     And A has:
    //       | material.ambient | 1.0            |
    //       | material.pattern | test_pattern() |
    //     And B â† the second object in w
    //     And B has:
    //       | material.transparency     | 1.0 |
    //       | material.refractive_index | 1.5 |
    //     And r â† ray(point(0, 0, 0.1), vector(0, 1, 0))
    //     And xs â† intersections(-0.9899:A, -0.4899:B, 0.4899:B, 0.9899:A)
    //   When comps â† prepare_computations(xs[2], r, xs)
    //     And c â† refracted_color(w, comps, 5)
    //   Then c = color(0, 0.99888, 0.04725)

    // Scenario: shade_hit() with a transparent material
    //   Given w â† default_world()
    //     And floor â† plane() with:
    //       | transform                 | translation(0, -1, 0) |
    //       | material.transparency     | 0.5                   |
    //       | material.refractive_index | 1.5                   |
    //     And floor is added to w
    //     And ball â† sphere() with:
    //       | material.color     | (1, 0, 0)                  |
    //       | material.ambient   | 0.5                        |
    //       | transform          | translation(0, -3.5, -0.5) |
    //     And ball is added to w
    //     And r â† ray(point(0, 0, -3), vector(0, -âˆš2/2, âˆš2/2))
    //     And xs â† intersections(âˆš2:floor)
    //   When comps â† prepare_computations(xs[0], r, xs)
    //     And color â† shade_hit(w, comps, 5)
    //   Then color = color(0.93642, 0.68642, 0.68642)

    // Scenario: shade_hit() with a reflective, transparent material
    //   Given w â† default_world()
    //     And r â† ray(point(0, 0, -3), vector(0, -âˆš2/2, âˆš2/2))
    //     And floor â† plane() with:
    //       | transform                 | translation(0, -1, 0) |
    //       | material.reflective       | 0.5                   |
    //       | material.transparency     | 0.5                   |
    //       | material.refractive_index | 1.5                   |
    //     And floor is added to w
    //     And ball â† sphere() with:
    //       | material.color     | (1, 0, 0)                  |
    //       | material.ambient   | 0.5                        |
    //       | transform          | translation(0, -3.5, -0.5) |
    //     And ball is added to w
    //     And xs â† intersections(âˆš2:floor)
    //   When comps â† prepare_computations(xs[0], r, xs)
    //     And color â† shade_hit(w, comps, 5)
    //   Then color = color(0.93391, 0.69643, 0.69243)
}
