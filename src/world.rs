use crate::{
    color::{Color, BLACK, WHITE},
    intersection::*,
    light::Light,
    ray::Ray,
    shape::Shape,
    transformations::scaling,
    tuple::Tuple,
};

const RECURSION_LIMIT: usize = 5;

pub struct World {
    pub objects: Vec<Shape>,
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

        let mut s1 = Shape::sphere();
        s1.material.color = Color::new(0.8, 1.0, 0.6);
        s1.material.diffuse = 0.7;
        s1.material.specular = 0.2;

        let mut s2 = Shape::sphere();
        s2.transform = scaling(0.5, 0.5, 0.5);

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

    fn shade_hit(&self, comps: &Computations, remaining: usize) -> Color {
        let surface = comps.object.material.lighting(
            &self.light,
            &comps.object,
            &comps.point,
            &comps.eyev,
            &comps.normalv,
            self.is_shadowed(&comps.over_point),
        );

        let reflected = self.reflected_color(&comps, remaining - 1);
        let refracted = self.refracted_color(&comps, remaining - 1);

        surface + reflected + refracted
    }

    pub fn color_at(&self, ray: &Ray) -> Color {
        self.internal_color_at(&ray, RECURSION_LIMIT)
    }

    fn internal_color_at(&self, ray: &Ray, remaining: usize) -> Color {
        if remaining < 1 {
            return BLACK;
        }

        let xs = self.intersect(&ray);
        let hit = hit(&xs);

        match hit {
            None => BLACK,
            Some(intersection) => {
                self.shade_hit(&intersection.prepare_computations(&ray, &xs), remaining - 1)
            }
        }
    }

    fn is_shadowed(&self, point: &Tuple) -> bool {
        let point = *point;
        let vector = self.light.position - point;
        let distance = vector.magnitude();
        let direction = vector.normalize();

        let ray = Ray::new(point, direction);
        let intersections = self.intersect(&ray);

        let hit = hit(&intersections);
        match hit {
            Some(hit) => hit.t < distance,
            None => false,
        }
    }

    fn reflected_color(&self, comps: &Computations, remaining: usize) -> Color {
        if remaining < 1 {
            return BLACK;
        }
        if comps.object.material.reflective == 0.0 {
            return BLACK;
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflectv);
        let color = self.internal_color_at(&reflect_ray, remaining - 1);

        color * comps.object.material.reflective
    }

    fn refracted_color(&self, comps: &Computations, remaining: usize) -> Color {
        if remaining == 0 {
            return BLACK;
        }
        if comps.object.material.transparency == 0.0 {
            return BLACK;
        }

        // Snell's law: sin(tetha_2) / sin(tetha_1) = v2 / v1 = n2 / n1
        let n_ratio = comps.n1 / comps.n2; // Yes, inverted from Snell's law...
        let cos_i = comps.eyev.dot(&comps.normalv); // cos(tetha_i) is the same as the dot product of the two vectors.
        let sin2_t = n_ratio.powf(2.0) * (1.0 - cos_i.powf(2.0)); // Find sin(theta_t)^2 via trigonometric identity​
        if sin2_t > 1.0 {
            return BLACK;
        }

        // Find cos(theta_t) via trigonometric identity​
        let cos_t = (1.0 - sin2_t).sqrt();

        // Compute the direction of the refracted ray​
        let direction = comps.normalv * (n_ratio * cos_i - cos_t) - comps.eyev * n_ratio;

        // Create the refracted ray​
        let refract_ray = Ray::new(comps.under_point, direction);

        // Find the color of the refracted ray, making sure to multiply​
        // by the transparency value to account for any opacity​
        let color = self.internal_color_at(&refract_ray, remaining - 1)
            * comps.object.material.transparency;

        color
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assert_almost_eq, color::RED, transformations::translation};
    use crate::{pattern::Pattern, utils::*};

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

        let mut s1 = Shape::sphere();
        s1.material.color = Color::new(0.8, 1.0, 0.6);
        s1.material.diffuse = 0.7;
        s1.material.specular = 0.2;

        let mut s2 = Shape::sphere();
        s2.transform = scaling(0.5, 0.5, 0.5);

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
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        let c = w.shade_hit(&comps, RECURSION_LIMIT);

        assert_eq!(c, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shading_an_intersection_from_the_inside() {
        let mut w = World::default_world();
        w.light = Light::new(Tuple::point(0.0, 0.25, 0.0), WHITE);
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = &w.objects[1];
        let i = Intersection::new(0.5, shape);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        let c = w.shade_hit(&comps, RECURSION_LIMIT);

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

    #[test]
    fn there_is_no_shadow_when_nothing_is_collinear_with_point_and_light() {
        let w = World::default_world();
        let p = Tuple::point(0.0, 10.0, 0.0);

        assert!(!w.is_shadowed(&p));
    }
    #[test]
    fn the_shadow_when_an_object_is_between_the_point_and_the_light() {
        let w = World::default_world();
        let p = Tuple::point(10.0, -10.0, 10.0);

        assert!(w.is_shadowed(&p));
    }

    #[test]
    fn there_is_no_shadow_when_an_object_is_behind_the_light() {
        let w = World::default_world();
        let p = Tuple::point(-20.0, 20.0, -20.0);

        assert!(!w.is_shadowed(&p));
    }

    #[test]
    fn there_is_no_shadow_when_an_object_is_behind_the_point() {
        let w = World::default_world();
        let p = Tuple::point(-2.0, 2.0, -2.0);

        assert!(!w.is_shadowed(&p));
    }

    #[test]
    fn shade_hit_is_given_an_intersection_in_shadow() {
        let light = Light::new(Tuple::point(0.0, 0.0, -10.0), WHITE);
        let mut w = World::new(light);
        let s1 = Shape::sphere();
        w.objects.push(s1);
        let mut s2 = Shape::sphere();
        s2.transform = translation(0.0, 0.0, 10.0);
        w.objects.push(s2);

        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let i = Intersection::new(4.0, &w.objects[1]);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        let c = w.shade_hit(&comps, RECURSION_LIMIT);

        assert_eq!(c, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn the_reflected_color_for_a_nonreflective_material() {
        let mut w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        w.objects[1].material.ambient = 1.0; // make sure we have something to reflect
        let i = Intersection::new(1.0, &w.objects[1]);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        let color = w.reflected_color(&comps, RECURSION_LIMIT);
        assert_eq!(color, BLACK);
    }

    #[test]
    fn the_reflected_color_for_a_reflective_material() {
        let mut w = World::default_world();
        let mut shape = Shape::plane();
        shape.material.reflective = 0.5;
        shape.transform = translation(0.0, -1.0, 0.0);
        w.objects.push(shape.clone());
        let r = Ray::new(
            Tuple::point(0.0, 0.0, -3.0),
            Tuple::vector(0.0, -(2 as f64).sqrt() / 2.0, (2 as f64).sqrt() / 2.0),
        );
        let i = Intersection::new((2 as f64).sqrt(), &shape);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        let color = w.reflected_color(&comps, RECURSION_LIMIT);
        assert_eq!(color, Color::new(0.19033, 0.23791, 0.14274));
    }

    #[test]
    fn shade_hit_with_a_reflective_material() {
        let mut w = World::default_world();
        let mut shape = Shape::plane();
        shape.material.reflective = 0.5;
        shape.transform = translation(0.0, -1.0, 0.0);
        w.objects.push(shape.clone());
        let r = Ray::new(
            Tuple::point(0.0, 0.0, -3.0),
            Tuple::vector(0.0, -(2 as f64).sqrt() / 2.0, (2 as f64).sqrt() / 2.0),
        );
        let i = Intersection::new((2 as f64).sqrt(), &shape);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        let color = w.shade_hit(&comps, RECURSION_LIMIT);
        assert_eq!(color, Color::new(0.87675, 0.92434, 0.82918));
    }

    #[test]
    fn color_at_with_mutually_reflective_surfaces() {
        let mut w = World::new(Light::new(Tuple::point(0.0, 0.0, 0.0), WHITE));

        let mut lower = Shape::plane();
        lower.material.reflective = 1.0;
        lower.transform = translation(0.0, -1.0, 0.0);
        w.objects.push(lower);

        let mut upper = Shape::plane();
        upper.material.reflective = 1.0;
        upper.transform = translation(0.0, 1.0, 0.0);
        w.objects.push(upper);

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let _color = w.color_at(&r);
        assert!(true, "should terminate successfully");
    }
    #[test]
    fn the_reflected_color_at_the_maximum_recursive_depth() {
        let mut w = World::default_world();
        let mut shape = Shape::plane();
        shape.material.reflective = 0.5;
        shape.transform = translation(0.0, -1.0, 0.0);
        w.objects.push(shape.clone());
        let r = Ray::new(
            Tuple::point(0.0, 0.0, -3.0),
            Tuple::vector(0.0, -(2 as f64).sqrt() / 2.0, (2 as f64).sqrt() / 2.0),
        );
        let i = Intersection::new((2 as f64).sqrt(), &shape);
        let comps = i.prepare_computations(&r, &vec![i.clone()]);
        let color = w.reflected_color(&comps, 0);
        assert_eq!(color, BLACK);
    }

    #[test]
    fn the_refracted_color_with_an_opaque_surface() {
        let w = World::default_world();
        let shape = &w.objects[0];
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = vec![Intersection::new(4.0, shape), Intersection::new(6.0, shape)];
        let comps = xs[0].prepare_computations(&r, &xs);

        let c = w.refracted_color(&comps, RECURSION_LIMIT);

        assert_eq!(c, BLACK);
    }

    #[test]
    fn the_refracted_color_at_the_maximum_recursive_depth() {
        let mut w = World::default_world();
        {
            let shape = &mut w.objects[0];
            shape.material.transparency = 1.0;
            shape.material.refractive_index = 1.5;
        }
        let shape = &w.objects[0];
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = vec![Intersection::new(4.0, shape), Intersection::new(6.0, shape)];
        let comps = xs[0].prepare_computations(&r, &xs);

        let c = w.refracted_color(&comps, 0);

        assert_eq!(c, BLACK);
    }

    #[test]
    fn the_refracted_color_under_total_internal_reflection() {
        let mut w = World::default_world();
        {
            let shape = &mut w.objects[0];
            shape.material.transparency = 1.0;
            shape.material.refractive_index = 1.5;
        }
        let shape = &w.objects[0];
        let r = Ray::new(
            Tuple::point(0.0, 0.0, (2.0 as f64).sqrt() / 2.0),
            Tuple::vector(0.0, 1.0, 0.0),
        );
        let xs = vec![
            Intersection::new(-(2.0 as f64).sqrt() / 2.0, shape),
            Intersection::new((2.0 as f64).sqrt() / 2.0, shape),
        ];
        // NOTE: this time you're inside the sphere, so you need to look at the second intersection, xs[1], not xs[0]
        let comps = xs[1].prepare_computations(&r, &xs);

        let c = w.refracted_color(&comps, RECURSION_LIMIT);

        assert_eq!(c, BLACK);
    }

    #[test]
    fn the_refracted_color_with_a_refracted_ray() {
        let mut w = World::default_world();

        {
            let a = &mut w.objects[0];
            a.material.ambient = 1.0; // Fully ambient, so that it shows up regardless of lighting.
            a.material.pattern = Some(Pattern::test_pattern());
        }
        {
            let b = &mut w.objects[1];
            b.material.transparency = 1.0;
            b.material.refractive_index = 1.5;
        }

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.1), Tuple::vector(0.0, 1.0, 0.0));

        let a = &w.objects[0];
        let b = &w.objects[1];
        let xs = vec![
            Intersection::new(-0.9899, &a),
            Intersection::new(-0.4899, &b),
            Intersection::new(0.4899, &b),
            Intersection::new(0.9899, &a),
        ];
        let comps = xs[2].prepare_computations(&r, &xs);
        let c = w.refracted_color(&comps, RECURSION_LIMIT);

        assert_eq!(c, Color::new(0.0, 0.99888, 0.04721));
    }

    #[test]
    fn shade_hit_with_a_transparent_material() {
        let mut w = World::default_world();

        let mut floor = Shape::plane();
        floor.transform = translation(0.0, -1.0, 0.0);
        floor.material.transparency = 0.5;
        floor.material.refractive_index = 1.5;
        w.objects.push(floor);

        let mut ball = Shape::sphere();
        ball.material.color = RED;
        ball.material.ambient = 0.5;
        ball.transform = translation(0.0, -3.5, -0.5);
        w.objects.push(ball);

        let r = Ray::new(
            Tuple::point(0.0, 0.0, -3.0),
            Tuple::vector(0.0, -2.0_f64.sqrt() / 2.0, 2.0_f64.sqrt() / 2.0),
        );
        let xs = vec![Intersection::new(
            2.0_f64.sqrt(),
            &w.objects[w.objects.len() - 2], // floor...
        )];
        let comps = xs[0].prepare_computations(&r, &xs);
        let color = w.shade_hit(&comps, RECURSION_LIMIT);

        assert_eq!(color, Color::new(0.93642, 0.68642, 0.68642));
    }

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
