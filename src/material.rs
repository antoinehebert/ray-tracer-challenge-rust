use crate::{color::*, light::Light, pattern::Pattern, shape::Shape, tuple::Tuple};

#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    pub color: Color,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
    pub reflective: f64,
    pub pattern: Option<Pattern>,
}

impl Material {
    pub fn new() -> Self {
        Self {
            color: WHITE,
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
            reflective: 0.0,
            pattern: None,
        }
    }

    // TOOD: Feels like we should pass in Computations entirely here or maybe move this method on Object?
    pub fn lighting(
        &self,
        light: &Light,
        object: &Shape,
        point: &Tuple,
        eyev: &Tuple,
        normalv: &Tuple,
        in_shadow: bool,
    ) -> Color {
        // Combine the surface color with the light's color/intensity.
        let color = if let Some(pattern) = &self.pattern {
            pattern.color_at_shape(object, point)
        } else {
            self.color
        };

        let effective_color = color * light.intensity;

        // Direction to the light source.
        let lightv = (light.position - *point).normalize();

        let ambient = effective_color * self.ambient;

        let mut diffuse = BLACK;
        let mut specular = BLACK;
        if !in_shadow {
            let light_dot_normal = lightv.dot(&normalv);
            // A negative number means the light is on the other side of the surface.
            if light_dot_normal >= 0. {
                diffuse = effective_color * self.diffuse * light_dot_normal;

                let reflectv = (-lightv).reflect(&normalv);
                // A negative number means the light reflects away from the eye.
                let reflect_dot_eye = reflectv.dot(&eyev);

                if reflect_dot_eye > 0. {
                    let factor = reflect_dot_eye.powf(self.shininess);
                    specular = light.intensity * self.specular * factor;
                }
            }
        }

        ambient + diffuse + specular
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_material() {
        let m = Material::new();

        assert_eq!(m.color, WHITE);
        assert_eq!(m.ambient, 0.1);
        assert_eq!(m.diffuse, 0.9);
        assert_eq!(m.specular, 0.9);
        assert_eq!(m.shininess, 200.0);
    }

    #[test]
    fn reflectivity_for_the_default_material() {
        let m = Material::new();
        assert_eq!(m.reflective, 0.0);
    }

    // Scenario: Transparency and Refractive Index for the default material
    //   Given m â† material()
    //   Then m.transparency = 0.0
    //     And m.refractive_index = 1.0

    #[test]
    fn lighting_with_the_eye_between_the_light_and_the_surface() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., 0., -1.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 0., -10.), Color::new(1., 1., 1.));

        let result = m.lighting(&light, &Shape::sphere(), &position, &eyev, &normalv, false);
        assert_eq!(result, Color::new(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_with_the_eye_between_light_and_surface_eye_offset_45_degrees() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., (2. as f64).sqrt() / 2., -(2. as f64).sqrt() / 2.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 0., -10.), Color::new(1., 1., 1.));
        let result = m.lighting(&light, &Shape::sphere(), &position, &eyev, &normalv, false);
        assert_eq!(result, Color::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn lighting_with_eye_opposite_surface_light_offset_45_degrees() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., 0., -1.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 10., -10.), Color::new(1., 1., 1.));
        let result = m.lighting(&light, &Shape::sphere(), &position, &eyev, &normalv, false);
        assert_eq!(result, Color::new(0.7364, 0.7364, 0.7364));
    }

    #[test]
    fn lighting_with_eye_in_the_path_of_the_reflection_vector() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., -(2. as f64).sqrt() / 2., -(2. as f64).sqrt() / 2.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 10., -10.), Color::new(1., 1., 1.));
        let result = m.lighting(&light, &Shape::sphere(), &position, &eyev, &normalv, false);
        assert_eq!(result, Color::new(1.6364, 1.6364, 1.6364));
    }

    #[test]
    fn lighting_with_the_light_behind_the_surface() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., 0., -1.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 0., 10.), Color::new(1., 1., 1.));
        let result = m.lighting(&light, &Shape::sphere(), &position, &eyev, &normalv, false);
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_the_surface_in_shadow() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., 0., -1.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 0., -10.), Color::new(1., 1., 1.));
        let in_shadow = true;
        let result = m.lighting(
            &light,
            &Shape::sphere(),
            &position,
            &eyev,
            &normalv,
            in_shadow,
        );
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_a_pattern_applied() {
        let mut m = Material::new();
        m.pattern = Some(Pattern::stripe(WHITE, BLACK));
        m.ambient = 1.0;
        m.diffuse = 0.0;
        m.specular = 0.0;
        let eyev = Tuple::vector(0., 0., -1.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 0., -10.), WHITE);

        let c1 = m.lighting(
            &light,
            &Shape::sphere(),
            &Tuple::point(0.9, 0., 0.),
            &eyev,
            &normalv,
            false,
        );
        let c2 = m.lighting(
            &light,
            &Shape::sphere(),
            &Tuple::point(1.1, 0., 0.),
            &eyev,
            &normalv,
            false,
        );
        assert_eq!(c1, Color::new(1.0, 1.0, 1.0));
        assert_eq!(c2, Color::new(0.0, 0.0, 0.0));
    }
}
