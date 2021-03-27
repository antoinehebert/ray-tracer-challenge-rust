use crate::color::*;
use crate::light::*;
use crate::tuple::*;

#[derive(Debug, Copy, Clone)]
pub struct Material {
    pub color: Color,
    pub ambient: f32,
    pub diffuse: f32,
    pub specular: f32,
    pub shininess: f32,
}

impl Material {
    pub fn new() -> Self {
        Self {
            color: WHITE,
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
        }
    }

    fn lighting(&self, light: Light, point: Tuple, eyev: Tuple, normalv: Tuple) -> Color {
        // Combine the surface color with the light's color/intensity.
        let effective_color = self.color * light.intensity;

        // Direction to the light source.
        let lightv = (light.position - point).normalize();

        let ambient = effective_color * self.ambient;

        // A negative number means the light is on the other side of the surface.
        let light_dot_normal = lightv.dot(&normalv);
        let mut diffuse = BLACK;
        let mut specular = BLACK;
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

        ambient + diffuse + specular
    }
}

impl PartialEq for Material {
    fn eq(&self, other: &Self) -> bool {
        self.color == other.color
            && self.ambient == other.ambient
            && self.diffuse == other.diffuse
            && self.specular == other.specular
            && self.shininess == other.shininess
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

    // Scenario: Reflectivity for the default material
    //   Given m â† material()
    //   Then m.reflective = 0.0

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

        let result = m.lighting(light, position, eyev, normalv);
        assert_eq!(result, Color::new(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_with_the_eye_between_light_and_surface_eye_offset_45_degrees() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., (2. as f32).sqrt() / 2., -(2. as f32).sqrt() / 2.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 0., -10.), Color::new(1., 1., 1.));
        let result = m.lighting(light, position, eyev, normalv);
        assert_eq!(result, Color::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn lighting_with_eye_opposite_surface_light_offset_45_degrees() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., 0., -1.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 10., -10.), Color::new(1., 1., 1.));
        let result = m.lighting(light, position, eyev, normalv);
        assert_eq!(result, Color::new(0.7364, 0.7364, 0.7364));
    }

    #[test]
    fn lighting_with_eye_in_the_path_of_the_reflection_vector() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., -(2. as f32).sqrt() / 2., -(2. as f32).sqrt() / 2.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 10., -10.), Color::new(1., 1., 1.));
        let result = m.lighting(light, position, eyev, normalv);
        assert_eq!(result, Color::new(1.63638, 1.63638, 1.63638));
    }

    #[test]
    fn lighting_with_the_light_behind_the_surface() {
        let m = Material::new();
        let position = Tuple::point(0., 0., 0.);
        let eyev = Tuple::vector(0., 0., -1.);
        let normalv = Tuple::vector(0., 0., -1.);
        let light = Light::new(Tuple::point(0., 0., 10.), Color::new(1., 1., 1.));
        let result = m.lighting(light, position, eyev, normalv);
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    // Scenario: Lighting with the surface in shadow
    //   Given eyev â† vector(0., 0., -1.)
    //     And normalv â† vector(0., 0., -1.)
    //     And light â† point_light(point(0., 0., -10.), color(1., 1., 1.))
    //     And in_shadow â† true
    //   When result â† lighting(m, light, position, eyev, normalv, in_shadow)
    //   Then result = color(0..1., 0.1, 0.1)

    // Scenario: Lighting with a pattern applied
    //   Given m.pattern â† stripe_pattern(color(1., 1., 1.), color(0., 0., 0.))
    //     And m.ambient â† 1.
    //     And m.diffuse â† 0.
    //     And m.specular â† 0.
    //     And eyev â† vector(0., 0., -1.)
    //     And normalv â† vector(0., 0., -1.)
    //     And light â† point_light(point(0., 0., -10.), color(1., 1., 1.))
    //   When c1 â† lighting(m, light, point(0.9, 0., 0.), eyev, normalv, false)
    //     And c2 â† lighting(m, light, point(1.1, 0., 0.), eyev, normalv, false)
    //   Then c1 = color(1., 1., 1.)
    //     And c2 = color(0., 0., 0.)
}
