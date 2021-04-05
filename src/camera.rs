use crate::{canvas::Canvas, matrix::Matrix, ray::Ray};
use crate::{tuple::Tuple, world::World};

pub struct Camera {
    hsize: usize,
    vsize: usize,
    field_of_view: f32,
    pub transform: Matrix<4>,
    pixel_size: f32,
    half_width: f32,
    half_height: f32,
}

impl Camera {
    pub fn new(hsize: usize, vsize: usize, field_of_view: f32) -> Self {
        let mut result = Self {
            hsize,
            vsize,
            field_of_view,
            transform: Matrix::<4>::identity(),
            pixel_size: 0.0,
            half_width: 0.0,
            half_height: 0.0,
        };

        let half_view = (field_of_view / 2.0).tan();
        let aspect = hsize as f32 / vsize as f32;

        if aspect >= 1.0 {
            result.half_width = half_view;
            result.half_height = half_view / aspect;
        } else {
            result.half_width = half_view * aspect;
            result.half_height = half_view;
        }
        result.pixel_size = (result.half_width * 2.0) / hsize as f32;

        result
    }

    fn ray_for_pixel(&self, px: usize, py: usize) -> Ray {
        // Offset from the edge of the canvas to the pixel's center.
        let xoffset = (px as f32 + 0.5) * self.pixel_size;
        let yoffset = (py as f32 + 0.5) * self.pixel_size;

        // The untransformed coordinates of the pixel in world space ​(remember that the camera looks toward -z, so +x
        // is to the *left*.)​
        let world_x = self.half_width - xoffset;
        let world_y = self.half_height - yoffset;

        // Using the camera matrix, transform the canvas point and the origin,​ and then compute the ray's direction
        // vector​ (remember that the canvas is at z=-1.)​
        let pixel = self.transform.inverse().expect("Should be invertible")
            * Tuple::point(world_x, world_y, -1.0);
        let origin =
            self.transform.inverse().expect("Should be invertible") * Tuple::point(0.0, 0.0, 0.0);
        let direction = (pixel - origin).normalize();

        Ray::new(origin, direction)
    }

    pub fn render(&self, world: &World) -> Canvas {
        let mut image = Canvas::new(self.hsize, self.vsize);

        for y in 0..self.vsize {
            for x in 0..self.hsize {
                let ray = self.ray_for_pixel(x, y);
                let color = world.color_at(&ray);
                image.set_pixel(x, y, color);
            }
        }

        image
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::color::Color;
    use crate::transformations::*;
    use crate::utils::*;
    use crate::{assert_almost_eq, world::World};
    use std::f32::consts::PI;

    #[test]
    fn constructing_a_camera() {
        let hsize = 160;
        let vsize = 120;
        let field_of_view = PI / 2.0;
        let c = Camera::new(hsize, vsize, field_of_view);

        assert_eq!(c.hsize, 160);
        assert_eq!(c.vsize, 120);
        assert_eq!(c.field_of_view, PI / 2.0);
        assert_eq!(c.transform, Matrix::<4>::identity());
    }

    #[test]
    fn the_pixel_size_for_a_horizontal_canvas() {
        let c = Camera::new(200, 125, PI / 2.0);
        assert_almost_eq!(c.pixel_size, 0.01);
    }

    #[test]
    fn the_pixel_size_for_a_vertical_canvas() {
        let c = Camera::new(125, 200, PI / 2.0);
        assert_almost_eq!(c.pixel_size, 0.01);
    }

    #[test]
    fn constructing_a_ray_through_the_center_of_the_canvas() {
        let c = Camera::new(201, 101, PI / 2.0);
        let r = c.ray_for_pixel(100, 50);
        assert_eq!(r.origin, Tuple::point(0.0, 0.0, 0.0));
        assert_eq!(r.direction, Tuple::vector(0.0, 0.0, -1.0));
    }

    #[test]
    fn constructing_a_ray_through_a_corner_of_the_canvas() {
        let c = Camera::new(201, 101, PI / 2.0);
        let r = c.ray_for_pixel(0, 0);
        assert_eq!(r.origin, Tuple::point(0.0, 0.0, 0.0));
        assert_eq!(r.direction, Tuple::vector(0.66519, 0.33259, -0.66851));
    }

    #[test]
    fn constructing_a_ray_when_the_camera_is_transformed() {
        let mut c = Camera::new(201, 101, PI / 2.0);
        c.transform = rotation_y(PI / 4.0) * translation(0.0, -2.0, 5.0);
        let r = c.ray_for_pixel(100, 50);
        assert_eq!(r.origin, Tuple::point(0.0, 2.0, -5.0));
        assert_eq!(
            r.direction,
            Tuple::vector((2.0 as f32).sqrt() / 2.0, 0.0, -(2.0 as f32).sqrt() / 2.0)
        );
    }

    #[test]
    fn rendering_a_world_with_a_camera() {
        let w = World::default_world();
        let mut c = Camera::new(11, 11, PI / 2.0);
        let from = Tuple::point(0.0, 0.0, -5.0);
        let to = Tuple::point(0.0, 0.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        c.transform = view_transform(from, to, up);
        let image = c.render(&w);

        assert_eq!(image.get_pixel(5, 5), Color::new(0.38066, 0.47583, 0.2855));
    }
}
