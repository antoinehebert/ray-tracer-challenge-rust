mod canvas;
mod color;
mod intersection;
mod light;
mod material;
mod matrix;
mod ray;
mod sphere;
mod test_utils;
mod transformations;
mod tuple;
mod utils;

use crate::canvas::*;
use crate::color::*;
use crate::intersection::*;
use crate::ray::*;
use crate::sphere::*;
use crate::tuple::*;

fn main() {
    chapter_5_ray_sphere_intersections();
}

//
// Putting it together. Exercices at the end of chapters.
//
fn chapter_5_ray_sphere_intersections() {
    let ray_origin = Tuple::point(0., 0., -5.);
    let wall_z = 10.;
    let wall_size = 7.;
    let canvas_pixels = 100;
    let pixel_size = wall_size / canvas_pixels as f32;
    let half = wall_size / 2.;

    let mut canvas = Canvas::new(canvas_pixels, canvas_pixels);
    let color = RED;
    let shape = Sphere::new();

    for y in 0..canvas_pixels {
        let world_y = half - pixel_size * y as f32;
        for x in 0..canvas_pixels {
            let world_x = half - pixel_size * x as f32;

            let wall_position = Tuple::point(world_x, world_y, wall_z);

            let r = Ray::new(ray_origin, (wall_position - ray_origin).normalize());
            let xs = shape.intersect(&r);

            if let Some(_) = hit(&xs) {
                canvas.set_pixel(x, y, color);
            }
        }
    }

    println!("{}", canvas.to_ppm());
}
