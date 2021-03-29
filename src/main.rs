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

use canvas::Canvas;
use color::{Color, WHITE};
use intersection::hit;
use light::Light;
use material::Material;
use ray::Ray;
use sphere::Sphere;
use tuple::Tuple;

fn main() {
    chapter_6_ray_sphere_intersections();
}

//
// Putting it together. Exercises at the end of chapters.
//
fn chapter_6_ray_sphere_intersections() {
    let ray_origin = Tuple::point(0., 0., -5.);
    let wall_z = 10.;
    let wall_size = 7.;
    let canvas_pixels = 200;
    let pixel_size = wall_size / canvas_pixels as f32;
    let half = wall_size / 2.;

    let mut canvas = Canvas::new(canvas_pixels, canvas_pixels);
    let mut shape = Sphere::new();
    shape.material = Material::new();
    shape.material.color = Color::new(1.0, 0.2, 1.0);

    let light = Light::new(Tuple::point(-10.0, 10.0, -10.0), WHITE);

    for y in 0..canvas_pixels {
        let world_y = half - pixel_size * y as f32;
        for x in 0..canvas_pixels {
            let world_x = -half + pixel_size * x as f32;

            let wall_position = Tuple::point(world_x, world_y, wall_z);

            let ray = Ray::new(ray_origin, (wall_position - ray_origin).normalize());
            let xs = shape.intersect(&ray);

            if let Some(hit) = hit(&xs) {
                let point = ray.position(hit.t);
                let normal = hit.object.normal_at(point);
                let eye = -ray.direction;

                let color = hit.object.material.lighting(&light, &point, &eye, &normal);
                canvas.set_pixel(x, y, color);
            }
        }
    }

    println!("{}", canvas.to_ppm());
}
