mod camera;
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
mod world;

use std::f32::consts::PI;

use camera::Camera;
use canvas::Canvas;
use color::{Color, WHITE};
use intersection::hit;
use light::Light;
use material::Material;
use ray::Ray;
use sphere::Sphere;
use transformations::*;
use tuple::Tuple;
use world::World;

fn main() {
    chapter_7_making_a_scene();
}

//
// Putting it together. Exercises at the end of chapters.
//
fn chapter_7_making_a_scene() {
    // Floors and walls are extremely flatten spheres...
    let mut floor = Sphere::new();
    floor.transform(scaling(10.0, 0.01, 10.));
    floor.material.color = Color::new(1.0, 0.9, 0.9);
    floor.material.specular = 0.0;

    let mut left_wall = Sphere::new();
    left_wall.transform(
        translation(0.0, 0.0, 5.0)
            * rotation_y(-PI / 4.0)
            * rotation_x(PI / 2.0)
            * scaling(10.0, 0.01, 10.),
    );
    left_wall.material = floor.material;

    let mut right_wall = Sphere::new();
    right_wall.transform(
        translation(0.0, 0.0, 5.0)
            * rotation_y(PI / 4.0)
            * rotation_x(PI / 2.0)
            * scaling(10.0, 0.01, 10.),
    );
    right_wall.material = floor.material;

    let mut middle = Sphere::new();
    middle.transform(translation(-0.5, 1.0, 0.5));
    middle.material.color = Color::new(0.1, 1.0, 0.5);
    middle.material.diffuse = 0.7;
    middle.material.specular = 0.3;

    let mut right = Sphere::new();
    right.transform(translation(1.5, 0.5, -0.5) * scaling(0.5, 0.5, 0.5));
    right.material.color = Color::new(0.5, 1.0, 0.1);
    right.material.diffuse = 0.7;
    right.material.specular = 0.3;

    let mut left = Sphere::new();
    left.transform(translation(-1.5, 0.33, -0.75) * scaling(0.33, 0.33, 0.33));
    left.material.color = Color::new(1.0, 0.8, 0.1);
    left.material.diffuse = 0.7;
    left.material.specular = 0.3;

    let light = Light::new(Tuple::point(-10.0, 10.0, -10.0), WHITE);
    let mut world = World::new(light);
    world.objects.push(floor);
    world.objects.push(left_wall);
    world.objects.push(right_wall);
    world.objects.push(middle);
    world.objects.push(right);
    world.objects.push(left);

    let mut camera = Camera::new(100, 50, PI / 3.0);
    camera.transform = view_transform(
        Tuple::point(0.0, 1.5, -5.0),
        Tuple::point(0.0, 1.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    );

    let canvas = camera.render(&world);
    println!("{}", canvas.to_ppm());
}
