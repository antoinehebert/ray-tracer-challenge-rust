mod camera;
mod canvas;
mod color;
mod intersection;
mod light;
mod material;
mod matrix;
mod ray;
mod shape;
mod test_utils;
mod transformations;
mod tuple;
mod utils;
mod world;

use std::f64::consts::PI;

use camera::Camera;
use color::{Color, WHITE};
use light::Light;
use material::Material;
use shape::Shape;
use transformations::*;
use tuple::Tuple;
use world::World;

fn main() {
    chapter_8_shadows();
}

//
// Putting it together. Exercises at the end of chapters.
//
fn chapter_8_shadows() {
    // Floors and walls are extremely flatten spheres...
    let mut floor = Shape::sphere();
    floor.set_transform(scaling(10.0, 0.01, 10.));
    let mut material = Material::new();
    material.color = Color::new(1.0, 0.9, 0.9);
    material.specular = 0.0;
    floor.set_material(material);

    let mut left_wall = Shape::sphere();
    left_wall.set_transform(
        translation(0.0, 0.0, 5.0)
            * rotation_y(-PI / 4.0)
            * rotation_x(PI / 2.0)
            * scaling(10.0, 0.01, 10.),
    );
    left_wall.set_material(floor.material().clone());

    let mut right_wall = Shape::sphere();
    right_wall.set_transform(
        translation(0.0, 0.0, 5.0)
            * rotation_y(PI / 4.0)
            * rotation_x(PI / 2.0)
            * scaling(10.0, 0.01, 10.),
    );
    right_wall.set_material(floor.material().clone());

    let mut middle = Shape::sphere();
    middle.set_transform(translation(-0.5, 1.0, 0.5));
    let mut material = Material::new();
    material.color = Color::new(0.1, 1.0, 0.5);
    material.diffuse = 0.7;
    material.specular = 0.3;
    middle.set_material(material);

    let mut right = Shape::sphere();
    right.set_transform(translation(1.5, 0.5, -0.5) * scaling(0.5, 0.5, 0.5));
    let mut material = Material::new();
    material.color = Color::new(0.5, 1.0, 0.1);
    material.diffuse = 0.7;
    material.specular = 0.3;
    right.set_material(material);

    let mut left = Shape::sphere();
    left.set_transform(translation(-1.5, 0.33, -0.75) * scaling(0.33, 0.33, 0.33));
    let mut material = Material::new();
    material.color = Color::new(1.0, 0.8, 0.1);
    material.diffuse = 0.7;
    material.specular = 0.3;
    left.set_material(material);

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
