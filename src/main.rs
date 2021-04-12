mod camera;
mod canvas;
mod color;
mod intersection;
mod light;
mod material;
mod matrix;
mod pattern;
mod ray;
mod shape;
mod test_utils;
mod transformations;
mod tuple;
mod utils;
mod world;

use std::f64::consts::PI;

use camera::Camera;
use color::*;
use light::Light;
use pattern::Pattern;
use shape::Shape;
use transformations::*;
use tuple::Tuple;
use world::World;

fn main() {
    chapter_10_patterns();
}

const WIDTH: usize = 100;
const HEIGHT: usize = WIDTH / 2;

//
// Putting it together. Exercises at the end of chapters.
//
fn chapter_10_patterns() {
    let mut floor = Shape::plane();
    let mut floor_pattern = Pattern::stripe(BLUE, WHITE);
    floor_pattern.transform = rotation_y(PI / 4.0);
    floor.material.pattern = Some(floor_pattern);
    floor.transform = translation(0.0, 0.01, 0.0);

    let mut middle = Shape::sphere();
    middle.transform = translation(-0.5, 1.0, 0.5);
    middle.material.diffuse = 0.7;
    middle.material.specular = 0.3;
    let mut middle_pattern =
        Pattern::gradient(Color::new(1.0, 0.5, 0.0), Color::new(0.0, 1.0, 1.0));
    middle_pattern.transform = translation(-1.5, 0.0, 0.0) * scaling(2.5, 2.5, 2.5);
    middle.material.pattern = Some(middle_pattern);

    let mut right = Shape::sphere();
    right.transform = translation(1.5, 0.5, -0.5) * scaling(0.5, 0.5, 0.5);
    right.material.diffuse = 0.7;
    right.material.specular = 0.3;
    let mut right_pattern = Pattern::checkers(Color::new(0.5, 1.0, 0.1), RED);
    right_pattern.transform = scaling(0.25, 0.25, 0.25);
    right.material.pattern = Some(right_pattern);

    let mut left = Shape::sphere();
    left.transform = translation(-1.5, 0.33, -0.75) * scaling(0.33, 0.33, 0.33);
    left.material.color = Color::new(1.0, 0.8, 0.1);
    left.material.diffuse = 0.7;
    left.material.specular = 0.3;
    let mut left_pattern = Pattern::ring(GREEN + BLUE, BLUE + RED);
    left_pattern.transform = rotation_x(PI / 2.0) * scaling(0.15, 0.15, 0.15);
    left.material.pattern = Some(left_pattern);

    let light = Light::new(Tuple::point(-10.0, 10.0, -10.0), WHITE);
    let mut world = World::new(light);
    world.objects.push(floor);
    world.objects.push(middle);
    world.objects.push(right);
    world.objects.push(left);

    let mut camera = Camera::new(WIDTH, HEIGHT, PI / 3.0);
    camera.transform = view_transform(
        Tuple::point(0.0, 1.5, -5.0),
        Tuple::point(0.0, 1.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    );

    let canvas = camera.render(&world);
    println!("{}", canvas.to_ppm());
}
