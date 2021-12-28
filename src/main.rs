#![allow(dead_code)] // Allow dead code while experimenting...

// TODO: use clippy?

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

use std::env;

use std::fs::File;

use camera::Camera;
use color::*;
use light::Light;
use pattern::Pattern;
use shape::Shape;
use transformations::*;
use tuple::Tuple;
use world::World;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        println!("Expected a filename argument");
        return;
    }

    putting_it_together(&args[1]);
}

const WIDTH: usize = 400;
const HEIGHT: usize = WIDTH / 2;

//
// Putting it together. Exercises at the end of chapters.
//
fn putting_it_together(filename: &str) {
    // From: https://forum.raytracerchallenge.com/thread/6/tables-scene-description
    let mut camera = Camera::new(WIDTH, HEIGHT, 0.785);
    camera.transform = view_transform(
        Tuple::point(8.0, 6.0, -8.0),
        Tuple::point(0.0, 3.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    );

    let light = Light::new(Tuple::point(0.0, 6.9, -5.0), Color::new(1.0, 1.0, 0.9));

    let mut world = World::new(light);

    let floor_ceiling = Shape::cube();
    {
        let mut floor_ceiling = floor_ceiling.borrow_mut();
        floor_ceiling.transform = scaling(20.0, 7.0, 20.0) * translation(0.0, 1.0, 0.1);
        let mut pattern = Pattern::checkers(BLACK, Color::new(0.25, 0.25, 0.25));
        pattern.transform = scaling(0.07, 0.07, 0.07);
        floor_ceiling.material.pattern = Some(pattern);
        floor_ceiling.material.ambient = 0.25;
        floor_ceiling.material.diffuse = 0.7;
        floor_ceiling.material.specular = 0.9;
        floor_ceiling.material.shininess = 300.0;
        floor_ceiling.material.reflective = 0.1;
    }
    world.objects.push(floor_ceiling);

    let walls = Shape::cube();
    {
        let mut walls = walls.borrow_mut();
        walls.transform = scaling(10.0, 10.0, 10.0);
        let mut pattern = Pattern::checkers(
            Color::new(0.4863, 0.3765, 0.2941),
            Color::new(0.3725, 0.2902, 0.2275),
        );
        pattern.transform = scaling(0.05, 20.0, 0.05);
        walls.material.pattern = Some(pattern);
        walls.material.ambient = 0.1;
        walls.material.diffuse = 0.7;
        walls.material.specular = 0.9;
        walls.material.shininess = 300.0;
        walls.material.reflective = 0.1;
    }
    world.objects.push(walls);

    let table_top = Shape::cube();
    {
        let mut table_top = table_top.borrow_mut();
        table_top.transform = translation(0.0, 3.1, 0.0) * scaling(3.0, 0.1, 2.0);
        let mut pattern = Pattern::stripe(
            Color::new(0.5529, 0.4235, 0.3255),
            Color::new(0.6588, 0.5098, 0.4000),
        );
        pattern.transform = scaling(0.05, 0.05, 0.05) * rotation_y(0.1);
        table_top.material.pattern = Some(pattern);
        table_top.material.ambient = 0.1;
        table_top.material.diffuse = 0.7;
        table_top.material.specular = 0.9;
        table_top.material.shininess = 300.0;
        table_top.material.reflective = 0.2;
    }
    world.objects.push(table_top);

    let leg_1 = Shape::cube();
    {
        let mut leg_1 = leg_1.borrow_mut();
        leg_1.transform = translation(2.7, 1.5, -1.7) * scaling(0.1, 1.5, 0.1);
        leg_1.material.color = Color::new(0.5529, 0.4235, 0.3255);
        leg_1.material.ambient = 0.2;
        leg_1.material.diffuse = 0.7;
    }

    let leg_2 = Shape::cube();
    {
        let mut leg_2 = leg_2.borrow_mut();
        leg_2.transform = translation(2.7, 1.5, 1.7) * scaling(0.1, 1.5, 0.1);
        leg_2.material.color = Color::new(0.5529, 0.4235, 0.3255);
        leg_2.material.ambient = 0.2;
        leg_2.material.diffuse = 0.7;
    }

    let leg_3 = Shape::cube();
    {
        let mut leg_3 = leg_3.borrow_mut();
        leg_3.transform = translation(-2.7, 1.5, -1.7) * scaling(0.1, 1.5, 0.1);
        leg_3.material.color = Color::new(0.5529, 0.4235, 0.3255);
        leg_3.material.ambient = 0.2;
        leg_3.material.diffuse = 0.7;
    }

    let leg_4 = Shape::cube();
    {
        let mut leg_4 = leg_4.borrow_mut();
        leg_4.transform = translation(-2.7, 1.5, 1.7) * scaling(0.1, 1.5, 0.1);
        leg_4.material.color = Color::new(0.5529, 0.4235, 0.3255);
        leg_4.material.ambient = 0.2;
        leg_4.material.diffuse = 0.7;
    }

    world.objects.push(leg_1);
    world.objects.push(leg_2);
    world.objects.push(leg_3);
    world.objects.push(leg_4);

    let glass_cube = Shape::cube();
    {
        let mut glass_cube = glass_cube.borrow_mut();
        glass_cube.transform =
            translation(0.0, 3.45001, 0.0) * rotation_y(0.2) * scaling(0.25, 0.25, 0.25);
        glass_cube.material.color = Color::new(1.0, 1.0, 0.8);
        glass_cube.material.ambient = 0.0;
        glass_cube.material.diffuse = 0.3;
        glass_cube.material.specular = 0.9;
        glass_cube.material.shininess = 300.0;
        glass_cube.material.reflective = 0.1;
        glass_cube.material.transparency = 0.7;
        glass_cube.material.refractive_index = 1.5;
    }
    world.objects.push(glass_cube);

    let little_cube_1 = Shape::cube();
    {
        let mut little_cube_1 = little_cube_1.borrow_mut();
        little_cube_1.transform =
            translation(1.0, 3.35, -0.9) * rotation_y(-0.4) * scaling(0.15, 0.15, 0.15);
        little_cube_1.material.color = Color::new(1.0, 0.5, 0.5);
        little_cube_1.material.reflective = 0.6;
        little_cube_1.material.diffuse = 0.4;
    }
    world.objects.push(little_cube_1);

    let little_cube_2 = Shape::cube();
    {
        let mut little_cube_2 = little_cube_2.borrow_mut();
        little_cube_2.transform =
            translation(-1.5, 3.27, 0.3) * rotation_y(0.4) * scaling(0.15, 0.7, 0.15);
        little_cube_2.material.color = Color::new(1.0, 1.0, 0.5);
    }
    world.objects.push(little_cube_2);

    let little_cube_3 = Shape::cube();
    {
        let mut little_cube_3 = little_cube_3.borrow_mut();
        little_cube_3.transform =
            translation(0.0, 3.25, 1.0) * rotation_y(0.4) * scaling(0.2, 0.05, 0.05);
        little_cube_3.material.color = Color::new(0.5, 1.0, 0.5);
    }
    world.objects.push(little_cube_3);

    let little_cube_4 = Shape::cube();
    {
        let mut little_cube_4 = little_cube_4.borrow_mut();
        little_cube_4.transform =
            translation(-0.6, 3.4, -1.0) * rotation_y(0.8) * scaling(0.05, 0.2, 0.05);
        little_cube_4.material.color = Color::new(0.5, 0.5, 1.0);
    }
    world.objects.push(little_cube_4);

    let little_cube_5 = Shape::cube();
    {
        let mut little_cube_5 = little_cube_5.borrow_mut();
        little_cube_5.transform =
            translation(2.0, 3.4, 1.0) * rotation_y(0.8) * scaling(0.05, 0.2, 0.05);
        little_cube_5.material.color = Color::new(0.5, 1.0, 1.0);
    }
    world.objects.push(little_cube_5);

    let frame_1 = Shape::cube();
    {
        let mut frame_1 = frame_1.borrow_mut();
        frame_1.transform = translation(-10.0, 4.0, 1.0) * scaling(0.05, 1.0, 1.0);
        frame_1.material.color = Color::new(0.7098, 0.2471, 0.2196);
        frame_1.material.diffuse = 0.6;
    }
    world.objects.push(frame_1);

    let frame_2 = Shape::cube();
    {
        let mut frame_2 = frame_2.borrow_mut();
        frame_2.transform = translation(-10.0, 3.4, 2.7) * scaling(0.05, 0.4, 0.4);
        frame_2.material.color = Color::new(0.2667, 0.2706, 0.6902);
        frame_2.material.diffuse = 0.6;
    }
    world.objects.push(frame_2);

    let frame_3 = Shape::cube();
    {
        let mut frame_3 = frame_3.borrow_mut();
        frame_3.transform = translation(-10.0, 4.6, 2.7) * scaling(0.05, 0.4, 0.4);
        frame_3.material.color = Color::new(0.3098, 0.5961, 0.3098);
        frame_3.material.diffuse = 0.6;
    }
    world.objects.push(frame_3);

    let mirror_frame = Shape::cube();
    {
        let mut mirror_frame = mirror_frame.borrow_mut();
        mirror_frame.transform = translation(-2.0, 3.5, 9.95) * scaling(5.0, 1.5, 0.05);
        mirror_frame.material.color = Color::new(0.3882, 0.2627, 0.1882);
        mirror_frame.material.diffuse = 0.7;
    }
    world.objects.push(mirror_frame);

    let mirror = Shape::cube();
    {
        let mut mirror = mirror.borrow_mut();
        mirror.transform = translation(-2.0, 3.5, 9.95) * scaling(4.8, 1.4, 0.06);
        mirror.material.color = BLACK;
        mirror.material.diffuse = 0.0;
        mirror.material.ambient = 0.0;
        mirror.material.specular = 0.0;
        mirror.material.shininess = 300.0;
        mirror.material.reflective = 1.0;
    }
    world.objects.push(mirror);

    let canvas = camera.render(&world);

    let file = File::create(filename);

    match file {
        Ok(mut file) => canvas.to_ppm(&mut file),
        Err(msg) => println!("Can't open {}: {}", filename, msg),
    }
}
