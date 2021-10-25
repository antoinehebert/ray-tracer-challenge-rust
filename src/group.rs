use uuid::Uuid;

use crate::{
    intersection::Intersections,
    matrix::Matrix,
    ray::Ray,
    shape::{Shape, Shapeable},
    tuple::Tuple,
};

// TODO: Figure out if there's a nicer way of dealing with parents. For now we use an id just to make things simple
// since Rust makes it hard to have circular references.
//
// different struct.
#[derive(Debug, PartialEq, Clone)]
pub struct Group {
    pub transform: Matrix<4>,
    pub shapes: Vec<Shape>,
    parent_id: Option<Uuid>,
    id: Uuid,
}

impl Group {
    pub fn new() -> Self {
        Self {
            transform: Matrix::<4>::identity(),
            shapes: Vec::new(),
            parent_id: None,
            id: Uuid::new_v4(),
        }
    }

    pub fn add_child(&mut self, mut shape: Shape) {
        shape.parent_id = Some(self.id);
        self.shapes.push(shape);
    }
}

impl Shapeable for Group {
    fn parent(&self) -> Option<&Group> {
        None
    }

    fn intersect(&self, world_ray: &Ray) -> Intersections {
        // TODO: DRY: Cut and pasted from Shape, this is what #local_intersect should do in the book.
        let local_ray = world_ray.transform(
            self.transform
                .inverse()
                .expect("shape transform should be invertible"),
        );

        let mut results: Intersections = (&self.shapes)
            .into_iter()
            .flat_map(|shape| shape.intersect(&local_ray))
            .collect();

        results.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));

        results
    }

    fn normal_at(&self, world_point: &Tuple) -> Tuple {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use crate::transformations::{scaling, translation};

    use super::*;

    //
    //Feature: Groups
    //

    #[test]
    fn creating_a_new_group() {
        let g = Group::new();

        assert_eq!(g.transform, Matrix::<4>::identity());
        assert_eq!(g.shapes.len(), 0);
    }

    #[test]
    fn adding_a_child_to_a_group() {
        let mut g = Group::new();
        let s = Shape::sphere();
        g.add_child(s);

        assert_eq!(g.shapes.len(), 1);
        assert_eq!(g.shapes[0].parent_id, Some(g.id));
    }

    #[test]
    fn intersecting_a_ray_with_an_empty_group() {
        let g = Group::new();
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 0.0));

        let xs = g.intersect(&r);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn intersecting_a_ray_with_a_nonempty_group() {
        let mut g = Group::new();

        let s1 = Shape::sphere();
        let mut s2 = Shape::sphere();
        s2.transform = translation(0., 0., -3.);
        let mut s3 = Shape::sphere();
        s3.transform = translation(5., 0., 0.);

        g.add_child(s1);
        g.add_child(s2);
        g.add_child(s3);

        let r = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));

        let xs = g.intersect(&r);
        assert_eq!(xs.len(), 4);
        assert_eq!(xs[0].object, &g.shapes[1]); // s2
        assert_eq!(xs[1].object, &g.shapes[1]); // s2
        assert_eq!(xs[2].object, &g.shapes[0]); // s1
        assert_eq!(xs[3].object, &g.shapes[0]); // s1
    }

    #[test]
    fn intersecting_a_transformed_group() {
        let mut g = Group::new();
        g.transform = scaling(2., 2., 2.);
        let mut s = Shape::sphere();
        s.transform = translation(5., 0., 0.);
        g.add_child(s);
        let r = Ray::new(Tuple::point(10., 0., -10.), Tuple::vector(0., 0., 1.));

        let xs = g.intersect(&r);
        assert_eq!(xs.len(), 2);
    }
}
