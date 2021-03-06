use crate::shape::Shape;
use crate::shape::ShapeKind;
use crate::tuple::Tuple;

pub struct Bounds {
    pub min: Tuple,
    pub max: Tuple,
}

impl Bounds {
    pub fn new(shape: &Shape) -> Self {
        let min: Tuple;
        let max: Tuple;

        match &shape.kind {
            ShapeKind::Sphere | ShapeKind::Cube => {
                min = Tuple::point(-1., -1., -1.);
                max = Tuple::point(1., 1., 1.);
            }
            ShapeKind::Plane => {
                min = Tuple::point(-1., -1., 0.);
                max = Tuple::point(1., 1., 0.);
            }
            ShapeKind::Cylinder {
                minimum,
                maximum,
                capped,
            } => {
                if *capped {
                    min = Tuple::point(-1., *minimum, -1.);
                    max = Tuple::point(1., *maximum, 1.);
                } else {
                    min = Tuple::point(-1., -f64::INFINITY, -1.);
                    max = Tuple::point(1., f64::INFINITY, 1.);
                }
            }
            ShapeKind::Cone {
                minimum,
                maximum,
                capped,
            } => {
                if *capped {
                    min = Tuple::point(-1., *minimum, -1.);
                    max = Tuple::point(1., *maximum, 1.);
                } else {
                    min = Tuple::point(-1., -f64::INFINITY, -1.);
                    max = Tuple::point(1., f64::INFINITY, 1.);
                }
            }
            ShapeKind::Group { shapes } => {
                let mut out = Self {
                    min: Tuple::point(0., 0., 0.),
                    max: Tuple::point(0., 0., 0.),
                };

                for shape in shapes {
                    // TODO: Does this make sense now that we don't use parents for groups but push transformations down
                    // to children?

                    // First, to convert a point from object space to its parent space, multiply the point by the
                    // object’s transformation matrix.
                    let parent_space_bounds = Bounds::new(&shape);
                    let transformation = shape.transform();

                    // Second, when transforming an entire bounding box, first transform all eight of the cube’s
                    // corners, and then find a single bounding box that fits them all. If you can’t quite see why you’d
                    // need to transform all eight points, imagine rotating the box 45° around any axis, and then figure
                    // out what the new axis-aligned bounding box ought to look like.
                    let p1 = transformation
                        * Tuple::point(
                            parent_space_bounds.min.x,
                            parent_space_bounds.min.y,
                            parent_space_bounds.min.z,
                        );
                    let p2 = transformation
                        * Tuple::point(
                            parent_space_bounds.min.x,
                            parent_space_bounds.min.y,
                            parent_space_bounds.max.z,
                        );
                    let p3 = transformation
                        * Tuple::point(
                            parent_space_bounds.min.x,
                            parent_space_bounds.max.y,
                            parent_space_bounds.min.z,
                        );
                    let p4 = transformation
                        * Tuple::point(
                            parent_space_bounds.min.x,
                            parent_space_bounds.max.y,
                            parent_space_bounds.max.z,
                        );
                    let p5 = transformation
                        * Tuple::point(
                            parent_space_bounds.max.x,
                            parent_space_bounds.min.y,
                            parent_space_bounds.min.z,
                        );
                    let p6 = transformation
                        * Tuple::point(
                            parent_space_bounds.max.x,
                            parent_space_bounds.min.y,
                            parent_space_bounds.max.z,
                        );
                    let p7 = transformation
                        * Tuple::point(
                            parent_space_bounds.max.x,
                            parent_space_bounds.max.y,
                            parent_space_bounds.min.z,
                        );
                    let p8 = transformation * parent_space_bounds.max;

                    out.add(&p1);
                    out.add(&p2);
                    out.add(&p3);
                    out.add(&p4);
                    out.add(&p5);
                    out.add(&p6);
                    out.add(&p7);
                    out.add(&p8);
                }

                min = out.min;
                max = out.max;
            }
            ShapeKind::Triangle { p1, p2, p3, .. } => {
                let mut tmp_bounds = Self {
                    min: Tuple::point(0., 0., 0.),
                    max: Tuple::point(0., 0., 0.),
                };
                tmp_bounds.add(p1);
                tmp_bounds.add(p2);
                tmp_bounds.add(p3);

                min = tmp_bounds.min;
                max = tmp_bounds.max;
            }
        }
        Self { min, max }
    }

    fn add(&mut self, point: &Tuple) {
        assert!(point.is_point());

        self.min.x = f64::min(self.min.x, point.x);
        self.min.y = f64::min(self.min.y, point.y);
        self.min.z = f64::min(self.min.z, point.z);
        self.max.x = f64::max(self.max.x, point.x);
        self.max.y = f64::max(self.max.y, point.y);
        self.max.z = f64::max(self.max.z, point.z);
    }
}
