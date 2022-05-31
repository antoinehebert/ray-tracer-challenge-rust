use std::collections::HashMap;

use crate::{
    shape::{ChildShape, Shape},
    tuple::Tuple,
};

struct Parser {
    vertices: Vec<Tuple>,
    ignored_lines: usize,
    default_group: ChildShape,
    named_groups: HashMap<String, ChildShape>,
}

impl Parser {
    pub fn new() -> Self {
        Self {
            ignored_lines: 0,
            vertices: Vec::new(),
            default_group: Shape::group(), // TODO: Do we need this? Should we have a "" key in named_groups?
            named_groups: HashMap::new(),
        }
    }

    fn from_obj_file(text: &str) -> Self {
        let text = text.clone();
        let mut result = Self::new();

        let current_group = None;

        for s in text.lines() {
            let mut tokens = s.split_whitespace();
            if let Some(token) = tokens.next() {
                match token {
                    "v" => {
                        let x = tokens
                            .next()
                            .unwrap_or_else(|| panic!("vertex token to have a x in \"{}\"", s))
                            .parse::<f64>()
                            .unwrap_or_else(|_| panic!("vertex x should be an f64 in \"{}\"", s));
                        let y = tokens
                            .next()
                            .unwrap_or_else(|| panic!("vertex token to have a y in \"{}\"", s))
                            .parse::<f64>()
                            .unwrap_or_else(|_| panic!("vertex y should be an f64 in \"{}\"", s));
                        let z = tokens
                            .next()
                            .unwrap_or_else(|| panic!("vertex token to have a z in \"{}\"", s))
                            .parse::<f64>()
                            .unwrap_or_else(|_| panic!("vertex z should be an f64 in \"{}\"", s));
                        result.vertices.push(Tuple::point(x, y, z));
                    }
                    "f" => {
                        let v1 = tokens
                            .next()
                            .unwrap_or_else(|| panic!("face should have a v1 in \"{}\"", s))
                            .parse::<usize>()
                            .unwrap_or_else(|_| panic!("face v1 should be a usize in \"{}\"", s));

                        let mut v2 = tokens
                            .next()
                            .unwrap_or_else(|| panic!("face should have a v2 in \"{}\"", s))
                            .parse::<usize>()
                            .unwrap_or_else(|_| panic!("face v2 should be a usize in \"{}\"", s));

                        'face_loop: loop {
                            let v3 = match tokens.next() {
                                Some(v3_str) => v3_str.parse::<usize>().unwrap_or_else(|_| {
                                    panic!("face v3 should be a usize in \"{}\"", s)
                                }),
                                None => break 'face_loop,
                            };

                            let p1 = result.vertices(v1);
                            let p2 = result.vertices(v2);
                            let p3 = result.vertices(v3);

                            let mut triangle = Shape::triangle(p1, p2, p3);
                            if let Some(&g) = current_group {
                                Shape::add_child(g, &mut triangle);
                            } else {
                                Shape::add_child(&mut result.default_group, &mut triangle);
                            }

                            v2 = v3;
                        }
                    }
                    "g" => {
                        let current_group_name = tokens
                            .next()
                            .unwrap_or_else(|| panic!("group should have a name in \"{}\"", s));

                        let current_group = Some(Shape::group());

                        result.named_groups.insert(
                            current_group_name.to_string(),
                            current_group.unwrap().clone(),
                        );
                    }
                    _ => result.ignored_lines += 1,
                }
            }
        }

        result
    }

    fn vertices(&self, one_based_index: usize) -> Tuple {
        self.vertices[one_based_index - 1]
    }
}

#[cfg(test)]
mod tests {

    use crate::shape::ShapeKind;

    use super::*;

    #[test]
    fn ignoring_unrecognized_lines() {
        let gibberish = String::from(
            "
            There was a young lady named Bright
            who traveled much faster than light.
            She set out one day
            in a relative way,
            and came back the previous night.
            ",
        );
        let parser = Parser::from_obj_file(&gibberish);
        assert_eq!(parser.ignored_lines, 5);
    }

    #[test]
    fn vertex_records() {
        let file = "
        v -1 1 0
        v -1.0000 0.5000 0.0000
        v 1 0 0
        v 1 1 0
        ";

        let parser = Parser::from_obj_file(&file);

        assert_eq!(parser.vertices(1), Tuple::point(-1, 1, 0));
        assert_eq!(parser.vertices(2), Tuple::point(-1.0, 0.5, 0.0));
        assert_eq!(parser.vertices(3), Tuple::point(1, 0, 0));
        assert_eq!(parser.vertices(4), Tuple::point(1, 1, 0));
    }

    #[test]
    fn parsing_triangle_faces() {
        let file = "
        v -1 1 0
        v -1 0 0
        v 1 0 0
        v 1 1 0

        f 1 2 3
        f 1 3 4
        ";

        let parser = Parser::from_obj_file(&file);
        let g = parser.default_group.borrow();
        let t1 = g.shapes().unwrap()[0].borrow();
        let t2 = g.shapes().unwrap()[1].borrow();

        // TODO: verticles function that takes 1-based indexes?

        match &t1.kind {
            ShapeKind::Triangle { p1, p2, p3, .. } => {
                assert_eq!(p1, &parser.vertices(1));
                assert_eq!(p2, &parser.vertices(2));
                assert_eq!(p3, &parser.vertices(3));
            }
            _ => panic!("Not a triangle!"),
        }
        match &t2.kind {
            ShapeKind::Triangle { p1, p2, p3, .. } => {
                assert_eq!(p1, &parser.vertices(1));
                assert_eq!(p2, &parser.vertices(3));
                assert_eq!(p3, &parser.vertices(4));
            }
            _ => panic!("Not a triangle!"),
        }
    }

    #[test]
    fn triangulating_polygons() {
        let file = "
      v -1 1 0
      v -1 0 0
      v 1 0 0
      v 1 1 0
      v 0 2 0

      f 1 2 3 4 5
    ";

        let parser = Parser::from_obj_file(&file);
        let g = parser.default_group.borrow();
        let t1 = g.shapes().unwrap()[0].borrow();
        let t2 = g.shapes().unwrap()[1].borrow();
        let t3 = g.shapes().unwrap()[2].borrow();

        match &t1.kind {
            ShapeKind::Triangle { p1, p2, p3, .. } => {
                assert_eq!(p1, &parser.vertices(1));
                assert_eq!(p2, &parser.vertices(2));
                assert_eq!(p3, &parser.vertices(3));
            }
            _ => panic!("Not a triangle!"),
        }
        match &t2.kind {
            ShapeKind::Triangle { p1, p2, p3, .. } => {
                assert_eq!(p1, &parser.vertices(1));
                assert_eq!(p2, &parser.vertices(3));
                assert_eq!(p3, &parser.vertices(4));
            }
            _ => panic!("Not a triangle!"),
        }
        match &t3.kind {
            ShapeKind::Triangle { p1, p2, p3, .. } => {
                assert_eq!(p1, &parser.vertices(1));
                assert_eq!(p2, &parser.vertices(4));
                assert_eq!(p3, &parser.vertices(5));
            }
            _ => panic!("Not a triangle!"),
        }
    }

    #[test]
    fn triangles_in_groups() {
        let file = "
            v -1 1 0
            v -1 0 0
            v 1 0 0
            v 1 1 0

            g FirstGroup
            f 1 2 3
            g SecondGroup
            f 1 3 4
    ";

        let parser = Parser::from_obj_file(&file);
        let g1 = parser.named_groups["FirstGroup"].borrow();
        let g2 = parser.named_groups["FirstGroup"].borrow();

        let t1 = g1.shapes().unwrap()[0].borrow();
        let t2 = g2.shapes().unwrap()[0].borrow();

        match &t1.kind {
            ShapeKind::Triangle { p1, p2, p3, .. } => {
                assert_eq!(p1, &parser.vertices(1));
                assert_eq!(p2, &parser.vertices(2));
                assert_eq!(p3, &parser.vertices(3));
            }
            _ => panic!("Not a triangle!"),
        }
        match &t2.kind {
            ShapeKind::Triangle { p1, p2, p3, .. } => {
                assert_eq!(p1, &parser.vertices(1));
                assert_eq!(p2, &parser.vertices(3));
                assert_eq!(p3, &parser.vertices(4));
            }
            _ => panic!("Not a triangle!"),
        }
    }

    /*
    #[test]
    fn converting_an_obj_file_to_a_group() {
      Given file ← the file "triangles.obj"
        And parser ← parse_obj_file(file)
      When g ← obj_to_group(parser)
      Then g includes "FirstGroup" from parser
        And g includes "SecondGroup" from parser

    #[test]
    fn vertex_normal_records() {
      Given file ← a file containing:
        """
        vn 0 0 1
        vn 0.707 0 -0.707
        vn 1 2 3
        """
      When parser ← parse_obj_file(file)
      Then parser.normals[1] = vector(0, 0, 1)
        And parser.normals[2] = vector(0.707, 0, -0.707)
        And parser.normals[3] = vector(1, 2, 3)

    #[test]
    fn faces_with_normals() {
      Given file ← a file containing:
        """
        v 0 1 0
        v -1 0 0
        v 1 0 0

        vn -1 0 0
        vn 1 0 0
        vn 0 1 0

        f 1//3 2//1 3//2
        f 1/0/3 2/102/1 3/14/2
        """
      When parser ← parse_obj_file(file)
        And g ← parser.default_group
        And t1 ← first child of g
        And t2 ← second child of g
      Then t1.p1 = parser.vertices[1]
        And t1.p2 = parser.vertices[2]
        And t1.p3 = parser.vertices[3]
        And t1.n1 = parser.normals[3]
        And t1.n2 = parser.normals[1]
        And t1.n3 = parser.normals[2]
        And t2 = t1
        */
}
