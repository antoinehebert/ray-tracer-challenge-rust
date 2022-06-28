use std::{collections::HashMap, fs};

use crate::{shape::Shape, tuple::Tuple};

pub struct Parser {
    vertices: Vec<Tuple>,
    ignored_lines: usize,
    default_group: Shape,
    named_groups: HashMap<String, Shape>,
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

    // TODO: Make path relative to project root.
    pub fn from_obj_file(filename: &str) -> Self {
        let file_content = fs::read_to_string(filename)
            .expect(format!("something went wrong reading {filename}.").as_str());
        Self::from_obj_str(file_content.as_str())
    }

    fn from_obj_str(text: &str) -> Self {
        let text = text.clone();
        let mut result = Self::new();

        let mut current_group_name = None; // TODO: @Performance: Find a way to keep a reference to the current group instead?

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

                            let triangle = Shape::triangle(p1, p2, p3);
                            if let Some(group_name) = current_group_name {
                                result
                                    .named_groups
                                    .get_mut(group_name)
                                    .unwrap()
                                    .push_shape(triangle);
                            } else {
                                result.default_group.push_shape(triangle);
                            }

                            v2 = v3;
                        }
                    }
                    "g" => {
                        let group_name = tokens
                            .next()
                            .unwrap_or_else(|| panic!("group should have a name in \"{}\"", s));

                        result
                            .named_groups
                            .insert(group_name.to_string(), Shape::group());

                        current_group_name = Some(group_name);
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

    // Consumes the parser...
    pub fn obj_to_group(self) -> Shape {
        let mut g = Shape::group();
        g.push_shape(self.default_group);
        for (_name, shape) in self.named_groups {
            g.push_shape(shape);
        }

        g
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
        let parser = Parser::from_obj_str(&gibberish);
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

        let parser = Parser::from_obj_str(&file);

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

        let parser = Parser::from_obj_str(&file);
        let g = &parser.default_group;
        let t1 = &g.shapes().unwrap()[0];
        let t2 = &g.shapes().unwrap()[1];

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

        let parser = Parser::from_obj_str(&file);
        let g = &parser.default_group;
        let t1 = &g.shapes().unwrap()[0];
        let t2 = &g.shapes().unwrap()[1];
        let t3 = &g.shapes().unwrap()[2];

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
        let parser = Parser::from_obj_file("src/test/files/triangles.obj");
        let g1 = &parser.named_groups["FirstGroup"];
        let g2 = &parser.named_groups["SecondGroup"];

        let t1 = &g1.shapes().unwrap()[0];
        let t2 = &g2.shapes().unwrap()[0];

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
    fn converting_an_obj_file_to_a_group() {
        let parser = Parser::from_obj_file("src/test/files/triangles.obj");

        let default_group = parser.default_group.clone();
        let first_group = parser.named_groups["FirstGroup"].clone();
        let second_group = parser.named_groups["SecondGroup"].clone();

        let g = parser.obj_to_group();
        let shapes = g.shapes().unwrap();

        assert_eq!(shapes.len(), 3);
        assert!(shapes.contains(&default_group));
        assert!(shapes.contains(&first_group));
        assert!(shapes.contains(&second_group));
    }

    /*
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
