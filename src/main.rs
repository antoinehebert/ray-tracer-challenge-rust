#[derive(Debug)]
struct Tuple {
    x: f32,
    y: f32,
    z: f32,
    w: f32,
}

impl Tuple {
    fn new(x: f32, y: f32, z: f32, w: f32) -> Tuple {
        Tuple { x, y, z, w }
    }

    fn is_point(&self) -> bool {
        self.w == 1.0
    }
    fn is_vector(&self) -> bool {
        self.w == 0.0
    }

    fn point(x: f32, y: f32, z: f32) -> Tuple {
        Tuple::new(x, y, z, 1.0)
    }

    fn vector(x: f32, y: f32, z: f32) -> Tuple {
        Tuple::new(x, y, z, 0.0)
    }
}

impl PartialEq for Tuple {
    fn eq(&self, other: &Self) -> bool {
        is_equal(self.x, other.x)
            && is_equal(self.y, other.y)
            && is_equal(self.z, other.z)
            && is_equal(self.w, other.w)
    }
}

pub fn is_equal(a: f32, b: f32) -> bool {
    (a - b).abs() < f32::EPSILON
}

fn main() {
    println!("Hello, world!");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_tuple_with_w_eq_1_is_a_point() {
        let t = Tuple::new(4.3, -4.2, 3.1, 1.0);

        assert_eq!(t.x, 4.3);
        assert_eq!(t.y, -4.2);
        assert_eq!(t.z, 3.1);
        assert_eq!(t.w, 1.0);

        assert!(t.is_point());
        assert!(!t.is_vector());
    }

    #[test]
    fn a_tuple_with_w_eq_0_is_a_vector() {
        let t = Tuple::new(4.3, -4.2, 3.1, 0.0);

        assert_eq!(t.x, 4.3);
        assert_eq!(t.y, -4.2);
        assert_eq!(t.z, 3.1);
        assert_eq!(t.w, 0.0);

        assert!(!t.is_point());
        assert!(t.is_vector());
    }

    #[test]
    fn point_creates_tuples_with_w_eq_1() {
        let point = Tuple::point(4.0, -4.0, 3.0);
        assert_eq!(point, Tuple::new(4.0, -4.0, 3.0, 1.0));
    }

    #[test]
    fn vector_creates_tuples_with_w_0() {
        let vector = Tuple::vector(4.0, -4.0, 3.0);
        assert_eq!(vector, Tuple::new(4.0, -4.0, 3.0, 0.0));
    }
}
