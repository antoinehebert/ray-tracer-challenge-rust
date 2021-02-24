use crate::color::*;

struct Canvas {
    width: usize,
    height: usize,
    pixels: Vec<Color>,
}

impl Canvas {
    fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            pixels: vec![BLACK; width * height],
        }
    }

    fn get_pixel(&self, x: usize, y: usize) -> Color {
        self.pixels[x * self.width + y]
    }

    fn set_pixel(&mut self, x: usize, y: usize, color: Color) {
        self.pixels[x * self.width + y] = color;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creating_a_canvas() {
        let c = Canvas::new(10, 20);

        assert_eq!(c.width, 10);
        assert_eq!(c.height, 20);

        for x in 0..c.width {
            for y in 0..c.height {
                assert_eq!(c.get_pixel(x, y), BLACK);
            }
        }
    }

    #[test]
    fn writing_pixels_to_canvas() {
        let mut c = Canvas::new(10, 20);

        c.set_pixel(2, 3, RED);

        assert_eq!(c.get_pixel(2, 3), RED);
    }
}
