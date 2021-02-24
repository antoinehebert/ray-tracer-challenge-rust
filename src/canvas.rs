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
        self.pixels[x + y * self.width]
    }

    fn set_pixel(&mut self, x: usize, y: usize, color: Color) {
        self.pixels[x + y * self.width] = color;
    }

    fn to_ppm(&self) -> String {
        let mut result = String::from("P3\n"); // Version/Flavor of PPM
        result.push_str(&format!("{} {}\n", self.width, self.height));
        result.push_str("255\n"); // Max color value

        for y in 0..self.height {
            let mut len = 0;
            for x in 0..self.width {
                let color = self.get_pixel(x, y);

                let mut color_strs = vec![
                    Self::color_to_str(color.red),
                    Self::color_to_str(color.green),
                    Self::color_to_str(color.blue),
                ];
                for color_str in color_strs {
                    if len + color_str.len() + 1 > 70 {
                        result.push_str("\n");
                        len = 0;
                    }
                    if len > 0 {
                        result.push_str(" ");
                        len += 1;
                    }
                    result.push_str(&color_str);
                    len += color_str.len();
                }
            }
            result.push_str("\n");
        }

        result
    }

    // Move to Color.
    fn color_to_str(color: f32) -> String {
        ((color.clamp(0., 1.) * 255.).round() as i32).to_string()
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

    #[test]
    fn constructing_the_ppm_header() {
        let c = Canvas::new(5, 3);
        let ppm = c.to_ppm();

        let lines: Vec<_> = ppm.split("\n").collect();

        assert_eq!(lines[0], "P3");
        assert_eq!(lines[1], "5 3");
        assert_eq!(lines[2], "255");
    }

    #[test]
    fn constructing_the_ppm_pixel_data() {
        let mut c = Canvas::new(5, 3);
        let c1 = Color::new(1.5, 0., 0.);
        let c2 = Color::new(0., 0.5, 0.);
        let c3 = Color::new(-0.5, 0., 1.);

        c.set_pixel(0, 0, c1);
        c.set_pixel(2, 1, c2);
        c.set_pixel(4, 2, c3);

        let ppm = c.to_ppm();

        let lines: Vec<_> = ppm.split("\n").collect();

        assert_eq!(lines.len(), 7);
        assert_eq!(lines[3], "255 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
        assert_eq!(lines[4], "0 0 0 0 0 0 0 128 0 0 0 0 0 0 0");
        assert_eq!(lines[5], "0 0 0 0 0 0 0 0 0 0 0 0 0 0 255");
    }

    #[test]
    fn splitting_long_lines_in_ppm_files() {
        let mut c = Canvas::new(10, 2);
        let color = Color::new(1., 0.8, 0.6);
        for y in 0..c.height {
            for x in 0..c.width {
                c.set_pixel(x, y, color);
            }
        }

        let ppm = c.to_ppm();

        let lines: Vec<_> = ppm.split("\n").collect();

        assert_eq!(lines.len(), 8);
        assert_eq!(
            lines[3],
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines[4],
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
        assert_eq!(
            lines[5],
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines[6],
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
    }
    #[test]
    fn ppm_files_are_terminated_by_a_newline_character() {
        let c = Canvas::new(5, 3);
        let ppm = c.to_ppm();

        assert_eq!(
            ppm.chars().rev().nth(0).expect("ppm should not be empty"),
            '\n'
        );
    }
}
