use crate::utils::*;
use std::ops;

pub static BLACK: Color = Color {
    red: 0.,
    green: 0.,
    blue: 0.,
};

pub static RED: Color = Color {
    red: 1.,
    green: 0.,
    blue: 0.,
};

pub static GREEN: Color = Color {
    red: 0.,
    green: 1.,
    blue: 0.,
};

pub static BLUE: Color = Color {
    red: 0.,
    green: 0.,
    blue: 1.,
};

pub static WHITE: Color = Color {
    red: 1.,
    green: 1.,
    blue: 1.,
};

#[derive(Debug, Copy, Clone)]
pub struct Color {
    pub red: f64,
    pub green: f64,
    pub blue: f64,
}

impl Color {
    pub fn new(red: f64, green: f64, blue: f64) -> Self {
        Self { red, green, blue }
    }
}

impl PartialEq for Color {
    fn eq(&self, other: &Self) -> bool {
        is_almost_equal(self.red, other.red)
            && is_almost_equal(self.green, other.green)
            && is_almost_equal(self.blue, other.blue)
    }
}

impl ops::Add<Color> for Color {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(
            self.red + rhs.red,
            self.green + rhs.green,
            self.blue + rhs.blue,
        )
    }
}

impl ops::Sub<Color> for Color {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(
            self.red - rhs.red,
            self.green - rhs.green,
            self.blue - rhs.blue,
        )
    }
}

impl ops::Mul<f64> for Color {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.red * rhs, self.green * rhs, self.blue * rhs)
    }
}

// Used when blending colors.
impl ops::Mul<Color> for Color {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.red * rhs.red,
            self.green * rhs.green,
            self.blue * rhs.blue,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use crate::assert_almost_eq;

    #[test]
    fn colors_are_red_green_blue() {
        let color = Color::new(-0.5, 0.4, 1.7);

        assert_eq!(color.red, -0.5);
        assert_eq!(color.green, 0.4);
        assert_eq!(color.blue, 1.7);
    }

    #[test]
    fn adding_colors() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert_eq!(c1 + c2, Color::new(1.6, 0.7, 1.0));
    }

    #[test]
    fn subtracting_colors() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert_eq!(c1 - c2, Color::new(0.2, 0.5, 0.5));
    }

    #[test]
    fn multiplying_a_color_by_a_scalar() {
        let c = Color::new(0.2, 0.3, 0.4);
        assert_eq!(c * 2., Color::new(0.4, 0.6, 0.8));
    }

    #[test]
    fn multiplying_a_color_by_a_color() {
        let c1 = Color::new(1., 0.2, 0.4);
        let c2 = Color::new(0.9, 1.0, 0.1);

        assert_eq!(c1 * c2, Color::new(0.9, 0.2, 0.04));
    }
}
