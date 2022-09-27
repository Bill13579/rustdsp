use num_complex::{Complex};

use crate::{math::*};

#[derive(Debug, Clone)]
pub struct GainFilter {
    gain: f64
}
impl GainFilter {
    pub fn new(gain: f64) -> GainFilter {
        GainFilter { gain }
    }
    pub fn set_gain(&mut self, gain: f64) {
        self.gain = gain;
    }
}
impl Filter for GainFilter {
    fn clear(&mut self) {  }
    fn compute(&mut self, signal: f64) -> f64 {
        signal * self.gain
    }
}
impl FR for GainFilter {
    fn frequency_response(&self, _: num_complex::Complex<f64>) -> num_complex::Complex<f64> {
        Complex::new(self.gain, 0.0)
    }
}