use std::{f64::consts::PI};

use crate::{math::*};

#[derive(Debug, Clone)]
pub struct NormalizedBilinLaplace {
    T: f64,
    bilin: [Vec<Vec<f64>>; 2],
    freq: f64,
    filter: IIR
}
impl NormalizedBilinLaplace {
    pub fn new(T: f64, n: &[f64], d: &[f64], freq: f64) -> Result<NormalizedBilinLaplace, NonmatchingCoefficients> {
        let bilin = bilinear_transform(n, d)?;
        let filter = apply_k(&bilin, k_prewrap_and_scale_normalized(T, freq * 2.0 * PI));
        Ok(NormalizedBilinLaplace {
            T, bilin, freq, filter: IIR::new(filter)
        })
    }
    fn recalculate_filter_coefficients_from_bilin(&mut self) {
        let filter = apply_k(&self.bilin, k_prewrap_and_scale_normalized(self.T, self.freq * 2.0 * PI));
        self.filter.set_filter(filter);
    }
    pub fn set_sampling_period(&mut self, T: f64) {
        self.T = T;
        self.recalculate_filter_coefficients_from_bilin();
    }
    pub fn set_s_transfer(&mut self, n: &[f64], d: &[f64]) -> Result<(), NonmatchingCoefficients> {
        self.bilin = bilinear_transform(n, d)?;
        self.recalculate_filter_coefficients_from_bilin();
        Ok(())
    }
    pub fn set_frequency(&mut self, freq: f64) {
        self.freq = freq;
        self.recalculate_filter_coefficients_from_bilin();
    }
    pub fn filter(&self) -> &IIR {
        &self.filter
    }
}
impl Filter for NormalizedBilinLaplace {
    fn clear(&mut self) {
        self.filter.clear();
    }
    fn compute(&mut self, signal: f64) -> f64 {
        self.filter.compute(signal)
    }
}
impl FR for NormalizedBilinLaplace {
    fn frequency_response(&self, x: num_complex::Complex<f64>) -> num_complex::Complex<f64> {
        self.filter.frequency_response(x)
    }
}

