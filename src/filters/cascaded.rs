use num_complex::{Complex};
use rustfft::num_traits::One;

use crate::{math::*};

#[derive(Debug, Clone)]
pub struct Cascaded<T> where T: Filter + FR {
    filters: Vec<T>
}
impl<T> Cascaded<T> where T: Filter + FR {
    pub fn new(filters: Vec<T>) -> Cascaded<T> {
        Cascaded { filters }
    }
    pub fn filters(&self) -> &Vec<T> {
        &self.filters
    }
    pub fn filters_mut(&mut self) -> &mut Vec<T> {
        &mut self.filters
    }
}
impl<T> Filter for Cascaded<T> where T: Filter + FR {
    fn clear(&mut self) {
        for f in self.filters.iter_mut() {
            f.clear();
        }
    }
    fn compute(&mut self, mut signal: f64) -> f64 {
        for f in self.filters.iter_mut() {
            signal = f.compute(signal);
        }
        signal
    }
}
impl<T> FR for Cascaded<T> where T: Filter + FR {
    fn frequency_response(&self, x: num_complex::Complex<f64>) -> num_complex::Complex<f64> {
        let mut resp = Complex::<f64>::one();
        for f in self.filters.iter() {
            resp *= f.frequency_response(x);
        }
        resp
    }
}

pub trait ZTransformFilter: Filter + FR {  }
pub struct CascadedDyn<'a> {
    filters: &'a mut Vec<Box<dyn ZTransformFilter>>
}
impl<'a> CascadedDyn<'a> {
    pub fn new(filters: &'a mut Vec<Box<dyn ZTransformFilter>>) -> CascadedDyn<'a> {
        CascadedDyn { filters }
    }
    pub fn filters(&self) -> &Vec<Box<dyn ZTransformFilter>> {
        &self.filters
    }
    pub fn filters_mut(&mut self) -> &mut Vec<Box<dyn ZTransformFilter>> {
        &mut self.filters
    }
}
impl<'a> Filter for CascadedDyn<'a> {
    fn clear(&mut self) {
        for f in self.filters.iter_mut() {
            f.clear();
        }
    }
    fn compute(&mut self, mut signal: f64) -> f64 {
        for f in self.filters.iter_mut() {
            signal = f.compute(signal);
        }
        signal
    }
}
impl<'a> FR for CascadedDyn<'a> {
    fn frequency_response(&self, mut x: num_complex::Complex<f64>) -> num_complex::Complex<f64> {
        for f in self.filters.iter() {
            x = f.frequency_response(x);
        }
        x
    }
}

