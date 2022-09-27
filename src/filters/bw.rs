use std::{f64::consts::PI};

use lazy_static::lazy_static;
use num_complex::{Complex};
use serde::{Deserialize, Serialize};

use crate::{math::*};
use super::core::{NormalizedBilinLaplace};
use super::cascaded::*;

#[derive(Debug, Clone)]
pub struct FilterOrderUnsupportedError(f64);
impl std::fmt::Display for FilterOrderUnsupportedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "filter order {} is unsupported!", self.0)
    }
}
impl std::error::Error for FilterOrderUnsupportedError { }

#[derive(Debug, Clone, PartialEq)]
pub struct ButterworthFilterOrder {
    base: usize,
    alpha: f64
}
impl ButterworthFilterOrder {
    pub fn new(base: usize, alpha: f64) -> ButterworthFilterOrder {
        ButterworthFilterOrder { base, alpha }
    }
    pub fn from_f64(val: f64) -> Result<ButterworthFilterOrder, FilterOrderUnsupportedError> {
        if val < 1.0 {
            return Err(FilterOrderUnsupportedError(val));
        }
        let base = val.floor();
        let alpha = val - base;
        Ok(ButterworthFilterOrder { base: base as usize, alpha })
    }
}

#[derive(Debug, Clone)]
pub struct FilterFactors {
    numerator: Vec<Fac>,
    denominator: Vec<Fac>
}
impl FilterFactors {
    fn new(numerator: Vec<Fac>, denominator: Vec<Fac>) -> Result<FilterFactors, NonmatchingCoefficients> {
        if numerator.len() != denominator.len() { return Err(NonmatchingCoefficients {  }); }
        Ok(FilterFactors { numerator, denominator })
    }
    // fn from_exact_iter<T: IntoIterator<Item = (Fac, Fac)> + ExactSizeIterator>(iter: T) -> FilterFactors {
    //     let mut numerator = Vec::with_capacity(iter.len());
    //     let mut denominator = Vec::with_capacity(iter.len());
    //     for (n, d) in iter {
    //         numerator.push(n);
    //         denominator.push(d);
    //     }
    //     Self::new(numerator, denominator)
    // }
    // fn from_exact_results_iter<T: IntoIterator<Item = Result<(Fac, Fac), Box<dyn std::error::Error>>> + ExactSizeIterator>(iter: T) -> Result<FilterFactors, Box<dyn std::error::Error>> {
    //     let mut numerator = Vec::with_capacity(iter.len());
    //     let mut denominator = Vec::with_capacity(iter.len());
    //     for val in iter {
    //         let (n, d) = val?;
    //         numerator.push(n);
    //         denominator.push(d);
    //     }
    //     Ok(Self::new(numerator, denominator))
    // }
    fn len(&self) -> usize {
        self.numerator.len() //self.denominator.len() works too
    }
    fn iter(&self) -> impl Iterator<Item = (&Fac, &Fac)> + ExactSizeIterator {
        self.numerator.iter().zip(self.denominator.iter())
    }
    fn iter_mut(&mut self) -> impl Iterator<Item = (&mut Fac, &mut Fac)> + ExactSizeIterator {
        self.numerator.iter_mut().zip(self.denominator.iter_mut())
    }
    fn into_iter(self) -> impl Iterator<Item = (Fac, Fac)> + ExactSizeIterator {
        self.numerator.into_iter().zip(self.denominator.into_iter())
    }
}
#[derive(Debug, Clone)]
pub struct Butterworth {
    T: f64,
    freq: f64,
    order: ButterworthFilterOrder,
    factors: FilterFactors,
    parameters: ButterworthParameters,
    filter: Cascaded<NormalizedBilinLaplace>,
    global_gain: f64,
}

#[derive(Debug, Clone)]
pub enum BiquadTransformation {
    Chain(Vec<BiquadTransformation>),
    MultiplyOut(Box<BiquadTransformation>, Box<BiquadTransformation>),
    Q(f64),
    InvertCutoff,
    ShiftCutoff(f64),
    ShiftCutoffInv(f64),
    Invert,
    Noop,
}
impl BiquadTransformation {
    fn process(&self, factors: &mut FilterFactors) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            BiquadTransformation::Chain(chain) => {
                for transform in chain {
                    transform.process(factors)?;
                }
                Ok(())
            },
            BiquadTransformation::MultiplyOut(t1, t2) => {
                let mut factors2 = factors.clone();
                t1.process(factors)?;
                t2.process(&mut factors2)?;
                for ((f1_n, f1_d), (f2_n, f2_d)) in factors.iter_mut().zip(factors2.into_iter()) {
                    let mul_n = f1_n.multiply(&f2_n);
                    let mul_d = f1_d.multiply(&f2_d);
                    *f1_n = mul_n;
                    *f1_d = mul_d;
                }
                Ok(())
            },
            _ => {
                for biquad in factors.iter_mut() {
                    self.process_biquad(biquad)?;
                }
                Ok(())
            }
        }
    }
    fn process_biquad(&self, biquad: (&mut Fac, &mut Fac)) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            BiquadTransformation::Invert => {
                std::mem::swap(biquad.0, biquad.1);
                Ok(())
            },
            BiquadTransformation::Q(_) |
            BiquadTransformation::InvertCutoff |
            BiquadTransformation::ShiftCutoff(_) |
            BiquadTransformation::ShiftCutoffInv(_) |
            BiquadTransformation::Noop => {
                self.process_one(biquad.0)?;
                self.process_one(biquad.1)?;
                Ok(())
            },
            _ => return Err(Box::new(NotImplementedError(String::from("transformation unsupported through `process_biquad`"))))
        }
    }
    fn process_one(&self, factor: &mut Fac) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            BiquadTransformation::Q(q) => {
                let mut factor_out = match factor {
                    Fac::Linear(_) => factor.multiply(factor),
                    Fac::Quadratic(_) => factor.clone(),
                    _ => return Err(Box::new(NotImplementedError(String::from("transformation type 'Q' is not implemented for factors of degree higher than 2 (quadratic)"))))
                };
                factor_out.get_mut()[1] *= q;
                *factor = factor_out;
                Ok(())
            },
            BiquadTransformation::InvertCutoff => {
                factor.get_mut().reverse();
                Ok(())
            },
            BiquadTransformation::ShiftCutoff(shift) => {
                for (i, coefficient) in factor.get_mut().iter_mut().enumerate() {
                    *coefficient *= shift.powi(i as i32);
                }
                Ok(())
            },
            BiquadTransformation::ShiftCutoffInv(shift) => {
                for (i, coefficient) in (0..factor.get().len()).rev().zip(factor.get_mut().iter_mut()) {
                    *coefficient *= shift.powi(i as i32);
                }
                Ok(())
            },
            BiquadTransformation::Noop => Ok(()),
            _ => return Err(Box::new(NotImplementedError(String::from("transformation unsupported through `process_one`"))))
        }
    }
}
#[derive(Debug, Clone)]
pub struct ButterworthParameters {
    /// Quality factor
    /// 
    /// The Q-factor measures the distance of the poles to the jw axis.
    /// The uppermost and bottommost poles (nearest to the target frequency) affect the resonance directly, 
    /// which is why this Q-factor is applied to those poles only. 
    /// 
    /// ## How's it applied?
    /// α = 1 / Q
    /// Where α in a 2nd-order Butterworth low-pass filter's denominator is commonly sqrt(2). Yeah, it's that middle term.
    /// It is easy to derive what α should be for the separate biquads; it's just the Butterworth factors' z^1 term.
    /// Thus it can be seen that the Q-factor actually changes for different orders of the filter.
    /// However the Q-factor is supposed to be a thing that you can change to modify the shape of the filter!
    /// How can we turn it into something that won't change with the order?
    /// 
    /// Let α be the z^1 term of an arbitrary biquad.
    /// 
    /// The Q factor is inversely proportional to α. Divide the α by sqrt(2),
    /// 
    /// x = α / sqrt(2)
    /// 
    /// and now x contains the scale factor of that α in relation to the base case (sqrt(2) and the 
    /// 2nd-order butterworth low-pass).
    /// Now use that scale factor to scale the α_u accordingly(this is a different α, calculated 
    /// from the user-provided Q-factor which assumes that the filter is just a simple 
    /// 2nd-order butterworth low-pass):
    /// 
    /// α_u = 1 / Q_u
    /// 
    /// α = (α / sqrt(2)) * α_u
    /// α = α * (α_u / sqrt(2))
    /// α *= α_u / sqrt(2)
    /// 
    /// This is only done for the pole pairs nearest to the target frequency since changing the other poles 
    /// will completely change the frequency response of the filter.
    pub transform: BiquadTransformation,
}
#[derive(Serialize, Deserialize, Debug, Clone)]
struct ButterworthFractionalCoefficients {
    base: f64,
    alpha: f64,
    b: Vec<Fac>,
    a: Vec<Fac>
}
impl Butterworth {
    fn get_factors(order: &ButterworthFilterOrder) -> FilterFactors {
        trait OrderMapDiscrete {
            fn map_discrete(&self, approximants_per_base: usize) -> Option<usize>;
        }
        impl OrderMapDiscrete for ButterworthFilterOrder {
            fn map_discrete(&self, approximants_per_base: usize) -> Option<usize> {
                let index: f64 = ((self.base - 1) * approximants_per_base) as f64 + (self.alpha / (1.0 / (approximants_per_base as f64 + 1.0)) - 1.0).round();
                if (index + 1.0).abs() < f64::EPSILON {
                    None
                } else {
                    Some(index as usize)
                }
            }
        }
        lazy_static! {
            static ref COEFFICIENTS_STR: &'static str = include_str!("fractional_coefficients.json");
            static ref COEFFICIENTS: Vec<ButterworthFractionalCoefficients> = serde_json::from_str(&COEFFICIENTS_STR).expect("failed to load fractional coefficients!");
        }
        let denominator_factors;
        let numerator_factors;
        if let Some(approximant_index) = order.map_discrete(99) { //TODO: Don't hardcode 99
            let coefficients = &COEFFICIENTS[approximant_index];
            // println!("{:?}", coefficients);
            numerator_factors = coefficients.b.clone();
            denominator_factors = coefficients.a.clone();
            // println!("{:?} {:?}", numerator_factors, denominator_factors);
        } else {
            denominator_factors = Self::polynomial_factors(order.base);
            numerator_factors = denominator_factors.iter().map(|x| match x {
                Fac::Linear(_) => Fac::Linear([0.0, 1.0]),
                Fac::Quadratic(_) => Fac::Quadratic([0.0, 0.0, 1.0]),
                _ => panic!("this is highly unexpected")
            }).collect();
        }
        let factors = FilterFactors::new(numerator_factors, denominator_factors).unwrap();//NOTE: This should never error since the numerator factors are explicitly matched by the above function
        factors
    }
    pub fn new(T: f64, freq: f64, order: ButterworthFilterOrder, parameters: ButterworthParameters, global_gain: f64) -> Butterworth {
        let factors = Self::get_factors(&order);
        let filter = Self::apply_q(T, freq, &factors, &parameters).unwrap();//NOTE: Will never error since the return of polynomial_factors is always of degree 2
        Butterworth { T, freq, order, factors, parameters, filter, global_gain }
    }
    fn polynomial_factors(order: usize) -> Vec<Fac> {
        let mut factors = Vec::with_capacity(order / 2 + 1);
        if order % 2 == 1 {
            factors.push(Fac::Linear([1.0, 1.0]));
        }
        for k in 1..order/2+1 {
            factors.push(Fac::Quadratic([1.0, -2.0 * ((2.0 * k as f64 + order as f64 - 1.0) * PI / (2.0 * order as f64)).cos(), 1.0]));
        }
        factors
    }
    fn apply_q(T: f64, freq: f64, factors: &FilterFactors, parameters: &ButterworthParameters) -> Result<Cascaded<NormalizedBilinLaplace>, Box<dyn std::error::Error>> {
        let mut factors = factors.clone();
        parameters.transform.process(&mut factors)?;
        let mut filters = Vec::with_capacity(factors.len());
        for (n, d) in factors.iter() {
            filters.push(NormalizedBilinLaplace::new(T, n.get(), d.get(), freq)?);
        }
        Ok(Cascaded::new(filters))
    }
    fn apply_q_inplace(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let mut factors = self.factors.clone();
        self.parameters.transform.process(&mut factors)?;
        let filters = self.filter.filters_mut();
        let mut last_i = 0;
        for (i, (n, d)) in factors.iter().enumerate() {
            if i < filters.len() {
                filters[i].set_s_transfer(n.get(), d.get())?;
            } else {
                filters.push(NormalizedBilinLaplace::new(self.T, n.get(), d.get(), self.freq)?);
            }
            last_i = i;
        }
        filters.truncate(last_i + 1);
        Ok(())
    }
    pub fn set_sampling_period(&mut self, T: f64) -> Result<(), Box<dyn std::error::Error>> {
        self.T = T;
        self.apply_q_inplace()
    }
    pub fn set_frequency(&mut self, freq: f64) -> Result<(), Box<dyn std::error::Error>> {
        self.freq = freq;
        self.apply_q_inplace()
    }
    pub fn set_order(&mut self, order: ButterworthFilterOrder) -> Result<(), Box<dyn std::error::Error>> {
        if order != self.order { //TODO: Float equality
            self.order = order;
            self.factors = Self::get_factors(&self.order);
            self.apply_q_inplace()
        } else {
            Ok(())
        }
    }
    pub fn set_parameters(&mut self, parameters: ButterworthParameters) -> Result<(), Box<dyn std::error::Error>> {
        self.parameters = parameters;
        self.apply_q_inplace()
    }
    pub fn set_global_gain(&mut self, global_gain: f64) {
        self.global_gain = global_gain;
    }
    pub fn filter(&self) -> &Cascaded<NormalizedBilinLaplace> {
        &self.filter
    }
}
impl Filter for Butterworth {
    fn clear(&mut self) {
        self.filter.clear();
    }
    fn compute(&mut self, signal: f64) -> f64 {
        self.filter.compute(signal) * self.global_gain
    }
}
impl FR for Butterworth {
    fn frequency_response(&self, x: Complex<f64>) -> Complex<f64> {
        self.filter.frequency_response(x) * self.global_gain
    }
}

