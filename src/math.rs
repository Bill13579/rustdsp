use std::{fmt, f64::consts::PI, iter, sync::{Arc, Mutex}, collections::{VecDeque}, time::Instant, ops::{AddAssign, SubAssign, MulAssign, DivAssign, RemAssign}};
use num_complex::{Complex, ComplexFloat};
use num_traits::{NumCast, Float, Signed};
use rustfft::{FftPlanner, num_traits::{Zero}};
use serde::{Serialize, Deserialize};

use crate::dtype::{ChunkedBuffer, RingBuffer};

#[derive(Debug, Clone)]
pub struct IsComplexError;
impl fmt::Display for IsComplexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "number is complex! cannot convert into real without losing data!")
    }
}
impl std::error::Error for IsComplexError { }

pub trait IntoReal<T> where T: PartialOrd<f64> + Signed + Clone {
    fn is_real(&self) -> bool;
    fn to_real(&self) -> Option<T>;
    fn try_into_real(self) -> Result<T, IsComplexError>;
}
impl<T> IntoReal<T> for Complex<T> where T: PartialOrd<f64> + Signed + Clone {
    fn is_real(&self) -> bool {
        //TODO: Is this enough?
        self.im.abs() <= 1e-10
    }
    fn to_real(&self) -> Option<T> {
        self.clone().try_into_real().ok()
    }
    fn try_into_real(self) -> Result<T, IsComplexError> {
        if self.is_real() {
            Ok(self.re)
        } else {
            Err(IsComplexError {})
        }
    }
}

#[derive(Debug, Clone)]
pub struct NotImplementedError(pub String);
impl fmt::Display for NotImplementedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", &self.0)
    }
}
impl std::error::Error for NotImplementedError { }
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(untagged)]
pub enum Fac {
    Linear([f64; 2]),
    Quadratic([f64; 3]),
    Cubic([f64; 4]),
    Arbitrary(Vec<f64>),
}
impl Fac {
    pub fn try_from_root(root: &Complex<f64>) -> Result<Fac, Box<dyn std::error::Error>> {
        Ok(Self::from_root(root.to_real().ok_or_else(|| IsComplexError {})?))
    }
    pub fn from_root(root: f64) -> Fac {
        let a = 1.0;
        let b = -root;
        Fac::Linear([a, b])
    }
    pub fn from_conjugates(r1: &Complex<f64>, r2: &Complex<f64>) -> Result<Fac, Box<dyn std::error::Error>> {
        let a = 1.0;
        let b = (-(r1 + r2)).try_into_real()?;
        let c = (r1 * r2).try_into_real()?;
        Ok(Fac::Quadratic([a, b, c]))
    }
    pub fn from_vec(vec: Vec<f64>) -> Fac {
        if vec.len() == 2 { Fac::Linear(vec.try_into().unwrap()) } // These should never error
        else if vec.len() == 3 { Fac::Quadratic(vec.try_into().unwrap()) }
        else if vec.len() == 4 { Fac::Cubic(vec.try_into().unwrap()) }
        else { Fac::Arbitrary(vec) }
    }
    pub fn multiply(&self, other: &Fac) -> Fac {
        let p1 = self.get();
        let p2 = other.get();
        let mut result = vec![0.0; p1.len() + p2.len() - 1];
        for i in 0..p1.len() {
            for j in 0..p2.len() {
                result[i+j] += p1[i] * p2[j];
            }
        }
        Self::from_vec(result)
    }
    pub fn get(&self) -> &[f64] {
        match self {
            Self::Linear(a) => a,
            Self::Quadratic(a) => a,
            Self::Cubic(a) => a,
            Self::Arbitrary(a) => a,
        }
    }
    pub fn get_mut(&mut self) -> &mut [f64] {
        match self {
            Self::Linear(a) => a,
            Self::Quadratic(a) => a,
            Self::Cubic(a) => a,
            Self::Arbitrary(a) => a,
        }
    }
    pub fn evaluate(&self, x: f64) -> f64 {
        let coefficients = self.get();
        coefficients.iter().zip((0..coefficients.len()).rev()).fold(0.0, |accum, (coefficient, power)| accum + coefficient * x.powi(power as i32))
    }
    pub fn roots(&self) -> Result<Vec<Complex<f64>>, NotImplementedError> {
        match self {
            Self::Linear(linear) => {
                Ok(vec![Complex::new(-linear[1] / linear[0], 0.0)])
            },
            Self::Quadratic(quadratic) => {
                let discriminant_sqrt = Complex::new(quadratic[1].powi(2) - 4.0 * quadratic[0] * quadratic[2], 0.0).powf(1.0 / 2.0);
                Ok(vec![
                    (-quadratic[1] + discriminant_sqrt) / (2.0 * quadratic[0]),
                    (-quadratic[1] - discriminant_sqrt) / (2.0 * quadratic[0])
                ])
            },
            Self::Cubic(cubic) => {
                let p = cubic[2] / cubic[0] - cubic[1].powi(2) / (3.0 * cubic[0].powi(2));
                let q = (2.0 * cubic[1].powi(3)) / (27.0 * cubic[0].powi(3)) - (cubic[1] * cubic[2]) / (3.0 * cubic[0].powi(2)) + cubic[3] / cubic[0];
                let t1 = -q / 2.0;
                let D = t1.powi(2) + (p / 3.0).powi(3);
                let t2 = Complex::new(D, 0.0).root(2);
                let u = (t1 + t2[0]).root_verbose(3, false);
                let v = (t1 + t2[1]).root_verbose(3, true);
                let shift = -cubic[1] / (3.0 * cubic[0]);

                let mut result: [Option<(f64, (Complex<f64>, Complex<f64>))>; 3] = [None; 3];
                for u in u.into_iter() {
                    for v in v.iter() {
                        let score = (u * v + p / 3.0).norm();
                        if result[0] == None || score < result[0].unwrap().0 {
                            result[2] = result[1];
                            result[1] = result[0];
                            result[0] = Some((score, (u, v.clone())));
                        } else if result[1] == None || score < result[1].unwrap().0 {
                            result[2] = result[1];
                            result[1] = Some((score, (u, v.clone())));
                        } else if result[2] == None || score < result[2].unwrap().0 {
                            result[2] = Some((score, (u, v.clone())));
                        }
                    }
                }

                Ok(result.into_iter().flatten().map(|(_, (u, v))| u + v + shift).collect())
            },
            Self::Arbitrary(_) => {
                Err(NotImplementedError(String::from("root-finding for polynomials of arbitrary degrees is not implemented!")))
            }
        }
    }
}

pub trait NthRoot<T> where T: Float + From<f64> + AddAssign + SubAssign + MulAssign + DivAssign + RemAssign {
    fn root(&self, nth: usize) -> Vec<Complex<T>>;
    fn root_verbose(&self, nth: usize, flip_angle_direction: bool) -> Vec<Complex<T>>;
    fn roots_of_unity(n: usize, inverse: bool) -> Vec<Complex<T>>;
}
impl<T> NthRoot<T> for Complex<T> where T: Float + From<f64> + AddAssign + SubAssign + MulAssign + DivAssign + RemAssign {
    fn root(&self, nth: usize) -> Vec<Complex<T>> {
        if self.im >= T::zero() {
            self.root_verbose(nth, false)
        } else {
            self.root_verbose(nth, true)
        }
    }
    fn root_verbose(&self, nth: usize, flip_angle_direction: bool) -> Vec<Complex<T>> {
        let mut unity = Self::roots_of_unity(nth, flip_angle_direction);
        let cast_usize = |value: usize| <T as NumCast>::from(value).unwrap();
        let cast_float = |value: f64| <T as From<f64>>::from(value);
        let base = Complex::from_polar(self.norm().powf(cast_float(1.0 / nth as f64)), self.arg() / cast_usize(nth));
        for root_of_unity in &mut unity {
            *root_of_unity *= base;
        }
        unity
    }
    fn roots_of_unity(n: usize, inverse: bool) -> Vec<Complex<T>> {
        let mut roots_of_unity: Vec<Complex<T>> = Vec::with_capacity(n);
        let cast_usize = |value: usize| <T as NumCast>::from(value).unwrap();
        let cast_float = |value: f64| <T as From<f64>>::from(value);
        for i in 0..n {
            let root;
            if inverse {
                root = Complex::from_polar(T::one(), -(cast_usize(i) * cast_float(2.0 * PI)) / cast_usize(n));
            } else {
                root = Complex::from_polar(T::one(), (cast_usize(i) * cast_float(2.0 * PI)) / cast_usize(n));
            }
            roots_of_unity.push(root);
        }
        roots_of_unity
    }
}

#[derive(Debug, Clone)]
pub struct NonmatchingCoefficients;
impl fmt::Display for NonmatchingCoefficients {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "transfer function coefficients must be matched in the numerator and the denominator!")
    }
}
impl std::error::Error for NonmatchingCoefficients {  }

pub fn pascal(row: i64) -> Vec<i64> {
    let mut current = 1;
    let mut vals = Vec::with_capacity(row as usize + 1);
    vals.push(current);
    for i in 1..(row+1) {
        current = (row+1-i) * current / i;
        vals.push(current);
    }
    vals
}

pub fn decibels(gain: f64) -> f64 {
    20.0 * gain.log10()
}
pub fn gain(decibels: f64) -> f64 {
    10.0_f64.powf(decibels / 20.0)
}

pub fn frequency(T: f64) -> f64 {
    1.0 / T
}
// m is zero indexed
pub fn fftfreq(m: usize, M: usize, T: Option<f64>) -> f64 {
    if let Some(T) = T {
        m as f64 / (M as f64 * T)
    } else {
        m as f64 / M as f64
    }
}
pub fn fftindex(f: f64, M: usize, T: Option<f64>) -> usize {
    if let Some(T) = T {
        (f * M as f64 * T) as usize
    } else {
        (f * M as f64) as usize
    }
}
// Radians per second or Radians per sample (T is sampling interval in seconds)
pub fn omega(f: f64, T: Option<f64>) -> f64 {
    if let Some(T) = T {
        2.0 * PI * f * T
    } else {
        2.0 * PI * f
    }
}
pub fn z(r: f64, omega: f64) -> Complex<f64> {
    Complex::from_polar(r, omega)
}

// Laplace to Z-Transform
fn bilinear_binomial(n: i64, p: i64) -> Vec<i64> {
    let mut vals = vec![0; n as usize + p as usize + 1];
    let n_pascal = pascal(n);
    let p_pascal = pascal(p);
    for (n_i, n_c) in n_pascal.iter().enumerate() {
        let n_c = if n_i % 2 != 0 { -n_c } else { *n_c };
        for (p_i, p_c) in p_pascal.iter().enumerate() {
            let c = n_c * p_c;
            let e = n_i + p_i;
            vals[e] += c;
        }
    }
    vals
}
fn partial_bilinear_transform(cf: &[f64]) -> Vec<Vec<f64>> {
    let mut vals = vec![vec![1.0; cf.len()]; cf.len()];
    for p in 0..cf.len() {
        let n = cf.len() - p - 1;
        for (e, term) in bilinear_binomial(n as i64, p as i64).into_iter().enumerate() {
            vals[e][p] *= term as f64 * cf[p];
        }
    }
    vals
}
pub fn bilinear_transform(n: &[f64], d: &[f64]) -> Result<[Vec<Vec<f64>>; 2], NonmatchingCoefficients> {
    if n.len() != d.len() { return Err(NonmatchingCoefficients {}); }
    Ok([partial_bilinear_transform(n), partial_bilinear_transform(d)])
}

pub fn to_analog_frequency(digital: f64, T: f64) -> f64 { //untested
    (2.0 / T) * (digital * T / 2.0).tan()
}
pub fn to_digital_frequency(analog: f64, T: f64) -> f64 { //untested
    (2.0 / T) * (analog * T / 2.0).atan()
}

fn partially_apply_k(bt: &Vec<Vec<f64>>, K: f64) -> Vec<f64> {
    let mut coefficients = Vec::with_capacity(bt.len());
    for components in bt.iter() {
        let mut sum = 0.0;
        for (i, component) in components.iter().enumerate() {
            let exponent = components.len() - 1 - i;
            sum += K.powi(exponent as i32) * component;
        }
        coefficients.push(sum);
    }
    coefficients
}
pub fn apply_k(filter: &[Vec<Vec<f64>>; 2], K: f64) -> [Vec<f64>; 2] {
    [partially_apply_k(&filter[0], K), partially_apply_k(&filter[1], K)]
}
pub fn k_standard(T: f64) -> f64 {
    2.0 / T
}
pub fn k_standard_and_scale_normalized(T: f64, critical_frequency: f64) -> f64 {
    2.0 / (T * critical_frequency)
}
pub fn k_prewrap(T: f64, critical_frequency: f64) -> f64 {
    critical_frequency / (critical_frequency * T / 2.0).tan()
}
pub fn k_prewrap_and_scale_normalized(T: f64, critical_frequency: f64) -> f64 {
    1.0 / (critical_frequency * T / 2.0).tan()
}

pub trait Filter {
    fn clear(&mut self);
    fn compute(&mut self, signal: f64) -> f64;
}

pub struct FFTConvolution<T> where T: FR + Clone {
    x: RingBuffer<Complex<f64>>,
    out: RingBuffer<f64>,
    last: VecDeque<Complex<f64>>,
    window_size: usize,
    T: f64,
    fr: T,
    fr_cache: Vec<Complex<f64>>,
    fft_planner: Arc<Mutex<FftPlanner<f64>>>,
}

impl<T> FFTConvolution<T> where T: FR + Clone {
    pub fn new(window_size: usize, T: f64, fr: T) -> FFTConvolution<T> {
        let mut conv = FFTConvolution {
            x: RingBuffer::new(window_size),
            out: RingBuffer::new(Self::padded_window_size(window_size)).initialize(0.0),
            last: VecDeque::from_iter(iter::repeat(Complex::zero()).take(window_size)),
            window_size,
            T,
            fr,
            fr_cache: Vec::new(),
            fft_planner: Arc::new(Mutex::new(FftPlanner::new())),
        };
        conv.refill_cache();
        conv
    }
    fn refill_cache(&mut self) {
        let window_size = self.window_size;
        self.fr_cache = Vec::with_capacity(window_size * 2);
        let sample_rate = 1.0 / self.T;
        for i in 0..window_size {
            let mut freq = fftfreq(i, window_size, Some(self.T));
            if freq >= sample_rate / 2.0 {
                freq = sample_rate - freq;
            }
            let sample = self.fr.frequency_response(z(1.0, omega(freq, Some(self.T))));
            self.fr_cache.push(Complex::new(sample.norm(), 0.0));
        }
        {
            let ifft = self.fft_planner.lock().unwrap().plan_fft_inverse(window_size);
            ifft.process(&mut self.fr_cache);
        }
        let norm_factor = self.fr_cache.len() as f64;
        for val in &mut self.fr_cache {
            *val /= norm_factor;
        }
        let swap_dist = window_size / 2;
        for i in 0..swap_dist {
            self.fr_cache.swap(i, i + swap_dist);
        }
        for _ in 0..window_size {
            self.fr_cache.push(Complex::zero());
        }
        {
            let fft = self.fft_planner.lock().unwrap().plan_fft_forward(window_size * 2);
            fft.process(&mut self.fr_cache);
        }
    }
    pub fn set_window_size(&mut self, window_size: usize) {
        self.x.to_capacity_front(Some(window_size));
        self.out.to_capacity_front(Some(window_size));
        self.out.fill_front(0.0);
        self.last = VecDeque::from_iter(iter::repeat(Complex::zero()).take(self.out.capacity()));
        self.window_size = window_size;
        self.refill_cache();
    }
    pub fn set_sampling_period(&mut self, T: f64) {
        self.T = T;
        self.refill_cache();
    }
    pub fn set_frequency_response(&mut self, fr: T) {
        self.fr = fr;
        self.refill_cache();
    }
    pub fn output_buffer(&self) -> &RingBuffer<f64> {
        &self.out
    }
    pub fn last_mut(&mut self) -> &mut VecDeque<Complex<f64>> {
        &mut self.last
    }
    fn padded_window_size(window_size: usize) -> usize {
        (window_size.next_power_of_two() + 1).next_power_of_two()
    }
}
impl<T> Filter for FFTConvolution<T> where T: FR + Clone {
    fn clear(&mut self) {
        self.x.clear();
        self.out.initialize_again(0.0);
        self.last = VecDeque::from_iter(iter::repeat(Complex::zero()).take(self.out.capacity()));
        self.fr_cache.clear();
    }
    fn compute(&mut self, signal: f64) -> f64 {
        let buffered_signal = self.out.pop_front().unwrap();
        self.out.push_back(0.0);

        if let Some(chunk) = self.x.buffer_back(Complex::new(signal, 0.0)) {
            let window_size = chunk.len();
            // let window_size = self.last.len() + chunk.len();
            // let mut buffer: Vec<Complex<f64>> = self.last.iter().copied().chain(chunk.iter().copied()).collect();
            // {
            //     let fft = self.fft_planner.lock().unwrap().plan_fft_forward(window_size);
            //     fft.process(&mut buffer);
            // }
            // // let fr = self.fr.frequency_responses(chunk.len(), window_size - chunk.len(), &mut self.fft_planner.lock().unwrap());
            // for (i, val) in buffer.iter_mut().enumerate() {
            //     // *val *= fr[i].norm(); //TODO: Strip phase?
            //     *val *= self.fr_cache[i].norm(); //TODO: Strip phase?
            // }
            // {
            //     let ifft = self.fft_planner.lock().unwrap().plan_fft_inverse(window_size);
            //     ifft.process(&mut buffer);
            // }
            // for (out_ref, buf_val) in self.out.inner_mut().iter_mut().zip(buffer.iter().skip(self.last.len())) {
            //     *out_ref += buf_val.re / window_size as f64; //TODO: Magnitude or Real part?
            // }
            // self.last = chunk;
            let padded_window_size = Self::padded_window_size(window_size);
            let mut buffer: Vec<Complex<f64>> = chunk.into_iter().chain(iter::repeat(Complex::zero()).take(padded_window_size - window_size)).collect();
            {
                let fft = self.fft_planner.lock().unwrap().plan_fft_forward(padded_window_size);
                fft.process(&mut buffer);
            }
            for (i, val) in buffer.iter_mut().enumerate() {
                *val *= self.fr_cache[i]; //TODO: Strip phase?
            }
            {
                let ifft = self.fft_planner.lock().unwrap().plan_fft_inverse(padded_window_size);
                ifft.process(&mut buffer);
            }
            for (out_ref, buf_val) in self.out.inner_mut().iter_mut().zip(buffer.into_iter()).take(padded_window_size) {
                *out_ref += buf_val.re / padded_window_size as f64; //TODO: Magnitude or Real part?
            }
        }
        
        buffered_signal
    }
}

#[derive(Debug, Clone)]
pub struct IIR {
    x: RingBuffer<f64>,
    y: RingBuffer<f64>,
    filter: [Vec<f64>; 2]
}

impl IIR {
    pub fn new(filter: [Vec<f64>; 2]) -> IIR {
        IIR {
            x: RingBuffer::new(filter[0].len()).initialize(0.0),
            y: RingBuffer::new(filter[1].len()).initialize(0.0),
            filter
        }
    }
    fn compute_with_coefficients(signal_buffer: &RingBuffer<f64>, coefficients: &Vec<f64>, normalization_factor: f64) -> f64 {
        let mut sum = 0.0;
        for (i, coefficient) in coefficients.iter().enumerate() {
            sum += signal_buffer.back_n(i).unwrap_or(&0.0) * (coefficient / normalization_factor);
        }
        sum
    }
    pub fn filter(&self) -> &[Vec<f64>; 2] {
        &self.filter
    }
    pub fn set_filter(&mut self, filter: [Vec<f64>; 2]) {
        self.x.to_capacity_front(Some(filter[0].len()));
        self.x.fill_front(0.0);
        self.y.to_capacity_front(Some(filter[1].len())); //TODO: Evaluate if it is better to just clear the Y buffer
        self.y.fill_front(0.0);
        self.filter = filter;
    }
}
impl Filter for IIR {
    fn clear(&mut self) {
        self.x.initialize_again(0.0);
        self.y.initialize_again(0.0);
    }
    fn compute(&mut self, signal: f64) -> f64 {
        self.x.push_back(signal);
        let x = Self::compute_with_coefficients(&self.x, &self.filter[0], self.filter[1][0]);
        self.y.push_back(0.0);
        let y = Self::compute_with_coefficients(&self.y, &self.filter[1], self.filter[1][0]);
        let new_y = x - y;
        if let Some(new_y_ref) = self.y.back_mut() {
            *new_y_ref = new_y;
        }
        new_y
    }
}

pub trait FR {
    fn frequency_response(&self, x: Complex<f64>) -> Complex<f64>;
    // fn frequency_responses(&self, L: usize, pad: usize, fft_planner: &mut FftPlanner<f64>) -> Vec<Complex<f64>>;
}

impl FR for IIR {
    fn frequency_response(&self, z: Complex<f64>) -> Complex<f64> {
        let num = self.filter[0].iter().enumerate().fold(Complex::<f64>::zero(), |accum, (i, coefficient)| accum + coefficient / z.powi(i as i32));
        let den = self.filter[1].iter().enumerate().fold(Complex::<f64>::zero(), |accum, (i, coefficient)| accum + coefficient / z.powi(i as i32));
        num / den
    }
    // fn frequency_responses(&self, L: usize, pad: usize, fft_planner: &mut FftPlanner<f64>) -> Vec<Complex<f64>> {
    //     let mut filter_clone = self.clone();
    //     let mut filter_impulse_response = Vec::with_capacity(L + pad);
    //     for signal in iter::repeat(1.0).take(1).chain(iter::repeat(0.0).take(L-1)) {
    //         filter_impulse_response.push(Complex::new(filter_clone.compute(signal), 0.0));
    //     }
    //     for _ in 0..pad {
    //         filter_impulse_response.push(Complex::zero());
    //     }
    //     let fft = fft_planner.plan_fft_forward(L + pad);
    //     fft.process(&mut filter_impulse_response);
    //     filter_impulse_response
    // }
}

