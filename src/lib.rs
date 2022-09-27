mod dtype;
pub mod math;
pub mod filters;
pub mod prelude;

use std::{fmt, collections::{VecDeque, vec_deque}, iter};
use crate::math::{IIR, Filter};
use crate::dtype::{ChunkedBuffer, RingBuffer};

fn polynomial_square(polynomial: &[f64]) -> Vec<f64> {
    let mut product = vec![0.0; (polynomial.len() - 1) * 2 + 1];
    for i in 0..polynomial.len() {
        for j in 0..polynomial.len() {
            product[j+i] += polynomial[i] * polynomial[j];
        }
    }
    product
}

// pub trait TIIR {
//     fn get_truncation_filter(&self, L: usize) -> IIR;
//     fn get_truncation_filter2(&self, L: usize) -> IIR;
// }

// impl TIIR for IIR {
//     fn get_truncation_filter(&self, L: usize) -> IIR {
//         let time_shifted: Vec<f64> = self.filter()[0].iter().copied().chain(iter::repeat(0.0).take(L)).collect();
//         let (_, remainder) = synthetic_division(&time_shifted, &self.filter()[1]);
//         IIR::new([remainder, self.filter()[1].clone()])
//     }
//     fn get_truncation_filter2(&self, L: usize) -> IIR {
//         let mut filter_clone = self.clone();
//         let mut filter_impulse_response = Vec::with_capacity(L+1);
//         for signal in iter::repeat(1.0).take(1).chain(iter::repeat(0.0).take(L)) {
//             filter_impulse_response.push(filter_clone.compute(signal));
//         }
//         let K = self.filter()[0].len() - 1;
//         let mut parameters = Vec::with_capacity(K);
//         for i in 0..K {
//             let mut sum = 0.0;
//             for j in 0..(K-i) {
//                 sum += -(self.filter()[1][i+j+1] / self.filter()[1][0]) * filter_impulse_response[filter_impulse_response.len()-2-j];//TODO: Minus by 1?
//             }
//             parameters.push(sum);
//         }
//         IIR::new([parameters, self.filter()[1].clone()])
//     }
// }

// pub struct PowellChau {
//     top: IIR,
//     bot: IIR,
//     tf: IIR,
//     buffer: RingBuffer<f64>,
//     buffer_truncation: RingBuffer<f64>,
//     out_buffer: RingBuffer<f64>,
//     buffer_iter: Option<vec_deque::IntoIter<f64>>,
//     buffer_truncation_iter: Option<vec_deque::IntoIter<f64>>,
//     out_buffer_iter: Option<vec_deque::IntoIter<f64>>,
//     switch: bool,
//     reverse_time: bool,
// }

// impl PowellChau {
//     pub fn new(f1: IIR, f2: IIR, tf: IIR, buffer_size: usize, reverse_time: bool) -> Self {
//         Self {
//             top: f1,
//             bot: f2,
//             tf: tf,
//             buffer: RingBuffer::new(buffer_size),
//             buffer_truncation: RingBuffer::new(buffer_size),
//             out_buffer: RingBuffer::new(buffer_size),
//             buffer_iter: None,
//             buffer_truncation_iter: None,
//             out_buffer_iter: None,
//             switch: true,
//             reverse_time
//         }
//     }
//     pub fn buffer(&mut self, signal: f64) -> f64 {
//         if self.reverse_time {
//             let mut buffered_signal = 0.0;

//             if let Some(iter) = self.out_buffer_iter.as_mut() {
//                 let next = iter.next();
//                 if let Some(signal) = next {
//                     buffered_signal = signal;
//                 }
//             }

//             let internal_buffer_out = self.internal_buffer(signal);
//             if let Some(new_chunk) = self.out_buffer.buffer_front(internal_buffer_out) { // Switch to buffer_back for normal processing
//                 self.out_buffer_iter = Some(new_chunk.into_iter());
//             }

//             buffered_signal
//         } else {
//             self.internal_buffer(signal)
//         }
//     }
//     pub fn internal_buffer(&mut self, signal: f64) -> f64 {
//         let active;
//         let passive;
//         if self.switch {
//             active = &mut self.top;
//             passive = &mut self.bot;
//         } else {
//             active = &mut self.bot;
//             passive = &mut self.top;
//         }

//         let mut leading = 0.0;
//         let mut trailing = 0.0;
//         let mut truncation = 0.0;
//         let mut possible_new_chunk_truncation = None;
//         if let Some(iter) = self.buffer_iter.as_mut() {
//             let next = iter.next();
//             if let Some(signal) = next {
//                 leading = active.compute(signal);
//                 trailing = passive.compute(0.0);
                
//                 //TODO: Evaluate if this should run if its not reverse time as well
//                 //TODO: Evaluate impact of this on initialization when there is no trailing signal
//                 // if self.reverse_time {
//                 //     truncation = self.tf.compute(signal);
//                 // }
//                 possible_new_chunk_truncation = self.buffer_truncation.buffer_back(self.tf.compute(signal));
//             }
//         }
//         if let Some(iter) = self.buffer_truncation_iter.as_mut() {
//             let next = iter.next();
//             if let Some(signal) = next {
//                 truncation = signal;
//             }
//         }
//         let sum = leading + trailing - truncation;

//         let possible_new_chunk;
//         if self.reverse_time {
//             possible_new_chunk = self.buffer.buffer_front(signal);
//         } else {
//             possible_new_chunk = self.buffer.buffer_back(signal);
//         }
//         if let Some(new_chunk) = possible_new_chunk { // Switch to buffer_back for normal processing
//             self.buffer_iter = Some(new_chunk.into_iter());
//             self.switch = !self.switch;
//             passive.clear();
//         }
//         if let Some(new_chunk_truncation) = possible_new_chunk_truncation {
//             self.buffer_truncation_iter = Some(new_chunk_truncation.into_iter());
//             self.tf.clear();
//         }

//         sum
//     }
// }

fn synthetic_division(dividend: &[f64], divisor: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let quotient_terms = dividend.len() - divisor.len() + 1;
    let mut quotient = dividend[..quotient_terms].to_vec();
    let mut remainder = dividend[quotient_terms..].to_vec();

    let normalization_factor = divisor[0];
    for i in 0..quotient_terms {
        quotient[i] /= normalization_factor;
        for (j, val) in divisor.iter().enumerate().skip(1) { // j starts at 1
            if i+j < quotient_terms {
                quotient[i+j] += -val * quotient[i];
            } else {
                remainder[i+j-quotient_terms] += -val * quotient[i];
            }
        }
    }

    (quotient, remainder)
}

