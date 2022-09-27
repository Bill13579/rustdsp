mod dtype;
mod math;
mod prelude;
use rustdsp::filters::bw::ButterworthFilterOrder;
use num_complex::Complex;
use num_complex::ComplexFloat;
use num_traits::Float;
use num_traits::NumCast;
use num_traits::One;

use crate::dtype::*;
use crate::math::*;
use crate::prelude::*;

use std::ops::AddAssign;
use std::ops::DivAssign;
use std::ops::MulAssign;
use std::ops::RemAssign;
use std::ops::SubAssign;
use std::{iter, f64::consts::PI};

fn main() {
    
}

