use super::cascaded::Cascaded;
use super::bw::{Butterworth, ButterworthParameters, BiquadTransformation};
use super::gain::GainFilter;

#[derive(Debug, Clone)]
pub enum FilterType {
    Bell,
    Highpass,
    Lowpass,
    Notch,
    Bandpass,
    Highshelf,
    Lowshelf,
    Allpass,
    Gain,
}
#[derive(Debug, Clone)]
pub struct Slope(f64);
impl Slope {
    pub fn new(db_per_oct: f64) -> Slope {
        Slope(db_per_oct)
    }
    pub fn slope(&self) -> f64 {
        self.0
    }
    pub fn order(&self) -> f64 {
        // order * 6 = db_slope
        self.0 / 6.0
    }
}
#[derive(Debug, Clone)]
pub struct Universal {
    ty: FilterType,
    freq: f64, // Frequency => in Hz
    q: f64, // Q-factor => sqrt(2)/2 in base case, between 0 and 1 for percentages
    db_gain: f64, // Gain in Decibels (dB)
    slope: Slope, // Also corresponds to filter order
    filter: Cascaded<Butterworth>,
    global_gain: GainFilter,
}
impl Universal {
    pub fn generate_parameters_low_pass(q: f64) -> ButterworthParameters {
        ButterworthParameters {
            // transform: BiquadTransformation::Noop
            transform: BiquadTransformation::Q(1.0 / (2.0_f64.sqrt() * q))
        }
    }
    // pub fn generate_parameters_low_pass(q: f64) -> ButterworthParameters {
    //     ButterworthParameters {
    //         scale_zeroes: [0.0, 0.0, 1.0],
    //         base_zeros_on_poles: false,
    //         scale_zero_of_odd_binomial: [0.0, 1.0],
    //         Q: q,
    //         bandpass_normalize: false,
    //         denom_normalize: [1.0, 1.0, 1.0],
    //         scale_poles: [1.0, 1.0, 1.0],
    //         scale_pole_of_odd_binomial: [1.0, 1.0],
    //     }
    // }
    // pub fn generate_parameters_high_shelf(q: f64, A: f64, order: usize) -> ButterworthParameters {
    //     ButterworthParameters {
    //         scale_zeroes: [A.powf(1.0 / order as f64), A.powf(1.0 / (2.0 * order as f64)), 1.0],
    //         base_zeros_on_poles: true,
    //         scale_zero_of_odd_binomial: [A.powf(1.0 / (2.0 * order as f64)), 1.0],
    //         Q: q,
    //         bandpass_normalize: false,
    //         denom_normalize: [1.0, 1.0, 1.0],
    //         scale_poles: [A.powf(-1.0 / order as f64), A.powf(-1.0 / (2.0 * order as f64)), 1.0],
    //         scale_pole_of_odd_binomial: [A.powf(-1.0 / (2.0 * order as f64)), 1.0],
    //     }
    // }
    // pub fn generate_parameters_low_shelf(q: f64, A: f64, order: usize) -> ButterworthParameters {
    //     ButterworthParameters {
    //         scale_zeroes: [1.0, A.powf(1.0 / (2.0 * order as f64)), A.powf(1.0 / order as f64)],
    //         base_zeros_on_poles: true,
    //         scale_zero_of_odd_binomial: [1.0, A.powf(1.0 / (2.0 * order as f64))],
    //         Q: q,
    //         bandpass_normalize: false,
    //         denom_normalize: [1.0, 1.0, 1.0],
    //         scale_poles: [1.0, A.powf(-1.0 / (2.0 * order as f64)), A.powf(-1.0 / order as f64)],
    //         scale_pole_of_odd_binomial: [1.0, A.powf(-1.0 / (2.0 * order as f64))],
    //     }
    // }
    // pub fn generate_parameters_low_shelf2(q: f64, A: f64, order: usize) -> ButterworthParameters {
    //     ButterworthParameters {
    //         scale_zeroes: [1.0, A.powf(1.0 / order as f64), A.powf(2.0 / order as f64)],
    //         base_zeros_on_poles: true,
    //         scale_zero_of_odd_binomial: [1.0, A.powf(1.0 / order as f64)],
    //         Q: q,
    //         bandpass_normalize: false,
    //         denom_normalize: [1.0, 1.0, 1.0],
    //         scale_poles: [1.0, 1.0, 1.0],
    //         scale_pole_of_odd_binomial: [1.0, 1.0],
    //     }
    // }
    // //These work, but are just unnecessary
    // // pub fn low_shelf2_gain_at_frequency(filter_resonance_freq: f64, freq: f64, A: f64, order: usize) -> f64 {
    // //     let omega_2m = (freq / filter_resonance_freq).powf(2.0 * order as f64);
    // //     ((omega_2m + A.powi(2)) / (omega_2m + 1.0)).sqrt()
    // // }
    // // pub fn low_shelf2_gain_to_achieve_response_at_frequency(filter_resonance_freq: f64, freq: f64, mut A: f64, order: usize) -> f64 {
    // //     A = A.powi(2);
    // //     let omega_2m = (freq / filter_resonance_freq).powf(2.0 * order as f64);
    // //     (omega_2m * (A - 1.0) + A).sqrt()
    // // }
    // pub fn generate_parameters_band_shelf(q: f64, db_gain: f64, order: usize) -> (ButterworthParameters, ButterworthParameters) {
    //     let low_shelf_1 = Self::generate_parameters_low_shelf2(q, gain(db_gain), order); // Fixed Q since resonance at edges of bandpass makes little sense
    //     let low_shelf_2 = Self::generate_parameters_low_shelf2(q, gain(-db_gain), order);
    //     (low_shelf_1, low_shelf_2)
    // }
    // /// Band pass filters MUST HAVE an EVEN order
    // pub fn generate_parameters_band_pass(q: f64) -> ButterworthParameters {
    //     ButterworthParameters {
    //         scale_zeroes: [0.0, 1.0, 0.0],
    //         base_zeros_on_poles: false,
    //         scale_zero_of_odd_binomial: [1.0, 0.0],
    //         Q: q,
    //         bandpass_normalize: true,
    //         denom_normalize: [1.0, 1.0, 1.0],
    //         scale_poles: [1.0, 1.0, 1.0],
    //         scale_pole_of_odd_binomial: [1.0, 1.0],
    //     }
    // }
    // pub fn generate_parameters_band_pass2(q: f64) -> (ButterworthParameters, ButterworthParameters) {
    //     let high_pass = Self::generate_parameters_high_pass(q);
    //     let low_pass = Self::generate_parameters_low_pass(q);
    //     (low_pass, high_pass)
    // }
    // /// Notch filters MUST HAVE an EVEN order
    // pub fn generate_parameters_notch(q: f64) -> ButterworthParameters {
    //     ButterworthParameters {
    //         scale_zeroes: [1.0, 0.0, 1.0],
    //         base_zeros_on_poles: false,
    //         scale_zero_of_odd_binomial: [1.0, 0.0],
    //         Q: q,
    //         bandpass_normalize: false,
    //         denom_normalize: [1.0, 1.0, 1.0],
    //         scale_poles: [1.0, 1.0, 1.0],
    //         scale_pole_of_odd_binomial: [1.0, 1.0],
    //     }
    // }
    // pub fn generate_parameters_all_pass(q: f64) -> ButterworthParameters {
    //     ButterworthParameters {
    //         scale_zeroes: [1.0, -1.0, 1.0],
    //         base_zeros_on_poles: true,
    //         scale_zero_of_odd_binomial: [1.0, -1.0],
    //         Q: q,
    //         bandpass_normalize: false,
    //         denom_normalize: [1.0, 1.0, 1.0],
    //         scale_poles: [1.0, 1.0, 1.0],
    //         scale_pole_of_odd_binomial: [1.0, 1.0],
    //     }
    // }
    // pub fn generate_parameters_bell(q: f64, mut A: f64) -> ButterworthParameters {
    //     A = A.sqrt();
    //     ButterworthParameters {
    //         scale_zeroes: [1.0, A, 1.0],
    //         base_zeros_on_poles: true,
    //         scale_zero_of_odd_binomial: [1.0, -1.0],
    //         Q: q,
    //         bandpass_normalize: false,
    //         denom_normalize: [1.0, 1.0 / A, 1.0],
    //         scale_poles: [1.0, 1.0, 1.0],
    //         scale_pole_of_odd_binomial: [1.0, 1.0],
    //     }
    // }
    /// Q's usual definition can be found in `ButterworthParameters`, but for filters made up of two cascaded filters creating 
    /// a band-pass or a band-shelf has a different definition of Q. Here Q is 0.0~1.0 with 0 representing the filters being right 
    /// next to each other and 1 representing the filters being sqrt(2) log distance away from the resonance frequency.
    /// 
    /// This particular method converts the Q into a logarithmic distance value.
    pub fn q_definition_band__q_to_log_dist(q: f64) -> f64 {
        // sqrt(2) log distance is the max distance allowed for the filters to be away from the resonance frequency
        2.0_f64.sqrt() * q
    }
    /// # Usage
    /// 
    /// ```rust
    /// // Example of creating filters from the two frequencies
    /// if upper_filter_frequency > 0.0 && upper_filter_frequency < (1.0 / T / 2.0) {
    ///     Butterworth::new(T, upper_filter_frequency, order, params, 1.0);
    /// } else {
    ///     // Gain filter
    /// }
    /// if lower_filter_frequency > 0.0 && lower_filter_frequency < (1.0 / T / 2.0) {
    ///     Butterworth::new(T, lower_filter_frequency, order, params2, 1.0);
    /// }
    /// ```
    pub fn q_definition_band__edge_frequencies(freq: f64, q: f64) -> (f64, f64) {
        let log_dist = Self::q_definition_band__q_to_log_dist(q);

        let log_freq = freq.log10(); // Convert the frequency in Hz into the logarithmic scale
        let upper_filter_frequency = 10.0_f64.powf(log_freq + log_dist); // Add/Subtract the log distance from the log frequency and turn it back into Hz
        let lower_filter_frequency = 10.0_f64.powf(log_freq - log_dist);

        (upper_filter_frequency, lower_filter_frequency)
    }
    // pub fn new(ty: FilterType, freq: f64, q: f64, db_gain: f64, slope: SlopeType) -> Universal {

    // }
}


// fn generate_filter(T: f64, freq: f64, Q: f64) -> IIR {
//     let transform = bilinear_transform(&[0.0, 0.0, 1.0], &[1.0, 1.0 / Q, 1.0]).unwrap();
//     let f = apply_k(&transform, k_prewrap_and_scale_normalized(T, freq * 2.0 * PI));
//     IIR::new(f)
// }

// fn generate_filter2(T: f64, freq: f64, Q: f64) -> IIR {
//     let transform = bilinear_transform(&[1.0, 0.0, 0.0], &[1.0, 1.0 / Q, 1.0]).unwrap();
//     let f = apply_k(&transform, k_prewrap_and_scale_normalized(T, freq * 2.0 * PI));
//     IIR::new(f)
// }

