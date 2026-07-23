use crate::common::ib_to_bf_prec;
use crate::common::{reset_prec_bits, set_prec_bits};
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::{IBig, UBig};
use rand::{rngs::StdRng, SeedableRng};
use std::str::FromStr;

/// Configuration data for Diophantine equation solving.
///
/// Contains timeout parameters and a random number generator for
/// probabilistic algorithms used in Diophantine equation solving.
#[derive(Debug)]
pub struct DiophantineData {
    pub diophantine_timeout: u128,
    pub factoring_timeout: u128,
    pub rng: StdRng,
}

/// Main configuration structure for the GridSynth algorithm.
///
/// Contains all parameters needed to run the gridsynth gate synthesis algorithm,
/// including the target rotation angle, precision requirements, timeout settings,
/// and various flags controlling the algorithm's behavior.
#[derive(Debug)]
pub struct GridSynthConfig {
    pub theta: FBig<HalfEven>,
    pub epsilon: FBig<HalfEven>,
    pub verbose: bool,
    pub measure_time: bool,
    pub diophantine_data: DiophantineData,
    pub up_to_phase: bool,
    pub compute_error: bool,
}

impl GridSynthConfig {
    /// Configures whether to compute and verify the approximation error.
    ///
    /// # Arguments
    ///
    /// * `compute_error` - If true, computes the error between the synthesized
    ///   gate sequence and the target unitary after synthesis completes
    ///
    /// # Returns
    ///
    /// A new `GridSynthConfig` with the `compute_error` flag updated
    pub fn with_compute_error(self, compute_error: bool) -> Self {
        Self {
            compute_error,
            ..self
        }
    }
}

/// The result of running the GridSynth algorithm.
///
/// Contains the synthesized gate sequence, global phase information,
/// and optional error metrics if error computation was enabled.
pub struct GridSynthResult {
    /// List of gates.
    pub gates: String,

    /// The global phase factor.
    pub global_phase: bool,

    /// If error is computed, stores the error.
    pub error: Option<f64>,

    /// If error is computed, stores whether approximation is correct.
    pub is_correct: Option<bool>,
}

/// Parses a decimal string with optional scientific notation into a rational number.
///
/// # Arguments
///
/// * `input` - A string representing a decimal number, optionally in scientific notation
///   (e.g., "3.14", "-2.5e-3", "1.23E+4")
///
/// # Returns
///
/// An `Option` containing a tuple `(numerator, denominator)` representing the rational
/// number as `numerator / denominator`. Returns `None` if the input is malformed.
///
/// # Examples
///
/// - "3.14" → (314, 100)
/// - "2.5e-3" → (25, 10000)
/// - "-1.5E2" → (-150, 1)
pub fn parse_decimal_with_exponent(input: &str) -> Option<(IBig, IBig)> {
    let input = input.trim();
    let (sign, body) = if let Some(s) = input.strip_prefix('-') {
        (-1, s)
    } else if let Some(s) = input.strip_prefix('+') {
        (1, s)
    } else {
        (1, input)
    };

    let (base_str, exp_str) = match body.split_once(['e', 'E']) {
        Some((b, e)) => (b, e),
        None => (body, "0"),
    };

    let mut parts = base_str.split('.');
    let int_part = match parts.next() {
        Some(part) => part,
        _ => return None,
    };
    let frac_part: &str = parts.next().unwrap_or_default();
    if parts.next().is_some() {
        return None;
    }
    let digits = format!("{}{}", int_part, frac_part);
    let decimal_digits = frac_part.len() as i32;

    let exponent: i32 = exp_str.parse().ok()?;
    let scale = exponent - decimal_digits;

    let mut numerator = IBig::from_str(&digits).ok()? * sign;
    let mut denominator = IBig::from(1);

    match scale.cmp(&0) {
        std::cmp::Ordering::Greater => {
            numerator *= IBig::from(10u8).pow(scale as usize);
        }
        std::cmp::Ordering::Less => {
            denominator = IBig::from(10u8).pow((-scale) as usize);
        }
        std::cmp::Ordering::Equal => {}
    }

    Some((numerator, denominator))
}

/// Creates a default GridSynth configuration from basic parameters.
///
/// This is a convenience function for calling GridSynth from other Rust packages
/// without needing to manually construct all configuration details.
///
/// # Arguments
///
/// * `theta` - The target rotation angle in radians
/// * `epsilon` - The desired precision/error tolerance
/// * `seed` - Random seed for the RNG used throughout the synthesis
/// * `verbose` - Enable verbose logging output
/// * `up_to_phase` - If true, synthesize up to a global phase factor
///
/// # Returns
///
/// A fully configured `GridSynthConfig` with:
/// - Default timeout values (200ms for Diophantine, 50ms for factoring)
/// - Automatically calculated precision based on epsilon
/// - Minimum precision of 16 bits
/// - Error computation disabled by default
///
/// # Precision Handling
///
/// The function automatically resets and recalculates the global precision
/// based on the epsilon value to ensure consistent behavior across multiple calls.
pub fn config_from_theta_epsilon(
    theta: f64,
    epsilon: f64,
    seed: u64,
    verbose: bool,
    up_to_phase: bool,
) -> GridSynthConfig {
    let (theta_num, theta_den) = parse_decimal_with_exponent(&theta.to_string()).unwrap();

    // The desired floating precision is initialized in module common.
    // It has the initial value the first time this function is called.
    // But, on subsequent calls, the precision has changed.
    // We reset it so that the precison is the same at the beginning of every synthesis.
    reset_prec_bits();
    let theta = ib_to_bf_prec(theta_num) / ib_to_bf_prec(theta_den);
    let (epsilon_num, epsilon_den) = parse_decimal_with_exponent(&epsilon.to_string()).unwrap();
    // The magic number 12 safely overapproximates the bits of precision.
    let calculated_prec_bits =
        12 * (epsilon_den.ilog(&UBig::from(10u8)) - epsilon_num.ilog(&UBig::from(10u8)));
    let prec_bits: usize = calculated_prec_bits;
    // Using precision that is too low can cause errors. For example stack overflow and sigabrt.
    // We don't actually need our target precision tied to working precision.
    let prec_bits = if prec_bits < 16 { 16 } else { prec_bits };
    set_prec_bits(prec_bits);
    let epsilon = ib_to_bf_prec(epsilon_num) / ib_to_bf_prec(epsilon_den);
    let diophantine_timeout = 200u128;
    let factoring_timeout = 50u128;
    let time = false;

    let rng: StdRng = SeedableRng::seed_from_u64(seed);
    let diophantine_data = DiophantineData {
        diophantine_timeout,
        factoring_timeout,
        rng,
    };

    GridSynthConfig {
        theta,
        epsilon,
        verbose,
        measure_time: time,
        diophantine_data,
        up_to_phase,
        compute_error: false,
    }
}
