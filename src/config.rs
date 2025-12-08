use crate::common::ib_to_bf_prec;
use crate::common::{reset_prec_bits, set_prec_bits};
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::{IBig, UBig};
use rand::{rngs::StdRng, SeedableRng};
use std::str::FromStr;

#[derive(Debug)]
pub struct DiophantineData {
    pub diophantine_timeout: u128,
    pub factoring_timeout: u128,
    pub rng: StdRng,
}

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
    /// Turns on or off checking solutions at the end of the run
    pub fn with_compute_error(self, compute_error: bool) -> Self {
        Self {
            compute_error,
            ..self
        }
    }
}

/// The result of running the gridsynth algorithm
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

/// Creates the default config to easily call the code from other rust packages.
/// `seed` is used to set single RNG that is used through the call to `gridsynth`.
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
