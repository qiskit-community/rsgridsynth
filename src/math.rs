// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::common::{fb_with_prec, get_prec_bits, ib_to_bf_prec};
use dashu_float::round::mode;
use dashu_float::Context;
use dashu_float::{round::mode::HalfEven, FBig};
use dashu_int::ops::SquareRoot;
use dashu_int::IBig;
use std::cmp::Ordering;

/// Computes √2 with the current global precision.
///
/// # Returns
///
/// A high-precision floating-point representation of √2.
pub fn sqrt2() -> FBig<HalfEven> {
    let ctx: Context<mode::HalfEven> = Context::<mode::HalfEven>::new(get_prec_bits());
    let x = ib_to_bf_prec(IBig::from(2));
    let a: FBig<HalfEven> = ctx.sqrt(x.repr()).value();
    a
}

/// Computes the square root of a high-precision floating-point number.
///
/// # Arguments
///
/// * `x` - The number to take the square root of
///
/// # Returns
///
/// The square root of `x` with the current global precision.
///
/// # Note
///
/// This function may allocate memory for precision adjustment.
pub fn sqrt_fbig(x: &FBig<HalfEven>) -> FBig<HalfEven> {
    let ctx: Context<mode::HalfEven> = Context::<mode::HalfEven>::new(get_prec_bits());
    let x = fb_with_prec(x.clone());
    let sx: FBig<HalfEven> = ctx.sqrt(x.repr()).value();
    sx
}

/// Counts the number of trailing zeros in the binary representation of an integer.
///
/// # Arguments
///
/// * `n` - The integer to analyze
///
/// # Returns
///
/// The number of trailing zeros, or 0 if `n` is zero.
pub fn ntz(n: &IBig) -> i64 {
    match n.trailing_zeros() {
        Some(k) => k as i64,
        None => 0,
    }
}

/// Computes the natural logarithm of a high-precision floating-point number.
///
/// # Arguments
///
/// * `x` - The number to take the logarithm of
///
/// # Returns
///
/// ln(x) with the current global precision.
pub fn log(x: FBig<HalfEven>) -> FBig<HalfEven> {
    let ctx: Context<mode::HalfEven> = Context::<mode::HalfEven>::new(get_prec_bits());
    ctx.ln(x.repr()).value()
}

/// Returns the sign of a floating-point number.
///
/// # Arguments
///
/// * `x` - The number to check
///
/// # Returns
///
/// - `1` if x > 0
/// - `-1` if x < 0
/// - `0` if x = 0 or comparison fails
pub fn sign(x: FBig<HalfEven>) -> i8 {
    match x.partial_cmp(&ib_to_bf_prec(IBig::ZERO)) {
        Some(Ordering::Greater) => 1,
        Some(Ordering::Less) => -1,
        _ => 0,
    }
}

/// Computes the floor of the square root of a floating-point number.
///
/// # Arguments
///
/// * `x` - A non-negative floating-point number
///
/// # Returns
///
/// The largest integer n such that n² ≤ x.
///
/// # Panics
///
/// Panics if `x` is negative.
pub fn floorsqrt(x: FBig<HalfEven>) -> IBig {
    assert!(
        x >= ib_to_bf_prec(IBig::ZERO),
        "Negative input to floorsqrt"
    );
    binary_search_sqrt(x)
}

/// Computes the integer square root using binary search.
///
/// # Arguments
///
/// * `x` - A non-negative floating-point number
///
/// # Returns
///
/// The largest integer n such that n² ≤ x.
fn binary_search_sqrt(x: FBig<HalfEven>) -> IBig {
    let mut ok = IBig::ZERO;
    let mut ng: IBig =
        <FBig<HalfEven> as TryInto<IBig>>::try_into(x.ceil()).unwrap() + IBig::from(1);
    while &ng - &ok > IBig::ONE {
        let mid = &ok + (&ng - &ok) / IBig::from(2);
        let mid_sqr = ib_to_bf_prec(&mid * &mid);
        if mid_sqr <= x {
            ok = mid;
        } else {
            ng = mid;
        }
    }
    ok
}

/// Performs integer division with rounding to nearest.
///
/// # Arguments
///
/// * `x` - The dividend
/// * `y` - The divisor
///
/// # Returns
///
/// The quotient x/y rounded to the nearest integer.
///
/// # Rounding
///
/// Rounds towards positive infinity for positive quotients,
/// towards negative infinity for negative quotients.
pub fn rounddiv(x: IBig, y: &IBig) -> IBig {
    let adjustment = y / 2;
    if (x >= IBig::ZERO && y > &IBig::ZERO) || (x <= IBig::ZERO && y < &IBig::ZERO) {
        (x + adjustment) / y
    } else {
        (x - adjustment) / y
    }
}

/// Computes (√2)^k efficiently.
///
/// # Arguments
///
/// * `k` - The exponent
///
/// # Returns
///
/// (√2)^k computed as 2^(k/2) * √2^(k mod 2).
pub fn pow_sqrt2(k: i64) -> FBig<HalfEven> {
    let k_div_2 = k >> 1;
    let k_mod_2 = k & 1;
    let base =
        ib_to_bf_prec(IBig::ONE) << k_div_2.try_into().expect("Shift amount must fit in usize");
    if k_mod_2 != 0 {
        base * sqrt2()
    } else {
        base
    }
}

/// Computes the floor of log_y(x) and the remainder.
///
/// # Arguments
///
/// * `x` - The value (must be positive)
/// * `y` - The base (must be positive)
///
/// # Returns
///
/// A tuple `(n, r)` where:
/// - `n` is floor(log_y(x))
/// - `r` is the remainder such that x = y^n * r
///
/// # Panics
///
/// Panics if `x` is not positive.
pub fn floorlog(x: FBig<HalfEven>, y: FBig<HalfEven>) -> (IBig, FBig<HalfEven>) {
    assert!(x > ib_to_bf_prec(IBig::ZERO), "math domain error");
    let m = compute_precision_requirement(&x, &y);
    let pow_y = compute_powers(&y, m);
    compute_logarithm(x, y, pow_y)
}

/// Computes the precision requirement for logarithm calculation.
///
/// # Arguments
///
/// * `x` - The value
/// * `y` - The base
///
/// # Returns
///
/// The number of squaring operations needed for the algorithm.
fn compute_precision_requirement(x: &FBig<HalfEven>, y: &FBig<HalfEven>) -> usize {
    let mut tmp = y.clone();
    let mut m = 0;
    while x >= &tmp || x * &tmp < ib_to_bf_prec(IBig::ONE) {
        tmp = &tmp * &tmp;
        m += 1;
    }
    m
}

/// Computes successive squares of y: [y^(2^(m-1)), y^(2^(m-2)), ..., y^2, y].
///
/// # Arguments
///
/// * `y` - The base value
/// * `m` - The number of powers to compute
///
/// # Returns
///
/// A vector of successive squares in descending order of exponent.
fn compute_powers(y: &FBig<HalfEven>, m: usize) -> Vec<FBig<HalfEven>> {
    let mut pow_y = Vec::with_capacity(m);
    pow_y.push(y.clone());
    for i in 1..m {
        let last = &pow_y[i - 1];
        pow_y.push(last * last);
    }
    pow_y.reverse();
    pow_y
}

/// Computes the logarithm using precomputed powers.
///
/// # Arguments
///
/// * `x` - The value
/// * `y` - The base
/// * `pow_y` - Precomputed successive squares of y
///
/// # Returns
///
/// A tuple `(n, r)` where n is floor(log_y(x)) and r is the remainder.
fn compute_logarithm(
    x: FBig<HalfEven>,
    y: FBig<HalfEven>,
    pow_y: Vec<FBig<HalfEven>>,
) -> (IBig, FBig<HalfEven>) {
    let (mut n, mut r) = if x >= ib_to_bf_prec(IBig::ONE) {
        (IBig::ZERO, x)
    } else {
        (IBig::NEG_ONE, x * pow_y.iter().fold(y, |acc, p| acc * p))
    };

    for p in pow_y {
        n <<= 1;
        if r > p {
            r /= p;
            n += 1;
        }
    }
    (n, r)
}

/// Solves the quadratic equation ax^2 + bx + c = 0 for real roots.
///
/// # Arguments
///
/// - `a`: Coefficient of x^2, assumed to be non-zero for a valid quadratic equation.
/// - `b`: Coefficient of x.
/// - `c`: Constant term.
///
/// # Returns
///
/// An `Option` containing a tuple of two roots if they exist; otherwise, `None` if the roots are not real (i.e., the discriminant is negative).
pub fn solve_quadratic(
    a: &FBig<HalfEven>,
    b: &FBig<HalfEven>,
    c: &FBig<HalfEven>,
) -> Option<(FBig<HalfEven>, FBig<HalfEven>)> {
    let zero = ib_to_bf_prec(IBig::ZERO);
    let two = ib_to_bf_prec(IBig::from(2));
    let four = ib_to_bf_prec(IBig::from(4));

    // Handle degenerate and easy cases
    if a == &zero {
        // Linear: b x + c = 0
        if b == &zero {
            return None;
        }
        let x = -c.clone() / b.clone();
        return Some((x.clone(), x));
    }
    if c == &zero {
        // Roots are 0 and -b/a; avoids 0/0 in the stable branch
        return Some((zero.clone(), -b.clone() / a.clone()));
    }

    // Discriminant
    let disc = b.clone() * b.clone() - four.clone() * a.clone() * c.clone();
    if disc < zero {
        return None;
    }
    let sqrt_d = disc.sqrt();

    // Stable computation:
    // s = -b - sign(b) * sqrt(d)
    // x1 = s / (2a), x2 = (2c) / s
    let s = if b >= &zero {
        -b.clone() - sqrt_d.clone()
    } else {
        -b.clone() + sqrt_d.clone()
    };

    // s == 0 can only happen when c == 0 (handled above), so safe here
    let x1 = s.clone() / (two.clone() * a.clone());
    let x2 = (two * c.clone()) / s;
    if b >= &zero {
        Some((x1, x2))
    } else {
        Some((x2, x1))
    }
}
