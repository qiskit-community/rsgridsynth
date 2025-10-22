// Copyright (c) IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use dashu_base::RemEuclid;
use dashu_float::{round::mode::HalfEven, FBig};
use dashu_int::IBig;
use once_cell::sync::Lazy;
use std::sync::atomic::{AtomicUsize, Ordering};

pub static PREC_BITS: AtomicUsize = AtomicUsize::new(1000);

pub fn set_prec_bits(bits: usize) {
    PREC_BITS.store(bits, Ordering::Relaxed);
}

pub fn get_prec_bits() -> usize {
    PREC_BITS.load(Ordering::Relaxed)
}

fn compute_pi() -> FBig<HalfEven> {
    let pi_str = "314159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127";
    let decimals = "100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    ib_to_bf_prec(IBig::from_str_radix(pi_str, 10).unwrap())
        / ib_to_bf_prec(IBig::from_str_radix(decimals, 10).unwrap())
}

pub static PI: Lazy<FBig<HalfEven>> = Lazy::new(compute_pi);
pub static TAU: Lazy<FBig<HalfEven>> = Lazy::new(|| 2 * PI.clone());

fn reduce_to_pi_range(mut x: FBig<HalfEven>) -> FBig<HalfEven> {
    let tau = TAU.clone();
    x = x.rem_euclid(tau.clone());
    if x > PI.clone() {
        x -= tau;
    }
    x
}

pub fn cos_fbig(x: &FBig<HalfEven>) -> FBig<HalfEven> {
    let t = reduce_to_pi_range(x.clone());

    let mut term = ib_to_bf_prec(IBig::ONE);
    let mut sum = term.clone();
    let t2 = &t * &t;

    for i in 1..get_prec_bits() {
        let denom = IBig::from((2 * i - 1) * (2 * i));
        term = -term * &t2 / denom;
        sum += &term;

        if term == ib_to_bf_prec(IBig::ZERO) {
            break;
        }
    }
    sum
}

pub fn sin_fbig(x: &FBig<HalfEven>) -> FBig<HalfEven> {
    let t = reduce_to_pi_range(x.clone());

    let mut term = t.clone();
    let mut sum = term.clone();
    let t2 = &t * &t;

    for i in 1..get_prec_bits() {
        let denom = IBig::from((2 * i) * (2 * i + 1));
        term = -term * &t2 / denom;
        sum += &term;

        if term == ib_to_bf_prec(IBig::ZERO) {
            break;
        }
    }
    sum
}

pub fn ib_to_bf_prec(x: IBig) -> FBig<HalfEven> {
    FBig::from(x).with_precision(get_prec_bits()).value()
}

pub fn fb_with_prec(x: FBig<HalfEven>) -> FBig<HalfEven> {
    x.with_precision(get_prec_bits()).value()
}

#[cfg(test)]
mod tests {
    use super::*;
    use dashu_float::round::mode::HalfEven;
    use dashu_float::FBig;
    use dashu_int::ops::Abs;
    use rand::Rng;
    use std::f64::consts::PI as PI_F64;

    fn to_fbig(x: f64) -> FBig<HalfEven> {
        FBig::<HalfEven>::try_from(x)
            .unwrap()
            .with_precision(get_prec_bits())
            .value()
    }

    fn approx_eq(a: &FBig<HalfEven>, b: &FBig<HalfEven>, tol_bits: usize) -> bool {
        let diff = (a - b).abs();
        let tol = ib_to_bf_prec(IBig::ONE)
            .with_precision(get_prec_bits())
            .value()
            / FBig::from(1u64 << tol_bits)
                .with_precision(get_prec_bits())
                .value();
        diff <= tol
    }

    #[test]
    fn test_sin_fbig_random() {
        let mut rng = rand::rng();
        for _ in 0..100 {
            let x_f64 = rng.random_range(-10.0 * PI_F64..=10.0 * PI_F64);
            let x = to_fbig(x_f64);
            let expected = to_fbig(x_f64.sin());
            let result = sin_fbig(&x);
            assert!(
                approx_eq(&result, &expected, 50),
                "sin({}) = {}, expected {}, diff = {}",
                x_f64,
                result,
                expected,
                (&result - &expected).abs()
            );
        }
    }

    #[test]
    fn test_cos_fbig_random() {
        let mut rng = rand::rng();
        for _ in 0..100 {
            let x_f64 = rng.random_range(-10.0 * PI_F64..=10.0 * PI_F64);
            let x = to_fbig(x_f64);
            let expected = to_fbig(x_f64.cos());
            let result = cos_fbig(&x);
            assert!(
                approx_eq(&result, &expected, 50),
                "cos({}) = {}, expected {}, diff = {}",
                x_f64,
                result,
                expected,
                (&result - &expected).abs()
            );
        }
    }
}
