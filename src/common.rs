// Copyright (c) IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use dashu_base::RemEuclid;
use dashu_float::{round::mode::HalfEven, FBig};
use dashu_int::IBig;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

const PREC_BITS_INITIAL: usize = 1000;
pub static PREC_BITS: AtomicUsize = AtomicUsize::new(PREC_BITS_INITIAL);

// Reset precision to the initial value
pub fn reset_prec_bits() {
    PREC_BITS.store(PREC_BITS_INITIAL, Ordering::Relaxed);
}

pub fn set_prec_bits(bits: usize) {
    PREC_BITS.store(bits, Ordering::Relaxed);
}

pub fn get_prec_bits() -> usize {
    PREC_BITS.load(Ordering::Relaxed)
}

fn compute_pi_at_current_prec() -> FBig<HalfEven> {
    let pi_str = "314159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127";
    let decimals = "100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    ib_to_bf_prec(IBig::from_str_radix(pi_str, 10).unwrap())
        / ib_to_bf_prec(IBig::from_str_radix(decimals, 10).unwrap())
}

struct PiCache {
    bits: usize,
    value: FBig<HalfEven>,
}

static PI_CACHE: Mutex<Option<PiCache>> = Mutex::new(None);

/// `pi` at the crate's current working precision (`get_prec_bits()`), cached but refreshed
/// whenever the requested precision no longer matches what's cached.
///
/// `PREC_BITS` changes across calls into this crate -- e.g. every `config_from_theta_epsilon`
/// call resizes it to fit the requested epsilon -- so a value computed and cached exactly
/// once (as a `Lazy`/`OnceCell` constant would do) stays pinned at whatever precision was
/// active on its *first* use for the rest of the process's lifetime. Since `cos_fbig`/
/// `sin_fbig` (via `reduce_to_pi_range`) are on this crate's hot path, that silently degrades
/// every later, higher-precision call in the same process to the first call's precision --
/// exactly the kind of silent-precision-loss failure mode this crate's own design already
/// worries about elsewhere (see the module-level docs on `PREC_BITS`), just from an
/// unexpected source: a process that synthesizes at a small epsilon early on and a much
/// tighter one later (e.g. a batch of rotations at increasing accuracy) would get silently
/// wrong `cos_fbig`/`sin_fbig` results for every call after the first.
pub fn pi() -> FBig<HalfEven> {
    let bits = get_prec_bits();
    {
        let cache = PI_CACHE.lock().unwrap();
        if let Some(c) = cache.as_ref() {
            if c.bits == bits {
                return c.value.clone();
            }
        }
    }
    let value = compute_pi_at_current_prec();
    *PI_CACHE.lock().unwrap() = Some(PiCache {
        bits,
        value: value.clone(),
    });
    value
}

pub fn tau() -> FBig<HalfEven> {
    2 * pi()
}

fn reduce_to_pi_range(mut x: FBig<HalfEven>) -> FBig<HalfEven> {
    let tau = tau();
    x = x.rem_euclid(tau.clone());
    if x > pi() {
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
    use serial_test::serial;
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
    #[serial]
    fn test_sin_fbig_random() {
        // This test asserts against a fixed bit-tolerance regardless of the crate's current
        // global working precision (`PREC_BITS`), so it must not run concurrently with (via
        // `#[serial]`) or inherit a stale reduced precision left behind by (via
        // `reset_prec_bits()`) another test that mutates it -- e.g. a low-epsilon synthesis
        // test can drive precision down to as little as 16 bits, which would make this test
        // fail spuriously despite `sin_fbig` being correct at that precision.
        //
        // The tolerance itself must also stay safely below the *reference*'s own precision:
        // `expected` is `x_f64.sin()`, an f64 (~53 bits of mantissa) computed by the
        // platform's libm, which is typically only guaranteed correctly-rounded to within
        // ~1 ULP, not exactly-rounded -- so comparing against it at 50 bits (only ~3 bits of
        // margin below f64's own precision) occasionally fails from ordinary libm rounding
        // noise over enough random trials, with no connection to PREC_BITS at all. 40 bits
        // leaves a comfortable margin while still being a far stricter tolerance (~1 part in
        // 10^12) than any practical use of `sin_fbig` needs.
        reset_prec_bits();
        let mut rng = rand::rng();
        for _ in 0..100 {
            let x_f64 = rng.random_range(-10.0 * PI_F64..=10.0 * PI_F64);
            let x = to_fbig(x_f64);
            let expected = to_fbig(x_f64.sin());
            let result = sin_fbig(&x);
            assert!(
                approx_eq(&result, &expected, 40),
                "sin({}) = {}, expected {}, diff = {}",
                x_f64,
                result,
                expected,
                (&result - &expected).abs()
            );
        }
    }

    #[test]
    #[serial]
    fn test_cos_fbig_random() {
        // See test_sin_fbig_random's comment: same reset/serialization and tolerance-margin
        // reasoning applies here.
        reset_prec_bits();
        let mut rng = rand::rng();
        for _ in 0..100 {
            let x_f64 = rng.random_range(-10.0 * PI_F64..=10.0 * PI_F64);
            let x = to_fbig(x_f64);
            let expected = to_fbig(x_f64.cos());
            let result = cos_fbig(&x);
            assert!(
                approx_eq(&result, &expected, 40),
                "cos({}) = {}, expected {}, diff = {}",
                x_f64,
                result,
                expected,
                (&result - &expected).abs()
            );
        }
    }
}
