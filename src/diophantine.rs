// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::ring::{DOmega, DRootTwo, ZOmega, ZRootTwo};
use dashu_base::{BitTest, Gcd, RemEuclid};
use dashu_int::ops::Abs;
use dashu_int::{IBig, UBig};
use log::warn;
use once_cell::sync::Lazy;
use rand::Rng;
use std::sync::{LazyLock, Mutex};
use std::time::Instant;
use std::{
    collections::{HashMap, VecDeque},
    ops::Mul,
};

static PRIMALITY_CACHE: LazyLock<Mutex<HashMap<IBig, bool>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

type SqrtCacheType = LazyLock<Mutex<HashMap<(IBig, IBig), Option<IBig>>>>;
static SQRT_CACHE: SqrtCacheType = LazyLock::new(|| Mutex::new(HashMap::new()));

// WARN: Distribution of random_ubig is different from the rng.random_range
fn random_ubig<R>(bits: usize, rng: &mut R) -> UBig
where
    R: Rng + ?Sized,
{
    let mut bytes = vec![0u8; (bits + 7) / 8];
    rng.fill_bytes(&mut bytes);
    let mut n = UBig::from_le_bytes(&bytes);
    n |= UBig::ONE << (bits - 1);
    n &= (UBig::ONE << bits) - UBig::ONE;
    n
}

fn is_small_prime(n: u64) -> bool {
    match n {
        0 | 1 => false,
        2 => true,
        _ if n % 2 == 0 => false,
        3..=16 => [3, 5, 7, 11, 13].contains(&n),
        _ => {
            let witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
            miller_rabin_deterministic(n, &witnesses)
        }
    }
}

fn miller_rabin_deterministic(n: u64, witnesses: &[u64]) -> bool {
    let mut d = n - 1;
    let mut r = 0;
    while d % 2 == 0 {
        d /= 2;
        r += 1;
    }

    'outer: for &a in witnesses {
        if a >= n {
            continue;
        }
        let mut x = modpow_u64(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }

        for _ in 0..r - 1 {
            x = ((x as u128 * x as u128) % n as u128) as u64;
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

fn modpow_u64(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut result = 1;
    base %= modulus;
    while exp > 0 {
        if exp % 2 == 1 {
            result = ((result as u128 * base as u128) % modulus as u128) as u64;
        }
        base = ((base as u128 * base as u128) % modulus as u128) as u64;
        exp /= 2;
    }
    result
}

fn modpow(mut base: IBig, exp: &IBig, modulus: &IBig) -> IBig {
    let mut result = IBig::ONE;
    base %= modulus;
    let mut tmp = exp.clone();
    while tmp > IBig::ZERO {
        if &tmp % IBig::from(2) == IBig::ONE {
            result = (result * &base) % modulus;
        }
        base = (&base * &base) % modulus;
        tmp /= IBig::from(2);
    }
    result
}

pub fn is_prime(mut n: IBig, iterations: usize) -> bool {
    if n < IBig::ZERO {
        n = -n;
    }

    let n_key = n.clone();

    if let Ok(cache) = PRIMALITY_CACHE.lock() {
        if let Some(&result) = cache.get(&n_key) {
            return result;
        }
    }

    static TWO: Lazy<IBig> = Lazy::new(|| IBig::from(2));
    static ONE: Lazy<IBig> = Lazy::new(|| IBig::ONE);

    if n <= *ONE {
        return false;
    }
    if n == *TWO {
        return true;
    }
    if &n % 2 == 0 {
        return false;
    }

    if n.bit_len() <= 64 {
        if let Ok(n_u64) = TryInto::<u64>::try_into(&n) {
            let result = is_small_prime(n_u64);
            if let Ok(mut cache) = PRIMALITY_CACHE.lock() {
                cache.insert(n_key, result);
            }
            return result;
        }
    }

    let reduced_iterations = std::cmp::min(iterations, 4);
    let mut d = &n - 1;
    let mut r = 0;
    while &d % 2 == 0 {
        d /= 2;
        r += 1;
    }

    let mut rng = rand::rng();
    let bits = n.bit_len();
    let n_minus1 = &n - 1;

    'outer: for _ in 0..reduced_iterations {
        let a = IBig::from(random_ubig(bits, &mut rng));
        let mut x = modpow(a, &d, &n);
        if x == *ONE || x == n_minus1 {
            continue;
        }
        for _ in 0..r - 1 {
            x = (&x * &x) % &n;
            if x == n_minus1 {
                continue 'outer;
            }
        }
        if let Ok(mut cache) = PRIMALITY_CACHE.lock() {
            cache.insert(n_key, false);
        }
        return false;
    }

    if let Ok(mut cache) = PRIMALITY_CACHE.lock() {
        cache.insert(n_key, true);
    }
    true
}

pub fn sqrt_negative_one(p: &IBig, trials: usize) -> Option<IBig> {
    let cache_key = (p.clone(), IBig::NEG_ONE);
    if let Ok(cache) = SQRT_CACHE.lock() {
        if let Some(result) = cache.get(&cache_key) {
            return result.clone();
        }
    }

    let mut rng = rand::rng();
    let bits = p.bit_len();
    let reduced_trials = std::cmp::min(trials, 8); // Reduce trials for performance

    for _ in 0..reduced_trials {
        let b = IBig::from(random_ubig(bits, &mut rng));
        let exp = (p - IBig::ONE) >> 2;
        let h = modpow(b, &exp, p);
        let r = (&h * &h) % p;
        if r == p - IBig::ONE {
            if let Ok(mut cache) = SQRT_CACHE.lock() {
                cache.insert(cache_key, Some(h.clone()));
            }
            return Some(h);
        } else if r != IBig::ONE {
            if let Ok(mut cache) = SQRT_CACHE.lock() {
                cache.insert(cache_key, None);
            }
            return None;
        }
    }

    if let Ok(mut cache) = SQRT_CACHE.lock() {
        cache.insert(cache_key, None);
    }
    None
}

fn root_mod(x: IBig, p: &IBig, trials: usize) -> Option<IBig> {
    let x: IBig = x.rem_euclid(p).into();

    let cache_key = (x.clone(), p.clone());
    if let Ok(cache) = SQRT_CACHE.lock() {
        if let Some(result) = cache.get(&cache_key) {
            return result.clone();
        }
    }

    if p == &IBig::from(2) {
        return Some(x);
    }
    if x == IBig::ZERO {
        return Some(IBig::ZERO);
    }
    if p % IBig::from(2) == IBig::ZERO {
        return None;
    }
    if modpow(x.clone(), &((p - IBig::ONE) / IBig::from(2)), p) != IBig::ONE {
        if let Ok(mut cache) = SQRT_CACHE.lock() {
            cache.insert(cache_key, None);
        }
        return None;
    }

    let mut rng = rand::rng();
    let bits = p.bit_len();
    let reduced_trials = std::cmp::min(trials, 8); // Reduce trials for performance

    for _ in 0..reduced_trials {
        let b = IBig::from(random_ubig(bits, &mut rng));
        let r = modpow(b.clone(), &(p - IBig::ONE), p);
        if r != IBig::ONE {
            if let Ok(mut cache) = SQRT_CACHE.lock() {
                cache.insert(cache_key, None);
            }
            return None;
        }
        let base = (b.clone() * b.clone() + p - x.clone()) % p;
        if modpow(base.clone(), &((p - IBig::ONE) / IBig::from(2)), p) != IBig::ONE {
            let mut f = Fp2::new_with_modulus(b, IBig::ONE, base, p);
            let result = Some(f.pow((p + 1) / 2).a);
            if let Ok(mut cache) = SQRT_CACHE.lock() {
                cache.insert(cache_key, result.clone());
            }
            return result;
        }
    }

    if let Ok(mut cache) = SQRT_CACHE.lock() {
        cache.insert(cache_key, None);
    }
    None
}

#[derive(Clone)]
pub struct Fp2 {
    a: IBig,
    b: IBig,
    base: IBig,
    p: IBig,
}

impl Fp2 {
    pub fn new_with_modulus(a: IBig, b: IBig, base: IBig, p: &IBig) -> Self {
        Self {
            a: a.rem_euclid(p).into(),
            b: b.rem_euclid(p).into(),
            base,
            p: p.clone(),
        }
    }

    pub fn pow(&mut self, mut exp: IBig) -> Self {
        let mut result = Fp2::new_with_modulus(IBig::ONE, IBig::ZERO, self.base.clone(), &self.p);
        while exp > IBig::ZERO {
            if &exp & IBig::ONE == IBig::ONE {
                result = result.mul(self.clone());
            }
            *self = self.clone().mul(self.clone());
            exp >>= 1;
        }
        result
    }
}

impl Mul for Fp2 {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let p = self.p;
        let base = self.base;
        let a =
            (self.a.clone() * other.a.clone() + self.b.clone() * other.b.clone() % &p * &base) % &p;
        let b = (self.a.clone() * other.b.clone() + self.b.clone() * other.a.clone()) % &p;
        Self::new_with_modulus(a, b, base, &p)
    }
}

fn decompose_relatively_prime(mut factors: Vec<(IBig, i32)>) -> (IBig, Vec<(IBig, i32)>) {
    let mut u = IBig::ONE;
    let mut facs: Vec<(IBig, i32)> = vec![];
    let mut stack: VecDeque<(IBig, i32)> = factors.drain(..).rev().collect();

    while let Some((b, k_b)) = &stack.pop_back() {
        let mut i = 0;
        while i < facs.len() {
            let (a, k_a) = &facs[i];
            if a == b || a == &-b {
                if a == &-b && k_b & 1 == 1 {
                    u = -u;
                }
                facs[i] = (a.clone(), k_a + k_b);
                break;
            }
            let g: IBig = a.clone().gcd(b).into();
            if g == IBig::ONE || g == IBig::from(-1) {
                i += 1;
            } else {
                let (u_a, mut facs_a) =
                    decompose_relatively_prime(vec![(a / &g, *k_a), (g.clone(), k_a + k_b)]);
                u *= u_a;
                facs[i] = facs_a.remove(0);
                facs.append(&mut facs_a);
                stack.push_back((b / &g, *k_b));
                break;
            }
        }
        if i == facs.len() {
            if b == &IBig::ONE || b == &IBig::NEG_ONE {
                if b == &IBig::NEG_ONE && k_b & IBig::ONE == IBig::ONE {
                    u = -u;
                }
            } else {
                facs.push((b.clone(), *k_b));
            }
        }
    }
    (u, facs)
}

static FACTOR_CACHE: LazyLock<Mutex<HashMap<IBig, Option<IBig>>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

fn find_factor(n: &IBig, timeout_ms: u128, _m: usize) -> Option<IBig> {
    if n.bit_len() > 1024 || n == &IBig::ZERO || n == &IBig::ONE {
        return None;
    }

    if let Ok(cache) = FACTOR_CACHE.try_lock() {
        if let Some(cached_result) = cache.get(n) {
            return cached_result.clone();
        }
    }

    if n % IBig::from(2) == IBig::ZERO && n > &IBig::from(2) {
        return Some(IBig::from(2));
    }

    let small_primes = [
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    ];
    let start = Instant::now();

    for &p in &small_primes {
        if start.elapsed().as_millis() >= timeout_ms / 4 {
            break;
        }
        if n % IBig::from(p) == IBig::ZERO && n > &IBig::from(p) {
            let result = Some(IBig::from(p));
            if let Ok(mut cache) = FACTOR_CACHE.try_lock() {
                cache.insert(n.clone(), result.clone());
                if cache.len() > 5000 {
                    cache.clear();
                }
            }
            return result;
        }
    }

    if n.bit_len() <= 32 {
        let n_u64: u64 = n.try_into().unwrap_or(0);
        if n_u64 > 0 {
            let sqrt_n = std::cmp::min((n_u64 as f64).sqrt() as u64 + 1, 100000); // Limit search range
            for i in (101..=sqrt_n).step_by(2) {
                if start.elapsed().as_millis() >= timeout_ms {
                    break;
                }
                if n_u64 % i == 0 {
                    let result = Some(IBig::from(i));
                    if let Ok(mut cache) = FACTOR_CACHE.try_lock() {
                        cache.insert(n.clone(), result.clone());
                    }
                    return result;
                }
            }
        }
        if let Ok(mut cache) = FACTOR_CACHE.try_lock() {
            cache.insert(n.clone(), None);
        }
        return None;
    }

    let remaining_timeout = timeout_ms.saturating_sub(start.elapsed().as_millis());
    let result = pollard_rho(n, remaining_timeout);

    if let Ok(mut cache) = FACTOR_CACHE.try_lock() {
        cache.insert(n.clone(), result.clone());
        if cache.len() > 5000 {
            cache.clear();
        }
    }

    result
}

fn pollard_rho(n: &IBig, timeout_ms: u128) -> Option<IBig> {
    let start = Instant::now();
    let mut rng = rand::rng();
    for _ in 0..5 {
        if start.elapsed().as_millis() >= timeout_ms {
            return None;
        }

        let a = IBig::from(2 + (rng.random_range(0..100) as u64));
        let mut x = IBig::from(2);
        let mut y = IBig::from(2);
        let mut d = IBig::ONE;

        while d == IBig::ONE {
            if start.elapsed().as_millis() >= timeout_ms {
                return None;
            }
            x = (&x * &x + &a) % n;
            y = (&y * &y + &a) % n;
            y = (&y * &y + &a) % n;

            d = (&x - &y).abs().gcd(n).into();

            if d == *n {
                break;
            }
        }

        if d != IBig::ONE && d != *n {
            return Some(d);
        }
    }

    None
}

fn adj_decompose_int_prime(p: &IBig) -> Result<Option<ZOmega>, String> {
    let p = p.abs();
    if p == IBig::ZERO || p == IBig::ONE {
        return Ok(Some(ZOmega::from_int(p)));
    }
    if p == IBig::from(2) {
        return Ok(Some(ZOmega::new(
            IBig::NEG_ONE,
            IBig::ZERO,
            IBig::ONE,
            IBig::ZERO,
        )));
    }
    if is_prime(p.clone(), 4) {
        if &p & 0b11 == IBig::ONE {
            if let Some(h) = sqrt_negative_one(&p, 100) {
                let t = ZOmega::gcd(
                    h + ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
                    ZOmega::from_int(p.clone()),
                );
                let t_conj = t.conj();
                if t_conj.clone() * t.clone() == p.clone() || t_conj * &t == -p.clone() {
                    return Ok(Some(t));
                }
            }
            Ok(None)
        } else if &p & 0b111 == IBig::from(3) {
            if let Some(h) = root_mod(IBig::from(-2), &p, 100) {
                let t = ZOmega::gcd(
                    h + ZOmega::new(IBig::ONE, IBig::ZERO, IBig::ONE, IBig::ZERO),
                    ZOmega::from_int(p.clone()),
                );
                let t_conj = t.conj();
                if t_conj.clone() * t.clone() == p.clone() || t_conj * &t == -p.clone() {
                    return Ok(Some(t));
                }
            }
            Ok(None)
        } else if &p & 0b111 == IBig::from(7) {
            if root_mod(IBig::from(2), &p, 100).is_some() {
                Err("No solution in adj_decompose_int_prime".to_string())
            } else {
                Ok(None)
            }
        } else {
            Ok(None)
        }
    } else if &p & 0b111 == IBig::from(7) {
        if root_mod(IBig::from(2), &p, 100).is_some() {
            Err("No solution in adj_decompose_int_prime".to_string())
        } else {
            Ok(None)
        }
    } else {
        Ok(None)
    }
}

fn adj_decompose_int_prime_power(p: &IBig, k: i32) -> Result<Option<ZOmega>, String> {
    if k & 1 == 0 {
        Ok(Some(ZOmega::from_int(p.pow((k / 2) as usize))))
    } else {
        match adj_decompose_int_prime(p) {
            Ok(Some(t)) => Ok(Some(t.pow(k.try_into().unwrap()))),
            Ok(None) => Ok(None),
            Err(e) => Err(e),
        }
    }
}

fn adj_decompose_int(
    mut n: IBig,
    diophantine_timeout: u128,
    factoring_timeout: u128,
    start_time: Instant,
) -> Result<Option<ZOmega>, String> {
    n = n.abs();
    let mut facs = vec![(n, 1)];
    let mut t = ZOmega::from_int(IBig::ONE);

    while let Some((p, k)) = facs.pop() {
        if start_time.elapsed().as_millis() >= diophantine_timeout {
            return Ok(None);
        }

        if p.bit_len() > 2048 {
            return Ok(None);
        }

        match adj_decompose_int_prime_power(&p, k) {
            Ok(None) => {
                let individual_timeout = std::cmp::min(factoring_timeout, 20);
                if let Some(fac) = find_factor(&p, individual_timeout, 128) {
                    facs.push((p / fac.clone(), k));
                    facs.push((fac, k));
                    let (_, new_facs) = decompose_relatively_prime(facs);
                    facs = new_facs;
                } else {
                    return Ok(None);
                }
            }
            Ok(Some(t_p)) => {
                t = t * t_p;
            }
            Err(e) => {
                return Err(e);
            }
        }
    }
    Ok(Some(t))
}

fn adj_decompose_selfassociate(
    xi: ZRootTwo,
    diophantine_timeout: u128,
    factoring_timeout: u128,
    start_time: Instant,
) -> Result<Option<ZOmega>, String> {
    if xi == ZRootTwo::from_int(IBig::ZERO) {
        return Ok(Some(ZOmega::from_int(IBig::ZERO)));
    }
    let n: IBig = (&xi.a).gcd(&xi.b).into();
    let r = xi / n.clone();
    let t1 = adj_decompose_int(n, diophantine_timeout, factoring_timeout, start_time);
    let t2 = if r % ZRootTwo::new(IBig::ZERO, IBig::ONE) == ZRootTwo::from_int(IBig::ZERO) {
        ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ONE, IBig::ONE)
    } else {
        ZOmega::from_int(IBig::ONE)
    };
    match t1 {
        Ok(Some(t1_val)) => Ok(Some(t1_val * t2)),
        Ok(None) => Ok(None),
        Err(e) => Err(e),
    }
}

pub fn decompose_relatively_zomega_prime(
    partial_facs: Vec<(ZRootTwo, i32)>,
) -> (ZRootTwo, Vec<(ZRootTwo, i32)>) {
    let mut u = ZRootTwo::from_int(IBig::ONE);
    let mut stack = partial_facs.clone();
    stack.reverse();
    let mut facs = Vec::new();

    while let Some((b, k_b)) = stack.pop() {
        let mut i = 0;
        loop {
            if i >= facs.len() {
                if ZRootTwo::sim(b.clone(), ZRootTwo::from_int(IBig::ONE)) {
                    u = u * b.pow(&IBig::from(k_b));
                } else {
                    facs.push((b.clone(), k_b));
                }
                break;
            }
            let (a, k_a) = &facs[i];
            if ZRootTwo::sim(a.clone(), b.clone()) {
                u = u * (b.clone() / a.clone()).pow(&IBig::from(k_b));
                facs[i] = (a.clone(), k_a + k_b);
                break;
            } else {
                let g = ZRootTwo::gcd(a.clone(), b.clone());
                if ZRootTwo::sim(g.clone(), ZRootTwo::from_int(IBig::ONE)) {
                    i += 1;
                    continue;
                } else {
                    let partial = vec![(a.clone() / g.clone(), *k_a), (g.clone(), k_a + k_b)];
                    let (u_a, mut facs_a) = decompose_relatively_zomega_prime(partial);
                    u = u * u_a;
                    facs[i] = facs_a.remove(0);
                    facs.append(&mut facs_a);
                    stack.push((b / g, k_b));
                    break;
                }
            }
        }
    }
    (u, facs)
}

pub fn adj_decompose_zomega_prime(eta: ZRootTwo) -> Result<Option<ZOmega>, String> {
    let mut p = eta.norm();
    if p < IBig::ZERO {
        p = -p;
    }
    if p == IBig::ZERO || p == IBig::ONE {
        return Ok(Some(ZOmega::from_int(p)));
    } else if p == IBig::from(2) {
        return Ok(Some(ZOmega::new(
            IBig::NEG_ONE,
            IBig::ZERO,
            IBig::ONE,
            IBig::ZERO,
        )));
    }

    if is_prime(p.clone(), 4) {
        let check: i32 = (&p & IBig::from(0b111)).try_into().unwrap();
        match check {
            1 => Ok(sqrt_negative_one(&p, 100).and_then(|h| {
                let t = ZOmega::gcd(
                    h + ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
                    ZOmega::from_zroottwo(&eta),
                );
                let a = t.clone().conj() * &t;
                let a_zrt = ZRootTwo::from_zomega(a);
                if a_zrt.clone() % eta.clone() == ZRootTwo::from_int(IBig::ZERO)
                    && eta % a_zrt == ZRootTwo::from_int(IBig::ZERO)
                {
                    Some(t)
                } else {
                    None
                }
            })),
            3 => Ok(root_mod(IBig::from(-2), &p, 100).and_then(|h| {
                let t = ZOmega::gcd(
                    h + ZOmega::new(IBig::ONE, IBig::ZERO, IBig::ONE, IBig::ZERO),
                    ZOmega::from_zroottwo(&eta),
                );
                let a = t.clone().conj() * &t;
                let a_zrt = ZRootTwo::from_zomega(a);
                if a_zrt.clone() % eta.clone() == ZRootTwo::from_int(IBig::ZERO)
                    && eta % a_zrt == ZRootTwo::from_int(IBig::ZERO)
                {
                    Some(t)
                } else {
                    None
                }
            })),
            7 => {
                if root_mod(IBig::from(2), &p, 100).is_some() {
                    Err("No solution in adj_decompose_zomega_prime".to_string())
                } else {
                    Ok(None)
                }
            }
            _ => Ok(None),
        }
    } else if &p & 0b111 == IBig::from(7) {
        if root_mod(IBig::from(2), &p, 100).is_some() {
            Err("No solution in adj_decompose_zomega_prime".to_string())
        } else {
            Ok(None)
        }
    } else {
        Ok(None)
    }
}

pub fn adj_decompose_zomega_prime_power(eta: ZRootTwo, k: i32) -> Result<Option<ZOmega>, String> {
    if k & 1 == 0 {
        Ok(Some(ZOmega::from_zroottwo(&eta.pow(&IBig::from(k / 2)))))
    } else {
        match adj_decompose_zomega_prime(eta) {
            Ok(Some(t)) => Ok(Some(t.pow(k.try_into().unwrap()))),
            Ok(None) => Ok(None),
            Err(e) => Err(e),
        }
    }
}

pub fn adj_decompose_selfcoprime(
    xi: ZRootTwo,
    diophantine_timeout: u128,
    factoring_timeout: u128,
    start_time: Instant,
) -> Result<Option<ZOmega>, String> {
    let mut facs: Vec<(ZRootTwo, i32)> = vec![(xi.clone(), 1)];
    let mut t = ZOmega::from_int(IBig::ONE);
    while let Some((eta, k)) = facs.pop() {
        if start_time.elapsed().as_millis() >= diophantine_timeout {
            return Ok(None);
        }

        match adj_decompose_zomega_prime_power(eta.clone(), k) {
            Ok(None) => {
                let mut n = eta.norm();
                if n < IBig::ZERO {
                    n = -n;
                }
                let individual_timeout = std::cmp::min(factoring_timeout, 15); // Max 15ms per factor
                if let Some(fac_n) = find_factor(&n, individual_timeout, 128) {
                    let fac = ZRootTwo::gcd(xi.clone(), ZRootTwo::from_int(fac_n));
                    facs.push((eta / fac.clone(), k));
                    facs.push((fac, k));
                    let (_, new_facs) = decompose_relatively_zomega_prime(facs);
                    facs = new_facs;
                } else {
                    return Ok(None);
                }
            }
            Ok(Some(t_eta)) => {
                t = t * t_eta;
            }
            Err(e) => {
                return Err(e);
            }
        }
    }
    Ok(Some(t))
}

fn adj_decompose(
    xi: ZRootTwo,
    diophantine_timeout: u128,
    factoring_timeout: u128,
    start_time: Instant,
) -> Result<Option<ZOmega>, String> {
    if xi == ZRootTwo::from_int(IBig::ZERO) {
        return Ok(Some(ZOmega::from_int(IBig::ZERO)));
    }
    let d = ZRootTwo::gcd(xi.clone(), xi.conj_sq2());
    let eta = xi / d.clone();
    let t1 = adj_decompose_selfassociate(d, diophantine_timeout, factoring_timeout, start_time);
    match t1 {
        Ok(Some(t1_val)) => {
            let t2 =
                adj_decompose_selfcoprime(eta, diophantine_timeout, factoring_timeout, start_time);
            match t2 {
                Ok(Some(t2_val)) => Ok(Some(t1_val * t2_val)),
                Ok(None) => Ok(None),
                Err(e) => Err(e),
            }
        }
        Ok(None) => Ok(None),
        Err(e) => Err(e),
    }
}
fn diophantine(
    xi: &ZRootTwo,
    diophantine_timeout: u128,
    factoring_timeout: u128,
    start_time: Instant,
) -> Result<Option<ZOmega>, String> {
    if xi == &ZRootTwo::from_int(IBig::ZERO) {
        return Ok(Some(ZOmega::from_int(IBig::ZERO)));
    } else if xi < &ZRootTwo::from_int(IBig::ZERO)
        || xi.clone().conj_sq2() < ZRootTwo::from_int(IBig::ZERO)
    {
        return Err("No solution".to_string());
    }
    let t = adj_decompose(
        xi.clone(),
        diophantine_timeout,
        factoring_timeout,
        start_time,
    );
    match t {
        Ok(Some(t)) => {
            let xi_associate = ZRootTwo::from_zomega(t.conj() * &t);
            let u = xi.clone() / xi_associate;
            match u.sqrt() {
                Some(v) => Ok(Some(ZOmega::from_zroottwo(&v) * t)),
                None => {
                    warn!("Cannot find square root of u");
                    Ok(None)
                }
            }
        }
        Ok(None) => Ok(None),
        Err(e) => Err(e),
    }
}

type DiophantineCacheType = LazyLock<Mutex<HashMap<(IBig, IBig), Option<DOmega>>>>;
static DIOPHANTINE_CACHE: DiophantineCacheType = LazyLock::new(|| Mutex::new(HashMap::new()));

pub fn diophantine_dyadic(
    xi: DRootTwo,
    diophantine_timeout: u128,
    factoring_timeout: u128,
) -> Option<DOmega> {
    let cache_key = (xi.alpha.a.clone(), xi.alpha.b.clone());
    if let Ok(cache) = DIOPHANTINE_CACHE.try_lock() {
        if let Some(cached_result) = cache.get(&cache_key) {
            return cached_result.clone();
        }
    }

    let k_div_2 = xi.k >> 1;
    let k_mod_2 = xi.k & 1;

    let alpha = if k_mod_2 == 1 {
        xi.alpha * ZRootTwo::new(IBig::ONE, IBig::ONE)
    } else {
        xi.alpha
    };
    let alpha_bits = std::cmp::max(alpha.a.bit_len(), alpha.b.bit_len());
    let (optimized_diophantine_timeout, optimized_factoring_timeout) = if alpha_bits > 1000 {
        (
            std::cmp::min(diophantine_timeout, 50),
            std::cmp::min(factoring_timeout, 10),
        )
    } else {
        (
            std::cmp::min(diophantine_timeout, 15),
            std::cmp::min(factoring_timeout, 3),
        )
    };

    let start_time = Instant::now();
    let t = diophantine(
        &alpha,
        optimized_diophantine_timeout,
        optimized_factoring_timeout,
        start_time,
    );
    let result = match t {
        Err(_) => None,
        Ok(None) => None,
        Ok(Some(mut t)) => {
            if k_mod_2 == 1 {
                t = t * ZOmega::new(IBig::ZERO, IBig::NEG_ONE, IBig::ONE, IBig::ZERO);
            }
            Some(DOmega::new(t, k_div_2 + k_mod_2))
        }
    };

    if let Ok(mut cache) = DIOPHANTINE_CACHE.try_lock() {
        cache.insert(cache_key, result.clone());
        if cache.len() > 10000 {
            cache.clear();
        }
    }

    result
}
