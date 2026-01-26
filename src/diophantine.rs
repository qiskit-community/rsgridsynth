// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::config::DiophantineData;
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

/// Generates a random unsigned big integer with a specified number of bits.
///
/// # Arguments
///
/// * `bits` - The number of bits for the generated number
/// * `rng` - A random number generator
///
/// # Returns
///
/// A random `UBig` with the specified bit length, with the most significant bit set to 1.
///
/// # Warning
///
/// The distribution of this function differs from `rng.random_range`.
fn random_ubig<R>(bits: usize, rng: &mut R) -> UBig
where
    R: Rng + ?Sized,
{
    let mut bytes = vec![0u8; bits.div_ceil(8)];
    rng.fill_bytes(&mut bytes);
    let mut n = UBig::from_le_bytes(&bytes);
    n |= UBig::ONE << (bits - 1);
    n &= (UBig::ONE << bits) - UBig::ONE;
    n
}

/// Checks if a small 64-bit number is prime using deterministic Miller-Rabin test.
///
/// # Arguments
///
/// * `n` - The number to test for primality
///
/// # Returns
///
/// `true` if `n` is prime, `false` otherwise.
///
/// # Implementation
///
/// Uses a deterministic Miller-Rabin test with a fixed set of witnesses
/// that guarantees correctness for all 64-bit integers.
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

/// Performs a deterministic Miller-Rabin primality test.
///
/// # Arguments
///
/// * `n` - The number to test for primality
/// * `witnesses` - Array of witness values for the test
///
/// # Returns
///
/// `true` if `n` passes all witness tests (likely prime), `false` if composite.
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

/// Computes modular exponentiation for 64-bit integers: (base^exp) mod modulus.
///
/// # Arguments
///
/// * `base` - The base value
/// * `exp` - The exponent
/// * `modulus` - The modulus
///
/// # Returns
///
/// The result of (base^exp) mod modulus.
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

/// Computes modular exponentiation for arbitrary-precision integers: (base^exp) mod modulus.
///
/// # Arguments
///
/// * `base` - The base value
/// * `exp` - The exponent
/// * `modulus` - The modulus
///
/// # Returns
///
/// The result of (base^exp) mod modulus using binary exponentiation.
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

/// Tests whether an arbitrary-precision integer is prime using probabilistic Miller-Rabin.
///
/// # Arguments
///
/// * `n` - The number to test for primality
/// * `iterations` - Maximum number of Miller-Rabin iterations (capped at 4 for performance)
/// * `rng` - Random number generator for selecting test witnesses
///
/// # Returns
///
/// `true` if `n` is likely prime, `false` if definitely composite.
///
/// # Caching
///
/// Results are cached in `PRIMALITY_CACHE` for performance. The function is thread-safe.
///
/// # Performance
///
/// For numbers ≤ 64 bits, uses deterministic testing. For larger numbers,
/// uses probabilistic Miller-Rabin with reduced iterations for efficiency.
pub fn is_prime<R: Rng + ?Sized>(mut n: IBig, iterations: usize, rng: &mut R) -> bool {
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

    let bits = n.bit_len();
    let n_minus1 = &n - 1;

    'outer: for _ in 0..reduced_iterations {
        let a = IBig::from(random_ubig(bits, rng));
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

/// Finds a square root of -1 modulo a prime p (if it exists).
///
/// # Arguments
///
/// * `p` - A prime modulus
/// * `trials` - Maximum number of random trials (capped at 8 for performance)
/// * `rng` - Random number generator
///
/// # Returns
///
/// `Some(h)` where h² ≡ -1 (mod p), or `None` if no solution exists or trials exhausted.
///
/// # Caching
///
/// Results are cached in `SQRT_CACHE` for performance.
pub fn sqrt_negative_one<R: Rng + ?Sized>(p: &IBig, trials: usize, rng: &mut R) -> Option<IBig> {
    let cache_key = (p.clone(), IBig::NEG_ONE);
    if let Ok(cache) = SQRT_CACHE.lock() {
        if let Some(result) = cache.get(&cache_key) {
            return result.clone();
        }
    }

    let bits = p.bit_len();
    let reduced_trials = std::cmp::min(trials, 8); // Reduce trials for performance

    for _ in 0..reduced_trials {
        let b = IBig::from(random_ubig(bits, rng));
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

/// Computes a square root of x modulo prime p using Tonelli-Shanks-like algorithm.
///
/// # Arguments
///
/// * `x` - The value to find the square root of
/// * `p` - A prime modulus
/// * `trials` - Maximum number of random trials (capped at 8 for performance)
/// * `rng` - Random number generator
///
/// # Returns
///
/// `Some(r)` where r² ≡ x (mod p), or `None` if x is not a quadratic residue mod p.
///
/// # Caching
///
/// Results are cached in `SQRT_CACHE` for performance.
fn root_mod<R: Rng + ?Sized>(x: IBig, p: &IBig, trials: usize, rng: &mut R) -> Option<IBig> {
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

    let bits = p.bit_len();
    let reduced_trials = std::cmp::min(trials, 8); // Reduce trials for performance

    for _ in 0..reduced_trials {
        let b = IBig::from(random_ubig(bits, rng));
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

/// Represents an element in the field extension F_p² = F_p[x]/(x² - base).
///
/// Elements are of the form a + b*√base where a, b ∈ F_p.
#[derive(Clone)]
pub struct Fp2 {
    a: IBig,
    b: IBig,
    base: IBig,
    p: IBig,
}

impl Fp2 {
    /// Creates a new element in F_p² with specified modulus.
    ///
    /// # Arguments
    ///
    /// * `a` - The rational part coefficient
    /// * `b` - The irrational part coefficient (multiplies √base)
    /// * `base` - The base for the field extension (x² = base)
    /// * `p` - The prime modulus
    ///
    /// # Returns
    ///
    /// A new `Fp2` element representing a + b*√base in F_p.
    pub fn new_with_modulus(a: IBig, b: IBig, base: IBig, p: &IBig) -> Self {
        Self {
            a: a.rem_euclid(p).into(),
            b: b.rem_euclid(p).into(),
            base,
            p: p.clone(),
        }
    }

    /// Computes self raised to the power of exp using binary exponentiation.
    ///
    /// # Arguments
    ///
    /// * `exp` - The exponent
    ///
    /// # Returns
    ///
    /// A new `Fp2` element representing self^exp.
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

/// Decomposes a factorization into relatively prime factors.
///
/// # Arguments
///
/// * `factors` - A vector of (base, exponent) pairs representing a factorization
///
/// # Returns
///
/// A tuple `(u, facs)` where:
/// - `u` is a unit factor
/// - `facs` is a vector of relatively prime factors with their exponents
///
/// # Algorithm
///
/// Processes factors to ensure all bases in the result are pairwise coprime,
/// splitting factors with common divisors as needed.
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

/// Attempts to find a non-trivial factor of n using trial division and Pollard's rho.
///
/// # Arguments
///
/// * `n` - The number to factor
/// * `timeout_ms` - Maximum time in milliseconds to spend factoring
/// * `_m` - Unused parameter (kept for API compatibility)
/// * `rng` - Random number generator for probabilistic algorithms
///
/// # Returns
///
/// `Some(factor)` if a non-trivial factor is found, `None` otherwise.
///
/// # Caching
///
/// Results are cached in `FACTOR_CACHE` for performance.
///
/// # Algorithm
///
/// 1. Checks for small prime factors using trial division
/// 2. For numbers ≤ 32 bits, uses optimized trial division
/// 3. For larger numbers, uses Pollard's rho algorithm
fn find_factor<R: Rng + ?Sized>(
    n: &IBig,
    timeout_ms: u128,
    _m: usize,
    rng: &mut R,
) -> Option<IBig> {
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
    let result = pollard_rho(n, remaining_timeout, rng);

    if let Ok(mut cache) = FACTOR_CACHE.try_lock() {
        cache.insert(n.clone(), result.clone());
        if cache.len() > 5000 {
            cache.clear();
        }
    }

    result
}

/// Pollard's rho algorithm for integer factorization.
///
/// # Arguments
///
/// * `n` - The number to factor
/// * `timeout_ms` - Maximum time in milliseconds to spend
/// * `rng` - Random number generator for selecting starting values
///
/// # Returns
///
/// `Some(factor)` if a non-trivial factor is found, `None` otherwise.
///
/// # Algorithm
///
/// Uses the cycle-finding algorithm with random polynomial functions
/// to find factors probabilistically.
fn pollard_rho<R: Rng + ?Sized>(n: &IBig, timeout_ms: u128, rng: &mut R) -> Option<IBig> {
    let start = Instant::now();
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

/// Decomposes a prime integer into Z[ω] if possible.
///
/// # Arguments
///
/// * `p` - A prime number
/// * `rng` - Random number generator
///
/// # Returns
///
/// - `Ok(Some(t))` if p can be decomposed as t*conj(t) in Z[ω]
/// - `Ok(None)` if p is inert in Z[ω]
/// - `Err(msg)` if an error occurs
///
/// # Algorithm
///
/// Uses different strategies based on p mod 8:
/// - p ≡ 1 (mod 4): Uses sqrt(-1) mod p
/// - p ≡ 3 (mod 8): Uses sqrt(-2) mod p
/// - p ≡ 7 (mod 8): Checks if 2 is a quadratic residue
fn adj_decompose_int_prime<R>(p: &IBig, rng: &mut R) -> Result<Option<ZOmega>, String>
where
    R: Rng + ?Sized,
{
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
    if is_prime(p.clone(), 4, rng) {
        if &p & 0b11 == IBig::ONE {
            if let Some(h) = sqrt_negative_one(&p, 100, rng) {
                let t = ZOmega::gcd(
                    h + ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
                    ZOmega::from_int(p.clone()),
                );
                let t_conj = t.conj();
                if t_conj * &t == p || t_conj * &t == -p {
                    return Ok(Some(t));
                }
            }
            Ok(None)
        } else if &p & 0b111 == IBig::from(3) {
            if let Some(h) = root_mod(IBig::from(-2), &p, 100, rng) {
                let t = ZOmega::gcd(
                    h + ZOmega::new(IBig::ONE, IBig::ZERO, IBig::ONE, IBig::ZERO),
                    ZOmega::from_int(p.clone()),
                );
                let t_conj = t.conj();
                if t_conj * &t == p || t_conj * &t == -p {
                    return Ok(Some(t));
                }
            }
            Ok(None)
        } else if &p & 0b111 == IBig::from(7) {
            if root_mod(IBig::from(2), &p, 100, rng).is_some() {
                Err("No solution in adj_decompose_int_prime".to_string())
            } else {
                Ok(None)
            }
        } else {
            Ok(None)
        }
    } else if &p & 0b111 == IBig::from(7) {
        if root_mod(IBig::from(2), &p, 100, rng).is_some() {
            Err("No solution in adj_decompose_int_prime".to_string())
        } else {
            Ok(None)
        }
    } else {
        Ok(None)
    }
}

/// Decomposes a prime power p^k into Z[ω] if possible.
///
/// # Arguments
///
/// * `p` - A prime number
/// * `k` - The exponent
/// * `rng` - Random number generator
///
/// # Returns
///
/// - `Ok(Some(t))` if p^k can be decomposed in Z[ω]
/// - `Ok(None)` if decomposition is not possible
/// - `Err(msg)` if an error occurs
///
/// # Algorithm
///
/// For even k, returns p^(k/2). For odd k, decomposes p and raises to power k.
fn adj_decompose_int_prime_power<R>(p: &IBig, k: i32, rng: &mut R) -> Result<Option<ZOmega>, String>
where
    R: Rng + ?Sized,
{
    if k & 1 == 0 {
        Ok(Some(ZOmega::from_int(p.pow((k / 2) as usize))))
    } else {
        match adj_decompose_int_prime(p, rng) {
            Ok(Some(t)) => Ok(Some(t.pow(k.try_into().unwrap()))),
            Ok(None) => Ok(None),
            Err(e) => Err(e),
        }
    }
}

/// Decomposes an integer into Z[ω] by factoring and decomposing each prime power.
///
/// # Arguments
///
/// * `n` - The integer to decompose
/// * `start_time` - Start time for timeout tracking
/// * `diophantine_data` - Configuration including timeout and RNG
///
/// # Returns
///
/// - `Ok(Some(t))` if n can be decomposed as t*conj(t) in Z[ω]
/// - `Ok(None)` if decomposition fails or times out
/// - `Err(msg)` if an error occurs
///
/// # Algorithm
///
/// Factors n and attempts to decompose each prime power factor,
/// respecting timeout constraints.
fn adj_decompose_int(
    mut n: IBig,
    start_time: Instant,
    diophantine_data: &mut DiophantineData,
) -> Result<Option<ZOmega>, String> {
    n = n.abs();
    let mut facs = vec![(n, 1)];
    let mut t = ZOmega::from_int(IBig::ONE);

    while let Some((p, k)) = facs.pop() {
        if start_time.elapsed().as_millis() >= diophantine_data.diophantine_timeout {
            return Ok(None);
        }

        if p.bit_len() > 2048 {
            return Ok(None);
        }

        match adj_decompose_int_prime_power(&p, k, &mut diophantine_data.rng) {
            Ok(None) => {
                let individual_timeout = std::cmp::min(diophantine_data.factoring_timeout, 20);
                if let Some(fac) =
                    find_factor(&p, individual_timeout, 128, &mut diophantine_data.rng)
                {
                    facs.push((p / fac.clone(), k));
                    facs.push((fac, k));
                    let (_, new_facs) = decompose_relatively_prime(facs);
                    facs = new_facs;
                } else {
                    return Ok(None);
                }
            }
            Ok(Some(t_p)) => {
                t = &t * &t_p;
            }
            Err(e) => {
                return Err(e);
            }
        }
    }
    Ok(Some(t))
}

/// Decomposes a self-associate element of Z[√2] into Z[ω].
///
/// # Arguments
///
/// * `xi` - A self-associate element (xi = conj(xi))
/// * `start_time` - Start time for timeout tracking
/// * `diophantine_data` - Configuration including timeout and RNG
///
/// # Returns
///
/// - `Ok(Some(t))` if xi can be decomposed in Z[ω]
/// - `Ok(None)` if decomposition fails or times out
/// - `Err(msg)` if an error occurs
fn adj_decompose_selfassociate(
    xi: ZRootTwo,
    start_time: Instant,
    diophantine_data: &mut DiophantineData,
) -> Result<Option<ZOmega>, String> {
    if xi == ZRootTwo::from_int(IBig::ZERO) {
        return Ok(Some(ZOmega::from_int(IBig::ZERO)));
    }
    let n: IBig = (&xi.a).gcd(&xi.b).into();
    let r = &xi / n.clone();
    let t1 = adj_decompose_int(n, start_time, diophantine_data);
    let t2 = if r % ZRootTwo::new(IBig::ZERO, IBig::ONE) == ZRootTwo::from_int(IBig::ZERO) {
        ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ONE, IBig::ONE)
    } else {
        ZOmega::from_int(IBig::ONE)
    };
    match t1 {
        Ok(Some(t1_val)) => Ok(Some(&t1_val * &t2)),
        Ok(None) => Ok(None),
        Err(e) => Err(e),
    }
}

/// Decomposes a factorization in Z[√2] into relatively prime factors.
///
/// # Arguments
///
/// * `partial_facs` - A vector of (base, exponent) pairs in Z[√2]
///
/// # Returns
///
/// A tuple `(u, facs)` where:
/// - `u` is a unit factor in Z[√2]
/// - `facs` is a vector of relatively prime factors with their exponents
///
/// # Algorithm
///
/// Similar to `decompose_relatively_prime` but works in the ring Z[√2],
/// using similarity and GCD operations specific to that ring.
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
                u = u * (&b / a).pow(&IBig::from(k_b));
                facs[i] = (a.clone(), k_a + k_b);
                break;
            } else {
                let g = ZRootTwo::gcd(a.clone(), b.clone());
                if ZRootTwo::sim(g.clone(), ZRootTwo::from_int(IBig::ONE)) {
                    i += 1;
                    continue;
                } else {
                    let partial = vec![(a / &g, *k_a), (g.clone(), k_a + k_b)];
                    let (u_a, mut facs_a) = decompose_relatively_zomega_prime(partial);
                    u = u * u_a;
                    facs[i] = facs_a.remove(0);
                    facs.append(&mut facs_a);
                    stack.push((&b / &g, k_b));
                    break;
                }
            }
        }
    }
    (u, facs)
}

/// Decomposes a prime element of Z[√2] into Z[ω] if possible.
///
/// # Arguments
///
/// * `eta` - A prime element in Z[√2]
/// * `rng` - Random number generator
///
/// # Returns
///
/// - `Ok(Some(t))` if eta can be decomposed in Z[ω]
/// - `Ok(None)` if eta is inert in Z[ω]
/// - `Err(msg)` if an error occurs
///
/// # Algorithm
///
/// Uses the norm of eta and applies similar techniques as `adj_decompose_int_prime`,
/// but verifies the result divides eta in Z[√2].
pub fn adj_decompose_zomega_prime<R: Rng + ?Sized>(
    eta: ZRootTwo,
    rng: &mut R,
) -> Result<Option<ZOmega>, String> {
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

    if is_prime(p.clone(), 4, rng) {
        let check: i32 = (&p & IBig::from(0b111)).try_into().unwrap();
        match check {
            1 => Ok(sqrt_negative_one(&p, 100, rng).and_then(|h| {
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
            3 => Ok(root_mod(IBig::from(-2), &p, 100, rng).and_then(|h| {
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
                if root_mod(IBig::from(2), &p, 100, rng).is_some() {
                    Err("No solution in adj_decompose_zomega_prime".to_string())
                } else {
                    Ok(None)
                }
            }
            _ => Ok(None),
        }
    } else if &p & 0b111 == IBig::from(7) {
        if root_mod(IBig::from(2), &p, 100, rng).is_some() {
            Err("No solution in adj_decompose_zomega_prime".to_string())
        } else {
            Ok(None)
        }
    } else {
        Ok(None)
    }
}

/// Decomposes a prime power in Z[√2] into Z[ω] if possible.
///
/// # Arguments
///
/// * `eta` - A prime element in Z[√2]
/// * `k` - The exponent
/// * `rng` - Random number generator
///
/// # Returns
///
/// - `Ok(Some(t))` if eta^k can be decomposed in Z[ω]
/// - `Ok(None)` if decomposition is not possible
/// - `Err(msg)` if an error occurs
pub fn adj_decompose_zomega_prime_power<R: Rng + ?Sized>(
    eta: ZRootTwo,
    k: i32,
    rng: &mut R,
) -> Result<Option<ZOmega>, String> {
    if k & 1 == 0 {
        Ok(Some(ZOmega::from_zroottwo(&eta.pow(&IBig::from(k / 2)))))
    } else {
        match adj_decompose_zomega_prime(eta, rng) {
            Ok(Some(t)) => Ok(Some(t.pow(k.try_into().unwrap()))),
            Ok(None) => Ok(None),
            Err(e) => Err(e),
        }
    }
}

/// Decomposes a self-coprime element of Z[√2] into Z[ω].
///
/// # Arguments
///
/// * `xi` - A self-coprime element (gcd(xi, conj(xi)) = 1)
/// * `start_time` - Start time for timeout tracking
/// * `diophantine_data` - Configuration including timeout and RNG
///
/// # Returns
///
/// - `Ok(Some(t))` if xi can be decomposed in Z[ω]
/// - `Ok(None)` if decomposition fails or times out
/// - `Err(msg)` if an error occurs
///
/// # Algorithm
///
/// Factors xi in Z[√2] and decomposes each prime power factor.
pub fn adj_decompose_selfcoprime(
    xi: ZRootTwo,
    start_time: Instant,
    diophantine_data: &mut DiophantineData,
) -> Result<Option<ZOmega>, String> {
    let mut facs: Vec<(ZRootTwo, i32)> = vec![(xi.clone(), 1)];
    let mut t = ZOmega::from_int(IBig::ONE);
    while let Some((eta, k)) = facs.pop() {
        if start_time.elapsed().as_millis() >= diophantine_data.diophantine_timeout {
            return Ok(None);
        }

        match adj_decompose_zomega_prime_power(eta.clone(), k, &mut diophantine_data.rng) {
            Ok(None) => {
                let mut n = eta.norm();
                if n < IBig::ZERO {
                    n = -n;
                }
                let individual_timeout = std::cmp::min(diophantine_data.factoring_timeout, 15); // Max 15ms per factor
                if let Some(fac_n) =
                    find_factor(&n, individual_timeout, 128, &mut diophantine_data.rng)
                {
                    //  &mut diophantine_data.rng) {
                    let fac = ZRootTwo::gcd(xi.clone(), ZRootTwo::from_int(fac_n));
                    facs.push((&eta / &fac, k));
                    facs.push((fac, k));
                    let (_, new_facs) = decompose_relatively_zomega_prime(facs);
                    facs = new_facs;
                } else {
                    return Ok(None);
                }
            }
            Ok(Some(t_eta)) => {
                t = &t * &t_eta;
            }
            Err(e) => {
                return Err(e);
            }
        }
    }
    Ok(Some(t))
}

/// Main decomposition function for elements of Z[√2] into Z[ω].
///
/// # Arguments
///
/// * `xi` - An element of Z[√2] to decompose
/// * `start_time` - Start time for timeout tracking
/// * `diophantine_data` - Configuration including timeout and RNG
///
/// # Returns
///
/// - `Ok(Some(t))` if xi can be decomposed in Z[ω]
/// - `Ok(None)` if decomposition fails or times out
/// - `Err(msg)` if an error occurs
///
/// # Algorithm
///
/// Splits xi into self-associate and self-coprime parts, then decomposes each separately.
fn adj_decompose(
    xi: ZRootTwo,
    start_time: Instant,
    diophantine_data: &mut DiophantineData,
) -> Result<Option<ZOmega>, String> {
    if xi == ZRootTwo::from_int(IBig::ZERO) {
        return Ok(Some(ZOmega::from_int(IBig::ZERO)));
    }
    let d = ZRootTwo::gcd(xi.clone(), xi.conj_sq2());
    let eta = &xi / &d;
    let t1 = adj_decompose_selfassociate(d, start_time, diophantine_data);
    match t1 {
        Ok(Some(t1_val)) => {
            let t2 = adj_decompose_selfcoprime(eta, start_time, diophantine_data);
            match t2 {
                Ok(Some(t2_val)) => Ok(Some(&t1_val * &t2_val)),
                Ok(None) => Ok(None),
                Err(e) => Err(e),
            }
        }
        Ok(None) => Ok(None),
        Err(e) => Err(e),
    }
}
/// Solves the Diophantine equation to find t ∈ Z[ω] such that t*conj(t) = xi.
///
/// # Arguments
///
/// * `xi` - A positive element of Z[√2]
/// * `start_time` - Start time for timeout tracking
/// * `diophantine_data` - Configuration including timeout and RNG
///
/// # Returns
///
/// - `Ok(Some(t))` where t*conj(t) = xi
/// - `Ok(None)` if no solution exists or computation times out
/// - `Err(msg)` if xi is negative or an error occurs
///
/// # Algorithm
///
/// 1. Decomposes xi using `adj_decompose`
/// 2. Computes the unit factor u = xi / (t*conj(t))
/// 3. Takes the square root of u and multiplies by t
fn diophantine(
    xi: &ZRootTwo,
    start_time: Instant,
    diophantine_data: &mut DiophantineData,
) -> Result<Option<ZOmega>, String> {
    if xi == &ZRootTwo::from_int(IBig::ZERO) {
        return Ok(Some(ZOmega::from_int(IBig::ZERO)));
    } else if xi < &ZRootTwo::from_int(IBig::ZERO)
        || xi.clone().conj_sq2() < ZRootTwo::from_int(IBig::ZERO)
    {
        return Err("No solution".to_string());
    }
    let t = adj_decompose(xi.clone(), start_time, diophantine_data);
    match t {
        Ok(Some(t)) => {
            let xi_associate = ZRootTwo::from_zomega(t.conj() * &t);
            let u = xi / &xi_associate;
            match u.sqrt() {
                Some(v) => Ok(Some(&ZOmega::from_zroottwo(&v) * &t)),
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

/// Clears all internal caches used by the Diophantine solver.
///
/// This includes:
/// - Primality test cache
/// - Square root modulo cache
/// - Factorization cache
/// - Diophantine solution cache
///
/// Useful for freeing memory or ensuring fresh computations in testing.
pub fn clear_caches() {
    if let Ok(mut cache) = PRIMALITY_CACHE.try_lock() {
        cache.clear();
    //        println!("cleared primality cache");
    } else {
        //        println!("Can't clear primality cache");
    }
    if let Ok(mut cache) = SQRT_CACHE.try_lock() {
        cache.clear();
        //        println!("cleared sqrt cache");
    }
    if let Ok(mut cache) = FACTOR_CACHE.try_lock() {
        cache.clear();
        //        println!("cleared factor cache");
    }
    if let Ok(mut cache) = DIOPHANTINE_CACHE.try_lock() {
        cache.clear();
        //        println!("cleared diophantine cache");
    }
}

/// Solves the Diophantine equation for dyadic elements (elements with denominator 2^k).
///
/// # Arguments
///
/// * `xi` - A dyadic element of D[√2] (rationals with denominator 2^k)
/// * `diophantine_data` - Configuration including timeout and RNG
///
/// # Returns
///
/// `Some(t)` where t ∈ D[ω] and t*conj(t) = xi, or `None` if no solution found.
///
/// # Caching
///
/// Results are cached in `DIOPHANTINE_CACHE` for performance.
///
/// # Algorithm
///
/// 1. Normalizes xi by handling the dyadic denominator
/// 2. Calls the main `diophantine` solver on the normalized value
/// 3. Adjusts the result back to dyadic form
/// 4. Applies timeout optimizations based on input size
pub(crate) fn diophantine_dyadic(
    xi: DRootTwo,
    diophantine_data: &mut DiophantineData,
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
            std::cmp::min(diophantine_data.diophantine_timeout, 50),
            std::cmp::min(diophantine_data.factoring_timeout, 10),
        )
    } else {
        (
            std::cmp::min(diophantine_data.diophantine_timeout, 15),
            std::cmp::min(diophantine_data.factoring_timeout, 3),
        )
    };

    // Ugh. I don't like mutating this. But I don't see an easy way around it.
    diophantine_data.diophantine_timeout = optimized_diophantine_timeout;
    diophantine_data.factoring_timeout = optimized_factoring_timeout;

    let start_time = Instant::now();
    let t = diophantine(&alpha, start_time, diophantine_data);
    let result = match t {
        Err(_) => None,
        Ok(None) => None,
        Ok(Some(mut t)) => {
            if k_mod_2 == 1 {
                t = &t * &ZOmega::new(IBig::ZERO, IBig::NEG_ONE, IBig::ONE, IBig::ZERO);
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
