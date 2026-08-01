// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::common::ib_to_bf_prec;
use crate::math::{floorsqrt, rounddiv, sign, sqrt2};
use crate::ring::ZOmega;
use dashu_base::Sign;
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;
use std::cmp::Ordering;
use std::fmt;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Rem, Sub};

/// Represents an element of Z[√2], the ring of integers extended by √2.
///
/// Elements are represented as a + b√2 where a, b ∈ Z.
/// This ring is used extensively in the GridSynth algorithm.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ZRootTwo {
    pub(crate) a: IBig,
    pub(crate) b: IBig,
}

pub const ONE: ZRootTwo = ZRootTwo {
    a: IBig::ONE,
    b: IBig::ZERO,
};

pub const DELTA_SQUARED: ZRootTwo = ZRootTwo {
    a: IBig::from_parts_const(Sign::Positive, 2),
    b: IBig::ONE,
};

// // Absolute value squared of root two conjugate of δ.
// pub const DELTA_SQUARED_M: ZRootTwo = ZRootTwo {
//     a: IBig::from_parts_const(Sign::Positive, 2),
//     b: IBig::NEG_ONE,
// };

// See Definition 3.5 on pg 3 of R+S for the definition of lambda (and delta)
// Careful! There is a different, unrelated defintion of lambda in
// the discussion in Definition 9.1 on page 19 of R+S. It is in fact
// the fixed phase factor used in the up-to-phase algorithm.
pub const LAMBDA: ZRootTwo = ZRootTwo {
    a: IBig::ONE,
    b: IBig::ONE,
};

impl ZRootTwo {
    /// Creates a new element of Z[√2].
    ///
    /// # Arguments
    ///
    /// * `a` - The rational part
    /// * `b` - The coefficient of √2
    ///
    /// # Returns
    ///
    /// An element representing a + b√2.
    pub fn new(a: IBig, b: IBig) -> Self {
        Self { a, b }
    }

    /// Creates an element from an integer.
    ///
    /// # Arguments
    ///
    /// * `x` - An integer value
    ///
    /// # Returns
    ///
    /// An element representing x in Z[√2].
    pub fn from_int(x: IBig) -> Self {
        Self {
            a: x,
            b: IBig::ZERO,
        }
    }

    pub fn from_borrowed_int(x: &IBig) -> Self {
        Self {
            a: x.clone(),
            b: IBig::ZERO,
        }
    }

    /// Returns the parity (least significant bit) of the rational part.
    ///
    /// # Returns
    ///
    /// 0 if a is even, 1 if a is odd.
    pub fn parity(&self) -> IBig {
        &self.a & IBig::ONE
    }

    /// Computes the norm of the element.
    ///
    /// # Returns
    ///
    /// The norm N(a + b√2) = a² - 2b².
    pub fn norm(&self) -> IBig {
        &self.a * &self.a - IBig::from(2) * &self.b * &self.b
    }

    /// Converts to a high-precision floating-point number.
    ///
    /// # Returns
    ///
    /// The real value a + b√2 as a floating-point number.
    pub fn to_real(&self) -> FBig<HalfEven> {
        &self.a + sqrt2() * &self.b
    }

    /// Computes the conjugate with respect to √2.
    ///
    /// # Returns
    ///
    /// The conjugate a - b√2.
    pub fn conj_sq2(&self) -> Self {
        Self {
            a: self.a.clone(),
            b: -self.b.clone(),
        }
    }

    /// Computes the multiplicative inverse if it exists.
    ///
    /// # Returns
    ///
    /// The inverse element.
    ///
    /// # Panics
    ///
    /// Panics if the norm is not ±1 (element is not a unit).
    pub fn inv(&self) -> Self {
        let norm = self.norm();
        if norm == IBig::ONE {
            self.conj_sq2()
        } else if norm == IBig::NEG_ONE {
            -self.conj_sq2()
        } else {
            panic!(
                "ZeroDivisionError: cannot invert ZRootTwo with norm {}",
                norm
            );
        }
    }

    /// Converts an element from Z[ω] to Z[√2] if possible.
    ///
    /// # Arguments
    ///
    /// * `x` - An element of Z[ω]
    ///
    /// # Returns
    ///
    /// The corresponding element in Z[√2].
    ///
    /// # Panics
    ///
    /// Panics if the element cannot be represented in Z[√2].
    pub fn from_zomega(x: ZOmega) -> Self {
        if x.b == IBig::ZERO && x.a == -&x.c {
            Self::new(x.d, x.c)
        } else {
            panic!("Cannot convert ZOmega to ZRootTwo: {:?}", x);
        }
    }

    /// Computes the square root if it exists in Z[√2].
    ///
    /// # Returns
    ///
    /// `Some(w)` where w² = self, or `None` if no such element exists.
    pub fn sqrt(&self) -> Option<Self> {
        let norm = ib_to_bf_prec(self.norm());
        if norm < ib_to_bf_prec(IBig::ZERO) || self.a < IBig::ZERO {
            return None;
        }
        let r = floorsqrt(norm);
        let a1 = floorsqrt(ib_to_bf_prec((&self.a + &r) / IBig::from(2)));
        let b1 = floorsqrt(ib_to_bf_prec((&self.a - &r) / IBig::from(4)));
        let a2 = floorsqrt(ib_to_bf_prec((&self.a - &r) / IBig::from(2)));
        let b2 = floorsqrt(ib_to_bf_prec((&self.a + &r) / IBig::from(4)));

        let (w1, w2) =
            if sign(ib_to_bf_prec(self.a.clone())) * sign(ib_to_bf_prec(self.b.clone())) >= 0 {
                (Self::new(a1, b1), Self::new(a2, b2))
            } else {
                (Self::new(a1, -b1), Self::new(a2, -b2))
            };

        if (*self) == (&w1 * &w1) {
            Some(w1)
        } else if (*self) == (&w2 * &w2) {
            Some(w2)
        } else {
            None
        }
    }

    /// Checks if two elements are similar (associates).
    ///
    /// # Arguments
    ///
    /// * `a` - First element
    /// * `b` - Second element
    ///
    /// # Returns
    ///
    /// `true` if a and b divide each other (are associates).
    pub fn sim(a: Self, b: Self) -> bool {
        a.clone() % b.clone() == ZRootTwo::from_int(IBig::ZERO)
            && b % a == ZRootTwo::from_int(IBig::ZERO)
    }

    /// Computes the extended GCD in Z[√2].
    ///
    /// # Arguments
    ///
    /// * `a` - First element
    /// * `b` - Second element
    ///
    /// # Returns
    ///
    /// A tuple `(x, z, gcd)` where gcd = gcd(a, b) and x·a + z·b = gcd.
    pub fn ext_gcd(mut a: Self, mut b: Self) -> (Self, Self, Self) {
        let mut x = ZRootTwo::from_int(IBig::ONE);
        let mut y = ZRootTwo::from_int(IBig::ZERO);
        let mut z = ZRootTwo::from_int(IBig::ZERO);
        let mut w = ZRootTwo::from_int(IBig::ONE);

        while b != ZRootTwo::from_int(IBig::ZERO) {
            let (q, r) = divmod(&a, &b);
            let tmp_y = y.clone();
            y = &x - &y * &q;
            x = tmp_y;

            let tmp_w = w.clone();
            w = &z - &w * &q;
            z = tmp_w;

            a = b;
            b = r;
        }
        (x, z, a)
    }

    /// Computes the greatest common divisor in Z[√2].
    ///
    /// # Arguments
    ///
    /// * `a` - First element
    /// * `b` - Second element
    ///
    /// # Returns
    ///
    /// The GCD of a and b.
    pub fn gcd(a: Self, b: Self) -> Self {
        let (_, _, g) = Self::ext_gcd(a, b);
        g
    }

    /// Raises the element to a power using binary exponentiation.
    ///
    /// # Arguments
    ///
    /// * `exp` - The exponent (can be negative for units)
    ///
    /// # Returns
    ///
    /// self^exp computed efficiently.
    pub fn pow(&self, exp: &IBig) -> Self {
        if *exp < IBig::ZERO {
            return self.inv().pow(&(-exp));
        }

        let mut base = self.clone();
        let mut result = ZRootTwo::from_int(IBig::ONE);
        let mut e = exp.clone();

        while e > IBig::ZERO {
            if &e & IBig::ONE == IBig::ONE {
                result = &result * &base;
            }
            base = &base * &base;
            e >>= 1;
        }

        result
    }
}

impl fmt::Display for ZRootTwo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{:+}√2", self.a, self.b)
    }
}

impl Neg for ZRootTwo {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self::new(-self.a, -self.b)
    }
}

impl Add for ZRootTwo {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::new(self.a + other.a, self.b + other.b)
    }
}

impl Add<&ZRootTwo> for &ZRootTwo {
    type Output = ZRootTwo;

    fn add(self, other: &ZRootTwo) -> ZRootTwo {
        ZRootTwo::new(&self.a + &other.a, &self.b + &other.b)
    }
}

impl Add<IBig> for ZRootTwo {
    type Output = Self;
    fn add(self, other: IBig) -> Self {
        self + ZRootTwo::from_int(other)
    }
}

impl Add<&IBig> for ZRootTwo {
    type Output = Self;
    fn add(self, other: &IBig) -> Self {
        self + ZRootTwo::from_borrowed_int(other)
    }
}

impl AddAssign<&ZRootTwo> for ZRootTwo {
    fn add_assign(&mut self, rhs: &ZRootTwo) {
        self.a += &rhs.a;
        self.b += &rhs.b;
    }
}

// impl Sub for ZRootTwo {
//     type Output = Self;
//     fn sub(self, other: Self) -> Self {
//         self + (-other)
//     }
// }

impl Sub<ZRootTwo> for &ZRootTwo {
    type Output = ZRootTwo;
    fn sub(self, other: ZRootTwo) -> ZRootTwo {
        self + &(-other)
    }
}

impl Sub<IBig> for &ZRootTwo {
    type Output = ZRootTwo;
    fn sub(self, other: IBig) -> ZRootTwo {
        self - ZRootTwo::from_int(other)
    }
}

impl Mul for ZRootTwo {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let a = &self.a * &other.a + 2 * &self.b * &other.b;
        let b = &self.a * &other.b + &self.b * &other.a;
        Self::new(a, b)
    }
}

impl Mul<&ZRootTwo> for ZRootTwo {
    type Output = ZRootTwo;

    fn mul(self, other: &ZRootTwo) -> ZRootTwo {
        let a = &self.a * &other.a + 2 * &self.b * &other.b;
        let b = &self.a * &other.b + &self.b * &other.a;
        ZRootTwo::new(a, b)
    }
}

impl Mul<&ZRootTwo> for &ZRootTwo {
    type Output = ZRootTwo;
    fn mul(self, other: &ZRootTwo) -> ZRootTwo {
        let a = &self.a * &other.a + 2 * &self.b * &other.b;
        let b = &self.a * &other.b + &self.b * &other.a;
        ZRootTwo::new(a, b)
    }
}

impl Mul<IBig> for ZRootTwo {
    type Output = Self;
    fn mul(self, other: IBig) -> Self {
        Self::new(self.a * other.clone(), self.b * other)
    }
}

impl Rem for ZRootTwo {
    type Output = Self;
    fn rem(self, other: Self) -> Self::Output {
        let (_, r) = divmod(&self, &other);
        r
    }
}

impl Div<&ZRootTwo> for &ZRootTwo {
    type Output = ZRootTwo;
    fn div(self, other: &ZRootTwo) -> ZRootTwo {
        let (q, _) = divmod(self, other);
        q
    }
}

impl Div<IBig> for &ZRootTwo {
    type Output = ZRootTwo;
    fn div(self, other: IBig) -> ZRootTwo {
        self / &ZRootTwo::from_int(other)
    }
}

impl Ord for ZRootTwo {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.b.cmp(&other.b) {
            Ordering::Less => {
                if self.a < other.a || (&self.a - &other.a).pow(2) < 2 * (&self.b - &other.b).pow(2)
                {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            }
            Ordering::Greater => {
                if self.a < other.a && (&self.a - &other.a).pow(2) > 2 * (&self.b - &other.b).pow(2)
                {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            }
            Ordering::Equal => self.a.cmp(&other.a),
        }
    }
}

impl PartialOrd for ZRootTwo {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Performs division with remainder in Z[√2].
///
/// # Arguments
///
/// * `x` - The dividend
/// * `y` - The divisor
///
/// # Returns
///
/// A tuple `(quotient, remainder)` such that x = y × quotient + remainder.
pub fn divmod(x: &ZRootTwo, y: &ZRootTwo) -> (ZRootTwo, ZRootTwo) {
    let p = x * &y.conj_sq2();
    let k = y.norm();
    let q = ZRootTwo::new(rounddiv(p.a, &k), rounddiv(p.b, &k));
    let r = x - y * &q;
    (q, r)
}
