// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::common::ib_to_bf_prec;
use crate::math::{floorsqrt, rounddiv, sign, sqrt2};
use crate::ring::ZOmega;
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;
use std::cmp::Ordering;
use std::fmt;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Rem, Sub};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ZRootTwo {
    pub(crate) a: IBig,
    pub(crate) b: IBig,
}

impl ZRootTwo {
    pub fn new(a: IBig, b: IBig) -> Self {
        Self { a, b }
    }

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

    pub fn parity(&self) -> IBig {
        &self.a & IBig::ONE
    }

    pub fn norm(&self) -> IBig {
        &self.a * &self.a - IBig::from(2) * &self.b * &self.b
    }

    pub fn to_real(&self) -> FBig<HalfEven> {
        &self.a + sqrt2() * &self.b
    }

    pub fn conj_sq2(&self) -> Self {
        Self {
            a: self.a.clone(),
            b: -self.b.clone(),
        }
    }

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

    pub fn from_zomega(x: ZOmega) -> Self {
        if x.b == IBig::ZERO && x.a == -&x.c {
            Self::new(x.d, x.c)
        } else {
            panic!("Cannot convert ZOmega to ZRootTwo: {:?}", x);
        }
    }

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

    pub fn sim(a: Self, b: Self) -> bool {
        a.clone() % b.clone() == ZRootTwo::from_int(IBig::ZERO)
            && b % a == ZRootTwo::from_int(IBig::ZERO)
    }

    pub fn ext_gcd(mut a: Self, mut b: Self) -> (Self, Self, Self) {
        let mut x = ZRootTwo::from_int(IBig::ONE);
        let mut y = ZRootTwo::from_int(IBig::ZERO);
        let mut z = ZRootTwo::from_int(IBig::ZERO);
        let mut w = ZRootTwo::from_int(IBig::ONE);

        while b != ZRootTwo::from_int(IBig::ZERO) {
            let (q, r) = divmod(&a, &b);
            let tmp_y = y.clone();
            y = x - &y * &q;
            x = tmp_y;

            let tmp_w = w.clone();
            w = z - &w * &q;
            z = tmp_w;

            a = b;
            b = r;
        }
        (x, z, a)
    }

    pub fn gcd(a: Self, b: Self) -> Self {
        let (_, _, g) = Self::ext_gcd(a, b);
        g
    }

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

    fn add(self, other: & ZRootTwo) -> ZRootTwo {
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
    fn add(self, other: & IBig) -> Self {
        self + ZRootTwo::from_borrowed_int(other)
    }
}


impl AddAssign<&ZRootTwo> for ZRootTwo {
    fn add_assign(&mut self, rhs: &ZRootTwo) {
        self.a += &rhs.a;
        self.b += &rhs.b;
    }
}

impl Sub for ZRootTwo {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl Sub<IBig> for ZRootTwo {
    type Output = Self;
    fn sub(self, other: IBig) -> Self {
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

impl Div for ZRootTwo {
    type Output = Self;
    fn div(self, other: Self) -> Self::Output {
        let (q, _) = divmod(&self, &other);
        q
    }
}

impl Div<IBig> for ZRootTwo {
    type Output = Self;
    fn div(self, other: IBig) -> Self::Output {
        self / ZRootTwo::from_int(other)
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

pub fn divmod(x: &ZRootTwo, y: &ZRootTwo) -> (ZRootTwo, ZRootTwo) {
    let p = x * &y.conj_sq2();
    let k = y.norm();
    let q = ZRootTwo::new(rounddiv(p.a, &k), rounddiv(p.b, &k));
    let r = x.clone() - y * &q;
    (q, r)
}
