// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;
use nalgebra::Complex;
use std::cell::OnceCell;
use std::cmp::{Ordering, PartialEq};
use std::fmt::{Debug, Display, Formatter, Result};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::math::{ntz, pow_sqrt2};
use crate::ring::{DRootTwo, ZOmega, ZRootTwo};

#[derive(Clone)]
pub struct DOmega {
    pub u: ZOmega,
    pub k: i64,
    // TODO: Naming caches for various computations is not correct because cloning may be expensive
    scale_cache: OnceCell<FBig<HalfEven>>,
    conj_cache: OnceCell<Box<Self>>,
    conj_sq2_cache: OnceCell<Box<Self>>,
    real_cache: OnceCell<FBig<HalfEven>>,
    imag_cache: OnceCell<FBig<HalfEven>>,
}

impl DOmega {
    pub fn new(u: ZOmega, k: i64) -> Self {
        Self {
            u,
            k,
            scale_cache: OnceCell::new(),
            conj_cache: OnceCell::new(),
            conj_sq2_cache: OnceCell::new(),
            real_cache: OnceCell::new(),
            imag_cache: OnceCell::new(),
        }
    }

    pub fn from_int(x: IBig) -> Self {
        Self::new(ZOmega::from_int(x), 0)
    }

    pub fn from_zroottwo(x: &ZRootTwo) -> Self {
        Self::new(ZOmega::from_zroottwo(x), 0)
    }

    pub fn from_zomega(x: ZOmega) -> Self {
        Self::new(x, 0)
    }

    pub fn from_droottwo(alpha: ZRootTwo, k: i64) -> Self {
        Self::new(ZOmega::from_zroottwo(&alpha), k)
    }

    pub fn from_droottwo_vector(x: &DRootTwo, y: &DRootTwo, k: i64) -> Self {
        let omega = ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO);
        let term0 = Self::new(ZOmega::from_zroottwo(&x.alpha), x.k);
        let term1 = Self::new(ZOmega::from_zroottwo(&y.alpha), y.k) * omega;
        (term0 + term1).renew_denomexp(k)
    }

    pub fn scale(&self) -> &FBig<HalfEven> {
        self.scale_cache.get_or_init(|| pow_sqrt2(self.k))
    }

    pub fn real(&self) -> &FBig<HalfEven> {
        self.real_cache.get_or_init(|| self.u.real() / self.scale())
    }

    pub fn imag(&self) -> &FBig<HalfEven> {
        self.imag_cache.get_or_init(|| self.u.imag() / self.scale())
    }

    pub fn to_complex(&self) -> Complex<&FBig<HalfEven>> {
        Complex::new(self.real(), self.imag())
    }

    pub fn conj(&self) -> &Self {
        self.conj_cache
            .get_or_init(|| Box::new(Self::new(self.u.conj().clone(), self.k)))
            .as_ref()
    }

    pub fn conj_sq2(&self) -> &Self {
        self.conj_sq2_cache
            .get_or_init(|| {
                if self.k & 1 == 1 {
                    Box::new(Self::new(-self.u.conj_sq2().clone(), self.k))
                } else {
                    Box::new(Self::new(self.u.conj_sq2().clone(), self.k))
                }
            })
            .as_ref()
    }

    pub fn renew_denomexp(&self, new_k: i64) -> Self {
        let d = new_k - self.k;
        let new_u = self.mul_by_sqrt2_power(d).u;
        Self::new(new_u, new_k)
    }

    pub fn reduce_denomexp(&self) -> Self {
        let zero = IBig::ZERO;
        let k_a = if self.u.a == zero {
            self.k
        } else {
            ntz(&self.u.a)
        };
        let k_b = if self.u.b == zero {
            self.k
        } else {
            ntz(&self.u.b)
        };
        let k_c = if self.u.c == zero {
            self.k
        } else {
            ntz(&self.u.c)
        };
        let k_d = if self.u.d == zero {
            self.k
        } else {
            ntz(&self.u.d)
        };

        let reduce_k = k_a.min(k_b).min(k_c).min(k_d);

        let mut new_k = self.k - (reduce_k << 1);
        let k_usize: usize = reduce_k.try_into().expect("k must fit in i128");

        let bit: IBig = (IBig::ONE << (k_usize + 1)) - IBig::ONE;

        let ca_add = &self.u.c + &self.u.a;
        let bd_add = &self.u.b + &self.u.d;

        if (&ca_add & &bit).is_zero() && (&bd_add & &bit).is_zero() {
            new_k -= 1;
        }

        self.renew_denomexp(new_k.max(0))
    }

    pub fn mul_by_inv_sqrt2(&self) -> Self {
        let (a, b, c, d) = (
            self.u.a.clone(),
            self.u.b.clone(),
            self.u.c.clone(),
            self.u.d.clone(),
        );

        if ((&b + &d) & 1) == IBig::ZERO && ((&c + &a) & 1) == IBig::ZERO {
            let new_u = ZOmega::new(
                (&b - &d) >> 1,
                (&c + &a) >> 1,
                (&b + &d) >> 1,
                (&c - &a) >> 1,
            );
            Self::new(new_u, self.k)
        } else {
            panic!("Invalid coefficients for inverse sqrt2 multiplication")
        }
    }

    pub fn mul_by_sqrt2_power(&self, d: i64) -> Self {
        if d < 0 {
            if d == -1 {
                return self.mul_by_inv_sqrt2();
            }

            let d_neg = -d;
            let d_div_2 = &d_neg >> 1;
            let d_mod_2 = &d_neg & 1;
            let d_div_2_usize: usize = d_div_2.try_into().expect("k must fit in i128");

            if d_mod_2 == 0 {
                let bit = (IBig::ONE << d_div_2_usize) - IBig::ONE;

                if (&self.u.a & &bit).is_zero()
                    && (&self.u.b & &bit).is_zero()
                    && (&self.u.c & &bit).is_zero()
                    && (&self.u.d & bit).is_zero()
                {
                    let new_u = ZOmega::new(
                        &self.u.a >> d_div_2_usize,
                        &self.u.b >> d_div_2_usize,
                        &self.u.c >> d_div_2_usize,
                        &self.u.d >> d_div_2_usize,
                    );
                    return Self::new(new_u, self.k);
                }
            } else {
                let shift = d_div_2_usize + 1;
                let bit = (IBig::ONE << shift) - IBig::ONE;

                let bd_sub = &self.u.b - &self.u.d;
                let ca_add = &self.u.c + &self.u.a;
                let bd_add = &self.u.b + &self.u.d;
                let ca_sub = &self.u.c - &self.u.a;

                if (&bd_sub & &bit).is_zero()
                    && (&ca_add & &bit).is_zero()
                    && (&bd_add & &bit).is_zero()
                    && (&ca_sub & bit).is_zero()
                {
                    let new_u = ZOmega::new(
                        bd_sub >> shift,
                        ca_add >> shift,
                        bd_add >> shift,
                        ca_sub >> shift,
                    );
                    return Self::new(new_u, self.k);
                }
            }

            panic!("Invalid coefficients for sqrt2^-d");
        } else {
            let d_div_2 = &d >> 1;
            let d_mod_2 = &d & 1;
            let d_div_2_usize: usize = d_div_2.try_into().expect("k must fit in i128");

            let mut new_u = self.u.clone() * (IBig::ONE << d_div_2_usize);
            if d_mod_2 == 1 {
                new_u = new_u * ZOmega::new(IBig::NEG_ONE, IBig::ZERO, IBig::ONE, IBig::ZERO);
            }

            Self::new(new_u, self.k)
        }
    }

    pub fn mul_by_omega(&self) -> Self {
        Self::new(self.u.mul_by_omega(), self.k)
    }

    pub fn mul_by_omega_inv(&self) -> Self {
        Self::new(self.u.mul_by_omega_inv(), self.k)
    }

    pub fn mul_by_omega_power(&self, n: usize) -> Self {
        Self::new(self.u.mul_by_omega_power(n), self.k)
    }

    pub fn residue(&self) -> u8 {
        self.u.residue()
    }

    pub fn eq_general(&self, other: &DOmega) -> bool {
        match self.k.cmp(&other.k) {
            Ordering::Less => self.renew_denomexp(other.k) == *other,
            Ordering::Greater => *self == other.renew_denomexp(self.k),
            Ordering::Equal => self.u == other.u && self.k == other.k,
        }
    }
}

impl PartialEq for DOmega {
    fn eq(&self, other: &Self) -> bool {
        self.eq_general(other)
    }
}

impl PartialEq<IBig> for DOmega {
    fn eq(&self, other: &IBig) -> bool {
        self.eq_general(&DOmega::from_int(other.clone()))
    }
}

impl PartialEq<ZRootTwo> for DOmega {
    fn eq(&self, other: &ZRootTwo) -> bool {
        self.eq_general(&DOmega::from_zroottwo(other))
    }
}

impl PartialEq<ZOmega> for DOmega {
    fn eq(&self, other: &ZOmega) -> bool {
        self.eq_general(&DOmega::from_zomega(other.clone()))
    }
}

impl Add<DOmega> for DOmega {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        match self.k.cmp(&rhs.k) {
            std::cmp::Ordering::Less => self.renew_denomexp(rhs.k) + rhs,
            std::cmp::Ordering::Greater => {
                let k = self.k;
                self + rhs.renew_denomexp(k)
            }
            std::cmp::Ordering::Equal => Self::new(self.u + rhs.u, self.k),
        }
    }
}

impl Add<IBig> for DOmega {
    type Output = Self;
    fn add(self, rhs: IBig) -> Self {
        self + DOmega::from_int(rhs)
    }
}

impl AddAssign<DOmega> for DOmega {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

impl Sub<DOmega> for DOmega {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl Sub<IBig> for DOmega {
    type Output = Self;
    fn sub(self, rhs: IBig) -> Self {
        self - DOmega::from_int(rhs)
    }
}

impl SubAssign<DOmega> for DOmega {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs;
    }
}

impl Neg for DOmega {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.u, self.k)
    }
}

impl Mul<i32> for DOmega {
    type Output = Self;
    fn mul(self, other: i32) -> Self {
        Self::new(self.u * IBig::from(other), self.k)
    }
}

impl Mul<DOmega> for i32 {
    type Output = DOmega;
    fn mul(self, other: DOmega) -> DOmega {
        DOmega::new(other.u * IBig::from(self), other.k)
    }
}

impl Mul<IBig> for DOmega {
    type Output = Self;
    fn mul(self, other: IBig) -> Self {
        Self::new(self.u * other, self.k)
    }
}

impl Mul<DOmega> for DOmega {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::new(self.u * rhs.u, self.k + rhs.k)
    }
}

impl Mul<&DOmega> for &DOmega {
    type Output = DOmega;

    fn mul(self, rhs: &DOmega) -> DOmega {
        DOmega::new(&self.u * &rhs.u, self.k + rhs.k)
    }
}

impl Mul<ZOmega> for DOmega {
    type Output = Self;
    fn mul(self, other: ZOmega) -> Self {
        self * DOmega::from_zomega(other)
    }
}

impl MulAssign<DOmega> for DOmega {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
    }
}

impl From<IBig> for DOmega {
    fn from(value: IBig) -> Self {
        Self::from_int(value)
    }
}

impl From<ZRootTwo> for DOmega {
    fn from(value: ZRootTwo) -> Self {
        Self::from_zroottwo(&value)
    }
}

impl From<ZOmega> for DOmega {
    fn from(value: ZOmega) -> Self {
        Self::from_zomega(value)
    }
}

impl Debug for DOmega {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "DOmega({:?}, k={})", self.u, self.k)
    }
}

impl Display for DOmega {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{} / √2^{}", self.u, self.k)
    }
}
