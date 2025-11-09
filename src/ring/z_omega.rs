// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;

use crate::common::ib_to_bf_prec;
use crate::math::sqrt2;
use crate::ring::ZRootTwo;
use std::cell::OnceCell;
use std::fmt::Debug;
use std::fmt::{Display, Formatter, Result};
use std::mem;
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

#[derive(Clone, Eq)]
pub struct ZOmega {
    pub a: IBig,
    pub b: IBig,
    pub c: IBig,
    pub d: IBig,
    coef: OnceCell<[IBig; 4]>,
    // TODO: Naming caches for various computations is not correct because cloning may be expensive
    conj_cache: OnceCell<Box<Self>>,
    conj_sq2_cache: OnceCell<Box<Self>>,
    norm_cache: OnceCell<IBig>,
    residue_cache: OnceCell<u8>,
}

impl ZOmega {
    pub fn new(a: IBig, b: IBig, c: IBig, d: IBig) -> Self {
        Self {
            a,
            b,
            c,
            d,
            coef: OnceCell::new(),
            conj_cache: OnceCell::new(),
            conj_sq2_cache: OnceCell::new(),
            norm_cache: OnceCell::new(),
            residue_cache: OnceCell::new(),
        }
    }

    pub fn coef(&self) -> &[IBig; 4] {
        self.coef.get_or_init(|| {
            [
                self.d.clone(),
                self.c.clone(),
                self.b.clone(),
                self.a.clone(),
            ]
        })
    }

    pub fn conj(&self) -> &Self {
        self.conj_cache
            .get_or_init(|| {
                Box::new(Self::new(
                    -self.c.clone(),
                    -self.b.clone(),
                    -self.a.clone(),
                    self.d.clone(),
                ))
            })
            .as_ref()
    }

    pub fn conj_sq2(&self) -> &Self {
        self.conj_sq2_cache
            .get_or_init(|| {
                Box::new(Self::new(
                    -self.a.clone(),
                    self.b.clone(),
                    -self.c.clone(),
                    self.d.clone(),
                ))
            })
            .as_ref()
    }

    pub fn residue(&self) -> u8 {
        *self.residue_cache.get_or_init(|| {
            (u8::try_from(&(&self.a & 1)).unwrap() << 3)
                | (u8::try_from(&(&self.b & 1)).unwrap() << 2)
                | (u8::try_from(&(&self.c & 1)).unwrap() << 1)
                | u8::try_from(&(&self.d & 1)).unwrap()
        })
    }

    pub fn norm(&self) -> &IBig {
        self.norm_cache.get_or_init(|| {
            let (a, b, c, d) = (&self.a, &self.b, &self.c, &self.d);
            let s1 = a * a + b * b + c * c + d * d;
            let s2 = a * b + b * c + c * d - d * a;
            &s1 * &s1 - 2 * &s2 * &s2
        })
    }

    pub fn from_int(x: IBig) -> Self {
        Self::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, x)
    }

    pub fn from_zroottwo(x: &ZRootTwo) -> Self {
        Self::new(-x.b.clone(), IBig::ZERO, x.b.clone(), x.a.clone())
    }

    pub fn real(&self) -> FBig<HalfEven> {
        ib_to_bf_prec(self.d.clone()) + sqrt2() * (&self.c - &self.a) / 2
    }

    pub fn imag(&self) -> FBig<HalfEven> {
        ib_to_bf_prec(self.b.clone()) + sqrt2() * (&self.c + &self.a) / 2
    }

    pub fn mul_by_omega(&self) -> Self {
        Self::new(
            self.b.clone(),
            self.c.clone(),
            self.d.clone(),
            -self.a.clone(),
        )
    }

    pub fn mul_by_omega_inv(&self) -> Self {
        Self::new(
            -self.d.clone(),
            self.a.clone(),
            self.b.clone(),
            self.c.clone(),
        )
    }

    pub fn mul_by_omega_power(&self, n: usize) -> Self {
        let n = n % 8;
        if n & 0b100 != 0 {
            return -self.mul_by_omega_power(n & 0b11);
        }
        let coef = self.coef();
        let mut new_coef = [IBig::ZERO; 4];
        for i in 0..n {
            new_coef[i] = -coef[i.wrapping_sub(n) % 4].clone();
        }
        new_coef[n..4].clone_from_slice(&coef[..(4 - n)]);
        Self::new(
            mem::take(&mut new_coef[3]),
            mem::take(&mut new_coef[2]),
            mem::take(&mut new_coef[1]),
            mem::take(&mut new_coef[0]),
        )
    }

    pub fn divmod(&self, other: &Self) -> (Self, Self) {
        let p0 = self * other.conj();
        let p1 = other.conj().conj_sq2() * other.conj_sq2();
        let p = &p0 * &p1;
        let k = other.norm();
        let q = ZOmega::new(
            rounddiv(&p.a, k),
            rounddiv(&p.b, k),
            rounddiv(&p.c, k),
            rounddiv(&p.d, k),
        );
        let r = self - &(other * &q);
        (q, r)
    }

    pub fn ext_gcd(mut a: Self, mut b: Self) -> (Self, Self, Self) {
        let mut x = ZOmega::from_int(IBig::ONE);
        let mut y = ZOmega::from_int(IBig::ZERO);
        let mut z = ZOmega::from_int(IBig::ZERO);
        let mut w = ZOmega::from_int(IBig::ONE);
        while b != ZOmega::from_int(IBig::ZERO) {
            let (q, r) = a.divmod(&b);
            let yq = &y * &q;
            let wq = &w * &q;

            let x_new = std::mem::replace(&mut x, y);
            let y_new = &x - &yq;

            let z_new = std::mem::replace(&mut z, w);
            let w_new = &z - &wq;
            a = b;
            b = r;
            x = x_new;
            y = y_new;
            z = z_new;
            w = w_new;
        }
        (x, z, a)
    }

    pub fn gcd(a: Self, b: Self) -> Self {
        Self::ext_gcd(a, b).2
    }

    pub fn pow(self, n: u32) -> Self {
        match n {
            0 => Self::from_int(IBig::ONE),
            1 => self,
            _ => {
                let mut result = Self::from_int(IBig::ONE);
                let mut base = self;
                let mut exp = n;

                while exp > 0 {
                    if exp & 1 == 1 {
                        result = &result * &base;
                    }
                    base = &base * &base;
                    exp >>= 1;
                }
                result
            }
        }
    }
}

// impl Add for ZOmega {
//     type Output = Self;
//     fn add(self, rhs: Self) -> Self {
//         Self::new(
//             self.a + rhs.a,
//             self.b + rhs.b,
//             self.c + rhs.c,
//             self.d + rhs.d,
//         )
//     }
// }

impl Add for &ZOmega {
    type Output = ZOmega;
    fn add(self, rhs: Self) -> ZOmega {
        ZOmega::new(
            &self.a + &rhs.a,
            &self.b + &rhs.b,
            &self.c + &rhs.c,
            &self.d + &rhs.d,
        )
    }
}

// impl Add<&ZOmega> for &IBig {
//     type Output = ZOmega;
//     fn add(self, rhs: &ZOmega) -> ZOmega {
//         &ZOmega::from_int(self) + &rhs
//     }
// }

impl Add<ZOmega> for IBig {
    type Output = ZOmega;
    fn add(self, rhs: ZOmega) -> ZOmega {
        &ZOmega::from_int(self) + &rhs
    }
}

impl Sub<&ZOmega> for &ZOmega {
    type Output = ZOmega;

    fn sub(self, rhs: &ZOmega) -> ZOmega {
        ZOmega::new(
            &self.a - &rhs.a,
            &self.b - &rhs.b,
            &self.c - &rhs.c,
            &self.d - &rhs.d,
        )
    }
}

// Never used apparently
// impl Sub<IBig> for ZOmega {
//     type Output = Self;
//     fn sub(self, rhs: IBig) -> Self {
//         self - Self::from_int(rhs)
//     }
// }

// Never used apparently
// impl Sub<&ZOmega> for IBig {
//     type Output = ZOmega;
//     fn sub(self, rhs: &ZOmega) -> ZOmega {
//         &ZOmega::from_int(self) - rhs
//     }
// }

impl Neg for ZOmega {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.a, -self.b, -self.c, -self.d)
    }
}

impl Mul<&ZOmega> for &ZOmega {
    type Output = ZOmega;

    fn mul(self, rhs: &ZOmega) -> ZOmega {
        let lhs_coef = self.coef();
        let rhs_coef = rhs.coef();

        let mut new_coef = [IBig::ZERO; 4];

        for i in 0..4 {
            for j in 0..4 {
                let product = &rhs_coef[i] * &lhs_coef[j];
                if i + j < 4 {
                    new_coef[i + j] += &product;
                } else {
                    new_coef[i + j - 4] -= &product;
                }
            }
        }

        ZOmega::new(
            std::mem::take(&mut new_coef[3]),
            std::mem::take(&mut new_coef[2]),
            std::mem::take(&mut new_coef[1]),
            std::mem::take(&mut new_coef[0]),
        )
    }
}

impl Mul<IBig> for ZOmega {
    type Output = Self;
    fn mul(self, rhs: IBig) -> Self {
        Self::new(self.a * &rhs, self.b * &rhs, self.c * &rhs, self.d * &rhs)
    }
}

impl Mul<ZOmega> for IBig {
    type Output = ZOmega;
    fn mul(self, rhs: ZOmega) -> ZOmega {
        rhs * self
    }
}

impl Div for ZOmega {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.divmod(&rhs).0
    }
}

impl Div<IBig> for ZOmega {
    type Output = Self;
    fn div(self, rhs: IBig) -> Self {
        self / ZOmega::from_int(rhs)
    }
}

impl Rem for ZOmega {
    type Output = Self;
    fn rem(self, rhs: Self) -> Self {
        self.divmod(&rhs).1
    }
}

impl Rem<ZRootTwo> for ZOmega {
    type Output = Self;
    fn rem(self, rhs: ZRootTwo) -> Self {
        self % ZOmega::from_zroottwo(&rhs)
    }
}

impl Rem<IBig> for ZOmega {
    type Output = Self;
    fn rem(self, rhs: IBig) -> Self {
        self % ZOmega::from_int(rhs)
    }
}

impl PartialEq for ZOmega {
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b && self.c == other.c && self.d == other.d
    }
}

impl PartialEq<IBig> for ZOmega {
    fn eq(&self, other: &IBig) -> bool {
        self == &ZOmega::from_int(other.clone())
    }
}

impl PartialEq<ZRootTwo> for ZOmega {
    fn eq(&self, other: &ZRootTwo) -> bool {
        self == &ZOmega::from_zroottwo(other)
    }
}

impl Debug for ZOmega {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "ZOmega({}, {}, {}, {})", self.a, self.b, self.c, self.d)
    }
}

impl Display for ZOmega {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}ω^3+{}ω^2+{}ω+{}", self.a, self.b, self.c, self.d)
    }
}

fn rounddiv(a: &IBig, b: &IBig) -> IBig {
    let ab = a * b;
    let bias = if ab >= IBig::ZERO { b } else { &-b };
    (a * 2 + bias) / (b * 2)
}
