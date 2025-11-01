// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use dashu_base::Sign;
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;

use crate::math::{ntz, pow_sqrt2};
use crate::ring::z_root_two::ZRootTwo;
use crate::ring::DOmega;
use std::cmp::Ordering;
use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DRootTwo {
    pub(crate) alpha: ZRootTwo,
    pub(crate) k: i64,
}

pub const ONE: DRootTwo = DRootTwo {
    alpha: ZRootTwo {
        a: IBig::ONE,
        b: IBig::ZERO,
    },
    k: 0,
};

pub const DELTA_SQUARED: DRootTwo = DRootTwo {
    alpha: ZRootTwo {
        a: IBig::from_parts_const(Sign::Positive, 2),
        b: IBig::ONE,
    },
    k: 0,
};

pub const DELTA_SQUARED_M: DRootTwo = DRootTwo {
    alpha: ZRootTwo {
        a: IBig::from_parts_const(Sign::Positive, 2),
        b: IBig::NEG_ONE,
    },
    k: 0,
};

impl DRootTwo {
    pub fn new(alpha: ZRootTwo, k: i64) -> Self {
        Self { alpha, k }
    }

    pub fn from_int(x: IBig) -> Self {
        Self::new(ZRootTwo::from_int(x), 0)
    }

    pub fn from_zroottwo(x: ZRootTwo) -> Self {
        Self::new(x, 0)
    }

    pub fn from_domega(x: DOmega) -> Self {
        Self::new(ZRootTwo::from_zomega(x.u), x.k)
    }

    pub fn parity(&self) -> IBig {
        self.alpha.parity()
    }

    pub fn scale(&self) -> FBig<HalfEven> {
        pow_sqrt2(self.k)
    }

    pub fn squared_scale(&self) -> IBig {
        let k: usize = self.k.try_into().expect("k must fit in i128");
        IBig::ONE << k
    }

    pub fn to_real(&self) -> FBig<HalfEven> {
        self.alpha.to_real() / self.scale()
    }

    pub fn conj_sq2(&self) -> Self {
        if self.k & 1 == 1 {
            Self::new(-self.alpha.conj_sq2(), self.k)
        } else {
            Self::new(self.alpha.conj_sq2(), self.k)
        }
    }

    pub fn renew_denomexp(&self, new_k: i64) -> Self {
        let d = new_k - self.k;
        let new_alpha = self.mul_by_sqrt2_power(d).alpha;
        Self::new(new_alpha, new_k)
    }

    pub fn reduce_denomexp(&self) -> Self {
        let k_a = if self.alpha.a == IBig::ZERO {
            self.k
        } else {
            ntz(&self.alpha.a)
        };
        let k_b = if self.alpha.b == IBig::ZERO {
            self.k
        } else {
            ntz(&self.alpha.b)
        };
        let new_k = if k_a <= k_b {
            self.k - 2 * k_a
        } else {
            self.k - 2 * k_b - 1
        };
        self.renew_denomexp(new_k.max(0))
    }

    pub fn mul_by_inv_sqrt2(&self) -> Self {
        if &self.alpha.a & 1 == IBig::ZERO {
            let new_alpha = ZRootTwo::new(self.alpha.b.clone(), self.alpha.a.clone() >> 1);
            Self::new(new_alpha, self.k)
        } else {
            panic!("ValueError")
        }
    }

    pub fn mul_by_sqrt2_power(&self, d: i64) -> Self {
        if d < 0 {
            if d == -1 {
                return self.mul_by_inv_sqrt2();
            }
            let d_div_2 = -d >> 1;
            let d_mod_2 = -d & 1;
            let d_div_2_usize: usize = d_div_2.try_into().expect("k must fit in i128");
            if d_mod_2 == 0 {
                let bit = (IBig::ONE << d_div_2_usize) - IBig::ONE;
                if &self.alpha.a & &bit == IBig::ZERO && &self.alpha.b & &bit == IBig::ZERO {
                    let new_alpha = ZRootTwo::new(
                        &self.alpha.a >> d_div_2_usize,
                        &self.alpha.b >> d_div_2_usize,
                    );
                    return Self::new(new_alpha, self.k);
                }
            } else {
                let bit = (IBig::ONE << d_div_2_usize) - IBig::ONE;
                let bit2 = (IBig::ONE << (d_div_2_usize + 1)) - IBig::ONE;
                if &self.alpha.a & &bit2 == IBig::ZERO && &self.alpha.b & &bit == IBig::ZERO {
                    let new_alpha = ZRootTwo::new(
                        &self.alpha.b >> d_div_2_usize,
                        &self.alpha.a >> (d_div_2_usize + 1),
                    );
                    return Self::new(new_alpha, self.k);
                }
            }
            panic!("ValueError")
        } else {
            let d_div_2 = &d >> 1;
            let d_mod_2 = &d & 1;
            let d_div_2_usize: usize = d_div_2.try_into().expect("k must fit in i128");
            let mut new_alpha = &self.alpha * &ZRootTwo::from_int(IBig::ONE << d_div_2_usize);
            if d_mod_2 == 1 {
                new_alpha = new_alpha * ZRootTwo::new(IBig::ZERO, IBig::ONE);
            }
            Self::new(new_alpha, self.k)
        }
    }

    pub fn mul_by_sqrt2_power_renewing_denomexp(&self, k: i64) -> Self {
        if k > self.k {
            panic!("ValueError")
        }
        Self::new(self.alpha.clone(), self.k - k)
    }

    pub fn power_of_inv_sqrt2(k: i64) -> Self {
        if k < 0 {
            panic!("ValueError")
        }
        Self::new(ZRootTwo::new(IBig::ONE, IBig::ZERO), k)
    }
}

impl Ord for DRootTwo {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.k.cmp(&other.k) {
            Ordering::Less => {
                let lhs = self.renew_denomexp(other.k);
                lhs.cmp(other)
            }
            Ordering::Greater => {
                let rhs = other.renew_denomexp(self.k);
                self.cmp(&rhs)
            }
            Ordering::Equal => self.alpha.cmp(&other.alpha),
        }
    }
}

impl PartialOrd for DRootTwo {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Add for DRootTwo {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        (&self).add(&rhs)
    }
}

impl<'a> Add<&'a DRootTwo> for &DRootTwo {
    type Output = DRootTwo;

    fn add(self, rhs: &'a DRootTwo) -> DRootTwo {
        match self.k.cmp(&rhs.k) {
            Ordering::Less => {
                let lhs_up = self.renew_denomexp(rhs.k);
                DRootTwo::new(lhs_up.alpha + rhs.alpha.clone(), rhs.k)
            }
            Ordering::Greater => {
                let rhs_up = rhs.renew_denomexp(self.k);
                DRootTwo::new(self.alpha.clone() + rhs_up.alpha, self.k)
            }
            Ordering::Equal => DRootTwo::new(self.alpha.clone() + rhs.alpha.clone(), self.k),
        }
    }
}

impl Sub for DRootTwo {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<'b> std::ops::Sub<&'b DRootTwo> for &DRootTwo {
    type Output = DRootTwo;

    fn sub(self, rhs: &'b DRootTwo) -> DRootTwo {
        self + &(-rhs.clone())
    }
}

impl Mul for DRootTwo {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(self.alpha * rhs.alpha, self.k + rhs.k)
    }
}

impl Neg for DRootTwo {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self::new(-self.alpha, self.k)
    }
}

impl Neg for &DRootTwo {
    type Output = DRootTwo;
    fn neg(self) -> DRootTwo {
        DRootTwo::new(-self.alpha.clone(), self.k)
    }
}
