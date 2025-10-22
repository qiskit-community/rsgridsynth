// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::{
    region::Ellipse,
    ring::{DOmega, ZOmega},
};
use dashu_float::{round::mode::HalfEven, FBig};
use dashu_int::IBig;
use nalgebra::{Matrix2, Vector2};
use num_traits::Pow;
use once_cell::unsync::OnceCell;
use std::{
    fmt::{Display, Formatter, Result},
    ops::{Mul, MulAssign},
};

#[derive(Debug, Clone)]
pub struct EllipsePair {
    pub a: Ellipse,
    pub b: Ellipse,
}

impl EllipsePair {
    pub fn new(a: Ellipse, b: Ellipse) -> Self {
        Self { a, b }
    }

    pub fn skew(&self) -> FBig<HalfEven> {
        self.a.skew() + self.b.skew()
    }

    pub fn bias(&self) -> FBig<HalfEven> {
        self.b.bias() / self.a.bias()
    }
}

impl Mul<EllipsePair> for GridOp {
    type Output = EllipsePair;
    fn mul(self, rhs: EllipsePair) -> EllipsePair {
        EllipsePair::new(&self * rhs.a, self.conj_sq2() * rhs.b)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GridOp {
    pub u0: ZOmega,
    pub u1: ZOmega,
    det_vec_cache: OnceCell<ZOmega>,
    adj_cache: OnceCell<Box<GridOp>>,
    conj_sq2_cache: OnceCell<Box<GridOp>>,
    inv_cache: OnceCell<Option<Box<GridOp>>>,
}

impl GridOp {
    pub fn new(u0: ZOmega, u1: ZOmega) -> Self {
        Self {
            u0,
            u1,
            det_vec_cache: OnceCell::new(),
            adj_cache: OnceCell::new(),
            conj_sq2_cache: OnceCell::new(),
            inv_cache: OnceCell::new(),
        }
    }

    pub fn a0(&self) -> &IBig {
        &self.u0.a
    }
    pub fn b0(&self) -> &IBig {
        &self.u0.b
    }
    pub fn c0(&self) -> &IBig {
        &self.u0.c
    }
    pub fn d0(&self) -> &IBig {
        &self.u0.d
    }
    pub fn a1(&self) -> &IBig {
        &self.u1.a
    }
    pub fn b1(&self) -> &IBig {
        &self.u1.b
    }
    pub fn c1(&self) -> &IBig {
        &self.u1.c
    }
    pub fn d1(&self) -> &IBig {
        &self.u1.d
    }

    pub fn det_vec(&self) -> ZOmega {
        self.det_vec_cache
            .get_or_init(|| self.u0.conj() * &self.u1)
            .clone()
    }

    pub fn is_special(&self) -> bool {
        let v = self.det_vec();
        v.a + v.c == IBig::ZERO && (v.b == IBig::ONE || v.b == IBig::NEG_ONE)
    }

    pub fn to_mat(&self) -> [[FBig<HalfEven>; 2]; 2] {
        [
            [self.u0.real(), self.u1.real()],
            [self.u0.imag(), self.u1.imag()],
        ]
    }

    pub fn inv(&self) -> Option<Self> {
        self.inv_cache
            .get_or_init(|| {
                if !self.is_special() {
                    return None;
                }
                let new_c0 = (self.c1() + self.a1() - self.c0() - self.a0()) / 2;
                let new_a0 = (-self.c1() - self.a1() - self.c0() - self.a0()) / 2;
                let mut new_u0 = ZOmega::new(new_a0, -self.b0(), new_c0, self.b1().clone());

                let new_c1 = (-self.c1() + self.a1() + self.c0() - self.a0()) / 2;
                let new_a1 = (self.c1() - self.a1() + self.c0() - self.a0()) / 2;
                let mut new_u1 = ZOmega::new(new_a1, self.d0().clone(), new_c1, -self.d1());

                if self.det_vec().b == IBig::NEG_ONE {
                    new_u0 = -new_u0;
                    new_u1 = -new_u1;
                }
                Some(Box::new(GridOp::new(new_u0, new_u1)))
            })
            .clone()
            .map(|b| *b)
    }

    pub fn adj(&self) -> &Self {
        self.adj_cache
            .get_or_init(|| {
                let new_c0 = (self.c1() - self.a1() + self.c0() - self.a0()) / 2;
                let new_a0 = (self.c1() - self.a1() - self.c0() + self.a0()) / 2;
                let new_u0 = ZOmega::new(new_a0, self.d0().clone(), new_c0, self.d1().clone());

                let new_c1 = (self.c1() + self.a1() + self.c0() + self.a0()) / 2;
                let new_a1 = (self.c1() + self.a1() - self.c0() - self.a0()) / 2;
                let new_u1 = ZOmega::new(new_a1, self.b1().clone(), new_c1, self.b0().clone());

                Box::new(GridOp::new(new_u0, new_u1))
            })
            .as_ref()
    }

    pub fn conj_sq2(&self) -> &Self {
        self.conj_sq2_cache
            .get_or_init(|| {
                Box::new(GridOp::new(
                    self.u0.conj_sq2().clone(),
                    self.u1.conj_sq2().clone(),
                ))
            })
            .as_ref()
    }
}

impl Display for GridOp {
    fn fmt(&self, f: &mut Formatter) -> Result {
        write!(
            f,
            "[[{}{:+}/√2, {}{:+}/√2],\n [{:+}{:+}/√2, {:+}{:+}/√2]]",
            self.d0(),
            self.c0() - self.a0(),
            self.d1(),
            self.c1() - self.a1(),
            self.b0(),
            self.c0() + self.a0(),
            self.b1(),
            self.c1() + self.a1()
        )
    }
}

impl Mul<GridOp> for GridOp {
    type Output = GridOp;
    fn mul(self, rhs: GridOp) -> Self::Output {
        GridOp::new(self.clone() * rhs.u0, self * rhs.u1)
    }
}

impl Mul<ZOmega> for GridOp {
    type Output = ZOmega;
    fn mul(self, other: ZOmega) -> ZOmega {
        let new_d = self.d0() * other.d.clone()
            + self.d1() * other.b.clone()
            + (self.c1() - self.a1() + self.c0() - self.a0()) / IBig::from(2) * other.c.clone()
            + (self.c1() - self.a1() - self.c0() + self.a0()) / IBig::from(2) * other.a.clone();
        let new_c = self.c0() * other.d.clone()
            + self.c1() * other.b.clone()
            + (self.b1() + self.d1() + self.b0() + self.d0()) / IBig::from(2) * other.c.clone()
            + (self.b1() + self.d1() - self.b0() - self.d0()) / IBig::from(2) * other.a.clone();
        let new_b = self.b0() * other.d.clone()
            + self.b1() * other.b.clone()
            + (self.c1() + self.a1() + self.c0() + self.a0()) / IBig::from(2) * other.c.clone()
            + (self.c1() + self.a1() - self.c0() - self.a0()) / IBig::from(2) * other.a.clone();
        let new_a = self.a0() * other.d.clone()
            + self.a1() * other.b.clone()
            + (self.b1() - self.d1() + self.b0() - self.d0()) / IBig::from(2) * other.c.clone()
            + (self.b1() - self.d1() - self.b0() + self.d0()) / IBig::from(2) * other.a.clone();

        ZOmega::new(new_a, new_b, new_c, new_d)
    }
}

impl Mul<DOmega> for GridOp {
    type Output = DOmega;
    fn mul(self, rhs: DOmega) -> DOmega {
        DOmega::new(self * rhs.u, rhs.k)
    }
}

impl MulAssign<GridOp> for GridOp {
    fn mul_assign(&mut self, rhs: GridOp) {
        *self = self.clone() * rhs;
    }
}

impl Pow<IBig> for GridOp {
    type Output = GridOp;
    fn pow(self, mut exp: IBig) -> GridOp {
        if exp < IBig::ZERO {
            return self.inv().unwrap().pow(-exp);
        }
        let mut result = GridOp::new(
            ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
            ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
        );
        let mut base = self;
        while exp.clone() > IBig::ZERO {
            if exp.clone() & IBig::ONE == IBig::ONE {
                result *= base.clone();
            }
            base *= base.clone();
            exp >>= 1;
        }
        result
    }
}

impl Pow<i32> for GridOp {
    type Output = GridOp;
    fn pow(self, mut exp: i32) -> GridOp {
        if exp < 0 {
            return self.inv().unwrap().pow(-exp);
        }
        let mut result = GridOp::new(
            ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
            ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
        );
        let mut base = self;
        while exp > 0 {
            if exp & 1 == 1 {
                result *= base.clone();
            }
            base *= base.clone();
            exp >>= 1;
        }
        result
    }
}

impl Mul<Ellipse> for &GridOp {
    type Output = Ellipse;
    fn mul(self, rhs: Ellipse) -> Ellipse {
        let inv = self.inv().expect("GridOp must be special").to_mat();
        let m00 = &inv[0][0];
        let m01 = &inv[0][1];
        let m10 = &inv[1][0];
        let m11 = &inv[1][1];

        let m00_sq = m00 * m00;
        let term1 = rhs.a() * &m00_sq;
        let m00_m10 = m00 * m10;
        let two_b = 2 * rhs.b();
        let term2 = two_b * m00_m10;
        let m10_sq = m10 * m10;
        let term3 = rhs.d() * m10_sq;
        let sum12_a = term1 + term2;
        let a = sum12_a + term3;

        let m00_m01 = m00 * m01;
        let term4 = rhs.a() * m00_m01;
        let m00_m11 = m00 * m11;
        let m01_m10 = m01 * m10;
        let sum_cross = m00_m11 + m01_m10;
        let term5 = rhs.b() * sum_cross;
        let m10_m11 = m10 * m11;
        let term6 = rhs.d() * m10_m11;
        let sum45_b = term4 + term5;
        let b = sum45_b + term6;

        let m01_sq = m01 * m01;
        let term7 = rhs.a() * m01_sq;
        let m11_m01 = m11 * m01;
        let two_b2 = 2 * rhs.b();
        let term8 = two_b2 * m11_m01;
        let m11_sq = m11 * m11;
        let term9 = rhs.d() * m11_sq;
        let sum78_d = term7 + term8;
        let d = sum78_d + term9;

        let new_d = Matrix2::new(a, b.clone(), b, d);

        let mat = self.to_mat();
        let mat00 = &mat[0][0];
        let mat01 = &mat[0][1];
        let mat10 = &mat[1][0];
        let mat11 = &mat[1][1];

        let px_term1 = mat00 * rhs.px();
        let px_term2 = mat01 * rhs.py();
        let px = px_term1 + px_term2;

        let py_term1 = mat10 * rhs.px();
        let py_term2 = mat11 * rhs.py();
        let py = py_term1 + py_term2;
        let new_p = Vector2::new(px, py);

        Ellipse::new(new_d, new_p)
    }
}
