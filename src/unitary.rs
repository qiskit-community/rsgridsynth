// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::ring::DOmega;
use dashu_int::IBig;
use once_cell::sync::OnceCell;
use std::fmt::{Debug, Display, Formatter, Result};

#[derive(Clone, PartialEq)]
pub struct DOmegaUnitary {
    pub(crate) z: DOmega,
    pub(crate) w: DOmega,
    pub(crate) n: u8,
    to_matrix_cache: OnceCell<[[DOmega; 2]; 2]>,
}

impl DOmegaUnitary {
    pub fn new(mut z: DOmega, mut w: DOmega, n: usize, k: Option<i64>) -> Self {
        let n = (n & 0b111) as u8;
        match k {
            Some(k_val) => {
                z = z.renew_denomexp(k_val);
                w = w.renew_denomexp(k_val);
            }
            None => match z.k > w.k {
                true => w = w.renew_denomexp(z.k),
                false if z.k < w.k => z = z.renew_denomexp(w.k),
                _ => {}
            },
        }
        Self {
            z,
            w,
            n,
            to_matrix_cache: OnceCell::new(),
        }
    }

    pub fn k(&self) -> i64 {
        self.w.k
    }

    pub fn to_matrix(&self) -> &[[DOmega; 2]; 2] {
        self.to_matrix_cache.get_or_init(|| {
            [
                [
                    self.z.clone(),
                    -self.w.conj().mul_by_omega_power(self.n.into()),
                ],
                [
                    self.w.clone(),
                    self.z.conj().mul_by_omega_power(self.n.into()),
                ],
            ]
        })
    }

    pub fn mul_by_t_from_left(&self) -> Self {
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega(),
            self.n as usize + 1,
            None,
        )
    }

    pub fn mul_by_t_inv_from_left(&self) -> Self {
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega_inv(),
            self.n as usize - 1,
            None,
        )
    }

    pub fn mul_by_t_power_from_left(&self, m: i32) -> Self {
        let m = m & 0b111;
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega_power(m.try_into().unwrap()),
            (self.n as i32 + m) as usize,
            None,
        )
    }

    pub fn mul_by_s_from_left(&self) -> Self {
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega_power(2),
            self.n as usize + 2,
            None,
        )
    }

    pub fn mul_by_s_power_from_left(&self, m: i32) -> Self {
        let m = m & 0b11;
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega_power(((m << 1) as u8).into()),
            (self.n as i32 + (m << 1)) as usize,
            None,
        )
    }

    pub fn mul_by_h_from_left(&self) -> Self {
        let new_z = (self.z.clone() + self.w.clone()).mul_by_inv_sqrt2();
        let new_w = (self.z.clone() - self.w.clone()).mul_by_inv_sqrt2();
        Self::new(new_z, new_w, self.n as usize + 4, None)
    }

    pub fn mul_by_h_and_t_power_from_left(&self, m: i32) -> Self {
        self.mul_by_t_power_from_left(m).mul_by_h_from_left()
    }

    pub fn mul_by_x_from_left(&self) -> Self {
        Self::new(self.w.clone(), self.z.clone(), self.n as usize + 4, None)
    }

    pub fn mul_by_w_from_left(&self) -> Self {
        Self::new(
            self.z.mul_by_omega(),
            self.w.mul_by_omega(),
            self.n as usize + 2,
            None,
        )
    }

    pub fn mul_by_w_power_from_left(&self, m: i32) -> Self {
        let m = m & 0b111;
        Self::new(
            self.z.mul_by_omega_power((m as u8).into()),
            self.w.mul_by_omega_power((m as u8).into()),
            self.n as usize + (m << 1) as usize,
            None,
        )
    }

    pub fn renew_denomexp(&self, new_k: i64) -> Self {
        Self::new(self.z.clone(), self.w.clone(), self.n as usize, Some(new_k))
    }

    pub fn reduce_denomexp(&self) -> Self {
        Self::new(
            self.z.reduce_denomexp(),
            self.w.reduce_denomexp(),
            self.n as usize,
            None,
        )
    }

    pub fn identity() -> Self {
        Self::new(
            DOmega::from_int(IBig::ONE),
            DOmega::from_int(IBig::ZERO),
            0,
            None,
        )
    }

    pub fn from_gates(gates: &str) -> Self {
        let mut unitary = Self::identity();
        for g in gates.chars().rev() {
            unitary = match g {
                'H' => unitary.renew_denomexp(unitary.k() + 1).mul_by_h_from_left(),
                'T' => unitary.mul_by_t_from_left(),
                'S' => unitary.mul_by_s_from_left(),
                'X' => unitary.mul_by_x_from_left(),
                'W' => unitary.mul_by_w_from_left(),
                _ => panic!("Unsupported gate: {}", g),
            };
        }
        unitary.reduce_denomexp()
    }
}

impl Display for DOmegaUnitary {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{:?}", self.to_matrix())
    }
}

impl Debug for DOmegaUnitary {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(
            f,
            "DOmegaUnitary({:?}, {:?}, n: {})",
            self.z, self.w, self.n
        )
    }
}
