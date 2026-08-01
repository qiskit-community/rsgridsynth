// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::ring::DOmega;
use dashu_float::{round::mode::HalfEven, FBig};
use dashu_int::IBig;
use nalgebra::Matrix2;
use num::Complex;
use once_cell::sync::OnceCell;
use std::fmt::{Debug, Display, Formatter, Result};

/// Represents a unitary matrix over D[ω] (dyadic cyclotomic integers).
///
/// A unitary matrix of the form:
/// ```text
/// [ z,  -conj(w) * ω^n ]
/// [ w,   conj(z) * ω^n ]
/// ```
/// where z, w ∈ D[ω] and ω = exp(iπ/4).
///
/// This representation is used in the GridSynth algorithm for Clifford+T synthesis.
#[derive(Clone, PartialEq)]
pub struct DOmegaUnitary {
    pub(crate) z: DOmega,
    pub(crate) w: DOmega,
    pub(crate) n: u8,
    to_matrix_cache: OnceCell<[[DOmega; 2]; 2]>,
}

impl DOmegaUnitary {
    /// Creates a new DOmega unitary matrix.
    ///
    /// # Arguments
    ///
    /// * `z` - The (0,0) entry of the matrix
    /// * `w` - The (1,0) entry of the matrix
    /// * `n` - The power of ω (reduced mod 8)
    /// * `k` - Optional denominator exponent; if provided, both z and w are adjusted to this exponent
    ///
    /// # Returns
    ///
    /// A new `DOmegaUnitary` with normalized denominator exponents.
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

    /// Returns the denominator exponent k (where denominator is 2^k).
    ///
    /// # Returns
    ///
    /// The denominator exponent of the matrix entries.
    pub fn k(&self) -> i64 {
        self.w.k
    }

    /// Converts the unitary to a 2×2 matrix representation.
    ///
    /// # Returns
    ///
    /// A reference to the 2×2 matrix with entries in D[ω], cached for performance.
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

    /// Converts the unitary to a complex matrix with high-precision floating-point entries.
    ///
    /// # Returns
    ///
    /// A 2×2 nalgebra matrix with complex entries, useful for numerical verification.
    pub fn to_complex_matrix(&self) -> Matrix2<Complex<FBig<HalfEven>>> {
        let mat = self.to_matrix();
        Matrix2::new(
            Complex::new(mat[0][0].real().clone(), mat[0][0].imag().clone()),
            Complex::new(mat[0][1].real().clone(), mat[0][1].imag().clone()),
            Complex::new(mat[1][0].real().clone(), mat[1][0].imag().clone()),
            Complex::new(mat[1][1].real().clone(), mat[1][1].imag().clone()),
        )
    }

    /// Multiplies the unitary by the T gate from the left.
    ///
    /// # Returns
    ///
    /// A new unitary representing T × self.
    pub fn mul_by_t_from_left(&self) -> Self {
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega(),
            self.n as usize + 1,
            None,
        )
    }

    /// Multiplies the unitary by T† (T-inverse) from the left.
    ///
    /// # Returns
    ///
    /// A new unitary representing T† × self.
    pub fn mul_by_t_inv_from_left(&self) -> Self {
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega_inv(),
            self.n as usize - 1,
            None,
        )
    }

    /// Multiplies the unitary by T^m from the left.
    ///
    /// # Arguments
    ///
    /// * `m` - The power of T (reduced mod 8)
    ///
    /// # Returns
    ///
    /// A new unitary representing T^m × self.
    pub fn mul_by_t_power_from_left(&self, m: i32) -> Self {
        let m = m & 0b111;
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega_power(m.try_into().unwrap()),
            (self.n as i32 + m) as usize,
            None,
        )
    }

    /// Multiplies the unitary by the S gate from the left.
    ///
    /// # Returns
    ///
    /// A new unitary representing S × self.
    pub fn mul_by_s_from_left(&self) -> Self {
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega_power(2),
            self.n as usize + 2,
            None,
        )
    }

    /// Multiplies the unitary by S^m from the left.
    ///
    /// # Arguments
    ///
    /// * `m` - The power of S (reduced mod 4)
    ///
    /// # Returns
    ///
    /// A new unitary representing S^m × self.
    pub fn mul_by_s_power_from_left(&self, m: i32) -> Self {
        let m = m & 0b11;
        Self::new(
            self.z.clone(),
            self.w.mul_by_omega_power(((m << 1) as u8).into()),
            (self.n as i32 + (m << 1)) as usize,
            None,
        )
    }

    /// Multiplies the unitary by the Hadamard gate from the left.
    ///
    /// # Returns
    ///
    /// A new unitary representing H × self.
    pub fn mul_by_h_from_left(&self) -> Self {
        let new_z = (self.z.clone() + self.w.clone()).mul_by_inv_sqrt2();
        let new_w = (self.z.clone() - self.w.clone()).mul_by_inv_sqrt2();
        Self::new(new_z, new_w, self.n as usize + 4, None)
    }

    /// Multiplies the unitary by T^m followed by H from the left.
    ///
    /// # Arguments
    ///
    /// * `m` - The power of T
    ///
    /// # Returns
    ///
    /// A new unitary representing H × T^m × self.
    pub fn mul_by_h_and_t_power_from_left(&self, m: i32) -> Self {
        self.mul_by_t_power_from_left(m).mul_by_h_from_left()
    }

    /// Multiplies the unitary by the Pauli-X gate from the left.
    ///
    /// # Returns
    ///
    /// A new unitary representing X × self.
    pub fn mul_by_x_from_left(&self) -> Self {
        Self::new(self.w.clone(), self.z.clone(), self.n as usize + 4, None)
    }

    /// Multiplies the unitary by the W gate (ω gate) from the left.
    ///
    /// # Returns
    ///
    /// A new unitary representing W × self, where W = diag(ω, ω).
    pub fn mul_by_w_from_left(&self) -> Self {
        Self::new(
            self.z.mul_by_omega(),
            self.w.mul_by_omega(),
            self.n as usize + 2,
            None,
        )
    }

    /// Multiplies the unitary by W^m from the left.
    ///
    /// # Arguments
    ///
    /// * `m` - The power of W (reduced mod 8)
    ///
    /// # Returns
    ///
    /// A new unitary representing W^m × self.
    pub fn mul_by_w_power_from_left(&self, m: i32) -> Self {
        let m = m & 0b111;
        Self::new(
            self.z.mul_by_omega_power((m as u8).into()),
            self.w.mul_by_omega_power((m as u8).into()),
            self.n as usize + (m << 1) as usize,
            None,
        )
    }

    /// Adjusts the denominator exponent to a new value.
    ///
    /// # Arguments
    ///
    /// * `new_k` - The new denominator exponent
    ///
    /// # Returns
    ///
    /// A new unitary with the same matrix but adjusted denominator.
    pub fn renew_denomexp(&self, new_k: i64) -> Self {
        Self::new(self.z.clone(), self.w.clone(), self.n as usize, Some(new_k))
    }

    /// Reduces the denominator exponent to its minimal value.
    ///
    /// # Returns
    ///
    /// A new unitary with the same matrix but minimal denominator exponent.
    pub fn reduce_denomexp(&self) -> Self {
        Self::new(
            self.z.reduce_denomexp(),
            self.w.reduce_denomexp(),
            self.n as usize,
            None,
        )
    }

    /// Creates the identity unitary matrix.
    ///
    /// # Returns
    ///
    /// A `DOmegaUnitary` representing the 2×2 identity matrix.
    pub fn identity() -> Self {
        Self::new(
            DOmega::from_int(IBig::ONE),
            DOmega::from_int(IBig::ZERO),
            0,
            None,
        )
    }

    /// Constructs a unitary from a gate sequence string.
    ///
    /// # Arguments
    ///
    /// * `gates` - A string of gate characters (I, H, T, S, X, W)
    ///
    /// # Returns
    ///
    /// A `DOmegaUnitary` representing the product of the gates.
    ///
    /// # Panics
    ///
    /// Panics if an unsupported gate character is encountered.
    pub fn from_gates(gates: &str) -> Self {
        let mut unitary = Self::identity();
        for g in gates.chars().rev() {
            unitary = match g {
                'I' => unitary,
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
