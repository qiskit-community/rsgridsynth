// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use std::fmt::{Debug, Display, Formatter, Result};
use std::ops::Mul;

/// Represents the axis for Clifford gate decomposition.
///
/// Used to classify Clifford gates into cosets based on their action.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    I = 0,
    H = 1,
    SH = 2,
}

/// Represents a syllable in the normal form decomposition.
///
/// Syllables are the basic building blocks for representing T-gates
/// and their conjugates by Clifford gates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Syllable {
    I = 0,
    T = 1,
    HT = 2,
    SHT = 3,
}

/// Represents a Clifford gate in a compact form.
///
/// A Clifford gate is represented as E^a X^b S^c ω^d where:
/// - E is a specific Clifford generator
/// - X is the Pauli-X gate
/// - S is the phase gate
/// - ω is the 8th root of unity (exp(iπ/4))
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Clifford {
    a: u8,
    b: u8,
    c: u8,
    d: u8,
}

impl Clifford {
    /// Creates a new Clifford gate with the given exponents.
    ///
    /// # Arguments
    ///
    /// * `a` - Exponent for E (reduced mod 3)
    /// * `b` - Exponent for X (reduced mod 2)
    /// * `c` - Exponent for S (reduced mod 4)
    /// * `d` - Exponent for ω (reduced mod 8)
    ///
    /// # Returns
    ///
    /// A normalized Clifford gate.
    pub fn new(mut a: i32, mut b: i32, mut c: i32, mut d: i32) -> Self {
        a = a.rem_euclid(3);
        b &= 1;
        c &= 0b11;
        d &= 0b111;
        Self {
            a: a as u8,
            b: b as u8,
            c: c as u8,
            d: d as u8,
        }
    }

    /// Computes the inverse of this Clifford gate.
    ///
    /// # Returns
    ///
    /// The inverse Clifford gate such that self * self.inv() = I.
    pub fn inv(&self) -> Self {
        let (a, b, c, d) = CINV_TABLE[((self.a << 3) | (self.b << 2) | self.c) as usize];
        Clifford::new(a as i32, b as i32, c as i32, d as i32 - self.d as i32)
    }

    /// Decomposes the Clifford into an axis and a coset representative.
    ///
    /// # Returns
    ///
    /// A tuple `(axis, clifford)` where the axis indicates the coset
    /// and clifford is the representative within that coset.
    pub fn decompose_coset(&self) -> (Axis, Self) {
        match self.a {
            0 => (Axis::I, *self),
            1 => (Axis::H, CLIFFORD_H.inv() * *self),
            2 => (Axis::SH, (CLIFFORD_S * CLIFFORD_H).inv() * *self),
            _ => unreachable!(),
        }
    }

    /// Decomposes the Clifford for T-gate conjugation.
    ///
    /// # Returns
    ///
    /// A tuple `(axis, clifford)` used in the normal form algorithm
    /// to handle T-gate conjugation efficiently.
    pub fn decompose_tconj(&self) -> (Axis, Self) {
        let (axis, c, d) = TCONJ_TABLE[((self.a << 1) | self.b) as usize];
        (
            axis,
            Clifford::new(
                0,
                self.b as i32,
                self.c as i32 + c as i32,
                self.d as i32 + d as i32,
            ),
        )
    }

    /// Converts the Clifford gate to a gate sequence string.
    ///
    /// # Returns
    ///
    /// A string representation using gates H, X, S, W, or "I" for identity.
    /// Converts the normal form back to a gate sequence string.
    ///
    /// # Returns
    ///
    /// A string representation of the circuit using syllables and Clifford gates,
    /// or "I" for the identity.
    pub fn to_gates(&self) -> String {
        let (axis, c) = self.decompose_coset();
        let mut gates = match axis {
            Axis::I => String::new(),
            _ => format!("{:?}", axis),
        };
        gates += &"X".repeat(c.b as usize);
        gates += &"S".repeat(c.c as usize);
        gates += &"W".repeat(c.d as usize);
        if gates.is_empty() {
            "I".to_string()
        } else {
            gates
        }
    }

    // pub fn from_str(g: &str) -> Self {
    //     match g {
    //         "H" => CLIFFORD_H,
    //         "S" => CLIFFORD_S,
    //         "X" => CLIFFORD_X,
    //         "W" => CLIFFORD_W,
    //         _ => panic!("Invalid gate"),
    //     }
    // }
}

impl Display for Clifford {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "E^{} X^{} S^{} ω^{}", self.a, self.b, self.c, self.d)
    }
}

impl Mul for Clifford {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let (a1, b1, c1, d1) = CONJ3_TABLE[((rhs.a << 3) | (self.b << 2) | self.c) as usize];
        let (c2, d2) = CONJ2_TABLE[((c1 << 1) | rhs.b) as usize];
        Clifford::new(
            self.a as i32 + a1 as i32,
            b1 as i32 + rhs.b as i32,
            c2 as i32 + rhs.c as i32,
            d1 as i32 + d2 as i32 + self.d as i32 + rhs.d as i32,
        )
    }
}

/// Represents a quantum circuit in normal form.
///
/// The normal form is a canonical representation of Clifford+T circuits
/// that facilitates optimization and analysis. It consists of a sequence
/// of syllables (T-gates and their Clifford conjugates) followed by a
/// final Clifford gate.
#[derive(Clone)]
pub struct NormalForm {
    syllables: Vec<Syllable>,
    c: Clifford,
}

impl NormalForm {
    /// Creates a new empty normal form (identity circuit).
    ///
    /// # Returns
    ///
    /// A `NormalForm` representing the identity operation.
    pub fn new() -> Self {
        Self {
            syllables: vec![],
            c: CLIFFORD_I,
        }
    }

    /// Appends a gate to the normal form, maintaining the canonical structure.
    ///
    /// # Arguments
    ///
    /// * `g` - The gate to append ("H", "S", "X", "W", or "T")
    ///
    /// # Panics
    ///
    /// Panics if an unsupported gate is provided.
    fn append_gate(&mut self, g: &str) {
        match g {
            "H" => self.c = self.c * CLIFFORD_H,
            "S" => self.c = self.c * CLIFFORD_S,
            "X" => self.c = self.c * CLIFFORD_X,
            "W" => self.c = self.c * CLIFFORD_W,
            "T" => {
                let (axis, new_c) = self.c.decompose_tconj();
                match axis {
                    Axis::I => {
                        if let Some(last) = self.syllables.last_mut() {
                            match last {
                                Syllable::T => {
                                    self.syllables.pop();
                                    self.c = CLIFFORD_S * new_c;
                                    return;
                                }
                                Syllable::HT => {
                                    self.syllables.pop();
                                    self.c = (CLIFFORD_H * CLIFFORD_S) * new_c;
                                    return;
                                }
                                Syllable::SHT => {
                                    self.syllables.pop();
                                    self.c = (CLIFFORD_H * CLIFFORD_S * CLIFFORD_H) * new_c;
                                    return;
                                }
                                _ => {}
                            }
                        }
                        self.syllables.push(Syllable::T);
                        self.c = new_c;
                    }
                    Axis::H => {
                        self.syllables.push(Syllable::HT);
                        self.c = new_c;
                    }
                    Axis::SH => {
                        self.syllables.push(Syllable::SHT);
                        self.c = new_c;
                    }
                }
            }
            _ => panic!("Unsupported gate"),
        }
    }

    /// Constructs a normal form from a gate sequence string.
    ///
    /// # Arguments
    ///
    /// * `gates` - A string of gate characters (H, S, X, W, T)
    ///
    /// # Returns
    ///
    /// A `NormalForm` representing the circuit in canonical form.
    pub fn from_gates(gates: &str) -> Self {
        let mut nf = Self::new();
        for ch in gates.chars() {
            nf.append_gate(&ch.to_string());
        }
        nf
    }

    pub fn to_gates(&self) -> String {
        let mut gates = String::new();
        for s in &self.syllables {
            if *s != Syllable::I {
                gates += &format!("{:?}", s);
            }
        }
        gates += &self.c.to_gates();
        if gates.is_empty() {
            "I".to_string()
        } else {
            gates
        }
    }
}

impl Default for NormalForm {
    fn default() -> Self {
        Self::new()
    }
}

impl Display for NormalForm {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(
            f,
            "NormalForm: {} | {}",
            self.syllables
                .iter()
                .map(|s| format!("{:?}", s))
                .collect::<Vec<_>>()
                .join(" "),
            self.c
        )
    }
}

impl Debug for NormalForm {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "NormalForm({:?}, {:?})", self.syllables, self.c)
    }
}

// Predefined Cliffords
const CLIFFORD_I: Clifford = Clifford {
    a: 0,
    b: 0,
    c: 0,
    d: 0,
};
const CLIFFORD_X: Clifford = Clifford {
    a: 0,
    b: 1,
    c: 0,
    d: 0,
};
const CLIFFORD_S: Clifford = Clifford {
    a: 0,
    b: 0,
    c: 1,
    d: 0,
};
const CLIFFORD_W: Clifford = Clifford {
    a: 0,
    b: 0,
    c: 0,
    d: 1,
};
const CLIFFORD_H: Clifford = Clifford {
    a: 1,
    b: 0,
    c: 1,
    d: 5,
};
// const CLIFFORD_SH: Clifford = Clifford { a: 1, b: 0, c: 2, d: 5 };
// const CLIFFORD_HS: Clifford = Clifford { a: 1, b: 0, c: 2, d: 1 };
// const CLIFFORD_SHS: Clifford = Clifford { a: 1, b: 0, c: 0, d: 1 };

// Lookup tables
const CONJ2_TABLE: [(u8, u8); 8] = [
    (0, 0),
    (0, 0),
    (1, 0),
    (3, 2),
    (2, 0),
    (2, 4),
    (3, 0),
    (1, 6),
];
const CONJ3_TABLE: [(u8, u8, u8, u8); 24] = [
    (0, 0, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 2, 0),
    (0, 0, 3, 0),
    (0, 1, 0, 0),
    (0, 1, 1, 0),
    (0, 1, 2, 0),
    (0, 1, 3, 0),
    (1, 0, 0, 0),
    (2, 0, 3, 6),
    (1, 1, 2, 2),
    (2, 1, 3, 6),
    (1, 0, 2, 0),
    (2, 1, 1, 0),
    (1, 1, 0, 6),
    (2, 0, 1, 4),
    (2, 0, 0, 0),
    (1, 1, 3, 4),
    (2, 1, 0, 0),
    (1, 0, 1, 2),
    (2, 1, 2, 2),
    (1, 1, 1, 0),
    (2, 0, 2, 6),
    (1, 0, 3, 2),
];
const CINV_TABLE: [(u8, u8, u8, u8); 24] = [
    (0, 0, 0, 0),
    (0, 0, 3, 0),
    (0, 0, 2, 0),
    (0, 0, 1, 0),
    (0, 1, 0, 0),
    (0, 1, 1, 6),
    (0, 1, 2, 4),
    (0, 1, 3, 2),
    (2, 0, 0, 0),
    (1, 0, 1, 2),
    (2, 1, 0, 0),
    (1, 1, 3, 4),
    (2, 1, 1, 2),
    (1, 1, 1, 6),
    (2, 0, 2, 2),
    (1, 0, 3, 4),
    (1, 0, 0, 0),
    (2, 1, 3, 6),
    (1, 1, 2, 2),
    (2, 0, 3, 6),
    (1, 0, 2, 0),
    (2, 1, 1, 6),
    (1, 1, 0, 2),
    (2, 0, 1, 6),
];
const TCONJ_TABLE: [(Axis, u8, u8); 6] = [
    (Axis::I, 0, 0),
    (Axis::I, 1, 7),
    (Axis::H, 3, 3),
    (Axis::H, 2, 0),
    (Axis::SH, 0, 5),
    (Axis::SH, 1, 4),
];
