// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use std::fmt::{Debug, Display, Formatter, Result};
use std::ops::Mul;

use crate::gate::{Gate, GateSeq};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    I = 0,
    H = 1,
    SH = 2,
}

impl Axis {
    /// Gate expansion of this coset representative, in left-to-right (application) order.
    ///
    /// This is the explicit replacement for relying on `format!("{:?}", axis)` spelling a valid
    /// gate string (`Axis::SH`'s `Debug` output happens to be `"SH"`) -- a coupling that would
    /// silently corrupt output if a variant were ever renamed.
    pub(crate) const fn gates(self) -> &'static [Gate] {
        match self {
            Axis::I => &[],
            Axis::H => &[Gate::H],
            Axis::SH => &[Gate::S, Gate::H],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Syllable {
    I = 0,
    T = 1,
    HT = 2,
    SHT = 3,
}

impl Syllable {
    /// Gate expansion of this syllable. Every non-`I` variant contains exactly one `Gate::T` --
    /// the invariant `GateSeq::t_count`/`NormalForm::t_count` rely on agreeing with each other.
    /// See [`Axis::gates`] for why this replaces `Debug`-formatting.
    pub(crate) const fn gates(self) -> &'static [Gate] {
        match self {
            Syllable::I => &[],
            Syllable::T => &[Gate::T],
            Syllable::HT => &[Gate::H, Gate::T],
            Syllable::SHT => &[Gate::S, Gate::H, Gate::T],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Clifford {
    a: u8,
    b: u8,
    c: u8,
    d: u8,
}

impl Clifford {
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

    pub fn inv(&self) -> Self {
        let (a, b, c, d) = CINV_TABLE[((self.a << 3) | (self.b << 2) | self.c) as usize];
        Clifford::new(a as i32, b as i32, c as i32, d as i32 - self.d as i32)
    }

    pub fn decompose_coset(&self) -> (Axis, Self) {
        match self.a {
            0 => (Axis::I, *self),
            1 => (Axis::H, CLIFFORD_H.inv() * *self),
            2 => (Axis::SH, (CLIFFORD_S * CLIFFORD_H).inv() * *self),
            _ => unreachable!(),
        }
    }

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

    pub fn to_gates(&self) -> GateSeq {
        let (axis, c) = self.decompose_coset();
        let mut gates =
            GateSeq::with_capacity(axis.gates().len() + c.b as usize + c.c as usize + c.d as usize);
        gates.extend(axis.gates());
        gates.push_n(Gate::X, c.b as usize);
        gates.push_n(Gate::S, c.c as usize);
        gates.push_n(Gate::W, c.d as usize);
        gates
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

#[derive(Clone)]
pub struct NormalForm {
    syllables: Vec<Syllable>,
    c: Clifford,
}

impl NormalForm {
    pub fn new() -> Self {
        Self {
            syllables: vec![],
            c: CLIFFORD_I,
        }
    }

    fn append_gate(&mut self, g: Gate) {
        match g {
            Gate::H => self.c = self.c * CLIFFORD_H,
            Gate::S => self.c = self.c * CLIFFORD_S,
            Gate::X => self.c = self.c * CLIFFORD_X,
            Gate::W => self.c = self.c * CLIFFORD_W,
            Gate::T => {
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
        }
    }

    pub fn from_gates(gates: &[Gate]) -> Self {
        let mut nf = Self::new();
        for &g in gates {
            nf.append_gate(g);
        }
        nf
    }

    pub fn to_gates(&self) -> GateSeq {
        let mut gates = GateSeq::with_capacity(3 * self.syllables.len() + 8);
        for s in &self.syllables {
            gates.extend(s.gates());
        }
        gates.extend(self.c.to_gates());
        gates
    }

    /// Number of T-gate applications in the Matsumoto-Amano normal form. Every syllable
    /// other than `I` (i.e. `T`, `HT`, `SHT`) is a Clifford conjugation of exactly one T,
    /// so this is simply the count of non-`I` syllables.
    pub fn t_count(&self) -> usize {
        self.syllables.iter().filter(|s| **s != Syllable::I).count()
    }

    /// The T-rotation syllables, in application order, matching the grammar `T?(HT|SHT)*`:
    /// a bare `Syllable::T` can only occur as the first element (`append_gate` merges any
    /// later `Axis::I` case into the preceding syllable instead of pushing a second one).
    pub fn syllables(&self) -> &[Syllable] {
        &self.syllables
    }

    /// The trailing Clifford correction (`Clifford` is `Copy`).
    pub fn clifford(&self) -> Clifford {
        self.c
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

#[cfg(test)]
mod tests {
    use super::*;

    fn manual_t_count(gates: &str) -> usize {
        gates.chars().filter(|&c| c == 'T').count()
    }

    fn seq(s: &str) -> GateSeq {
        s.parse().expect("test gate strings are valid gate words")
    }

    // Pins `Axis::gates`/`Syllable::gates` to the same spelling the old `format!("{:?}", _)`
    // trick produced, so the table transcription is checked independently of the larger
    // GateSeq migration.
    #[test]
    fn axis_and_syllable_gate_tables_match_debug_formatting() {
        for (axis, debug_str) in [(Axis::I, ""), (Axis::H, "H"), (Axis::SH, "SH")] {
            let via_table: String = axis.gates().iter().map(|g| g.as_char()).collect();
            assert_eq!(via_table, debug_str);
            if axis != Axis::I {
                assert_eq!(format!("{axis:?}"), debug_str);
            }
        }
        for (syllable, debug_str) in [
            (Syllable::I, ""),
            (Syllable::T, "T"),
            (Syllable::HT, "HT"),
            (Syllable::SHT, "SHT"),
        ] {
            let via_table: String = syllable.gates().iter().map(|g| g.as_char()).collect();
            assert_eq!(via_table, debug_str);
            if syllable != Syllable::I {
                assert_eq!(format!("{syllable:?}"), debug_str);
            }
        }
    }

    // A lone "T" does not trigger any Clifford-conjugation cancellation, so `t_count()`
    // should agree with a direct character count of the (un-normalized) input.
    #[test]
    fn t_count_matches_char_count_for_single_t() {
        let nf = NormalForm::from_gates(&seq("T"));
        assert_eq!(nf.t_count(), manual_t_count("T"));
    }

    // Repeated adjacent "T"s are NOT expected to preserve a raw T-count: e.g. T*T = S (a
    // Clifford, per T^2 = S), so `NormalForm::from_gates("TT")` legitimately normalizes to
    // zero T syllables. This is correct behavior of Matsumoto-Amano normal form, not a bug
    // in `t_count()` -- so here we instead verify `t_count()` is self-consistent with the
    // *normalized* output, i.e. equals a direct count of 'T' characters in
    // `nf.to_gates()`. `to_gates()` emits exactly one 'T' per non-`I` syllable (as "T",
    // "HT", or "SHT") and the trailing Clifford suffix never contains a 'T', so this must
    // hold for any normal form, and doubles as a lightweight check that `t_count()` and
    // `to_gates()` haven't drifted out of sync with each other.
    #[test]
    fn t_count_matches_normalized_char_count_for_repeated_t_strings() {
        for gates in ["", "T", "TT", "TTT", "TTTTTTTTTT"] {
            let nf = NormalForm::from_gates(&seq(gates));
            let normalized = nf.to_gates();
            assert_eq!(
                nf.t_count(),
                manual_t_count(&normalized.to_string()),
                "mismatch for input {gates:?}, normalized to {normalized}"
            );
        }
    }

    // Golden gate strings from tests/integration_test.rs. These are already the
    // (gate-only, i.e. no measurement/other non Clifford+T symbols) output of
    // `decompose_domega_unitary`, which itself already produces Matsumoto-Amano normal
    // form output. Re-normalizing normal-form output should be idempotent, so we expect
    // `t_count()` on these strings to equal a direct count of 'T' characters in them.
    const GOLDEN_GATES: &[&str] = &[
        "HTHTSHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTSHTHTSHTHTSHTXSSWWW",
        "HTSHTSHTSHTHTHTSHTHTHTSHTSHTHTHTHTHTHTSHTHTHTHTSHTSHTHTSHTHTHTHTHTHTHTHTSHTHTHTHTSHTHTSHTSHTSHTSHTHTSHTSHTHTSHTSHTHTSHTHTHTSHTSHTHTHTHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTHTSHTHTSHTSHTHTHTHTSHTHTHTSHTHTHTHTSHTHTSHTSHTHTSHTHTSHTHTHTHTHTHTHTHTSHTHTHTSHTSSSWW",
        "HTSHTHTHTSHTHTHTHTHTHTHTHTSHTHTSHTHTSHTSHTHTSHTHTHTHTSHTHTHTSHTHTHTHTSHTSHTHTSHTHTHTHTHTSHTSHTHTHTSHTHTSHTHTSHTHTHTHTSHTSHTHTHTSHTHTSHTSHTHTSHTSHTHTSHTSHTSHTSHTHTSHTHTHTHTSHTHTHTHTHTHTHTHTSHTHTSHTSHTHTHTHTSHTHTHTHTHTHTSHTSHTHTHTSHTHTHTSHTSHTSHTSSSWW",
        "SHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTHTSHTSHTSHTHTHTSHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTSHTHTHTHTHTHTSHTSHTHTHTSHTSHTHTHTHTSHTSHTHTSHTSHTSHTHTSHTHTHTSHTSHTHTSHTSHTSHTHTSHTHTHTSHTHTHTSHTSHTSHTSHTHTHTHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTHTHTHTSHTSHTHTSHTSHTSHTHTSHTSHTSHTHTSHTHTHTSHTHTHTSHTSHTHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTHTHTSHTSHTHTSHSWW",
        "TWWWWWWW",
    ];

    #[test]
    fn t_count_matches_char_count_for_golden_gate_strings() {
        for gates in GOLDEN_GATES {
            let nf = NormalForm::from_gates(&seq(gates));
            assert_eq!(
                nf.t_count(),
                manual_t_count(gates),
                "t_count mismatch for golden string {gates:?}"
            );
        }
    }

    /// The direct regression test for the fix: a trivial trailing Clifford correction
    /// contributes nothing (not the historical `"I"` sentinel).
    #[test]
    fn clifford_identity_to_gates_is_empty() {
        assert!(Clifford::new(0, 0, 0, 0).to_gates().is_empty());
    }

    /// Brute-force every input word over `{H, S, T, X, W}` up to length 5 (5^0 + .. + 5^5 =
    /// 19531 cases, fast) and check two invariants that together prove
    /// `NormalForm::to_gates` never produces output `NormalForm::from_gates` can't parse back
    /// (the class of bug this migration fixes, without needing to hand-construct the specific
    /// leaking case), and that the `t_count()`/`to_gates().t_count()` bridge invariant holds
    /// everywhere, not just on the hand-picked golden strings above.
    #[test]
    fn to_gates_is_idempotent_and_t_count_bridges_for_all_short_words() {
        fn recurse(prefix: &mut Vec<Gate>, depth: usize, max_depth: usize) {
            if depth > 0 {
                let seq = GateSeq::from(prefix.clone());
                let nf = NormalForm::from_gates(&seq);
                let rendered = nf.to_gates();

                // Idempotence: re-parsing and re-rendering the normal form's own output
                // reproduces it exactly -- in particular, `rendered` must contain no
                // unparseable embedded 'I', or this round-trip would silently drop gates
                // instead of matching.
                let round_tripped = NormalForm::from_gates(&rendered).to_gates();
                assert_eq!(
                    round_tripped, rendered,
                    "to_gates not idempotent for input {prefix:?} (rendered {rendered})"
                );

                // t_count bridge: structural count (non-`I` syllables) agrees with a count
                // of `Gate::T` in the rendered sequence.
                assert_eq!(
                    nf.t_count(),
                    rendered.t_count(),
                    "t_count bridge broken for input {prefix:?} (rendered {rendered})"
                );
            }
            if depth == max_depth {
                return;
            }
            for g in Gate::ALL {
                prefix.push(g);
                recurse(prefix, depth + 1, max_depth);
                prefix.pop();
            }
        }

        recurse(&mut Vec::new(), 0, 5);
    }
}
