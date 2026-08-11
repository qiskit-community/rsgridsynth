// Copyright (c) 2025 IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

//! Typed Clifford+T gate sequences, replacing the crate's historical `String` gate
//! representation. [`Gate`] is one of the five alphabet gates; [`GateSeq`] is a sequence of them.
//!
//! [`GateSeq::is_empty`] is the semantic predicate ("is this the identity?"); [`Display`] is the
//! wire format. The two disagree at the empty sequence by design: the wire format has no
//! representation for a zero-length word, so an empty (identity) sequence renders as `"I"`, and
//! `"I"` is the *only* position `Display` ever emits it — never mid-sequence. [`FromStr`] accepts
//! `'I'` anywhere as a no-op, for backward compatibility with strings produced before this type
//! existed (some of which legitimately contained an embedded, unparseable `'I'`; see
//! `NormalForm::to_gates`).

use std::fmt::{self, Display, Formatter, Write};
use std::iter::FromIterator;
use std::ops::Deref;
use std::str::FromStr;

/// A single Clifford+T gate. There is no `I` (identity) variant: identity is represented as the
/// structural absence of gates, i.e. an empty [`GateSeq`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Gate {
    H,
    S,
    T,
    X,
    W,
}

impl Gate {
    pub const ALL: [Gate; 5] = [Gate::H, Gate::S, Gate::T, Gate::X, Gate::W];

    /// The single source of truth for the gate alphabet's char mapping.
    pub const fn as_char(self) -> char {
        match self {
            Gate::H => 'H',
            Gate::S => 'S',
            Gate::T => 'T',
            Gate::X => 'X',
            Gate::W => 'W',
        }
    }

    /// Rejects `'I'`: identity-tolerance is a [`GateSeq`]-level concern (see [`FromStr`] on
    /// `GateSeq`), not a single-gate one.
    pub const fn from_char(c: char) -> Option<Gate> {
        match c {
            'H' => Some(Gate::H),
            'S' => Some(Gate::S),
            'T' => Some(Gate::T),
            'X' => Some(Gate::X),
            'W' => Some(Gate::W),
            _ => None,
        }
    }
}

impl Display for Gate {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_char(self.as_char())
    }
}

impl TryFrom<char> for Gate {
    type Error = ParseGateError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        Gate::from_char(c).ok_or(ParseGateError { found: c })
    }
}

/// Returned when a character outside the gate alphabet `{H, S, T, X, W}` (and, for [`GateSeq`]
/// parsing, `I`) is encountered while parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParseGateError {
    pub found: char,
}

impl Display for ParseGateError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "unsupported gate character: {:?}", self.found)
    }
}

impl std::error::Error for ParseGateError {}

/// A sequence of [`Gate`]s. The empty sequence represents the identity.
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
pub struct GateSeq(Vec<Gate>);

impl GateSeq {
    pub const fn new() -> Self {
        Self(Vec::new())
    }

    /// Alias for [`GateSeq::new`] that reads better at call sites building up the identity.
    pub fn identity() -> Self {
        Self::new()
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    pub fn from_vec(v: Vec<Gate>) -> Self {
        Self(v)
    }

    pub fn into_vec(self) -> Vec<Gate> {
        self.0
    }

    pub fn as_slice(&self) -> &[Gate] {
        &self.0
    }

    /// Number of `Gate::T` occurrences in this sequence.
    ///
    /// Bridge invariant (see `Axis::gates`/`Syllable::gates` in `normal_form.rs`): for any
    /// `nf: NormalForm`, `nf.to_gates().t_count() == nf.t_count()`, because every non-`I`
    /// syllable expands to exactly one `T` and `Clifford::to_gates` never emits a `T`.
    pub fn t_count(&self) -> usize {
        self.0.iter().filter(|&&g| g == Gate::T).count()
    }

    pub fn push(&mut self, g: Gate) {
        self.0.push(g);
    }

    /// Pushes `g` `n` times. The typed replacement for `"X".repeat(n)`-style string building.
    pub fn push_n(&mut self, g: Gate, n: usize) {
        self.0.extend(std::iter::repeat_n(g, n));
    }

    pub fn pop(&mut self) -> Option<Gate> {
        self.0.pop()
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    pub fn reserve(&mut self, additional: usize) {
        self.0.reserve(additional);
    }
}

impl Deref for GateSeq {
    type Target = [Gate];

    fn deref(&self) -> &[Gate] {
        &self.0
    }
}

impl Display for GateSeq {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            return f.write_str("I");
        }
        for g in &self.0 {
            f.write_char(g.as_char())?;
        }
        Ok(())
    }
}

impl FromStr for GateSeq {
    type Err = ParseGateError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.chars()
            .filter(|&c| c != 'I')
            .map(Gate::try_from)
            .collect()
    }
}

impl From<Vec<Gate>> for GateSeq {
    fn from(v: Vec<Gate>) -> Self {
        Self(v)
    }
}

impl From<GateSeq> for Vec<Gate> {
    fn from(seq: GateSeq) -> Self {
        seq.0
    }
}

impl From<Gate> for GateSeq {
    fn from(g: Gate) -> Self {
        Self(vec![g])
    }
}

impl<const N: usize> From<[Gate; N]> for GateSeq {
    fn from(gates: [Gate; N]) -> Self {
        Self(gates.to_vec())
    }
}

impl FromIterator<Gate> for GateSeq {
    fn from_iter<I: IntoIterator<Item = Gate>>(iter: I) -> Self {
        Self(Vec::from_iter(iter))
    }
}

impl Extend<Gate> for GateSeq {
    fn extend<I: IntoIterator<Item = Gate>>(&mut self, iter: I) {
        self.0.extend(iter);
    }
}

impl<'a> Extend<&'a Gate> for GateSeq {
    fn extend<I: IntoIterator<Item = &'a Gate>>(&mut self, iter: I) {
        self.0.extend(iter);
    }
}

impl IntoIterator for GateSeq {
    type Item = Gate;
    type IntoIter = std::vec::IntoIter<Gate>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a GateSeq {
    type Item = &'a Gate;
    type IntoIter = std::slice::Iter<'a, Gate>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GOLDEN_GATES: [&str; 6] = [
        "HTHTSHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTSHTHTSHTHTSHTXSSWWW",
        "HTSHTSHTSHTHTHTSHTHTHTSHTSHTHTHTHTHTHTSHTHTHTHTSHTSHTHTSHTHTHTHTHTHTHTHTSHTHTHTHTSHTHTSHTSHTSHTSHTHTSHTSHTHTSHTSHTHTSHTHTHTSHTSHTHTHTHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTHTSHTHTSHTSHTHTHTHTSHTHTHTSHTHTHTHTSHTHTSHTSHTHTSHTHTSHTHTHTHTHTHTHTHTSHTHTHTSHTSSSWW",
        "HTSHTHTHTSHTHTHTHTHTHTHTHTSHTHTSHTHTSHTSHTHTSHTHTHTHTSHTHTHTSHTHTHTHTSHTSHTHTSHTHTHTHTHTSHTSHTHTHTSHTHTSHTHTSHTHTHTHTSHTSHTHTHTSHTHTSHTSHTHTSHTSHTHTSHTSHTSHTSHTHTSHTHTHTHTSHTHTHTHTHTHTHTHTSHTHTSHTSHTHTHTHTSHTHTHTHTHTHTSHTSHTHTHTSHTHTHTSHTSHTSHTSSSWW",
        "SWWWWWWW",
        "I",
        "H",
    ];

    #[test]
    fn round_trips_historical_literals() {
        for s in GOLDEN_GATES {
            let seq: GateSeq = s.parse().expect("golden strings are valid gate words");
            assert_eq!(seq.to_string(), s, "round-trip mismatch for {s:?}");
        }
    }

    #[test]
    fn structural_round_trip_including_empty() {
        let seqs: Vec<GateSeq> = vec![
            GateSeq::new(),
            GateSeq::from([Gate::H]),
            GateSeq::from([Gate::S, Gate::H, Gate::T]),
            GateSeq::from([
                Gate::T,
                Gate::W,
                Gate::W,
                Gate::W,
                Gate::W,
                Gate::W,
                Gate::W,
                Gate::W,
            ]),
        ];
        for seq in seqs {
            let round_tripped: GateSeq = seq.to_string().parse().unwrap();
            assert_eq!(
                round_tripped, seq,
                "structural round-trip mismatch for {seq:?}"
            );
        }
    }

    #[test]
    fn empty_renders_as_identity_sentinel() {
        assert_eq!(GateSeq::new().to_string(), "I");
        assert_eq!(GateSeq::default(), GateSeq::identity());
    }

    #[test]
    fn identity_sentinel_parses_to_empty() {
        let seq: GateSeq = "I".parse().unwrap();
        assert_eq!(seq, GateSeq::new());
        assert!(seq.is_empty());
    }

    #[test]
    fn embedded_identity_sentinel_is_dropped_on_parse() {
        // The historical leak this type fixes: `NormalForm::to_gates` used to be able to emit a
        // mid-word "I" (e.g. "HTSHTI") that its own `from_gates` could not parse back. `GateSeq`
        // treats 'I' as a no-op wherever it appears, so this no longer panics.
        let seq: GateSeq = "HTSHTI".parse().unwrap();
        assert_eq!(seq.to_string(), "HTSHT");

        let seq: GateSeq = "IHITI".parse().unwrap();
        assert_eq!(seq, GateSeq::from([Gate::H, Gate::T]));
    }

    #[test]
    fn t_count_agrees_with_char_scan() {
        for s in GOLDEN_GATES {
            let seq: GateSeq = s.parse().unwrap();
            let manual = s.chars().filter(|&c| c == 'T').count();
            assert_eq!(seq.t_count(), manual, "t_count mismatch for {s:?}");
        }
    }

    #[test]
    fn parse_error_reports_offending_char() {
        let err = "HQT".parse::<GateSeq>().unwrap_err();
        assert_eq!(err.found, 'Q');
    }

    #[test]
    fn gate_char_round_trip() {
        for g in Gate::ALL {
            assert_eq!(Gate::from_char(g.as_char()), Some(g));
            assert_eq!(g.to_string(), g.as_char().to_string());
        }
        assert_eq!(Gate::from_char('I'), None);
    }
}
