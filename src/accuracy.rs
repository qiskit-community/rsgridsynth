// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

//! Shared, on-demand accuracy computation, usable by both the plain single-candidate
//! synthesis path (`crate::gridsynth`) and the "mixed diagonal"/"fallback"/"mixed fallback"
//! protocols (`crate::protocol`). Lives below both so neither has to depend on the other to
//! reach it.
//!
//! [`AchievedDiamondError`] is the uniform interface: every result type recomputes its
//! diamond-norm distance to the ideal target rotation on demand, straight from its own public
//! gate string(s) (decoded via `DOmegaUnitary::from_gates`), rather than caching it eagerly
//! during synthesis. This mirrors the crate's own internal formulas (`WFrame`/
//! `diagonal_diamond_distance`) -- it is a convenience for checking a result you already have,
//! not an independent derivation from a different method. `examples/pauli_transfer_verification.rs`
//! implements a genuinely independent (Pauli-transfer-matrix-based) check, for that purpose.

use crate::common::{cos_fbig, fb_with_prec, ib_to_bf_prec, sin_fbig};
use crate::gate::Gate;
use crate::math::sqrt_fbig;
use crate::ring::DOmega;
use crate::unitary::DOmegaUnitary;
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;
use num::Complex;

/// The rotated-frame projection used to measure rotation error against a target angle
/// `theta`. Given a candidate unitary's top-left entry `u`, defines
/// `w := u * e^{i theta/2}`, i.e. `u` rotated so that the target direction `e^{-i theta/2}`
/// maps to `1`. `Re(w)` is exactly the `cos_similarity` quantity computed inline (and
/// duplicated) inside [`crate::gridsynth::EpsilonRegion::inside`]/`intersect`; this type
/// factors that arithmetic out.
///
/// Derivation: with `z := z_x + i*z_y = e^{-i theta/2}` (as `EpsilonRegion::new` computes),
/// `e^{i theta/2} = conj(z) = z_x - i*z_y`. So
/// `w = (z_x - i*z_y) * (Re(u) + i*Im(u))
///    = [z_x*Re(u) + z_y*Im(u)] + i*[z_x*Im(u) - z_y*Re(u)]`.
/// Hence `Re(w) = z_x*Re(u) + z_y*Im(u)` and `Im(w) = z_x*Im(u) - z_y*Re(u)`.
#[derive(Debug, Clone)]
pub struct WFrame {
    z_x: FBig<HalfEven>,
    z_y: FBig<HalfEven>,
}

impl WFrame {
    /// Builds the rotated frame for target angle `theta`, computing `z_x = cos(-theta/2)`,
    /// `z_y = sin(-theta/2)` exactly as `EpsilonRegion::new` does.
    pub fn new(theta: &FBig<HalfEven>) -> Self {
        let two = fb_with_prec(FBig::try_from(2.0).unwrap());
        let theta_half = fb_with_prec(theta / &two);
        let neg_theta_half = -fb_with_prec(theta_half);
        let z_x: FBig<HalfEven> = fb_with_prec(cos_fbig(&neg_theta_half));
        let z_y: FBig<HalfEven> = fb_with_prec(sin_fbig(&neg_theta_half));
        Self { z_x, z_y }
    }

    /// Builds the same frame as [`WFrame::new`], but from the target direction's half-angle
    /// `(cos(-phi/2), sin(-phi/2))` directly, for a caller that already has that pair from
    /// algebraic angle-addition/half-angle identities (e.g. a fallback correction's residual
    /// angle) rather than a raw angle -- avoiding an `atan2`-style angle round-trip, mirroring
    /// [`crate::gridsynth::EpsilonRegion::from_target_direction`]. `(z_x, z_y)` must satisfy
    /// `z_x^2 + z_y^2 == 1`; not checked.
    ///
    /// Only used by `crate::protocol::*`, which `src/main.rs`'s separate module tree (the
    /// `cli`-gated binary target reuses these source files as its own crate root, without
    /// declaring `protocol`) never compiles -- hence `allow(dead_code)`, since that
    /// compilation would otherwise warn on this method despite it being used elsewhere.
    #[allow(dead_code)]
    pub(crate) fn from_target_direction(z_x: FBig<HalfEven>, z_y: FBig<HalfEven>) -> Self {
        Self { z_x, z_y }
    }

    /// `Re(w)` where `w = u * e^{i theta/2}`, for `u` already expressed as an `FBig` complex
    /// pair rather than a `DOmega` -- needed when `u` isn't (or can't be) an exact ring
    /// element, e.g. after rotating by an extra phase that isn't ring-representable (like the
    /// `e^{i pi/8}` shift `PhaseMode::Shifted` uses).
    pub fn re_w_fbig(&self, re: &FBig<HalfEven>, im: &FBig<HalfEven>) -> FBig<HalfEven> {
        let term1 = fb_with_prec(&self.z_x * re);
        let term2 = fb_with_prec(&self.z_y * im);
        fb_with_prec(&term1 + &term2)
    }

    /// `Re(w)` where `w = u * e^{i theta/2}`. Matches `EpsilonRegion`'s existing
    /// `cos_similarity` exactly (same formula, same operand order).
    pub fn re_w(&self, u: &DOmega) -> FBig<HalfEven> {
        self.re_w_fbig(u.real(), u.imag())
    }

    /// `Im(w)` where `w = u * e^{i theta/2}`.
    pub fn im_w(&self, u: &DOmega) -> FBig<HalfEven> {
        let term1 = fb_with_prec(&self.z_x * u.imag());
        let term2 = fb_with_prec(&self.z_y * u.real());
        fb_with_prec(&term1 - &term2)
    }
}

/// Exact diamond-norm distance between a target Z-rotation and its diagonal-unitary
/// approximation, given only the achieved `Re(w)` (`w = u * e^{i theta/2}`, see
/// [`WFrame::re_w`]): `||Z_phi - U||_diamond = 2*sqrt(1 - Re(w)^2)`.
pub fn diagonal_diamond_distance(re_w: &FBig<HalfEven>) -> FBig<HalfEven> {
    let one = ib_to_bf_prec(IBig::ONE);
    let re_w_sq = fb_with_prec(re_w * re_w);
    let one_minus_re_w_sq = fb_with_prec(&one - &re_w_sq);
    // Guard against tiny negative values from rounding error, matching the analogous
    // clamp in `gridsynth::compute_error`.
    let zero = ib_to_bf_prec(IBig::ZERO);
    let clamped = one_minus_re_w_sq.max(zero);
    let two = fb_with_prec(FBig::try_from(2.0).unwrap());
    fb_with_prec(&two * sqrt_fbig(&clamped))
}

/// Decodes `gates`, optionally rotating the decoded top-left entry by the extra `e^{i pi/8}`
/// global phase `PhaseMode::Shifted` uses (not ring-representable in `Z[omega]`, so this has
/// to happen in `FBig` on the decoded complex entry rather than in exact ring arithmetic), and
/// returns the diamond-norm distance to the ideal target rotation by `theta`.
pub(crate) fn gate_seq_diamond_error(
    theta: &FBig<HalfEven>,
    gates: &[Gate],
    extra_eighth_turn: bool,
) -> FBig<HalfEven> {
    let wframe = WFrame::new(theta);
    let u = DOmegaUnitary::from_gates(gates);

    let re_w = if extra_eighth_turn {
        let p = fb_with_prec(FBig::<HalfEven>::try_from(std::f64::consts::PI / 8.).unwrap());
        let phase = Complex::new(fb_with_prec(cos_fbig(&p)), fb_with_prec(sin_fbig(&p)));
        let z = Complex::new(u.z().real().clone(), u.z().imag().clone());
        let shifted = &z * &phase;
        wframe.re_w_fbig(&shifted.re, &shifted.im)
    } else {
        wframe.re_w(u.z())
    };

    diagonal_diamond_distance(&re_w)
}

/// Decodes `gates` back into the unitary it represents and returns the diamond-norm distance
/// between it and the ideal single-candidate diagonal target rotation by `theta`, i.e.
/// `diagonal_diamond_distance(Re(w))` for `w` derived from the *decoded* top-left entry
/// rather than from whatever `FBig` value the search happened to hold internally.
///
/// Exposed so callers (and tests) can check a synthesized diagonal candidate's accuracy
/// straight from its public gate string, without re-deriving `WFrame`/`Re(w)` inline -- and so
/// that check exercises the actual gate-encode/decode round trip (`decompose_domega_unitary`
/// followed by `DOmegaUnitary::from_gates`), not just the value computed mid-search.
pub fn achieved_diagonal_diamond_error(theta: &FBig<HalfEven>, gates: &[Gate]) -> FBig<HalfEven> {
    gate_seq_diamond_error(theta, gates, false)
}

/// Like [`achieved_diagonal_diamond_error`], but normalizes the decoded top-left entry's
/// *phase* first (`z / |z|`) before comparing it to `theta`, rather than feeding in the raw
/// (possibly non-unit-modulus) entry.
///
/// This matters specifically for the fallback family's "projective" candidates
/// (`crate::protocol::fallback`/`crate::protocol::mixed_fallback`): unlike an ordinary
/// gridsynth/mixed-diagonal candidate (whose `z` is close to unit modulus *because* the whole
/// unitary is close to the target in operator norm), a fallback projective candidate's `z` is
/// deliberately allowed `|z|^2 = q` well below `1` -- the magnitude deficit `1 - q` is the
/// probability of the *separate* failure branch, not angular error. Feeding the raw `z` into
/// `diagonal_diamond_distance` conflates that magnitude deficit with angular error, giving a
/// nonsensical, epsilon-independent `~2*sqrt(1-q)` answer instead of the actual (typically much
/// smaller) angular error.
pub fn achieved_phase_diamond_error(theta: &FBig<HalfEven>, gates: &[Gate]) -> FBig<HalfEven> {
    let u = DOmegaUnitary::from_gates(gates);
    let z = u.z();
    let norm_sq =
        fb_with_prec(fb_with_prec(z.real() * z.real()) + fb_with_prec(z.imag() * z.imag()));
    let norm = sqrt_fbig(&norm_sq);
    let re_n = fb_with_prec(z.real() / &norm);
    let im_n = fb_with_prec(z.imag() / &norm);

    let wframe = WFrame::new(theta);
    let re_w = wframe.re_w_fbig(&re_n, &im_n);
    diagonal_diamond_distance(&re_w)
}

/// Implemented by result types whose diamond-norm distance to an ideal Z-rotation by `theta`
/// can be recomputed on demand, directly from their own public gate-string data, rather than
/// being cached eagerly during synthesis.
pub trait AchievedDiamondError {
    /// Recomputes the diamond-norm distance between this result's synthesized channel and the
    /// ideal Z-rotation by `theta`.
    fn achieved_diamond_error(&self, theta: &FBig<HalfEven>) -> FBig<HalfEven>;
}
