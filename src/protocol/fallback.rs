// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

//! Stage 2: the projective/fallback annulus-sector region.
//!
//! Implements the plain "fallback" protocol of Bocharov-Roetteler-Svore
//! (arXiv:1409.3552), as adapted by Kliuchnikov, Lauter, Minko, Paetznick, Petit
//! (arXiv:2203.10064v2, Prop 3.9). Instead of requiring the synthesized candidate's
//! top-left entry `z` to lie in a small angular window around the target direction (as
//! [`crate::gridsynth::EpsilonRegion`] does), fallback relaxes the angular tolerance a lot
//! (a much shorter T-count) but only accepts a candidate whose magnitude `|z|^2` is at
//! least some threshold `q` close to 1. With probability `|z|^2` the projective step alone
//! is already close enough (the "success" branch); with probability `1 - |z|^2` a
//! classically-computed diagonal "correction" gate is needed (the "failure"/fallback
//! branch). Because failures are rare (`q` close to 1), the *expected* T-count is much
//! lower than paying for full angular precision unconditionally.
//!
//! Unlike "mixed fallback" (a later stage), plain fallback only needs a single solved
//! candidate satisfying the new [`SectorRegion`] -- no straddling pair, no new search-loop
//! logic. [`crate::gridsynth::search_for_solution`] (already generic over any `Region`) is
//! reused unchanged.

use crate::accuracy::{
    achieved_phase_diamond_error, diagonal_diamond_distance, AchievedDiamondError, WFrame,
};
use crate::common::{cos_fbig, fb_with_prec, ib_to_bf_prec, sin_fbig};
use crate::config::config_from_theta_epsilon;
use crate::gate::{Gate, GateSeq};
use crate::gridsynth::{
    search_for_solution, setup_regions_and_transform, EpsilonRegion, PhaseMode, UnitDisk,
};
use crate::math::{sign, sqrt2, sqrt_fbig};
use crate::protocol::mixing::diamond_to_spec_epsilon;
use crate::region::Ellipse;
use crate::ring::{DOmega, DRootTwo, ZRootTwo};
use crate::synthesis_of_clifford_t::decompose_domega_unitary;
use crate::tdgp::Region;
use crate::unitary::DOmegaUnitary;

use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;

use nalgebra::{Matrix2, Vector2};

/// `(cos(phi), sin(phi), |v|^2)` for `phi = Arg(v)`, where `v` is a projective candidate's
/// off-diagonal `w` entry. Degenerate case: `v = 0` exactly (the projective candidate is an
/// exact diagonal unitary, i.e. `|z|^2 = 1` and the success probability is already 1) makes
/// `Arg(v)` undefined and the naive `re_v / |v|` division panic (dashu-int's "divisor must
/// not be 0"). Since the correction branch is then never selected (its weight is
/// `1 - achieved_success_probability() = 0`), any convention for `phi` is equally valid;
/// `phi = 0` is used so the residual angle collapses to `theta` itself. The returned `|v|^2`
/// is likewise substituted with `1` (rather than the true `0`) so callers that scale a
/// correction-search epsilon budget by `1 / |v|^2` (irrelevant in this branch, since the
/// correction is never used) don't divide by zero either.
pub(crate) fn phase_cos_sin(v: &DOmega) -> (FBig<HalfEven>, FBig<HalfEven>, FBig<HalfEven>) {
    let re_v = v.real().clone();
    let im_v = v.imag().clone();
    let v_norm_sq = fb_with_prec(fb_with_prec(&re_v * &re_v) + fb_with_prec(&im_v * &im_v));
    if v_norm_sq.repr().is_zero() {
        return (
            ib_to_bf_prec(IBig::ONE),
            ib_to_bf_prec(IBig::ZERO),
            ib_to_bf_prec(IBig::ONE),
        );
    }
    let v_norm = sqrt_fbig(&v_norm_sq);
    let cos_phi = fb_with_prec(&re_v / &v_norm);
    let sin_phi = fb_with_prec(&im_v / &v_norm);
    (cos_phi, sin_phi, v_norm_sq)
}

fn to_fbig(x: f64) -> FBig<HalfEven> {
    FBig::<HalfEven>::try_from(x)
        .unwrap()
        .with_precision(crate::common::get_prec_bits())
        .value()
}

/// Local re-implementation of the 2x2 matrix product used by
/// [`crate::gridsynth::EpsilonRegion::new`]'s `D1 * D2 * D3` ellipse construction.
/// `nalgebra`'s own `Mul` for `Matrix2` requires its scalar type to implement
/// `num_traits`/`Scalar`/`ClosedAdd`/`ClosedMul`, which `FBig<HalfEven>` does not, so
/// `gridsynth.rs` implements this by hand -- but that helper is module-private there and
/// this module may not edit `gridsynth.rs` to export it, so it is duplicated here.
fn matrix_multiply_2x2(
    a: &Matrix2<FBig<HalfEven>>,
    b: &Matrix2<FBig<HalfEven>>,
) -> Matrix2<FBig<HalfEven>> {
    let mut result = Matrix2::from_element(ib_to_bf_prec(IBig::ZERO));
    for i in 0..2 {
        for j in 0..2 {
            let mut sum = ib_to_bf_prec(IBig::ZERO);
            for k in 0..2 {
                sum += &a[(i, k)] * &b[(k, j)];
            }
            result[(i, j)] = sum;
        }
    }
    result
}

/// `q = 1 - 2^{-m}`, represented exactly in `ℤ[1/√2]` (this crate's `DRootTwo`), rather than
/// as a raw `f64`.
///
/// The paper's success-probability threshold is `|z|^2 >= q` with `1 - q <= 0.01`. `0.01` is
/// not exactly representable in `ℤ[1/√2]`, so accepting an arbitrary `f64` `q` and comparing
/// in floating point would silently reintroduce the class of "silently misses valid
/// solutions at low precision" bug this crate has already been burned by once (see
/// `test_low_precision_bug` in `tests/integration_test.rs`). `m = 7` gives `1 - q =
/// 0.0078125 <= 0.01`, satisfying the paper's bound, and is a reasonable default for
/// callers that don't need a different `m`.
///
/// The numerator `2^m - 1` has `b = 0` (a plain rational, no `sqrt(2)` component needed),
/// and the denominator exponent `k = 2*m` is EVEN so that `sqrt(2)^k = 2^(k/2) = 2^m` is
/// itself a plain rational (for an ODD `k`, `sqrt(2)^k` would be irrational).
pub fn exact_q(m: u32) -> DRootTwo {
    let numerator = (IBig::from(1) << (m as usize)) - IBig::from(1); // 2^m - 1
    DRootTwo::new(ZRootTwo::new(numerator, IBig::ZERO), 2 * (m as i64))
}

/// The annulus-sector region of Prop 3.9 (arXiv:2203.10064v2): `q*scale <= |u|^2 <= scale`
/// AND `|sin(Arg(w))| <= sin_alpha`, where `w = u * e^{i theta/2}` is `u` rotated into the
/// frame where the target direction `e^{-i theta/2}` maps to the positive real axis (same
/// convention as [`crate::gridsynth::EpsilonRegion`]/[`crate::protocol::mixing::WFrame`]).
///
/// Parameterized by `sin_alpha` directly (rather than by `alpha` itself) since the angular
/// half-width is only ever needed through its sine, avoiding an inverse-trig call this
/// crate does not implement.
#[derive(Debug)]
pub struct SectorRegion {
    scale: ZRootTwo,
    /// `q * scale`: the exact lower-bound threshold on `|u|^2`.
    q_scaled: DRootTwo,
    sin_alpha: FBig<HalfEven>,
    z_x: FBig<HalfEven>,
    z_y: FBig<HalfEven>,
    ellipse: Ellipse,
}

impl SectorRegion {
    /// Builds the sector region for target angle `theta`, radial-magnitude threshold `q`
    /// (see [`exact_q`]), angular half-width `sin_alpha = sin(alpha)`, and the same `scale`
    /// convention `EpsilonRegion`/`UnitDisk` use (region area is scaled by `scale` in the
    /// "up to phase" case; plain fallback uses `scale = 1`).
    pub fn new(
        theta: &FBig<HalfEven>,
        q: DRootTwo,
        sin_alpha: FBig<HalfEven>,
        scale: ZRootTwo,
    ) -> Self {
        let two = fb_with_prec(FBig::try_from(2.0).unwrap());
        let theta_half = fb_with_prec(theta / &two);
        let neg_theta_half = -fb_with_prec(theta_half);
        let z_x: FBig<HalfEven> = fb_with_prec(cos_fbig(&neg_theta_half));
        let z_y: FBig<HalfEven> = fb_with_prec(sin_fbig(&neg_theta_half));

        let q_scaled = q * DRootTwo::from_zroottwo(scale.clone());

        let one = ib_to_bf_prec(IBig::ONE);
        let sin_sq = fb_with_prec(&sin_alpha * &sin_alpha);
        let cos_alpha = sqrt_fbig(&fb_with_prec(&one - &sin_sq));

        let sqrt_s = sqrt_fbig(&scale.to_real());
        let qs_real = q_scaled.to_real();
        let sqrt_qs = sqrt_fbig(&qs_real);

        // Box radial half-width / center: the box spans the radial interval
        // [sqrt(q*scale)*cos(alpha), sqrt(scale)] (the CHORD's x-value as the inner edge --
        // a valid over-approximation of the true arc-bounded inner edge, since every point
        // on the true inner arc within the angular wedge has x-coordinate >=
        // sqrt(q*scale)*cos(alpha)).
        let inner_x = fb_with_prec(&sqrt_qs * &cos_alpha);
        let a0 = fb_with_prec(fb_with_prec(&sqrt_s - &inner_x) / &two);
        let x_c = fb_with_prec(fb_with_prec(&inner_x + &sqrt_s) / &two);
        // Box tangential half-width, at the outer radius (the sector's widest point).
        let b0 = fb_with_prec(&sqrt_s * &sin_alpha);

        // Circumscribe the box with an ellipse using semi-axes sqrt(2)*A0, sqrt(2)*B0: for
        // |x|<=A0, |y|<=B0, (x/(sqrt(2)*A0))^2 + (y/(sqrt(2)*B0))^2 <= 1/2 + 1/2 = 1, so the
        // box (hence the true sector, which the box contains) is provably contained.
        let sqrt2_val = sqrt2();
        let a_axis = fb_with_prec(&sqrt2_val * &a0);
        let b_axis = fb_with_prec(&sqrt2_val * &b0);

        let zero: FBig<HalfEven> = ib_to_bf_prec(IBig::ZERO);
        let neg_z_y: FBig<HalfEven> = -fb_with_prec(z_y.clone());
        let d1: Matrix2<FBig<HalfEven>> =
            Matrix2::new(z_x.clone(), neg_z_y.clone(), z_y.clone(), z_x.clone());
        let inv_a2 = fb_with_prec(&one / fb_with_prec(&a_axis * &a_axis));
        let inv_b2 = fb_with_prec(&one / fb_with_prec(&b_axis * &b_axis));
        let d2: Matrix2<FBig<HalfEven>> = Matrix2::new(inv_a2, zero.clone(), zero.clone(), inv_b2);
        let d3: Matrix2<FBig<HalfEven>> =
            Matrix2::new(z_x.clone(), z_y.clone(), neg_z_y, z_x.clone());
        let m1 = matrix_multiply_2x2(&d1, &d2);
        let m = matrix_multiply_2x2(&m1, &d3);

        let px = fb_with_prec(&x_c * &z_x);
        let py = fb_with_prec(&x_c * &z_y);
        let p = Vector2::new(px, py);
        let ellipse = Ellipse::new(m, p);

        Self {
            scale,
            q_scaled,
            sin_alpha,
            z_x,
            z_y,
            ellipse,
        }
    }
}

/// Clips `(t0, t1)` to satisfy the linear half-plane constraint `t * slope >= rhs`, or
/// returns `None` if doing so empties the interval. Unlike
/// [`crate::gridsynth::EpsilonRegion::intersect`]'s single half-plane clip (where an
/// inverted `(t0, t1)` pair is rare enough to be mostly harmless), `SectorRegion::intersect`
/// applies three sequential clips, so checking for emptiness after each is required.
fn clip_ge(
    t0: FBig<HalfEven>,
    t1: FBig<HalfEven>,
    slope: &FBig<HalfEven>,
    rhs: &FBig<HalfEven>,
) -> Option<(FBig<HalfEven>, FBig<HalfEven>)> {
    let zero = ib_to_bf_prec(IBig::ZERO);
    if slope > &zero {
        let bound = fb_with_prec(rhs / slope);
        let new_t0 = if t0 > bound { t0 } else { bound };
        if new_t0 > t1 {
            None
        } else {
            Some((new_t0, t1))
        }
    } else if slope < &zero {
        let bound = fb_with_prec(rhs / slope);
        let new_t1 = if t1 < bound { t1 } else { bound };
        if t0 > new_t1 {
            None
        } else {
            Some((t0, new_t1))
        }
    } else if rhs > &zero {
        None
    } else {
        Some((t0, t1))
    }
}

/// Clips `(t0, t1)` to satisfy the linear half-plane constraint `t * slope <= rhs`, with the
/// same emptiness check as [`clip_ge`].
fn clip_le(
    t0: FBig<HalfEven>,
    t1: FBig<HalfEven>,
    slope: &FBig<HalfEven>,
    rhs: &FBig<HalfEven>,
) -> Option<(FBig<HalfEven>, FBig<HalfEven>)> {
    let zero = ib_to_bf_prec(IBig::ZERO);
    if slope > &zero {
        let bound = fb_with_prec(rhs / slope);
        let new_t1 = if t1 < bound { t1 } else { bound };
        if t0 > new_t1 {
            None
        } else {
            Some((t0, new_t1))
        }
    } else if slope < &zero {
        let bound = fb_with_prec(rhs / slope);
        let new_t0 = if t0 > bound { t0 } else { bound };
        if new_t0 > t1 {
            None
        } else {
            Some((new_t0, t1))
        }
    } else if rhs < &zero {
        None
    } else {
        Some((t0, t1))
    }
}

impl Region for SectorRegion {
    fn ellipse(&self) -> Ellipse {
        self.ellipse.clone()
    }

    fn inside(&self, u: &DOmega) -> bool {
        let norm = DRootTwo::from_domega(u.conj() * u);
        if norm > DRootTwo::from_zroottwo(self.scale.clone()) {
            return false;
        }
        if norm < self.q_scaled {
            return false;
        }

        // Im(w) = z_x*Im(u) - z_y*Re(u), matching `WFrame::im_w` exactly.
        let term1 = fb_with_prec(&self.z_x * u.imag());
        let term2 = fb_with_prec(&self.z_y * u.real());
        let im_w = fb_with_prec(&term1 - &term2);
        let im_w_sq = fb_with_prec(&im_w * &im_w);

        let norm_real = norm.to_real();
        let sin_alpha_sq = fb_with_prec(&self.sin_alpha * &self.sin_alpha);
        let rhs = fb_with_prec(&sin_alpha_sq * &norm_real);

        im_w_sq <= rhs
    }

    fn intersect(&self, u0: &DOmega, v: &DOmega) -> Option<(FBig<HalfEven>, FBig<HalfEven>)> {
        // Outer disc: |L(t)|^2 <= scale (exact quadratic, same as
        // EpsilonRegion/UnitDisk::intersect).
        let a = v.conj() * v;
        let b = 2 * (v.conj() * u0);
        let c = u0.conj() * u0 - DOmega::from_zroottwo(&self.scale);
        let (t0, t1) = crate::math::solve_quadratic(a.real(), b.real(), c.real())?;

        let re_w_u0 =
            fb_with_prec(fb_with_prec(&self.z_x * u0.real()) + fb_with_prec(&self.z_y * u0.imag()));
        let im_w_u0 =
            fb_with_prec(fb_with_prec(&self.z_x * u0.imag()) - fb_with_prec(&self.z_y * u0.real()));
        let re_w_v =
            fb_with_prec(fb_with_prec(&self.z_x * v.real()) + fb_with_prec(&self.z_y * v.imag()));
        let im_w_v =
            fb_with_prec(fb_with_prec(&self.z_x * v.imag()) - fb_with_prec(&self.z_y * v.real()));

        let one = ib_to_bf_prec(IBig::ONE);
        let sin_sq = fb_with_prec(&self.sin_alpha * &self.sin_alpha);
        let cos_alpha = sqrt_fbig(&fb_with_prec(&one - &sin_sq));
        let tan_alpha = fb_with_prec(&self.sin_alpha / &cos_alpha);

        // (a) Im(w) <= Re(w)*tan(alpha)  <=>  t*gv <= rhs_a, gv = Im(w_v) - tan*Re(w_v).
        let gv = fb_with_prec(&im_w_v - fb_with_prec(&tan_alpha * &re_w_v));
        let rhs_a = fb_with_prec(fb_with_prec(&tan_alpha * &re_w_u0) - &im_w_u0);
        let (t0, t1) = clip_le(t0, t1, &gv, &rhs_a)?;

        // (b) Im(w) >= -Re(w)*tan(alpha)  <=>  t*hv >= rhs_b, hv = Im(w_v) + tan*Re(w_v).
        let hv = fb_with_prec(&im_w_v + fb_with_prec(&tan_alpha * &re_w_v));
        let rhs_b = -fb_with_prec(fb_with_prec(&tan_alpha * &re_w_u0) + &im_w_u0);
        let (t0, t1) = clip_ge(t0, t1, &hv, &rhs_b)?;

        // (c) Re(w) >= sqrt(q_scaled)*cos(alpha) -- the chord replacing the inner arc,
        // which is what makes the sector's hull convex (excluding an inner disc is what
        // makes an annulus non-convex).
        let qs_real = self.q_scaled.to_real();
        let sqrt_qs = sqrt_fbig(&qs_real);
        let d_inner = fb_with_prec(&sqrt_qs * &cos_alpha);
        let rhs_c = fb_with_prec(&d_inner - &re_w_u0);
        clip_ge(t0, t1, &re_w_v, &rhs_c)
    }
}

/// `(cos(phi/2), sin(phi/2))` from `(cos(phi), sin(phi))`, via the half-angle formulas,
/// avoiding `atan2`/any inverse-trig call.
///
/// `cos(phi/2) = sqrt((1+cos(phi))/2)` (always the "+" root; valid for `phi` in `(-pi,
/// pi]`), and `sin(phi/2) = sign(sin(phi)) * sqrt((1-cos(phi))/2)` -- except when `phi = pi`
/// exactly (`sin(phi)` exactly zero and `cos(phi) < 0`), where `sign(sin(phi)) = 0` would
/// wrongly give `sin(phi/2) = 0` instead of the correct `sin(pi/2) = 1`; that degenerate
/// case is detected and handled directly.
pub(crate) fn half_angle_cos_sin(
    cos_phi: &FBig<HalfEven>,
    sin_phi: &FBig<HalfEven>,
) -> (FBig<HalfEven>, FBig<HalfEven>) {
    let zero = ib_to_bf_prec(IBig::ZERO);
    if *sin_phi == zero && *cos_phi < zero {
        // phi == pi exactly: cos(phi/2) = 0, sin(phi/2) = 1 (either sign is a valid,
        // self-consistent branch choice; +1 is chosen).
        return (zero, ib_to_bf_prec(IBig::ONE));
    }

    let one = ib_to_bf_prec(IBig::ONE);
    let two = to_fbig(2.0);
    // `cos_phi` is a ratio of low-precision `FBig` values (see callers), so it can round to
    // just outside `[-1, 1]`; guard against the resulting tiny negative `sqrt_fbig` input,
    // matching the analogous clamp in `gridsynth::compute_error`/`mixing::diagonal_diamond_distance`.
    let zero = ib_to_bf_prec(IBig::ZERO);
    let one_plus_cos = fb_with_prec(&one + cos_phi).max(zero.clone());
    let one_minus_cos = fb_with_prec(&one - cos_phi).max(zero);
    let cos_half = sqrt_fbig(&fb_with_prec(&one_plus_cos / &two));
    let sin_half_mag = sqrt_fbig(&fb_with_prec(&one_minus_cos / &two));

    let sin_half = if sign(sin_phi.clone()) < 0 {
        -sin_half_mag
    } else {
        sin_half_mag
    };

    (cos_half, sin_half)
}

/// The output of [`synth_fallback`]: gate sequences for the projective step and its (rare)
/// classical correction, plus the `q` threshold used. Call
/// [`FallbackResult::achieved_success_probability`] or
/// [`AchievedDiamondError::achieved_diamond_error`] to compute accuracy on demand.
#[derive(Debug, Clone)]
pub struct FallbackResult {
    /// Gate sequence for the projective/fallback step, applied unconditionally.
    pub projective_gates: GateSeq,
    /// Gate sequence for the classical correction, applied only on the "failure" branch
    /// (probability `1 - achieved_success_probability()`).
    pub correction_gates: GateSeq,
    /// The `q` threshold used to find the projective candidate.
    pub q: DRootTwo,
}

impl FallbackResult {
    /// Recomputes the achieved success probability directly from the returned
    /// `projective_gates` string (decoding it back into a unitary and taking its top-left
    /// entry's squared magnitude).
    pub fn achieved_success_probability(&self) -> FBig<HalfEven> {
        let u = DOmegaUnitary::from_gates(&self.projective_gates);
        let z = u.z();
        fb_with_prec(fb_with_prec(z.real() * z.real()) + fb_with_prec(z.imag() * z.imag()))
    }
}

/// Diamond-norm distance between `correction_gates` and the *residual* target rotation
/// `theta - Arg(v)` that it actually approximates -- `v` decoded from `projective_gates`'s own
/// off-diagonal (`w`) entry, per the same `atan2`-free half-angle algebra
/// [`synth_fallback`]/[`crate::protocol::mixed_fallback::build_side`] use to derive that
/// residual angle in the first place. NOT `theta` itself, and NOT composed with the
/// projective step: the correction is a standalone approximation of the residual angle (see
/// this module's doc comment on why the projective/failure split works the way it does), so
/// its own accuracy has to be measured against that residual, decoded here purely from the
/// public `projective_gates` string rather than any transient search state.
/// Builds the [`WFrame`] for the *residual* target rotation `theta - Arg(v)` that a
/// projective step's correction actually approximates -- `v` decoded from `projective_gates`'s
/// own off-diagonal (`w`) entry, via the same `atan2`-free half-angle algebra
/// [`synth_fallback`]/[`crate::protocol::mixed_fallback::build_side`] use to derive that
/// residual angle in the first place. Shared by [`residual_diamond_error`] (a single
/// correction gate string) and [`residual_diamond_error_mixed`] (a `MixedDiagonalResult`
/// correction, e.g. mixed fallback's).
pub(crate) fn residual_wframe(theta: &FBig<HalfEven>, projective_gates: &[Gate]) -> WFrame {
    let v = DOmegaUnitary::from_gates(projective_gates).w().clone();
    let (cos_phi, sin_phi, _) = phase_cos_sin(&v);
    let (cos_half_phi, sin_half_phi) = half_angle_cos_sin(&cos_phi, &sin_phi);

    let two = to_fbig(2.0);
    let neg_theta_half = -fb_with_prec(theta / &two);
    let z_x = fb_with_prec(cos_fbig(&neg_theta_half));
    let z_y = fb_with_prec(sin_fbig(&neg_theta_half));

    // cos(-theta_B/2) = cos(A+B) = Z_X*cos(phi/2) - Z_Y*sin(phi/2)
    // sin(-theta_B/2) = sin(A+B) = Z_Y*cos(phi/2) + Z_X*sin(phi/2)
    let cos_neg_theta_b_half =
        fb_with_prec(fb_with_prec(&z_x * &cos_half_phi) - fb_with_prec(&z_y * &sin_half_phi));
    let sin_neg_theta_b_half =
        fb_with_prec(fb_with_prec(&z_y * &cos_half_phi) + fb_with_prec(&z_x * &sin_half_phi));

    WFrame::from_target_direction(cos_neg_theta_b_half, sin_neg_theta_b_half)
}

/// Diamond-norm distance between `correction_gates` and the residual target rotation it
/// actually approximates (see [`residual_wframe`]) -- NOT `theta` itself, and NOT composed
/// with the projective step: the correction is a standalone approximation of the residual
/// angle.
pub(crate) fn residual_diamond_error(
    theta: &FBig<HalfEven>,
    projective_gates: &[Gate],
    correction_gates: &[Gate],
) -> FBig<HalfEven> {
    let wframe = residual_wframe(theta, projective_gates);
    let u = DOmegaUnitary::from_gates(correction_gates);
    let re_w = wframe.re_w(u.z());
    diagonal_diamond_distance(&re_w)
}

/// Like [`residual_diamond_error`], but for a correction that's itself a
/// [`MixedDiagonalResult`] (mixed fallback's twirled-branch corrections) rather than a single
/// gate string. Delegates to
/// [`MixedDiagonalResult::achieved_diamond_error_with_frame`] so the mixture's own quadratic
/// error cancellation (see that method's docs) is preserved, instead of naively
/// triangle-inequality-summing each individual branch's (much larger) distance to the
/// residual target.
pub(crate) fn residual_diamond_error_mixed(
    theta: &FBig<HalfEven>,
    projective_gates: &[Gate],
    correction: &crate::protocol::mixed_diagonal::MixedDiagonalResult,
) -> FBig<HalfEven> {
    let wframe = residual_wframe(theta, projective_gates);
    correction.achieved_diamond_error_with_frame(&wframe)
}

impl AchievedDiamondError for FallbackResult {
    /// Triangle-inequality upper bound on the *whole protocol's* diamond-norm distance to
    /// `theta`: `p_success * dist_phase(projective, theta) + (1 - p_success) *
    /// dist(correction, residual target)`, where `p_success` is
    /// [`FallbackResult::achieved_success_probability`].
    ///
    /// The success term compares the projective candidate's *phase* (`z/|z|`) to `theta`, not
    /// its raw (magnitude-deficient) top-left entry -- see
    /// [`crate::accuracy::achieved_phase_diamond_error`]'s docs for why that distinction
    /// matters for this family. The failure term compares `correction_gates` to the residual
    /// angle it actually approximates (see [`residual_diamond_error`]), not to `theta` and not
    /// composed with the projective step. This is a valid upper bound (diamond norm is a
    /// proper norm, and this is a probabilistic mixture of two unitary channels against a
    /// fixed target), not necessarily the tightest possible one.
    fn achieved_diamond_error(&self, theta: &FBig<HalfEven>) -> FBig<HalfEven> {
        let p_success = self.achieved_success_probability();
        let success_dist = achieved_phase_diamond_error(theta, &self.projective_gates);
        let failure_dist =
            residual_diamond_error(theta, &self.projective_gates, &self.correction_gates);

        let one = ib_to_bf_prec(IBig::ONE);
        let one_minus_p = fb_with_prec(&one - &p_success);
        fb_with_prec(
            fb_with_prec(&p_success * &success_dist) + fb_with_prec(&one_minus_p * &failure_dist),
        )
    }
}

/// Synthesizes a fallback (projective) approximation to `R_z(theta)` within diamond-norm
/// budget `epsilon_diamond`, following Bocharov-Roetteler-Svore's fallback protocol
/// (arXiv:1409.3552) as adapted by Kliuchnikov et al. (arXiv:2203.10064v2, Prop 3.9).
///
/// `sin_alpha` (the sine of the sector's angular half-width) is exposed as an explicit
/// parameter rather than hardcoded from `epsilon_diamond` -- for plain fallback, callers
/// typically want `sin_alpha` set to half of this crate's operator-norm-style epsilon (i.e.
/// `diamond_to_spec_epsilon(epsilon_diamond) / 2`, mirroring `EpsilonRegion`'s own `epsilon
/// / 2` factor), but leaving it as a parameter lets a later "mixed fallback" stage reuse
/// this same [`SectorRegion`]/`synth_fallback` machinery with a different `sin_alpha`.
///
/// Returns `None` if the projective or correction search exceeds its own internal bound
/// without finding a solution (propagated from `search_for_solution`), rather than
/// panicking -- a `SectorRegion` search failing within budget is a more "expected" outcome
/// than the plain single-candidate `EpsilonRegion` search failing.
pub fn synth_fallback(
    theta: f64,
    epsilon_diamond: f64,
    q: DRootTwo,
    sin_alpha: f64,
    seed: u64,
    verbose: bool,
) -> Option<FallbackResult> {
    let mut config = config_from_theta_epsilon(theta, epsilon_diamond, seed, verbose, false);

    let eps_diamond_fbig = config.epsilon.clone();
    let epsilon_spec = diamond_to_spec_epsilon(&eps_diamond_fbig);

    let sin_alpha_fbig = to_fbig(sin_alpha);
    let exact_scale = ZRootTwo::new(IBig::from(1), IBig::from(0));

    // Projective step: find a single candidate inside the annulus sector.
    let sector_region = SectorRegion::new(
        &config.theta,
        q.clone(),
        sin_alpha_fbig,
        exact_scale.clone(),
    );
    let unit_disk = UnitDisk::new(exact_scale.clone());
    let transform = setup_regions_and_transform(
        &sector_region,
        &unit_disk,
        config.verbose,
        config.measure_time,
    );
    let projective_unitary = search_for_solution(
        &sector_region,
        &unit_disk,
        &transform,
        &mut config,
        PhaseMode::Exact,
        None,
    )?;

    let v = projective_unitary.w().clone();
    let projective_gates = decompose_domega_unitary(projective_unitary);

    // Correction step: residual angle theta_B = theta - Arg(v), via the half-angle algebra
    // derived in this module's docs (avoids atan2).
    let (cos_phi, sin_phi, v_norm_sq) = phase_cos_sin(&v);
    let (cos_half_phi, sin_half_phi) = half_angle_cos_sin(&cos_phi, &sin_phi);

    let two = to_fbig(2.0);
    let neg_theta_half = -fb_with_prec(&config.theta / &two);
    let z_x = fb_with_prec(cos_fbig(&neg_theta_half));
    let z_y = fb_with_prec(sin_fbig(&neg_theta_half));

    // cos(-theta_B/2) = cos(A+B) = Z_X*cos(phi/2) - Z_Y*sin(phi/2)
    // sin(-theta_B/2) = sin(A+B) = Z_Y*cos(phi/2) + Z_X*sin(phi/2)
    let cos_neg_theta_b_half =
        fb_with_prec(fb_with_prec(&z_x * &cos_half_phi) - fb_with_prec(&z_y * &sin_half_phi));
    let sin_neg_theta_b_half =
        fb_with_prec(fb_with_prec(&z_y * &cos_half_phi) + fb_with_prec(&z_x * &sin_half_phi));

    let epsilon_for_correction = fb_with_prec(fb_with_prec(&epsilon_spec / &two) / &v_norm_sq);

    let correction_region = EpsilonRegion::from_target_direction(
        cos_neg_theta_b_half,
        sin_neg_theta_b_half,
        epsilon_for_correction,
        exact_scale.clone(),
    );
    let correction_unit_disk = UnitDisk::new(exact_scale.clone());
    let correction_transform = setup_regions_and_transform(
        &correction_region,
        &correction_unit_disk,
        config.verbose,
        config.measure_time,
    );
    let correction_unitary = search_for_solution(
        &correction_region,
        &correction_unit_disk,
        &correction_transform,
        &mut config,
        PhaseMode::Exact,
        None,
    )?;
    let correction_gates = decompose_domega_unitary(correction_unitary);

    Some(FallbackResult {
        projective_gates,
        correction_gates,
        q,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::reset_prec_bits;
    use crate::ring::ZOmega;
    use dashu_base::Abs;
    use serial_test::serial;
    use std::f64::consts::PI;

    fn approx_eq(a: &FBig<HalfEven>, b: &FBig<HalfEven>, tol_bits: usize) -> bool {
        let diff = (a - b).abs();
        let tol = ib_to_bf_prec(IBig::ONE) / ib_to_bf_prec(IBig::ONE << tol_bits);
        diff <= tol
    }

    // ---- Task 1: exact_q ----

    #[test]
    #[serial]
    fn exact_q_is_representable_and_correct() {
        reset_prec_bits();
        let q = exact_q(7);
        let real = q.to_real();
        let expected = to_fbig(1.0 - 1.0 / 128.0);
        assert!(
            approx_eq(&real, &expected, 200),
            "exact_q(7).to_real() = {real}, expected {expected}"
        );
        // Sanity: 1 - q <= 0.01, matching the paper's bound.
        let one_minus_q = fb_with_prec(ib_to_bf_prec(IBig::ONE) - &real);
        assert!(one_minus_q <= to_fbig(0.01));
    }

    // ---- Task 2: sector ellipse containment ----

    // Verifies the ellipse-construction derivation (Task 2 of the design doc) by sampling
    // points across the TRUE annulus-sector region (both the inner and outer arcs, and the
    // two angular rays at several radii) using plain f64 geometry, and confirming every
    // sampled point lies inside the constructed bounding ellipse. `theta = 0` is used so
    // that the region's rotated (z_x, z_y) = (1, 0) frame coincides with world coordinates,
    // which simplifies generating the sample points without changing what's being tested
    // (the region is rotation-covariant: theta only rotates the whole picture).
    #[test]
    #[serial]
    fn sector_ellipse_contains_true_region_samples() {
        reset_prec_bits();
        let alpha = 0.2_f64;
        let sin_alpha_f64 = alpha.sin();
        let q = exact_q(7);
        let scale = ZRootTwo::from_int(IBig::from(1));

        let theta = ib_to_bf_prec(IBig::ZERO);
        let region = SectorRegion::new(&theta, q, to_fbig(sin_alpha_f64), scale);
        let ellipse = region.ellipse();

        let sqrt_s = 1.0_f64;
        let q_real = 1.0_f64 - 1.0_f64 / 128.0_f64;
        let sqrt_qs = q_real.sqrt();

        let n_r = 6;
        let n_theta = 9;
        for i in 0..n_r {
            let t = i as f64 / (n_r - 1) as f64;
            let r = sqrt_qs + t * (sqrt_s - sqrt_qs);
            for j in 0..n_theta {
                let s = j as f64 / (n_theta - 1) as f64;
                let th = -alpha + s * 2.0 * alpha;
                let x = r * th.cos();
                let y = r * th.sin();
                let v = Vector2::new(to_fbig(x), to_fbig(y));
                assert!(
                    ellipse.inside(&v),
                    "sample point r={r}, theta={th} (x={x}, y={y}) is outside the bounding ellipse"
                );
            }
        }
    }

    // ---- Task 2: emptiness check in `intersect` ----

    // Hand-picked u0/v pair whose three sequential half-plane clips in `SectorRegion::intersect`
    // empty the running (t0, t1) interval: `u0` sits just inside the sector near one edge of
    // the angular wedge, and `v` is chosen so the line immediately leaves through the
    // *opposite* angular boundary before it can re-enter -- i.e. the (a)/(b) clips (or (b)/(c))
    // narrow the interval past emptiness. This directly exercises the "must return None,
    // not an inverted pair" requirement.
    #[test]
    #[serial]
    fn intersect_returns_none_on_emptied_interval() {
        reset_prec_bits();
        let alpha = 0.05_f64; // narrow wedge, easy to miss entirely with a generic direction
        let sin_alpha_f64 = alpha.sin();
        let q = exact_q(7);
        let scale = ZRootTwo::from_int(IBig::from(1));
        let theta = ib_to_bf_prec(IBig::ZERO);
        let region = SectorRegion::new(&theta, q, to_fbig(sin_alpha_f64), scale);

        // u0 well outside the sector (e.g. far along the negative real axis), v chosen so the
        // line through u0 in direction v never crosses the narrow wedge around the positive
        // real axis at all -- the outer-disc interval from the quadratic is nonempty, but the
        // angular clips must annihilate it.
        let u0 = DOmega::new(
            ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::from(-100)),
            0,
        );
        let v = DOmega::new(
            ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
            0,
        );

        assert_eq!(region.intersect(&u0, &v), None);
    }

    // ---- Task 3: half-angle correction algebra ----

    #[test]
    #[serial]
    fn half_angle_handles_phi_equals_pi() {
        reset_prec_bits();
        let cos_phi = to_fbig(-1.0);
        let sin_phi = to_fbig(0.0);
        let (cos_half, sin_half) = half_angle_cos_sin(&cos_phi, &sin_phi);
        assert!(approx_eq(&cos_half, &to_fbig(0.0), 200));
        assert!(approx_eq(&sin_half, &to_fbig(1.0), 200));
    }

    #[test]
    #[serial]
    fn half_angle_round_trips_for_generic_angles() {
        reset_prec_bits();
        for phi in [0.0_f64, 0.3, 1.0, 2.0, -0.7, -2.5, 3.0] {
            let cos_phi = to_fbig(phi.cos());
            let sin_phi = to_fbig(phi.sin());
            let (cos_half, sin_half) = half_angle_cos_sin(&cos_phi, &sin_phi);

            // Double-angle reconstruction.
            let reconstructed_cos = fb_with_prec(
                fb_with_prec(&cos_half * &cos_half) - fb_with_prec(&sin_half * &sin_half),
            );
            let reconstructed_sin =
                fb_with_prec(to_fbig(2.0) * fb_with_prec(&cos_half * &sin_half));

            assert!(
                approx_eq(&reconstructed_cos, &cos_phi, 40),
                "phi={phi}: reconstructed cos={reconstructed_cos}, expected {cos_phi}"
            );
            assert!(
                approx_eq(&reconstructed_sin, &sin_phi, 40),
                "phi={phi}: reconstructed sin={reconstructed_sin}, expected {sin_phi}"
            );
        }
    }

    // ---- Task 4: absolute oracle + success-probability guarantee ----

    #[test]
    #[serial]
    fn fallback_result_meets_success_probability_guarantee() {
        crate::clear_caches();
        let q = exact_q(7);
        for (theta, eps_diamond) in [
            (PI / 8.0, 1e-4),
            (PI / 3.0, 1e-6),
            (1.0_f64, 1e-4),
            (-0.5_f64, 1e-5),
        ] {
            crate::clear_caches();
            let sin_alpha = eps_diamond / 4.0;
            let result = synth_fallback(theta, eps_diamond, q.clone(), sin_alpha, 42, false)
                .expect("expected a solution within budget");

            // achieved_success_probability must meet the region's own guarantee.
            let q_real = result.q.to_real();
            let achieved = result.achieved_success_probability();
            assert!(
                achieved >= q_real,
                "achieved success probability {achieved} is below q floor {q_real}"
            );
        }
    }

    // ---- Task 4: slope-fit ----

    // Fits expected T-count (projective + failure-weighted correction) against
    // log2(1/epsilon_diamond), and reports the measured slope. The paper's plain fallback
    // protocol should achieve a slope near 1.03 (vs plain diagonal synthesis's ~3.02).
    #[test]
    #[serial]
    fn fallback_expected_cost_slope() {
        crate::clear_caches();
        let q = exact_q(7);
        let epsilons: [f64; 3] = [1e-4, 1e-6, 1e-8];
        let n_angles = 8;

        let mut xs = Vec::new();
        let mut ys = Vec::new();

        for &eps in &epsilons {
            let sin_alpha = eps / 4.0;
            for i in 0..n_angles {
                // Deterministic pseudo-random angles (avoids pulling in a new RNG
                // dependency for a test): spread across (0, 2*pi) with an irrational-ish
                // stride so consecutive angles don't land on special values.
                let theta = 0.37 + (i as f64) * 0.91 + eps.log10();
                let theta = theta % (2.0 * PI);

                crate::clear_caches();
                let Some(result) =
                    synth_fallback(theta, eps, q.clone(), sin_alpha, 7 + i as u64, false)
                else {
                    continue;
                };

                let projective_t = result.projective_gates.t_count();
                let correction_t = result.correction_gates.t_count();

                let p_f64 = match result.achieved_success_probability().to_f64() {
                    dashu_base::Approximation::Exact(v) => v,
                    dashu_base::Approximation::Inexact(v, _) => v,
                };
                let expected_cost = projective_t as f64 + (1.0 - p_f64) * correction_t as f64;

                xs.push((1.0 / eps).log2());
                ys.push(expected_cost);
            }
        }

        assert!(
            xs.len() >= 10,
            "too few successful data points ({}) to fit a slope",
            xs.len()
        );

        // Simple least-squares linear fit.
        let n = xs.len() as f64;
        let mean_x = xs.iter().sum::<f64>() / n;
        let mean_y = ys.iter().sum::<f64>() / n;
        let mut cov = 0.0;
        let mut var_x = 0.0;
        for i in 0..xs.len() {
            cov += (xs[i] - mean_x) * (ys[i] - mean_y);
            var_x += (xs[i] - mean_x) * (xs[i] - mean_x);
        }
        let slope = cov / var_x;

        eprintln!(
            "fallback_expected_cost_slope: measured slope = {slope:.4} over {} points \
             (expected near 1.03, vs plain diagonal's ~3.02)",
            xs.len()
        );

        // Sanity bound only -- this is a coarse fit over a handful of angles/epsilons, not
        // a tight statistical claim. The real reporting happens via the eprintln! above.
        assert!(
            slope > 0.0 && slope < 3.0,
            "measured slope {slope} is not even qualitatively better than plain diagonal \
             synthesis's ~3.02 slope"
        );
    }
}
