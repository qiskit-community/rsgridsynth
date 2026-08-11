// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

//! Standalone math helpers for the "mixed diagonal" / "mixed fallback" rotation-synthesis
//! protocols of Kliuchnikov et al. (arXiv:2203.10064v2).
//!
//! The core idea of those protocols: instead of finding one candidate `u` close to a
//! target `R_z(theta)`, find two candidates `u_lo`, `u_hi` that straddle the target (one
//! under-rotated, one over-rotated) and mix them with a classical probability `p` chosen
//! so that the first-order rotation error cancels exactly, leaving only a
//! quadratically-smaller residual. This module provides the arithmetic building blocks
//! (rotated-frame projection, the closed-form mixture weight, and diamond-norm error
//! evaluation) that a later stage will use to actually build those protocols. It does not
//! itself implement the mixed-diagonal or mixed-fallback region/synthesis logic.

use crate::common::{fb_with_prec, ib_to_bf_prec};
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::ops::Abs;
use dashu_int::IBig;

// `WFrame` and `diagonal_diamond_distance` live in `crate::accuracy` now (shared by both this
// protocol module and the plain single-candidate synthesis path in `crate::gridsynth`), but
// re-exported here since every existing internal call site in `protocol/*.rs` refers to them
// as `crate::protocol::mixing::{WFrame, diagonal_diamond_distance}`.
pub use crate::accuracy::{achieved_diagonal_diamond_error, diagonal_diamond_distance, WFrame};

/// The result of [`mixture_weight`]: the classical mixing probability `p` and the
/// resulting projective-step diamond-norm error achieved by mixing.
#[derive(Debug, Clone)]
pub struct MixtureWeight {
    /// Probability of using the "lo" (under-rotation) branch; `1 - p` is the probability of
    /// using the "hi" (over-rotation) branch.
    pub p: FBig<HalfEven>,
    /// The achieved projective-step diamond-norm error after mixing, `2*(p*Im(w_lo)^2 +
    /// (1-p)*Im(w_hi)^2)`.
    pub projective_diamond_error: FBig<HalfEven>,
}

/// Closed-form classical mixture weight from the "mixed diagonal" protocol's mixing
/// theorem (Kliuchnikov et al., arXiv:2203.10064v2). Given two straddling candidates for a
/// single diagonal Z-rotation, expressed only through `w := u * e^{i theta/2}` (see
/// [`WFrame`]), returns the probability `p` of using the "lo" branch that makes the
/// first-order rotation error cancel exactly, together with the resulting (quadratically
/// smaller) projective-step diamond-norm error.
///
/// Closed form (using `q*sin(2*delta) = 2*Re(w)*Im(w)` and `q*sin(delta)^2 = Im(w)^2`,
/// where `q = |u|^2` and `delta` is the argument of `w`):
/// ```text
/// p     = Re(w_hi)*Im(w_hi) / (Re(w_hi)*Im(w_hi) - Re(w_lo)*Im(w_lo))
/// error = 2 * (p*Im(w_lo)^2 + (1-p)*Im(w_hi)^2)
/// ```
///
/// # Precondition (caller's responsibility -- required, not enforced beyond a
/// `debug_assert`)
/// `w_lo` MUST be the under-rotation branch and `w_hi` the over-rotation branch, i.e.
/// `Im(w_lo) <= 0 <= Im(w_hi)` (equivalently `sin(delta_lo) <= 0 <= sin(delta_hi)`). This
/// module deliberately does NOT resolve which physical candidate is "lo" vs "hi": this
/// crate's `ZOmega` ring is documented inconsistently elsewhere as both `omega = e^{i
/// pi/4}` and `omega = e^{-i pi/4}`, which could flip that assignment. Callers (a future
/// stage) must determine the correct branch assignment for whatever sign convention the
/// rest of the crate settles on; this function only consumes the already-disambiguated
/// `(Re(w), Im(w))` pairs.
///
/// # Degenerate exact solutions
/// If `Im(w_lo) == 0` or `Im(w_hi) == 0`, that branch IS an exact solution (zero rotation
/// error): mixing is unnecessary, and evaluating the closed form directly would be `0/0`.
/// In these cases this returns all the weight on the exact branch (`p = 1` or `p = 0`) with
/// zero error, instead of dividing by zero.
///
/// Returns `None` only when neither `Im(w)` is zero but the closed form's denominator is
/// zero anyway -- a genuinely degenerate, ill-posed pair that callers should not construct
/// (rather than returning a bogus/NaN-equivalent `p`).
pub fn mixture_weight(
    w_lo: (&FBig<HalfEven>, &FBig<HalfEven>),
    w_hi: (&FBig<HalfEven>, &FBig<HalfEven>),
) -> Option<MixtureWeight> {
    let (re_lo, im_lo) = w_lo;
    let (re_hi, im_hi) = w_hi;

    let zero = ib_to_bf_prec(IBig::ZERO);
    debug_assert!(
        im_lo <= &zero,
        "mixture_weight precondition violated: Im(w_lo) must be <= 0 (w_lo must be the \
         under-rotation branch)"
    );
    debug_assert!(
        im_hi >= &zero,
        "mixture_weight precondition violated: Im(w_hi) must be >= 0 (w_hi must be the \
         over-rotation branch)"
    );

    // Degenerate case: one branch is already an exact solution (zero rotation error), so
    // mixing is unnecessary and the closed form below would divide 0/0.
    if im_lo.repr().is_zero() {
        return Some(MixtureWeight {
            p: ib_to_bf_prec(IBig::ONE),
            projective_diamond_error: zero,
        });
    }
    if im_hi.repr().is_zero() {
        return Some(MixtureWeight {
            p: zero.clone(),
            projective_diamond_error: zero,
        });
    }

    let hi_cross = fb_with_prec(re_hi * im_hi);
    let lo_cross = fb_with_prec(re_lo * im_lo);
    let denom = fb_with_prec(&hi_cross - &lo_cross);
    if denom.repr().is_zero() {
        // Neither branch is exact, yet the closed form is 0/0: a genuinely degenerate,
        // ill-posed pair. Refuse to guess rather than divide by zero.
        return None;
    }
    let p = fb_with_prec(&hi_cross / &denom);

    let one = ib_to_bf_prec(IBig::ONE);
    let one_minus_p = fb_with_prec(&one - &p);
    let im_lo_sq = fb_with_prec(im_lo * im_lo);
    let im_hi_sq = fb_with_prec(im_hi * im_hi);
    let lo_term = fb_with_prec(&p * &im_lo_sq);
    let hi_term = fb_with_prec(&one_minus_p * &im_hi_sq);
    let sum_terms = fb_with_prec(&lo_term + &hi_term);
    let two = fb_with_prec(FBig::try_from(2.0).unwrap());
    let projective_diamond_error = fb_with_prec(&two * &sum_terms);

    Some(MixtureWeight {
        p,
        projective_diamond_error,
    })
}

/// General Pauli-channel diamond-norm closed form: `||E1 - E2||_diamond = sum_P |q_P -
/// r_P|` over the four Pauli components `(I, X, Y, Z)`. Kept generic/simple for a future
/// stage's mixed-diagonal Pauli-channel error computation.
pub fn pauli_diamond_distance(a: &[FBig<HalfEven>; 4], b: &[FBig<HalfEven>; 4]) -> FBig<HalfEven> {
    let mut sum = ib_to_bf_prec(IBig::ZERO);
    for (x, y) in a.iter().zip(b.iter()) {
        let diff = fb_with_prec(x - y);
        sum = fb_with_prec(&sum + diff.abs());
    }
    sum
}

/// Converts a diamond-norm error budget `eps_diamond` into this crate's operator-norm-style
/// `epsilon` convention (the one already used by `EpsilonRegion`/`GridSynthConfig::epsilon`),
/// via the bridge `||U - V||_diamond <= 2*min(||U-V||, ||U+V||)`, i.e.
/// `eps_diamond ~= 2*eps_spec` to first order.
///
/// This is the single place this conversion happens. It is a first-order/linearized bridge
/// in general; it is exact only for the diagonal case, per [`diagonal_diamond_distance`]'s
/// closed form. A future stage implementing the mixed-diagonal/mixed-fallback protocols
/// must derive its working precision from the *tighter* of the requested diamond epsilon
/// and any per-branch spec epsilon it computes downstream -- that budgeting is out of scope
/// here.
pub fn diamond_to_spec_epsilon(eps_diamond: &FBig<HalfEven>) -> FBig<HalfEven> {
    let two = fb_with_prec(FBig::try_from(2.0).unwrap());
    fb_with_prec(eps_diamond / &two)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{cos_fbig, reset_prec_bits, sin_fbig};
    use crate::gridsynth::EpsilonRegion;
    use crate::math::sqrt_fbig;
    use crate::ring::{DOmega, ZOmega, ZRootTwo};
    use crate::tdgp::Region;
    use serial_test::serial;

    // NOTE: `crate::common::PREC_BITS` is a single process-global atomic (see
    // `src/common.rs`). Some tests elsewhere in this crate (e.g. anything that calls
    // `config_from_theta_epsilon`) permanently change it and never restore the default.
    // Every test below that does FBig arithmetic at meaningful precision therefore starts
    // by calling `reset_prec_bits()` to pin a known precision regardless of what leaked in
    // from a previously-run test, and is marked `#[serial]` to match this crate's
    // convention (see `tests/integration_test.rs`) for tests sensitive to that global state.

    fn to_fbig(x: f64) -> FBig<HalfEven> {
        FBig::<HalfEven>::try_from(x)
            .unwrap()
            .with_precision(crate::common::get_prec_bits())
            .value()
    }

    fn approx_eq(a: &FBig<HalfEven>, b: &FBig<HalfEven>, tol_bits: usize) -> bool {
        let diff = (a - b).abs();
        let tol = ib_to_bf_prec(IBig::ONE) / ib_to_bf_prec(IBig::ONE << tol_bits);
        diff <= tol
    }

    fn sample_domegas() -> Vec<DOmega> {
        vec![
            DOmega::from_int(IBig::ONE),
            DOmega::new(
                ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
                0,
            ),
            DOmega::new(ZOmega::new(IBig::ONE, IBig::ONE, IBig::ZERO, IBig::ONE), 0),
            DOmega::new(
                ZOmega::new(IBig::from(2), IBig::from(-1), IBig::from(3), IBig::ONE),
                0,
            ),
            DOmega::new(
                ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ONE, IBig::ZERO),
                2,
            ),
        ]
    }

    // B1 mandatory sanity check: w = u * e^{i theta/2} is a unit-modulus phase rotation of
    // u, so it must preserve the modulus exactly: Re(w)^2 + Im(w)^2 == Re(u)^2 + Im(u)^2.
    #[test]
    #[serial]
    fn wframe_preserves_modulus() {
        reset_prec_bits();
        for theta_f64 in [0.0_f64, 0.3, 1.0, std::f64::consts::PI / 2.0, -2.1, 5.9] {
            let theta = to_fbig(theta_f64);
            let frame = WFrame::new(&theta);
            for u in sample_domegas() {
                let re_w = frame.re_w(&u);
                let im_w = frame.im_w(&u);
                let lhs = fb_with_prec(fb_with_prec(&re_w * &re_w) + fb_with_prec(&im_w * &im_w));
                let rhs = fb_with_prec(
                    fb_with_prec(u.real() * u.real()) + fb_with_prec(u.imag() * u.imag()),
                );
                assert!(
                    approx_eq(&lhs, &rhs, 200),
                    "modulus not preserved for theta={theta_f64}: |w|^2={lhs}, |u|^2={rhs}"
                );
            }
        }
    }

    // Cross-check against `EpsilonRegion`'s own inline `cos_similarity` derivation: pick a
    // `scale` large enough that the norm-scale half of `EpsilonRegion::inside` is always
    // satisfied for our small sample points, so `inside(u)` reduces exactly to
    // `cos_similarity >= d`. We recompute `d` with the same formula as
    // `EpsilonRegion::new` and confirm `WFrame::re_w(u) >= d` agrees with `inside(u)`.
    #[test]
    #[serial]
    fn wframe_re_w_matches_epsilon_region_cos_similarity() {
        reset_prec_bits();
        let theta = to_fbig(0.7);
        let epsilon = to_fbig(0.5);
        let scale = ZRootTwo::from_int(IBig::from(1_000_000));

        let region = EpsilonRegion::new(theta.clone(), epsilon.clone(), scale.clone());
        let frame = WFrame::new(&theta);

        let one = ib_to_bf_prec(IBig::ONE);
        let four = to_fbig(4.0);
        let eps_sq = fb_with_prec(&epsilon * &epsilon);
        let half_eps_sq = fb_with_prec(&eps_sq / &four);
        let one_minus_half_eps_sq = fb_with_prec(&one - &half_eps_sq);
        let d = fb_with_prec(sqrt_fbig(&one_minus_half_eps_sq) * sqrt_fbig(&scale.to_real()));

        for u in sample_domegas() {
            let re_w = frame.re_w(&u);
            let expected_inside = re_w >= d;
            assert_eq!(
                region.inside(&u),
                expected_inside,
                "WFrame::re_w disagreed with EpsilonRegion::inside for u=({}, {})",
                u.real(),
                u.imag()
            );
        }
    }

    #[test]
    #[serial]
    fn mixture_weight_degenerate_lo_exact() {
        reset_prec_bits();
        let re_lo = to_fbig(1.0);
        let im_lo = to_fbig(0.0);
        let re_hi = to_fbig(0.98);
        let im_hi = to_fbig(0.02);

        let result = mixture_weight((&re_lo, &im_lo), (&re_hi, &im_hi))
            .expect("degenerate exact case must not return None");
        assert_eq!(result.p, ib_to_bf_prec(IBig::ONE));
        assert_eq!(result.projective_diamond_error, ib_to_bf_prec(IBig::ZERO));
    }

    #[test]
    #[serial]
    fn mixture_weight_degenerate_hi_exact() {
        reset_prec_bits();
        let re_lo = to_fbig(0.98);
        let im_lo = to_fbig(-0.02);
        let re_hi = to_fbig(1.0);
        let im_hi = to_fbig(0.0);

        let result = mixture_weight((&re_lo, &im_lo), (&re_hi, &im_hi))
            .expect("degenerate exact case must not return None");
        assert_eq!(result.p, ib_to_bf_prec(IBig::ZERO));
        assert_eq!(result.projective_diamond_error, ib_to_bf_prec(IBig::ZERO));
    }

    #[test]
    #[serial]
    fn mixture_weight_zero_denominator_returns_none() {
        reset_prec_bits();
        // Neither branch is exact (both Im(w) nonzero), but Re(w_hi)*Im(w_hi) ==
        // Re(w_lo)*Im(w_lo), so the closed form's denominator is exactly zero.
        let re_lo = to_fbig(2.0);
        let im_lo = to_fbig(-1.0);
        let re_hi = to_fbig(-2.0);
        let im_hi = to_fbig(1.0);

        assert!(mixture_weight((&re_lo, &im_lo), (&re_hi, &im_hi)).is_none());
    }

    // Convention-agnostic sanity check for the generic (non-degenerate) case, using plain
    // numbers rather than ring elements (see module docs: this crate's ZOmega sign
    // convention for `omega` is documented inconsistently, so this test deliberately does
    // not go through DOmega/ZOmega at all).
    //
    // Note on the comparison baseline: the mixed error `2*(p*Im(w_lo)^2 + (1-p)*Im(w_hi)^2)`
    // is *by construction* a convex combination of `2*Im(w_lo)^2` and `2*Im(w_hi)^2` (since
    // the closed-form `p` lies in [0, 1] for straddling branches), so it can never beat
    // *both* of those two quantities simultaneously -- a convex combination of two numbers
    // always lies between them. The real advantage of mixing is that it turns the
    // branches' actual (linear-order) unmixed diamond-norm error `2*sqrt(1-Re(w)^2) ~=
    // 2*|Im(w)|` into a quadratic-order residual `~= 2*Im(w)^2`, which is what we check
    // here using this module's own `diagonal_diamond_distance`.
    #[test]
    #[serial]
    fn mixture_weight_beats_each_branchs_own_unmixed_diamond_error() {
        reset_prec_bits();
        let re_lo = to_fbig(0.9999);
        let im_lo = to_fbig(-0.01);
        let re_hi = to_fbig(0.9999);
        let im_hi = to_fbig(0.02);

        let result = mixture_weight((&re_lo, &im_lo), (&re_hi, &im_hi))
            .expect("generic straddling pair must produce a mixture");

        assert!(
            result.p >= ib_to_bf_prec(IBig::ZERO) && result.p <= ib_to_bf_prec(IBig::ONE),
            "p={} is not a valid probability",
            result.p
        );

        let unmixed_lo = diagonal_diamond_distance(&re_lo);
        let unmixed_hi = diagonal_diamond_distance(&re_hi);
        assert!(
            result.projective_diamond_error < unmixed_lo,
            "mixed error {} should be less than the lo branch's own unmixed diamond error {}",
            result.projective_diamond_error,
            unmixed_lo
        );
        assert!(
            result.projective_diamond_error < unmixed_hi,
            "mixed error {} should be less than the hi branch's own unmixed diamond error {}",
            result.projective_diamond_error,
            unmixed_hi
        );
    }

    #[test]
    #[serial]
    fn diagonal_diamond_distance_zero_at_perfect_match() {
        reset_prec_bits();
        let re_w = ib_to_bf_prec(IBig::ONE);
        let dist = diagonal_diamond_distance(&re_w);
        assert!(approx_eq(&dist, &ib_to_bf_prec(IBig::ZERO), 200));
    }

    #[test]
    #[serial]
    fn diagonal_diamond_distance_matches_closed_form() {
        reset_prec_bits();
        // Re(w) = cos(delta) for delta = 0.05: distance should be 2*sin(0.05).
        let delta = to_fbig(0.05);
        let re_w = fb_with_prec(cos_fbig(&delta));
        let expected = fb_with_prec(to_fbig(2.0) * fb_with_prec(sin_fbig(&delta)));
        let dist = diagonal_diamond_distance(&re_w);
        assert!(
            approx_eq(&dist, &expected, 200),
            "dist={dist}, expected={expected}"
        );
    }

    #[test]
    #[serial]
    fn pauli_diamond_distance_sums_absolute_differences() {
        reset_prec_bits();
        let a = [to_fbig(1.0), to_fbig(0.0), to_fbig(0.0), to_fbig(0.0)];
        let b = [to_fbig(0.0), to_fbig(1.0), to_fbig(0.0), to_fbig(0.0)];
        let dist = pauli_diamond_distance(&a, &b);
        assert!(approx_eq(&dist, &to_fbig(2.0), 200));
    }

    #[test]
    #[serial]
    fn pauli_diamond_distance_zero_for_equal_vectors() {
        reset_prec_bits();
        let a = [to_fbig(0.3), to_fbig(-0.1), to_fbig(0.2), to_fbig(0.4)];
        let dist = pauli_diamond_distance(&a, &a.clone());
        assert!(approx_eq(&dist, &ib_to_bf_prec(IBig::ZERO), 200));
    }

    #[test]
    #[serial]
    fn diamond_to_spec_epsilon_halves() {
        reset_prec_bits();
        let eps_diamond = to_fbig(0.02);
        let eps_spec = diamond_to_spec_epsilon(&eps_diamond);
        assert!(approx_eq(&eps_spec, &to_fbig(0.01), 200));
    }
}
