// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

//! Stage 1: the mixed-diagonal region and the single-enumeration straddling-pair search.
//!
//! Implements the "mixed diagonal" rotation-synthesis protocol of Kliuchnikov, Lauter,
//! Minko, Paetznick, Petit, "Shorter quantum circuits via single-qubit gate approximation"
//! (Quantum 7, 1208 (2023), arXiv:2203.10064v2). Instead of finding one candidate close to
//! the target direction, this finds two candidates that straddle it (one under-rotated, one
//! over-rotated) and mixes them with a classical probability chosen so their first-order
//! rotation errors cancel exactly, at roughly half the T-count of the plain single-candidate
//! ("diagonal") protocol for the same diamond-norm accuracy.

use crate::accuracy::{diagonal_diamond_distance, AchievedDiamondError, WFrame};
use crate::common::{cos_fbig, fb_with_prec, get_prec_bits, ib_to_bf_prec, sin_fbig};
use crate::config::{config_from_theta_epsilon, GridSynthConfig};
use crate::diophantine::diophantine_dyadic;
use crate::gridsynth::{process_solution_candidate, setup_regions_and_transform, PhaseMode};
use crate::gridsynth::{UnitDisk, UprightTransform};
use crate::math::{solve_quadratic, sqrt_fbig};
use crate::protocol::mixing::{diamond_to_spec_epsilon, mixture_weight};
use crate::region::Ellipse;
use crate::ring::{DOmega, DRootTwo, ZRootTwo};
use crate::synthesis_of_clifford_t::decompose_domega_unitary;
use crate::tdgp::{solve_tdgp, Region};
use crate::unitary::DOmegaUnitary;

use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::ops::Abs;
use dashu_int::IBig;
use nalgebra::{Matrix2, Vector2};

/// Same 2x2 matrix product helper as `gridsynth::matrix_multiply_2x2` (that one is private to
/// its module and not reachable from here, so it is duplicated verbatim rather than plumbed
/// through a new `pub(crate)` export, per the "only touch `mixed_diagonal.rs`" constraint).
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

/// The relaxed-offset cap region used by the mixed-diagonal protocol's straddling-pair
/// search. Structurally identical to [`crate::gridsynth::EpsilonRegion`] (same ellipse
/// construction pattern, same disc+half-plane `inside`/`intersect` predicates), but built
/// from the *exact* (not asymptotic) closed forms for the offset, radial semi-axis, and
/// tangential semi-axis of a circular cap cut at `Re(w) = d`:
///
/// - `d = sqrt(s) * sqrt(1 - eps/2)` (offset of the cutting line from the origin, along the
///   target direction; note `eps/2`, not the plain-diagonal protocol's `eps^2/4`).
/// - `h = sqrt(s) - d` (radial semi-axis: exact, not the asymptotic `eps^2*sqrt(s)/8`).
/// - `c = sqrt(s * eps/2)` (tangential semi-axis: exact circle-chord half-length at `Re(w) =
///   d`, derived from `sqrt(s - d^2) = sqrt(s - s*(1-eps/2)) = sqrt(s*eps/2)`).
///
/// Both branches of the mixed protocol (under- and over-rotation) share this single region;
/// which branch a solved candidate belongs to is decided later by the sign of `Im(w)`, not
/// by the region itself -- see [`search_for_straddling_pair`].
#[derive(Debug)]
pub struct MixedDiagonalRegion {
    scale: ZRootTwo,
    d: FBig<HalfEven>,
    z_x: FBig<HalfEven>,
    z_y: FBig<HalfEven>,
    ellipse: Ellipse,
}

impl MixedDiagonalRegion {
    pub fn new(theta: &FBig<HalfEven>, epsilon: &FBig<HalfEven>, scale: ZRootTwo) -> Self {
        let two = fb_with_prec(FBig::try_from(2.0).unwrap());
        let theta_half = fb_with_prec(theta / &two);
        let neg_theta_half = -fb_with_prec(theta_half);
        let z_x: FBig<HalfEven> = fb_with_prec(cos_fbig(&neg_theta_half));
        let z_y: FBig<HalfEven> = fb_with_prec(sin_fbig(&neg_theta_half));
        Self::from_target_direction_impl(z_x, z_y, epsilon, scale)
    }

    /// Builds the same region as [`MixedDiagonalRegion::new`], but from the target
    /// direction's half-angle `(cos(-phi/2), sin(-phi/2))` directly, avoiding an `atan2`-style
    /// angle round-trip -- mirrors
    /// [`crate::gridsynth::EpsilonRegion::from_target_direction`]. Used by "mixed fallback"
    /// (a later stage) to build a mixed-diagonal *correction* region for a residual angle
    /// that only exists as an algebraically-derived `(cos, sin)` pair, not a raw `theta`.
    pub(crate) fn from_target_direction(
        z_x: FBig<HalfEven>,
        z_y: FBig<HalfEven>,
        epsilon: &FBig<HalfEven>,
        scale: ZRootTwo,
    ) -> Self {
        Self::from_target_direction_impl(z_x, z_y, epsilon, scale)
    }

    fn from_target_direction_impl(
        z_x: FBig<HalfEven>,
        z_y: FBig<HalfEven>,
        epsilon: &FBig<HalfEven>,
        scale: ZRootTwo,
    ) -> Self {
        let two = fb_with_prec(FBig::try_from(2.0).unwrap());
        let one = ib_to_bf_prec(IBig::ONE);
        let zero = ib_to_bf_prec(IBig::ZERO);
        let half_eps = fb_with_prec(epsilon / &two);
        // `epsilon` >= 2 (e.g. a derived, rescaled correction-step epsilon on a
        // near-degenerate candidate, as in `mixed_fallback::build_side`) is already past the
        // point where any point of the disk fails to qualify; the sane mathematical limit of
        // the formula below is `one_minus_half_eps = 0`, not a negative radicand. Clamp
        // rather than let that panic in `sqrt_fbig`, matching the analogous clamp in
        // `gridsynth::EpsilonRegion::from_target_direction_impl`.
        let one_minus_half_eps = fb_with_prec(&one - &half_eps).max(zero.clone());
        let scale_to_real = scale.to_real();

        // Exact offset, radial semi-axis, and tangential semi-axis -- see struct docs.
        let sqrt_s = sqrt_fbig(&scale_to_real);
        let sqrt_one_minus_half_eps = sqrt_fbig(&one_minus_half_eps);
        let d = fb_with_prec(&sqrt_s * &sqrt_one_minus_half_eps);
        let h = fb_with_prec(&sqrt_s - &d);
        let s_half_eps = fb_with_prec(&scale_to_real * &half_eps);
        let c = sqrt_fbig(&s_half_eps);

        let h_sq = fb_with_prec(&h * &h);
        let c_sq = fb_with_prec(&c * &c);
        let inv_h_sq = fb_with_prec(&one / &h_sq);
        let inv_c_sq = fb_with_prec(&one / &c_sq);

        // Same d1/d2/d3 matrix-product pattern as `EpsilonRegion::new`: d1/d3 rotate into and
        // out of the (radial, tangential) frame, d2 is the diagonal quadratic form in that
        // frame with the radial (thinner) direction first.
        let neg_z_y: FBig<HalfEven> = -fb_with_prec(z_y.clone());
        let d1: Matrix2<FBig<HalfEven>> =
            Matrix2::new(z_x.clone(), neg_z_y.clone(), z_y.clone(), z_x.clone());
        let d2: Matrix2<FBig<HalfEven>> =
            Matrix2::new(inv_h_sq, zero.clone(), zero.clone(), inv_c_sq);
        let d3: Matrix2<FBig<HalfEven>> =
            Matrix2::new(z_x.clone(), z_y.clone(), neg_z_y, z_x.clone());

        let px = fb_with_prec(&d * &z_x);
        let py = fb_with_prec(&d * &z_y);
        let p = Vector2::new(px, py);
        let m1: Matrix2<FBig<HalfEven>> = matrix_multiply_2x2(&d1, &d2);
        let m: Matrix2<FBig<HalfEven>> = matrix_multiply_2x2(&m1, &d3);
        let ellipse = Ellipse::new(m, p);

        Self {
            scale,
            d,
            z_x,
            z_y,
            ellipse,
        }
    }
}

impl Region for MixedDiagonalRegion {
    fn ellipse(&self) -> Ellipse {
        self.ellipse.clone()
    }

    fn inside(&self, u: &DOmega) -> bool {
        let cos_term1 = fb_with_prec(&self.z_x * u.real());
        let cos_term2 = fb_with_prec(&self.z_y * u.imag());
        let cos_similarity = fb_with_prec(&cos_term1 + &cos_term2);

        DRootTwo::from_domega(u.conj() * u) <= DRootTwo::from_zroottwo(self.scale.clone())
            && cos_similarity >= self.d
    }

    fn intersect(&self, u0: &DOmega, v: &DOmega) -> Option<(FBig<HalfEven>, FBig<HalfEven>)> {
        let a = v.conj() * v;
        let b = 2 * (v.conj() * u0);
        let c = u0.conj() * u0 - DOmega::from_zroottwo(&self.scale);
        let vz_term1 = fb_with_prec(&self.z_x * v.real());
        let vz_term2 = fb_with_prec(&self.z_y * v.imag());
        let vz = fb_with_prec(&vz_term1 + &vz_term2);

        let term1 = fb_with_prec(&self.z_x * u0.real());
        let term2 = fb_with_prec(&self.z_y * u0.imag());
        let temp_sub = fb_with_prec(&self.d - &term1);
        let rhs = fb_with_prec(&temp_sub - &term2);
        // t0 <= t1
        let (t0, t1) = solve_quadratic(a.real(), b.real(), c.real())?;
        let zero = fb_with_prec(ib_to_bf_prec(IBig::ZERO));

        if vz > zero {
            let t2 = fb_with_prec(&rhs / &vz);
            Some(if t0 > t2 { (t0, t1) } else { (t2, t1) })
        } else if vz < zero {
            let t2 = fb_with_prec(&rhs / &vz);
            Some(if t1 < t2 { (t0, t1) } else { (t0, t2) })
        } else if rhs <= zero {
            Some((t0, t1))
        } else {
            None
        }
    }
}

/// Result of [`search_for_straddling_pair`].
#[derive(Debug)]
pub(crate) enum StraddleOutcome {
    /// A candidate with `Im(w) == 0` was found: it IS the exact target direction, so no
    /// mixing is needed (or possible).
    Unmixed(DOmegaUnitary),
    /// `(lo, hi)`: `Im(w_lo) <= 0 <= Im(w_hi)`, i.e. under- and over-rotation respectively.
    /// `hi` is boxed purely to keep this enum's variants closer in size (clippy
    /// `large_enum_variant`); there is no semantic difference from an unboxed field.
    Mixed(DOmegaUnitary, Box<DOmegaUnitary>),
    /// Exceeded `max_k` without completing a pair. This should not happen for well-formed
    /// inputs; see [`crate::gridsynth::search_for_solution`] for why.
    NotFound,
}

/// Single-enumeration search for a straddling pair of candidates (one under-rotating, one
/// over-rotating the target direction), partitioned by the sign of `Im(w)` as the `k`-loop's
/// lazy candidate stream is consumed -- deliberately not two separate calls to
/// `search_for_solution`, which would duplicate the (expensive) `solve_tdgp` enumeration.
///
/// Generic over the region type (`A: Region + Debug`, matching the bound
/// `crate::gridsynth`'s own generic driver functions use) rather than hardcoded to
/// `MixedDiagonalRegion`, so a sibling region shape (e.g. the annulus sector used by "mixed
/// fallback") can reuse this same search unchanged -- only the region's `inside`/`intersect`
/// predicates differ; the bucket-by-sign-of-`Im(w)` logic here is region-agnostic.
///
/// `phase_tolerance` bounds how far off-angle (in `|Im(w)|`) an exact-ring-unitary candidate
/// (`|z| == 1`, see below) may be before the "no mixing needed" fast path is allowed to claim
/// it as the final answer -- callers should pass the same spec epsilon used to build `region`.
/// See the fast path's inline comment for why this check exists.
pub(crate) fn search_for_straddling_pair<A: Region + std::fmt::Debug>(
    region: &A,
    unit_disk: &UnitDisk,
    transformed: &UprightTransform,
    config: &mut GridSynthConfig,
    wframe: &WFrame,
    phase_tolerance: &FBig<HalfEven>,
) -> StraddleOutcome {
    // See `gridsynth::search_for_solution` for the rationale behind this bound: it is
    // effectively unreachable for well-formed inputs and exists only to fail loudly (rather
    // than hang) if a `Region` predicate is broken.
    let max_k = 4 * get_prec_bits() as i64;
    let zero = ib_to_bf_prec(IBig::ZERO);

    let mut lo: Option<DOmegaUnitary> = None;
    let mut hi: Option<DOmegaUnitary> = None;

    let mut k = 0;
    while k <= max_k {
        let solutions = solve_tdgp(
            region,
            unit_disk,
            &transformed.op_g,
            &transformed.bbox_a,
            &transformed.bbox_b,
            k,
            config.verbose,
        );

        if let Some(solutions) = solutions {
            for z in solutions {
                if (&z * z.conj()).residue() == 0 {
                    continue;
                }

                let xi = DRootTwo::from_int(IBig::ONE) - DRootTwo::from_domega(z.conj() * &z);

                // Floating-point `im_w` is used both to bucket non-exact candidates into the
                // under-/over-rotation slots below, and (see immediately below) to gate the
                // ring-exactness fast path. A landed-on-the-boundary `im == 0` for bucketing
                // purposes is just rounding noise, so it is folded into the `lo` (`<= 0`)
                // bucket, which still satisfies `mixture_weight`'s `im_lo <= 0 <= im_hi`
                // precondition.
                let im = wframe.im_w(&z);

                // Ring-exact zero-*synthesis*-error check: `xi == 1 - |z|^2` is computed
                // purely from exact ring arithmetic on `z`, so it is genuinely, bit-exactly
                // zero whenever `|z| == 1`, i.e. `z` is already a unitary and needs no `w`
                // correction. This is NOT, by itself, "z is the target direction": `z` only
                // has to satisfy the *search region's* (at loose `epsilon` values, sometimes
                // wide) containment check, not the tight target tolerance. So an exact-ring
                // point can land inside the region while still being measurably off-angle
                // (see issue #8) -- hence the additional `im.abs() <= *phase_tolerance` gate
                // below, using floating `im_w` (the only place that measures the *angle*, as
                // opposed to just the magnitude `xi`) to confirm the phase genuinely matches
                // before short-circuiting on it. `im_w` mixes two independently rounded
                // floating approximations of the same irrational value (`cos_fbig`/`sin_fbig`'s
                // Taylor series for `z_x`/`z_y` vs. `sqrt2()`'s Newton iteration inside
                // `u.real()`/`u.imag()`), so it is not bit-exact even for a genuine exact
                // match -- but that rounding noise is many orders of magnitude below
                // `phase_tolerance` (working precision scales with `-log(epsilon)`), while a
                // genuinely off-angle candidate that merely passed the region's wide
                // containment check sits well above it, so the threshold cleanly separates the
                // two cases without needing bit-exact equality.
                if xi.to_real() == zero && im.clone().abs() <= *phase_tolerance {
                    if let Some(w_val) =
                        diophantine_dyadic(xi.clone(), &mut config.diophantine_data)
                    {
                        return StraddleOutcome::Unmixed(process_solution_candidate(
                            z,
                            w_val,
                            PhaseMode::Exact,
                        ));
                    }
                    // Falls through to the ordinary lo/hi bucketing below in the (should not
                    // normally happen) case that the trivial xi=0 diophantine solve fails.
                }

                if im <= zero {
                    if lo.is_some() {
                        // Slot already filled: don't waste a factoring attempt on it.
                        continue;
                    }
                    if let Some(w_val) = diophantine_dyadic(xi, &mut config.diophantine_data) {
                        lo = Some(process_solution_candidate(z, w_val, PhaseMode::Exact));
                        if let Some(hi_val) = hi {
                            return StraddleOutcome::Mixed(lo.unwrap(), Box::new(hi_val));
                        }
                    }
                } else {
                    if hi.is_some() {
                        continue;
                    }
                    if let Some(w_val) = diophantine_dyadic(xi, &mut config.diophantine_data) {
                        hi = Some(process_solution_candidate(z, w_val, PhaseMode::Exact));
                        if let Some(lo_val) = lo {
                            return StraddleOutcome::Mixed(lo_val, Box::new(hi.unwrap()));
                        }
                    }
                }
            }
        }
        k += 1;
    }

    match (lo, hi) {
        (Some(l), Some(h)) => StraddleOutcome::Mixed(l, Box::new(h)),
        _ => StraddleOutcome::NotFound,
    }
}

/// The four `{I, S, Z, SZ}`-twirl variants of a solved `DOmegaUnitary`: conjugating
/// `U = [[z, ...], [w, ...]]` by a diagonal Clifford `diag(1, e^{i*beta})` for `beta in
/// {0, pi/2, pi, 3pi/2}` leaves `z` (the top-left entry, which is all that matters for the
/// diagonal-rotation approximation) exactly unchanged, and multiplies `w` by `e^{i*beta}`.
/// Since `omega^2` is a primitive 4th root of unity (regardless of this crate's `omega` sign
/// convention), `w.mul_by_omega_power(2*m)` for `m in 0..4` realizes all four phases on `w`
/// as an unordered set.
pub(crate) fn twirl_variants(u: &DOmegaUnitary) -> [DOmegaUnitary; 4] {
    std::array::from_fn(|m| {
        let twirled_w = u.w().mul_by_omega_power(2 * m);
        DOmegaUnitary::new(u.z().clone(), twirled_w, u.n() as usize, Some(u.k()))
    })
}

/// One weighted branch of a [`MixedDiagonalResult`]: a Clifford+T gate string and the
/// classical probability with which the mixed channel applies it.
#[derive(Debug, Clone)]
pub struct MixedDiagonalBranch {
    pub gates: String,
    pub weight: FBig<HalfEven>,
}

/// The output of [`synth_mixed_diagonal`]: a classical mixture of Clifford+T circuits
/// (`branches`, with weights summing to 1) implementing a probabilistic-channel
/// approximation of `R_z(theta)`. Call [`AchievedDiamondError::achieved_diamond_error`] to
/// compute the achieved projective-step diamond-norm error on demand.
#[derive(Debug, Clone)]
pub struct MixedDiagonalResult {
    pub branches: Vec<MixedDiagonalBranch>,
}

impl MixedDiagonalResult {
    /// Recomputes the achieved diamond-norm error to the target direction encoded by `wframe`,
    /// directly from the public `branches`, decoding each gate string back into a unitary.
    /// [`AchievedDiamondError::achieved_diamond_error`] is a thin wrapper over this that builds
    /// `wframe` from a raw `theta` via [`WFrame::new`] -- this lower-level entry point exists
    /// so a caller that already has a target direction as a `(cos, sin)` half-angle pair (e.g.
    /// a fallback correction's residual angle) can reuse the exact same mixture-aware
    /// computation without a lossy angle round-trip.
    ///
    /// For the single-branch (exact-angle) case this decodes that one candidate and computes
    /// its diamond distance directly. For the general (8-twirled-branch) case:
    /// [`twirl_variants`] only ever rotates a candidate's `w`, never its `z`, so every branch
    /// within the `lo` family (the first half) and within the `hi` family (the second half)
    /// shares the same `z` and hence the same `Re(w)`/`Im(w)`. This decodes one representative
    /// from each family and reapplies [`mixture_weight`]'s closed form -- naively
    /// triangle-inequality-summing the distance of all 8 individual branches (an earlier,
    /// wrong version of this crate's `mixed_fallback` correction accuracy check made exactly
    /// this mistake) would throw away the quadratic cancellation the mixture achieves, since
    /// each *individual* branch is only as close to the target as the (much looser)
    /// straddling-search tolerance, not the mixture's own tighter budget.
    pub(crate) fn achieved_diamond_error_with_frame(&self, wframe: &WFrame) -> FBig<HalfEven> {
        if self.branches.len() == 1 {
            let u = DOmegaUnitary::from_gates(&self.branches[0].gates);
            let re_w = wframe.re_w(u.z());
            return diagonal_diamond_distance(&re_w);
        }
        assert_eq!(
            self.branches.len() % 2,
            0,
            "a mixed result's branches must split evenly into lo/hi twirl families"
        );
        let half = self.branches.len() / 2;

        let lo_u = DOmegaUnitary::from_gates(&self.branches[0].gates);
        let hi_u = DOmegaUnitary::from_gates(&self.branches[half].gates);
        let re_lo = wframe.re_w(lo_u.z());
        let im_lo = wframe.im_w(lo_u.z());
        let re_hi = wframe.re_w(hi_u.z());
        let im_hi = wframe.im_w(hi_u.z());

        mixture_weight((&re_lo, &im_lo), (&re_hi, &im_hi))
            .expect("a real assembled Mixed result must yield a valid mixture")
            .projective_diamond_error
    }
}

impl AchievedDiamondError for MixedDiagonalResult {
    fn achieved_diamond_error(&self, theta: &FBig<HalfEven>) -> FBig<HalfEven> {
        self.achieved_diamond_error_with_frame(&WFrame::new(theta))
    }
}

/// Turns a [`StraddleOutcome`] into the final weighted branch list.
///
/// For the `Unmixed` (exact-angle) case, this emits a single branch of weight 1 rather than
/// 4 twirl variants of weight 1/4 each: twirling an exact solution (zero rotation error) is
/// harmless but has no error-cancellation purpose, so the single-branch form is the simpler,
/// equally-correct choice.
///
/// For the `Mixed` case, this emits all 8 twirl variants (4 of `lo` at weight `p/4` each, 4
/// of `hi` at weight `(1-p)/4` each), so the classical mixture is invariant under an
/// additional random {Z,S} twirl -- required for the mixture to implement a genuine
/// depolarizing-style probabilistic channel rather than leak phase information.
pub(crate) fn assemble_result(outcome: StraddleOutcome, wframe: &WFrame) -> MixedDiagonalResult {
    match outcome {
        StraddleOutcome::NotFound => panic!(
            "search_for_straddling_pair: exceeded max_k without finding a straddling pair \
             (or an exact solution) -- region predicate is likely incorrect"
        ),
        StraddleOutcome::Unmixed(u) => {
            let gates = decompose_domega_unitary(u);
            MixedDiagonalResult {
                branches: vec![MixedDiagonalBranch {
                    gates,
                    weight: ib_to_bf_prec(IBig::ONE),
                }],
            }
        }
        StraddleOutcome::Mixed(lo, hi) => {
            let hi = *hi;
            let re_lo = wframe.re_w(lo.z());
            let im_lo = wframe.im_w(lo.z());
            let re_hi = wframe.re_w(hi.z());
            let im_hi = wframe.im_w(hi.z());

            let mw = mixture_weight((&re_lo, &im_lo), (&re_hi, &im_hi)).expect(
                "mixture_weight returned None for a real solved straddling pair -- this \
                 indicates a genuine bug (search_for_straddling_pair only ever produces \
                 im_lo <= 0 <= im_hi pairs), not an expected degenerate input",
            );

            let one = ib_to_bf_prec(IBig::ONE);
            let four = fb_with_prec(FBig::try_from(4.0).unwrap());
            let one_minus_p = fb_with_prec(&one - &mw.p);
            let p_over_4 = fb_with_prec(&mw.p / &four);
            let one_minus_p_over_4 = fb_with_prec(&one_minus_p / &four);

            let mut branches = Vec::with_capacity(8);
            for variant in twirl_variants(&lo) {
                branches.push(MixedDiagonalBranch {
                    gates: decompose_domega_unitary(variant),
                    weight: p_over_4.clone(),
                });
            }
            for variant in twirl_variants(&hi) {
                branches.push(MixedDiagonalBranch {
                    gates: decompose_domega_unitary(variant),
                    weight: one_minus_p_over_4.clone(),
                });
            }

            MixedDiagonalResult { branches }
        }
    }
}

/// Synthesizes a mixed-diagonal probabilistic-channel approximation of `R_z(theta)` to
/// diamond-norm accuracy `epsilon_diamond`.
///
/// `epsilon_diamond` is converted to this crate's operator-norm-style `epsilon` convention
/// (via [`diamond_to_spec_epsilon`]) before building the search region. Only exact-phase
/// synthesis (`PhaseMode::Exact`) is implemented at this stage; `up_to_phase` mixing is out
/// of scope.
///
/// # Panics
/// Panics if the internal search exceeds its (very generous) bound on `k` without finding a
/// solution; see [`search_for_straddling_pair`]. Not expected to trigger for any well-formed
/// input.
pub fn synth_mixed_diagonal(
    theta: f64,
    epsilon_diamond: f64,
    seed: u64,
    verbose: bool,
) -> MixedDiagonalResult {
    // `config_from_theta_epsilon` is reused purely as a scaffold: it parses `theta` exactly,
    // and sizes working precision from the decimal magnitude of its `epsilon` argument, which
    // is close enough to the actual (post-conversion) spec epsilon for that purpose. The
    // *value* stored in `config.epsilon` here is still the diamond-norm epsilon; the spec
    // epsilon actually used to build the region is derived from it just below.
    let mut config = config_from_theta_epsilon(theta, epsilon_diamond, seed, verbose, false);
    let epsilon_spec = diamond_to_spec_epsilon(&config.epsilon);

    let scale = ZRootTwo::new(IBig::from(1), IBig::from(0));
    let region = MixedDiagonalRegion::new(&config.theta, &epsilon_spec, scale.clone());
    let unit_disk = UnitDisk::new(scale);
    let wframe = WFrame::new(&config.theta);

    let transformed =
        setup_regions_and_transform(&region, &unit_disk, config.verbose, config.measure_time);

    let outcome = search_for_straddling_pair(
        &region,
        &unit_disk,
        &transformed,
        &mut config,
        &wframe,
        &epsilon_spec,
    );

    assemble_result(outcome, &wframe)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::reset_prec_bits;
    use dashu_base::Approximation;
    use dashu_int::ops::Abs;
    use rand::{rngs::StdRng, Rng, SeedableRng};
    use serial_test::serial;
    use std::f64::consts::PI;

    fn to_fbig(x: f64) -> FBig<HalfEven> {
        FBig::<HalfEven>::try_from(x)
            .unwrap()
            .with_precision(get_prec_bits())
            .value()
    }

    fn fbig_to_f64(x: &FBig<HalfEven>) -> f64 {
        match x.to_f64() {
            Approximation::Inexact(v, _) => v,
            Approximation::Exact(v) => v,
        }
    }

    /// A tolerance (in bits) safely below the *actual* working precision configured by
    /// `config_from_theta_epsilon` for a given epsilon (which allocates working precision
    /// roughly proportional to epsilon's decimal digit count -- e.g. only ~60 bits for
    /// epsilon=1e-5, not this crate's default 1000-bit `PREC_BITS_INITIAL`). Using a fixed,
    /// generous tolerance like 200 bits regardless of `get_prec_bits()` would spuriously fail
    /// whenever the configured working precision is lower than that, even though the
    /// underlying arithmetic identity holds exactly up to genuine rounding.
    fn safe_tol_bits() -> usize {
        get_prec_bits().saturating_sub(30)
    }

    fn approx_eq(a: &FBig<HalfEven>, b: &FBig<HalfEven>, tol_bits: usize) -> bool {
        let diff = (a - b).abs();
        let tol = ib_to_bf_prec(IBig::ONE) / ib_to_bf_prec(IBig::ONE << tol_bits);
        diff <= tol
    }

    /// Builds a `(MixedDiagonalRegion, UnitDisk, UprightTransform, WFrame, GridSynthConfig)`
    /// tuple for a given `(theta, epsilon)` pair, mirroring `synth_mixed_diagonal`'s setup but
    /// exposing the intermediate pieces for tests that need to inspect the raw search.
    fn setup(
        theta_f64: f64,
        epsilon: f64,
        seed: u64,
    ) -> (
        MixedDiagonalRegion,
        UnitDisk,
        UprightTransform,
        WFrame,
        GridSynthConfig,
    ) {
        let config = config_from_theta_epsilon(theta_f64, epsilon, seed, false, false);
        let scale = ZRootTwo::new(IBig::from(1), IBig::from(0));
        let region = MixedDiagonalRegion::new(&config.theta, &config.epsilon, scale.clone());
        let unit_disk = UnitDisk::new(scale);
        let wframe = WFrame::new(&config.theta);
        let transformed =
            setup_regions_and_transform(&region, &unit_disk, config.verbose, config.measure_time);
        (region, unit_disk, transformed, wframe, config)
    }

    // ---- Task 1: ellipse containment ----

    // Samples points on the TRUE boundary of the circular cap (the two chord endpoints where
    // the line `Re(w) = d` meets the circle `|u| = sqrt(s)`, and the midpoint of the outer arc
    // between them) and confirms they lie inside the bounding ellipse. All three points are,
    // by construction (see module docs), exactly on the ellipse boundary in exact arithmetic:
    // the chord endpoints sit at (radial=0, tangential=c) in the (z_x,z_y)-aligned frame, and
    // the arc midpoint sits at (radial=h, tangential=0), each giving quadratic-form value
    // exactly 1. This is precisely the check that would catch an asymptotic-height bug (using
    // eps^2/8 instead of the exact `sqrt(s) - d`), since that would move the arc-midpoint
    // radial coordinate away from `h` and push its quadratic-form value above 1.
    #[test]
    #[serial]
    fn mixed_diagonal_region_boundary_points_inside_ellipse() {
        reset_prec_bits();
        let theta = to_fbig(0.7);
        let epsilon = to_fbig(0.3);
        let scale = ZRootTwo::from_int(IBig::from(1));

        let region = MixedDiagonalRegion::new(&theta, &epsilon, scale.clone());

        let two = to_fbig(2.0);
        let half_eps = fb_with_prec(&epsilon / &two);
        let scale_to_real = scale.to_real();
        let sqrt_s = sqrt_fbig(&scale_to_real);
        let c = sqrt_fbig(&fb_with_prec(&scale_to_real * &half_eps));

        let d = region.d.clone();
        let z_x = region.z_x.clone();
        let z_y = region.z_y.clone();

        // Chord endpoints: (z_x*d -/+ z_y*c, z_y*d +/- z_x*c).
        let ep_plus = Vector2::new(
            fb_with_prec(fb_with_prec(&z_x * &d) - fb_with_prec(&z_y * &c)),
            fb_with_prec(fb_with_prec(&z_y * &d) + fb_with_prec(&z_x * &c)),
        );
        let ep_minus = Vector2::new(
            fb_with_prec(fb_with_prec(&z_x * &d) + fb_with_prec(&z_y * &c)),
            fb_with_prec(fb_with_prec(&z_y * &d) - fb_with_prec(&z_x * &c)),
        );
        // Outer-arc midpoint: sqrt(s) * (z_x, z_y).
        let cap_top = Vector2::new(fb_with_prec(&sqrt_s * &z_x), fb_with_prec(&sqrt_s * &z_y));

        let ellipse = region.ellipse();

        // These three points are, in *exact* arithmetic, precisely ON the ellipse boundary
        // (quadratic-form value exactly 1 -- see the doc comment above): the chord endpoints
        // sit at (radial=0, tangential=c) in the (z_x,z_y)-aligned frame, and the arc midpoint
        // sits at (radial=h, tangential=0). At working precision, independently-rounded sqrt
        // evaluations (for `d`/`c`/`sqrt_s`, computed once inside the region constructor and
        // once more here) can round the value to a few ULP above 1 rather than landing exactly
        // on it, purely as floating-point noise unrelated to region correctness. To avoid that
        // false negative while still exercising the true boundary geometry, nudge each point
        // a minuscule (1e-9 relative) amount towards the ellipse center before testing
        // containment -- utterly negligible next to the kind of gross error (e.g. an
        // asymptotic-height formula off by a factor of ~2) this test is meant to catch.
        let center = ellipse.p.clone();
        let shrink = fb_with_prec(&ib_to_bf_prec(IBig::ONE) - &to_fbig(1e-9));
        let nudge_towards_center = |pt: &Vector2<FBig<HalfEven>>| -> Vector2<FBig<HalfEven>> {
            Vector2::new(
                fb_with_prec(
                    &center[0] + fb_with_prec(&shrink * fb_with_prec(&pt[0] - &center[0])),
                ),
                fb_with_prec(
                    &center[1] + fb_with_prec(&shrink * fb_with_prec(&pt[1] - &center[1])),
                ),
            )
        };
        let ep_plus = nudge_towards_center(&ep_plus);
        let ep_minus = nudge_towards_center(&ep_minus);
        let cap_top = nudge_towards_center(&cap_top);

        assert!(
            ellipse.inside(&ep_plus),
            "chord endpoint (+) not inside bounding ellipse"
        );
        assert!(
            ellipse.inside(&ep_minus),
            "chord endpoint (-) not inside bounding ellipse"
        );
        assert!(
            ellipse.inside(&cap_top),
            "outer-arc midpoint not inside bounding ellipse"
        );

        // Sanity: a point diametrically opposite the cap must NOT be inside.
        let far_side = Vector2::new(
            fb_with_prec(-(&sqrt_s * &z_x)),
            fb_with_prec(-(&sqrt_s * &z_y)),
        );
        assert!(
            !ellipse.inside(&far_side),
            "diametrically opposite point should not be inside"
        );
    }

    // ---- Task 2: single-enumeration straddling-pair search ----

    #[test]
    #[serial]
    fn straddling_search_finds_mixed_pair_for_generic_angles() {
        for k in [1i64, 3, 5, 7] {
            crate::clear_caches();
            let theta_f64 = k as f64 * PI / 32.0;
            let epsilon = 1e-6;
            let (region, unit_disk, transformed, wframe, mut config) =
                setup(theta_f64, epsilon, 42 + k as u64);

            let phase_tolerance = config.epsilon.clone();
            let outcome = search_for_straddling_pair(
                &region,
                &unit_disk,
                &transformed,
                &mut config,
                &wframe,
                &phase_tolerance,
            );
            match outcome {
                StraddleOutcome::Mixed(_, _) => {}
                other => panic!("expected Mixed for theta={theta_f64} (k={k}), got {other:?}"),
            }
        }
    }

    // NOTE on angle choice: the target direction `u = e^{-i*theta/2}` is exactly
    // representable in the ring Z[omega] (omega = e^{i*pi/4}) iff `theta` is a multiple of
    // pi/2 -- e.g. theta=pi/2 gives u=e^{-i*pi/4}=-omega^3 exactly, verified against this
    // crate's own `pi_over_two_test` fixture, which produces the identical trivial 8-gate
    // string "SWWWWWWW" regardless of epsilon or seed (proof of a zero-error, exact
    // solution). theta=pi/4 does NOT have this property (its target direction e^{-i*pi/8}
    // is not an omega power): `pi_over_4_exact_test`'s "exact" in its name refers to
    // `PhaseMode::Exact` (vs. `Shifted`), not to zero rotation error -- that test's long,
    // epsilon-dependent gate string is itself evidence there is no exact solution at
    // theta=pi/4. So this test (and `degenerate_angles_produce_unmixed_result` below) use
    // theta in {pi/2, pi} (both multiples of pi/2) rather than the {pi/2, pi/4} the original
    // task spec suggested -- pi/4 genuinely belongs in the "Mixed" bucket, not "Unmixed".
    #[test]
    #[serial]
    fn straddling_search_finds_unmixed_for_exact_angles() {
        for theta_f64 in [PI / 2.0, PI] {
            crate::clear_caches();
            let epsilon = 1e-6;
            let (region, unit_disk, transformed, wframe, mut config) = setup(theta_f64, epsilon, 7);

            let phase_tolerance = config.epsilon.clone();
            let outcome = search_for_straddling_pair(
                &region,
                &unit_disk,
                &transformed,
                &mut config,
                &wframe,
                &phase_tolerance,
            );
            match outcome {
                StraddleOutcome::Unmixed(_) => {}
                other => panic!("expected Unmixed for theta={theta_f64}, got {other:?}"),
            }
        }
    }

    // ---- Task 3: twirl construction ----

    #[test]
    #[serial]
    fn twirl_variants_preserve_top_left_entry_exactly() {
        crate::clear_caches();
        let theta_f64 = 3.0 * PI / 32.0;
        let epsilon = 1e-6;
        let (region, unit_disk, transformed, wframe, mut config) = setup(theta_f64, epsilon, 99);

        let phase_tolerance = config.epsilon.clone();
        let outcome = search_for_straddling_pair(
            &region,
            &unit_disk,
            &transformed,
            &mut config,
            &wframe,
            &phase_tolerance,
        );
        let (lo, hi) = match outcome {
            StraddleOutcome::Mixed(lo, hi) => (lo, *hi),
            other => panic!("expected Mixed, got {other:?}"),
        };

        for original in [lo, hi] {
            let orig_top_left = original.to_complex_matrix()[(0, 0)].clone();
            for variant in twirl_variants(&original) {
                let variant_top_left = variant.to_complex_matrix()[(0, 0)].clone();
                assert_eq!(
                    variant_top_left, orig_top_left,
                    "twirl variant changed the (0,0) entry"
                );

                let gates = decompose_domega_unitary(variant);
                let reconstructed_top_left =
                    DOmegaUnitary::from_gates(&gates).to_complex_matrix()[(0, 0)].clone();
                assert_eq!(
                    reconstructed_top_left, orig_top_left,
                    "re-multiplying the decomposed gate string changed the (0,0) entry"
                );
            }
        }
    }

    // ---- Task 4: public entry point and validation ----

    // Required test 2: "absolute oracle" -- re-verify, using REAL solved candidates from the
    // search (not synthetic numbers), that `mixture_weight`'s reported error matches its own
    // closed form, and that the mixture strictly improves on each branch's own unmixed
    // diamond-norm error (`crate::protocol::mixing`'s own tests already cover the formula in
    // isolation; this re-checks it end-to-end against this module's real search output).
    #[test]
    #[serial]
    fn absolute_oracle_end_to_end_error_matches_closed_form() {
        crate::clear_caches();
        let theta_f64 = 3.0 * PI / 32.0;
        let epsilon = 1e-5;
        let (region, unit_disk, transformed, wframe, mut config) = setup(theta_f64, epsilon, 55);

        let phase_tolerance = config.epsilon.clone();
        let outcome = search_for_straddling_pair(
            &region,
            &unit_disk,
            &transformed,
            &mut config,
            &wframe,
            &phase_tolerance,
        );
        let (lo, hi) = match outcome {
            StraddleOutcome::Mixed(lo, hi) => (lo, *hi),
            other => panic!("expected Mixed, got {other:?}"),
        };

        let re_lo = wframe.re_w(lo.z());
        let im_lo = wframe.im_w(lo.z());
        let re_hi = wframe.re_w(hi.z());
        let im_hi = wframe.im_w(hi.z());

        assert!(
            im_lo <= ib_to_bf_prec(IBig::ZERO),
            "lo branch must under-rotate"
        );
        assert!(
            im_hi >= ib_to_bf_prec(IBig::ZERO),
            "hi branch must over-rotate"
        );

        let mw = mixture_weight((&re_lo, &im_lo), (&re_hi, &im_hi))
            .expect("mixture_weight should succeed for a real straddling pair");

        // Cross-check the closed form directly: error == 2*(p*im_lo^2 + (1-p)*im_hi^2).
        let one = ib_to_bf_prec(IBig::ONE);
        let two = to_fbig(2.0);
        let one_minus_p = fb_with_prec(&one - &mw.p);
        let lo_term = fb_with_prec(&mw.p * fb_with_prec(&im_lo * &im_lo));
        let hi_term = fb_with_prec(&one_minus_p * fb_with_prec(&im_hi * &im_hi));
        let expected_error = fb_with_prec(&two * fb_with_prec(&lo_term + &hi_term));
        assert!(
            approx_eq(
                &mw.projective_diamond_error,
                &expected_error,
                safe_tol_bits()
            ),
            "mixed error {} != closed form {}",
            mw.projective_diamond_error,
            expected_error
        );

        // Mixing must beat each branch's own unmixed diamond error (quadratic vs linear-order
        // residual), using the real solved re_w values.
        let unmixed_lo = crate::protocol::mixing::diagonal_diamond_distance(&re_lo);
        let unmixed_hi = crate::protocol::mixing::diagonal_diamond_distance(&re_hi);
        assert!(mw.projective_diamond_error < unmixed_lo);
        assert!(mw.projective_diamond_error < unmixed_hi);

        // Round-trip both candidates through decompose/from_gates and confirm the (0,0) entry
        // used in the error computation above matches what actually gets synthesized.
        let gates_lo = decompose_domega_unitary(lo.clone());
        let reconstructed_lo = DOmegaUnitary::from_gates(&gates_lo).to_complex_matrix();
        assert_eq!(reconstructed_lo[(0, 0)], lo.to_complex_matrix()[(0, 0)]);
        let gates_hi = decompose_domega_unitary(hi.clone());
        let reconstructed_hi = DOmegaUnitary::from_gates(&gates_hi).to_complex_matrix();
        assert_eq!(reconstructed_hi[(0, 0)], hi.to_complex_matrix()[(0, 0)]);
    }

    // Required test 4: branch weights sum to 1.
    #[test]
    #[serial]
    fn branch_weights_sum_to_one() {
        for (theta_f64, epsilon) in [(3.0 * PI / 32.0, 1e-5), (PI / 2.0, 1e-5), (PI / 4.0, 1e-5)] {
            crate::clear_caches();
            let result = synth_mixed_diagonal(theta_f64, epsilon, 321, false);
            let mut total = ib_to_bf_prec(IBig::ZERO);
            for branch in &result.branches {
                total = fb_with_prec(&total + &branch.weight);
            }
            assert!(
                approx_eq(&total, &ib_to_bf_prec(IBig::ONE), safe_tol_bits()),
                "branch weights for theta={theta_f64} summed to {total}, not 1"
            );
        }
    }

    // Required test 5: degenerate angles produce a single-total-weight Unmixed result. See
    // the note above `straddling_search_finds_unmixed_for_exact_angles` for why pi/2 and pi
    // (both multiples of pi/2, hence ring-exact target directions) are used here rather than
    // the original task spec's {pi/2, pi/4} -- pi/4 has no exact solution and correctly goes
    // through the Mixed path (already exercised by `branch_weights_sum_to_one` above).
    #[test]
    #[serial]
    fn degenerate_angles_produce_unmixed_result() {
        for theta_f64 in [PI / 2.0, PI] {
            crate::clear_caches();
            let result = synth_mixed_diagonal(theta_f64, 1e-6, 654, false);
            assert_eq!(
                result.branches.len(),
                1,
                "expected the single-branch Unmixed form for theta={theta_f64}"
            );
            assert!(approx_eq(
                &result.branches[0].weight,
                &ib_to_bf_prec(IBig::ONE),
                safe_tol_bits()
            ));
            // `achieved_diamond_error` is computed via `WFrame::re_w`/`diagonal_diamond_distance`,
            // which mixes two independently-rounded floating trig approximations
            // (`cos_fbig`/`sin_fbig`) that do not cancel to bit-exact zero even when the true
            // rotation error is exactly zero (same caveat this module already documents for
            // `im_w`) -- so this must be an approximate, not exact, equality check.
            // `diagonal_diamond_distance`'s `sqrt` roughly halves the number of reliable bits
            // (a `re_w` deviation of `d` near 1 becomes an error of order `sqrt(d)`), so the
            // usual `safe_tol_bits()` is too tight here.
            let theta = to_fbig(theta_f64);
            assert!(approx_eq(
                &result.achieved_diamond_error(&theta),
                &ib_to_bf_prec(IBig::ZERO),
                safe_tol_bits() / 2
            ));
        }
    }

    // Required test 3: slope fit. Computes the branch-weighted T-count for a batch of random
    // (non-pi/4-multiple) angles at a few epsilons spanning a couple of orders of magnitude,
    // and estimates the cost-vs-log2(1/epsilon) slope. Reported in the test's own assertion
    // message (and by the harness driving this test) rather than hard-gated at exactly 1.52,
    // since a modest sample size is inherently noisy; the assertion instead guards against a
    // gross regression (e.g. a slope near the plain-diagonal baseline's ~3.02, which would
    // indicate the mixing isn't buying anything).
    //
    // T-count is computed via a plain `'T'`-char count rather than
    // `NormalForm::from_gates(..).t_count()`: `decompose_domega_unitary`'s output can, for a
    // branch whose trailing Clifford correction happens to be trivial, legitimately contain
    // an embedded literal `'I'` character (via `Clifford::to_gates`'s own "return \"I\" if
    // empty" convention leaking into `NormalForm::to_gates`'s concatenation of the syllable
    // prefix with the Clifford suffix) that `NormalForm::from_gates` cannot parse back
    // (`append_gate` only recognizes H/S/X/W/T) -- a pre-existing gap in `src/normal_form.rs`,
    // which is out of scope to fix here (only `mixed_diagonal.rs` is mine to edit). Since the
    // task spec explicitly allows either method, the char-count form sidesteps it entirely.
    #[test]
    #[serial]
    fn slope_fit_cost_vs_log2_inv_epsilon() {
        let mut rng = StdRng::seed_from_u64(2024);
        let num_angles = 24;
        let mut angles = Vec::with_capacity(num_angles);
        while angles.len() < num_angles {
            let theta: f64 = rng.random_range(0.0..(2.0 * PI));
            let nearest_pi4_multiple = (theta / (PI / 4.0)).round();
            if (theta - nearest_pi4_multiple * PI / 4.0).abs() > 1e-2 {
                angles.push(theta);
            }
        }

        let epsilons = [1e-4, 1e-6, 1e-8];
        let mut mean_cost = Vec::with_capacity(epsilons.len());

        for (eps_idx, &eps) in epsilons.iter().enumerate() {
            let mut total_cost = 0.0;
            for (i, &theta) in angles.iter().enumerate() {
                crate::clear_caches();
                let seed = 10_000 + (eps_idx * 1000 + i) as u64;
                let result = synth_mixed_diagonal(theta, eps, seed, false);
                let mut cost = 0.0;
                for branch in &result.branches {
                    let t_count = branch.gates.chars().filter(|&c| c == 'T').count() as f64;
                    cost += fbig_to_f64(&branch.weight) * t_count;
                }
                total_cost += cost;
            }
            mean_cost.push(total_cost / angles.len() as f64);
        }

        let log2_inv_eps: Vec<f64> = epsilons.iter().map(|e| (1.0_f64 / e).log2()).collect();
        let slope = (mean_cost[mean_cost.len() - 1] - mean_cost[0])
            / (log2_inv_eps[log2_inv_eps.len() - 1] - log2_inv_eps[0]);

        eprintln!(
            "mixed-diagonal slope fit: mean_cost={mean_cost:?} at log2(1/eps)={log2_inv_eps:?} \
             -> measured slope = {slope:.4} (target ~1.52, plain-diagonal baseline ~3.02)"
        );

        assert!(
            slope > 0.5 && slope < 2.5,
            "measured slope {slope:.4} is far from the expected ~1.52 (plain-diagonal baseline \
             is ~3.02); mean_cost={mean_cost:?}"
        );
    }
}
