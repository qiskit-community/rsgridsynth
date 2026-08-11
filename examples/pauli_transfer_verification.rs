// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

//! Independent numerical verification that the synthesized channels actually meet their
//! requested diamond-norm accuracy budget -- entirely independent of this crate's own
//! closed-form `mixture_weight`/`diagonal_diamond_distance` formulas (this check does not
//! call either; it only reads final gate strings and probabilities and recomputes everything
//! from scratch).
//!
//! ## Arbitrary precision throughout, not `f64`
//!
//! Every quantity in this file -- matrix entries, Pauli transfer matrices, angular
//! differences -- is computed as `FBig<HalfEven>` at the crate's own working precision
//! (`rsgridsynth::common::get_prec_bits()`, set by whichever `synth_*` call ran most
//! recently), never downcast to `f64` until the final printed digits. An earlier version of
//! this file converted synthesized matrices to `f64` immediately and computed
//! `2*sqrt(1-Re(w)^2)`: at epsilon ~1e-8 the true angular error can be ~1e-9, at which point
//! `Re(w)` rounds to exactly `1.0` in `f64` (whose resolution near 1.0 is ~2.2e-16), and
//! squaring then subtracting from 1 amplified that rounding noise into a *floor* of
//! `2*sqrt(f64::EPSILON) ~= 3e-8` regardless of the true (much smaller) error -- a spurious
//! "budget exceeded" false positive that had nothing to do with the actual synthesis. Working
//! entirely in `FBig` removes that floor (down to whatever precision the crate itself used to
//! synthesize the candidate), and every angular quantity below is derived algebraically from
//! ring/complex arithmetic rather than by taking `atan2`/`acos` of an early-rounded `f64` and
//! re-differencing it (this crate does not implement arbitrary-precision inverse trig at all,
//! so the derivations below route around ever needing one -- see `phase_sq_from_z` and
//! `residual_cos_sin`).
//!
//! ## Two different, non-interchangeable ways to measure "how far" a channel is from target
//!
//! **Single Z-rotation vs. target Z-rotation** (`single_rotation_diamond_distance`): exact,
//! closed form, `2*|Im(w)|` (equal to `2*sqrt(1-Re(w)^2)` since `w` is unit-modulus, but
//! numerically stable -- see above). Always valid for comparing one pure rotation to another.
//!
//! **Pauli transfer matrix** (`ptm`/`pauli_diamond_distance_from_target`): computes
//! `Sum_P |delta_q_P|` from a channel's Pauli-basis decomposition. This closed form is only
//! mathematically valid when the *difference* between the channel and the target is itself
//! diagonal in the Pauli basis -- true for mixed diagonal's output (the {Z,S} twirl exists
//! specifically to enforce this by symmetrizing the off-diagonal X/Y error), but NOT true in
//! general for an untwisted single rotation mixed with something else: a lone Z-rotation by a
//! small angle delta genuinely mixes X and Y in its Pauli transfer matrix (that mixing IS what
//! the twirl is for), so applying the Pauli-diagonal shortcut to a channel that was never
//! twirled computes a meaningless number. An earlier version of this file made exactly that
//! mistake for fallback/mixed fallback and reported achieved errors around 1.0 (nonsense) for
//! everything, subsequently reduced to a state that *looked* plausible but was still wrong,
//! because the shortcut simply does not apply to the untwisted projective branch.
//!
//! ## The correct decomposition for fallback / mixed fallback
//!
//! Per the paper's own circuit description, the "success" outcome of the ancilla+measurement
//! circuit (out of scope for this crate) is an idealized rotation by `Arg(z)` -- not the raw
//! synthesized unitary's full matrix, whose off-diagonal component the projection removes.
//! The "failure" outcome is an idealized rotation by `Arg(v)`, followed by the ordinary
//! (non-projected) correction gate `B`. Since `B` is applied as an ordinary unconditional
//! unitary after a *fixed* rotation, and diamond norm is invariant under composing both sides
//! of a difference with the same fixed unitary channel, `B`'s own accuracy is exactly its
//! distance to the *residual* target angle `theta - Arg(v)` -- which is a fair question to ask
//! via the Pauli-diagonal shortcut, because `B` (a mixed-diagonal result) genuinely IS
//! twirled internally.
//!
//! So the right check is a **triangle-inequality upper bound**, combining an exact formula for
//! the untwisted success branch with a Pauli-diagonal-valid check for the twirled correction:
//!   `‖q*Z_{Arg z} + (1-q)*(B o Z_{Arg v}) - Z_theta‖⋄`
//!   `<= q * ‖Z_{Arg z} - Z_theta‖⋄  +  (1-q) * ‖B - Z_{theta-Arg(v)}‖⋄`
//! Both right-hand terms are independently exact/valid; the sum is a genuine, conservative
//! upper bound on the true (harder to compute exactly without twirling) diamond distance. If
//! this bound is <= the requested budget, that is a rigorous pass; if the bound exceeds the
//! budget, it does not by itself prove a failure (the true distance could still be smaller),
//! but it does mean this check cannot certify a pass and warrants a tighter look.

use dashu_base::{Abs, Approximation};
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;
use num::Complex;
use rsgridsynth::clear_caches;
use rsgridsynth::common::{cos_fbig, fb_with_prec, ib_to_bf_prec, sin_fbig};
use rsgridsynth::math::{sign, sqrt_fbig};
use rsgridsynth::protocol::fallback::exact_q;
use rsgridsynth::protocol::{
    synth_fallback, synth_mixed_diagonal, synth_mixed_fallback, FallbackResult,
    MixedDiagonalResult, MixedFallbackResult, MixedFallbackSide,
};
use rsgridsynth::unitary::DOmegaUnitary;

type Fb = FBig<HalfEven>;

/// Converts an `f64` to `Fb` the same way `config_from_theta_epsilon` converts its `theta`/
/// `epsilon` parameters: by parsing the **decimal string** representation (`x.to_string()`)
/// as an exact fraction, via `rsgridsynth::config::parse_decimal_with_exponent`. This must
/// match exactly, not just approximately -- `FBig::try_from(f64)` instead would give the
/// f64's exact *binary* value (e.g. `4.2_f64` is really
/// `4.20000000000000017763568...`), which differs from the crate's own "4.2" target by
/// ~1.8e-16. That is far smaller than any epsilon this file tested before, but at
/// epsilon=1e-15 it is comparable to the budget itself, and shows up as a spurious ~35x
/// "budget exceeded" on `mixed_fallback` at theta=4.2 that had nothing to do with the
/// synthesized channel -- purely an artifact of this file targeting a very slightly
/// different angle than the crate actually solved for.
fn to_fbig(x: f64) -> Fb {
    let (num, den) = rsgridsynth::config::parse_decimal_with_exponent(&x.to_string()).unwrap();
    fdiv(&ib_to_bf_prec(num), &ib_to_bf_prec(den))
}

fn to_f64(x: &Fb) -> f64 {
    match x.to_f64() {
        Approximation::Exact(v) => v,
        Approximation::Inexact(v, _) => v,
    }
}

fn fzero() -> Fb {
    ib_to_bf_prec(IBig::ZERO)
}

fn fone() -> Fb {
    ib_to_bf_prec(IBig::ONE)
}

fn fadd(a: &Fb, b: &Fb) -> Fb {
    fb_with_prec(a + b)
}

fn fsub(a: &Fb, b: &Fb) -> Fb {
    fb_with_prec(a - b)
}

fn fmul(a: &Fb, b: &Fb) -> Fb {
    fb_with_prec(a * b)
}

fn fdiv(a: &Fb, b: &Fb) -> Fb {
    fb_with_prec(a / b)
}

fn fneg(a: &Fb) -> Fb {
    -fb_with_prec(a.clone())
}

/// Arbitrary-precision complex number. `+`/`-`/`*`/[`Complex::conj`]/[`Complex::scale`]/
/// [`Complex::norm_sqr`] all come from `num::Complex<Fb>`'s own generic impls (it only requires
/// `Fb: Clone + Num`, which `FBig` satisfies -- see issue #16); only the precision-aware zero/
/// real/magnitude helpers below are crate-specific, since `Complex::zero()`'s `Fb::zero()` would
/// build an unbounded-precision `FBig::ZERO` rather than one pinned to the working `PREC_BITS`
/// context (see `fzero`/`ib_to_bf_prec`).
type Cx = Complex<Fb>;

fn cx_zero() -> Cx {
    Cx::new(fzero(), fzero())
}

fn cx_real(x: Fb) -> Cx {
    Cx::new(x, fzero())
}

fn cx_abs(z: &Cx) -> Fb {
    sqrt_fbig(&z.norm_sqr())
}

type M2 = [[Cx; 2]; 2];

fn zero_m2() -> M2 {
    std::array::from_fn(|_| std::array::from_fn(|_| cx_zero()))
}

fn zero_mat4() -> [[Fb; 4]; 4] {
    std::array::from_fn(|_| std::array::from_fn(|_| fzero()))
}

fn matrix_from_gates(gates: &[rsgridsynth::gate::Gate]) -> M2 {
    let m = DOmegaUnitary::from_gates(gates).to_complex_matrix();
    [
        [
            Cx::new(m[(0, 0)].re.clone(), m[(0, 0)].im.clone()),
            Cx::new(m[(0, 1)].re.clone(), m[(0, 1)].im.clone()),
        ],
        [
            Cx::new(m[(1, 0)].re.clone(), m[(1, 0)].im.clone()),
            Cx::new(m[(1, 1)].re.clone(), m[(1, 1)].im.clone()),
        ],
    ]
}

fn mat_mul(a: &M2, b: &M2) -> M2 {
    let mut r = zero_m2();
    for i in 0..2 {
        for j in 0..2 {
            let mut s = cx_zero();
            for k in 0..2 {
                s = &s + &(&a[i][k] * &b[k][j]);
            }
            r[i][j] = s;
        }
    }
    r
}

fn mat_dagger(a: &M2) -> M2 {
    [
        [a[0][0].conj(), a[1][0].conj()],
        [a[0][1].conj(), a[1][1].conj()],
    ]
}

fn pauli_basis() -> [M2; 4] {
    let zero = cx_zero();
    let one = cx_real(fone());
    let neg_one = cx_real(fneg(&fone()));
    let i_val = Cx::new(fzero(), fone());
    let neg_i = Cx::new(fzero(), fneg(&fone()));
    let id = [[one.clone(), zero.clone()], [zero.clone(), one.clone()]];
    let x = [[zero.clone(), one.clone()], [one.clone(), zero.clone()]];
    let y = [[zero.clone(), neg_i], [i_val, zero.clone()]];
    let z = [[one, zero.clone()], [zero, neg_one]];
    [id, x, y, z]
}

/// Pauli transfer matrix `R[P][Q] = (1/2) Re(Tr[sigma_P * Lambda(sigma_Q)])` for
/// `Lambda(rho) = sum_k weight_k * U_k * rho * U_k^dagger`.
fn ptm(branches: &[(Fb, M2)]) -> [[Fb; 4]; 4] {
    let sigma = pauli_basis();
    let mut r = zero_mat4();
    let half = fdiv(&fone(), &ib_to_bf_prec(IBig::from(2)));
    for (q_idx, sq) in sigma.iter().enumerate() {
        let mut lambda_sq = zero_m2();
        for (weight, u) in branches {
            let u_dagger = mat_dagger(u);
            let term = mat_mul(&mat_mul(u, sq), &u_dagger);
            for i in 0..2 {
                for j in 0..2 {
                    lambda_sq[i][j] = &lambda_sq[i][j] + &term[i][j].scale(weight.clone());
                }
            }
        }
        for (p_idx, sp) in sigma.iter().enumerate() {
            let prod = mat_mul(sp, &lambda_sq);
            let trace = &prod[0][0] + &prod[1][1];
            r[p_idx][q_idx] = fmul(&half, &trace.re);
        }
    }
    r
}

/// Target channel `Z_theta(rho) = R_z(theta) rho R_z(theta)^dagger`,
/// `R_z(theta) = diag(e^{-i*theta/2}, e^{i*theta/2})` (this crate's convention), built directly
/// from the half-angle's cosine/sine so callers with only a residual half-angle (no plain
/// `theta` to hand to `cos_fbig`/`sin_fbig`, e.g. `theta - Arg(v)`) can reuse it -- see
/// `residual_cos_sin`/`half_angle_cos_sin`.
fn target_matrix_from_half(cos_half: &Fb, sin_half: &Fb) -> M2 {
    let zero = cx_zero();
    [
        [Cx::new(cos_half.clone(), fneg(sin_half)), zero.clone()],
        [zero, Cx::new(cos_half.clone(), sin_half.clone())],
    ]
}

fn target_matrix(theta: &Fb) -> M2 {
    let two = ib_to_bf_prec(IBig::from(2));
    let half = fdiv(theta, &two);
    target_matrix_from_half(&cos_fbig(&half), &sin_fbig(&half))
}

/// `delta_q_P` from the diagonal of a *difference* of two channels' PTMs (each channel's own
/// `q_I+q_X+q_Y+q_Z=1` normalization cancels exactly in a genuine difference, so only the
/// linear part applies -- verified in `debug_tests` below by checking a channel against
/// itself gives exactly 0, not the `1` a naive "reapply the affine +1 formula" mistake gives).
fn delta_q_from_ptm_diagonal_difference(diff_diag: &[Fb; 4]) -> [Fb; 4] {
    let four = ib_to_bf_prec(IBig::from(4));
    let xx = &diff_diag[1];
    let yy = &diff_diag[2];
    let zz = &diff_diag[3];
    let neg_xx = fneg(xx);
    let neg_yy = fneg(yy);

    let s_xyz = fadd(&fadd(xx, yy), zz);
    let s_x_yz = fsub(&fsub(xx, yy), zz);
    let s_mx_y_z = fsub(&fadd(&neg_xx, yy), zz);
    let s_mx_my_z = fadd(&fadd(&neg_xx, &neg_yy), zz);

    [
        fdiv(&s_xyz, &four),
        fdiv(&s_x_yz, &four),
        fdiv(&s_mx_y_z, &four),
        fdiv(&s_mx_my_z, &four),
    ]
}

/// Diamond distance of a channel to a `target` channel matrix, valid ONLY when the channel is
/// genuinely (or was constructed to be, e.g. via a {Z,S} twirl) Pauli-diagonal relative to
/// that target -- see the module docs. Also returns the max absolute off-diagonal entry of
/// the difference PTM as a diagnostic: large values here mean the Pauli-diagonal precondition
/// does not hold and the returned "distance" should not be trusted.
fn pauli_diamond_distance_from_target(branches: &[(Fb, M2)], target: M2) -> (Fb, Fb) {
    let actual = ptm(branches);
    let target_ptm = ptm(&[(fone(), target)]);

    let mut diff = zero_mat4();
    for i in 0..4 {
        for j in 0..4 {
            diff[i][j] = fsub(&actual[i][j], &target_ptm[i][j]);
        }
    }

    let mut max_off_diag = fzero();
    for (i, row) in diff.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            if i != j {
                let av = value.clone().abs();
                if av > max_off_diag {
                    max_off_diag = av;
                }
            }
        }
    }

    let diag = [
        diff[0][0].clone(),
        diff[1][1].clone(),
        diff[2][2].clone(),
        diff[3][3].clone(),
    ];
    let dq = delta_q_from_ptm_diagonal_difference(&diag);
    let mut diamond_distance = fzero();
    for v in &dq {
        diamond_distance = fadd(&diamond_distance, &v.clone().abs());
    }
    (diamond_distance, max_off_diag)
}

/// Convenience wrapper of [`pauli_diamond_distance_from_target`] for callers that have a plain
/// target angle (rather than an already-derived residual half-angle).
fn pauli_diamond_distance_from_branches(branches: &[(Fb, M2)], theta: &Fb) -> (Fb, Fb) {
    pauli_diamond_distance_from_target(branches, target_matrix(theta))
}

/// `(cos(phi/2), sin(phi/2))` from `(cos(phi), sin(phi))`, via the half-angle formulas,
/// avoiding `atan2`/any inverse-trig call -- this crate does not implement arbitrary-precision
/// inverse trig at all, so this (deliberately duplicated from
/// `rsgridsynth::protocol::fallback`'s private `half_angle_cos_sin`, which this example -- an
/// independent check -- may not import) is the only way to get a residual angle's own
/// cosine/sine without ever computing the angle itself as a number.
fn half_angle_cos_sin(cos_phi: &Fb, sin_phi: &Fb) -> (Fb, Fb) {
    let zero = fzero();
    if *sin_phi == zero && *cos_phi < zero {
        return (zero, fone());
    }
    let two = ib_to_bf_prec(IBig::from(2));
    let one_plus_cos = fadd(&fone(), cos_phi);
    let one_minus_cos = fsub(&fone(), cos_phi);
    let cos_half = sqrt_fbig(&fdiv(&one_plus_cos, &two));
    let sin_half_mag = sqrt_fbig(&fdiv(&one_minus_cos, &two));
    let sin_half = if sign(sin_phi.clone()) < 0 {
        fneg(&sin_half_mag)
    } else {
        sin_half_mag
    };
    (cos_half, sin_half)
}

/// `(cos(theta - Arg(v)), sin(theta - Arg(v)))`, via the angle-subtraction identities and
/// `(cos(Arg v), sin(Arg v)) = (Re(v), Im(v)) / |v|` -- no `atan2` needed, `Arg(v)` itself is
/// never materialized as a number.
fn residual_cos_sin(theta: &Fb, v: &Cx) -> (Fb, Fb) {
    let norm = sqrt_fbig(&v.norm_sqr());
    let inv_norm = fdiv(&fone(), &norm);
    let cos_argv = fmul(&v.re, &inv_norm);
    let sin_argv = fmul(&v.im, &inv_norm);
    let cos_theta = cos_fbig(theta);
    let sin_theta = sin_fbig(theta);
    let cos_res = fadd(&fmul(&cos_theta, &cos_argv), &fmul(&sin_theta, &sin_argv));
    let sin_res = fsub(&fmul(&sin_theta, &cos_argv), &fmul(&cos_theta, &sin_argv));
    (cos_res, sin_res)
}

/// Exact diamond distance between a single IDEALIZED rotation by `Arg(z)` (i.e. the unit-
/// modulus phase `z/|z|`, NOT `z` itself -- for fallback/mixed fallback `|z|^2 = q < 1`, and
/// the idealized success action is a pure phase rotation, with the magnitude deficit `1-q`
/// entirely accounted for by the separate failure-branch weight, not by this term) and the
/// target rotation by `theta`: `2*|Im(w)|` with `w = (z/|z|) * e^{i*theta/2}` -- equal to
/// `2*sqrt(1-Re(w)^2)` since `w` is unit-modulus, but numerically stable at tight epsilon (see
/// module docs). Always valid (no Pauli-diagonal precondition -- this is a direct
/// unitary-vs-unitary comparison, not a mixture). Feeding this the raw (non-unit-modulus) `z`
/// instead would conflate the magnitude deficit with angular error and give a nonsensical,
/// epsilon-independent ~`2*sqrt(1-q)` answer -- exactly the bug an earlier version of this file
/// had, caught by comparing against `plain_diagonal_epsilon_convention_calibration`'s
/// epsilon-scaling expectation.
fn single_rotation_diamond_distance(z: &Cx, theta: &Fb) -> Fb {
    let norm = sqrt_fbig(&z.norm_sqr());
    let phase = z.scale(fdiv(&fone(), &norm));
    let two = ib_to_bf_prec(IBig::from(2));
    let half = fdiv(theta, &two);
    let c = cos_fbig(&half);
    let s = sin_fbig(&half);
    // Im(w) = z_x*Im(u) - z_y*Re(u), z_x = cos(half), z_y = -sin(half) (WFrame's convention).
    let im_w = fadd(&fmul(&phase.im, &c), &fmul(&phase.re, &s));
    fmul(&ib_to_bf_prec(IBig::from(2)), &im_w.abs())
}

/// `(z/|z|)^2 = e^{2i*Arg(z)}`, computed as `z^2 / |z|^2` (a single division by the *squared*
/// norm, no square root) rather than normalizing `z` first and squaring -- one fewer
/// transcendental step, and avoids ever computing `Arg(z)` as a number (see module docs).
fn phase_sq_from_z(z: &Cx) -> Cx {
    let inv_norm_sq = fdiv(&fone(), &z.norm_sqr());
    (z * z).scale(inv_norm_sq)
}

/// Exact diamond distance of a weighted combination of same-axis (Z-)rotations minus a
/// target Z-rotation, `‖Sum_k w_k*Z_{psi_k} - Z_theta‖⋄`. Derivation: any Z-rotation acts on
/// a density matrix's diagonal entries as the identity and on its off-diagonal entry `rho_01`
/// by multiplying by a single complex phase `e^{-2i*psi}` -- so the *difference* channel acts
/// as zero on the diagonal (regardless of weights, since each term individually fixes it) and
/// as multiplication by a single complex number
/// `c = Sum_k w_k*(e^{-2i*psi_k} - e^{-i*theta})` on the off-diagonal. For a qubit map of
/// exactly this form (diagonal-preserving, off-diagonal scaled by one fixed complex number),
/// the diamond norm equals `|c|` exactly, achieved by the maximally-coherent probe state --
/// verified against the ALREADY-independently-checked single-rotation formula
/// `2*|sin(delta)|` as the single-branch special case in `debug_tests` below (both formulas
/// must agree there, and they do, in closed form: `|e^{-2i*delta}-1| = 2*|sin(delta)|`
/// exactly).
///
/// Takes each branch's `e^{-2i*psi_k}` directly (as `phase_sq_from_z` produces, see that
/// function's docs) rather than `psi_k` itself, so this never needs `Arg`/`atan2` either.
///
/// This is what lets mixed fallback's projective-mixing cancellation (the whole point of the
/// protocol) be verified independently of the crate's own `mixture_weight` closed form: a
/// naive triangle-inequality bound on the two branches SEPARATELY (an earlier version of this
/// check) is too loose to see the cancellation at all, and wrongly looks like a failure.
fn rotation_mixture_diamond_distance(weighted_phase_sq: &[(Fb, Cx)], theta: &Fb) -> Fb {
    let mut c = cx_zero();
    let target_phase = Cx::new(cos_fbig(theta), fneg(&sin_fbig(theta))); // e^{-i*theta}
    for (w, phase_sq) in weighted_phase_sq {
        let diff = phase_sq - &target_phase;
        c = &c + &diff.scale(w.clone());
    }
    cx_abs(&c)
}

fn mixed_diagonal_pauli_branches(result: &MixedDiagonalResult) -> Vec<(Fb, M2)> {
    result
        .branches
        .iter()
        .map(|b| (b.weight.clone(), matrix_from_gates(&b.gates)))
        .collect()
}

/// Triangle-inequality upper bound on the total diamond distance of a plain fallback result:
/// `q * single_rotation_distance(z, theta) + (1-q) * pauli_distance(B, theta - Arg(v))`.
/// Returns `(bound, correction_off_diagonal_diagnostic)`.
fn fallback_upper_bound(result: &FallbackResult, z: &Cx, v: &Cx, theta: &Fb) -> (Fb, Fb) {
    let q = result.achieved_success_probability();
    let success_term = single_rotation_diamond_distance(z, theta);

    let (cos_res, sin_res) = residual_cos_sin(theta, v);
    let (cos_half, sin_half) = half_angle_cos_sin(&cos_res, &sin_res);
    let target = target_matrix_from_half(&cos_half, &sin_half);

    let b_matrix = matrix_from_gates(&result.correction_gates);
    let (correction_term, off_diag) =
        pauli_diamond_distance_from_target(&[(fone(), b_matrix)], target);

    let one_minus_q = fsub(&fone(), &q);
    let bound = fadd(
        &fmul(&q, &success_term),
        &fmul(&one_minus_q, &correction_term),
    );
    (bound, off_diag)
}

/// Total diamond-distance bound for a mixed fallback result, correctly exploiting the
/// projective-mixing cancellation (see `rotation_mixture_diamond_distance`'s docs) instead of
/// a naive triangle-inequality split of the two sides' success terms, which is too loose to
/// see the cancellation at all and looks like a failure even for a correct implementation.
///
/// Groups the four underlying pieces (lo success, lo correction, hi success, hi correction)
/// as: `[p*q_lo*Z_{psi_lo} + (1-p)*q_hi*Z_{psi_hi} - Z_theta] + [p*(1-q_lo)*(B_lo-related)] +
/// [(1-p)*(1-q_hi)*(B_hi-related)]`. The first bracket has EXACTLY the same
/// "diagonal-preserving, single-complex-off-diagonal-scaling" structure
/// `rotation_mixture_diamond_distance` handles (the weights `p*q_lo`/`(1-p)*q_hi` need not
/// themselves sum to 1 -- the derivation there only used that each individual rotation fixes
/// the diagonal, which holds regardless of the weight). The other two brackets are each a
/// twirled mixed-diagonal correction (Pauli-diagonal-valid) composed with a fixed rotation
/// (removable by diamond-norm's unitary invariance), exactly as in `fallback_upper_bound`.
/// Combining the three brackets via the triangle inequality is a valid (if not perfectly
/// tight, since the theorem's own formula could in principle also correlate the correction
/// terms with the success terms) upper bound.
struct MixedFallbackSideInput<'a> {
    side: &'a MixedFallbackSide,
    z: Cx,
    v: Cx,
}

fn mixed_fallback_total_bound(
    lo: &MixedFallbackSideInput,
    hi: &MixedFallbackSideInput,
    p: &Fb,
    theta: &Fb,
) -> (Fb, Fb) {
    let q_lo = lo.side.achieved_success_probability();
    let q_hi = hi.side.achieved_success_probability();
    let one_minus_p = fsub(&fone(), p);

    let phase_sq_lo = phase_sq_from_z(&lo.z);
    let phase_sq_hi = phase_sq_from_z(&hi.z);
    let w_lo = fmul(p, &q_lo);
    let w_hi = fmul(&one_minus_p, &q_hi);
    let projective_term =
        rotation_mixture_diamond_distance(&[(w_lo, phase_sq_lo), (w_hi, phase_sq_hi)], theta);

    let (cos_res_lo, sin_res_lo) = residual_cos_sin(theta, &lo.v);
    let (cos_half_lo, sin_half_lo) = half_angle_cos_sin(&cos_res_lo, &sin_res_lo);
    let target_lo = target_matrix_from_half(&cos_half_lo, &sin_half_lo);
    let (correction_lo, off_diag_lo) = pauli_diamond_distance_from_target(
        &mixed_diagonal_pauli_branches(&lo.side.correction),
        target_lo,
    );

    let (cos_res_hi, sin_res_hi) = residual_cos_sin(theta, &hi.v);
    let (cos_half_hi, sin_half_hi) = half_angle_cos_sin(&cos_res_hi, &sin_res_hi);
    let target_hi = target_matrix_from_half(&cos_half_hi, &sin_half_hi);
    let (correction_hi, off_diag_hi) = pauli_diamond_distance_from_target(
        &mixed_diagonal_pauli_branches(&hi.side.correction),
        target_hi,
    );

    let one_minus_q_lo = fsub(&fone(), &q_lo);
    let one_minus_q_hi = fsub(&fone(), &q_hi);
    let term_lo = fmul(&fmul(p, &one_minus_q_lo), &correction_lo);
    let term_hi = fmul(&fmul(&one_minus_p, &one_minus_q_hi), &correction_hi);
    let bound = fadd(&fadd(&projective_term, &term_lo), &term_hi);

    let off_diag = if off_diag_lo > off_diag_hi {
        off_diag_lo
    } else {
        off_diag_hi
    };
    (bound, off_diag)
}

fn main() {
    let thetas = [
        0.1,
        3.0 * std::f64::consts::PI / 32.0,
        5.0 * std::f64::consts::PI / 32.0,
        7.0 * std::f64::consts::PI / 32.0,
        std::f64::consts::PI / 3.0,
        std::f64::consts::PI / 6.0,
        1.9,
        4.2,
    ];
    let epsilons = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-15];
    let q = exact_q(7);

    println!("protocol,theta,epsilon,epsilon_diamond_requested,achieved_diamond_distance_upper_bound,off_diagonal_diagnostic,meets_budget");

    for &eps in &epsilons {
        for (i, &theta_f64) in thetas.iter().enumerate() {
            let seed = 2000 + i as u64;

            clear_caches();
            let md = synth_mixed_diagonal(theta_f64, eps, seed, false);
            let theta = to_fbig(theta_f64);
            let eps_fb = to_fbig(eps);
            let branches = mixed_diagonal_pauli_branches(&md);
            let (dd, max_od) = pauli_diamond_distance_from_branches(&branches, &theta);
            println!(
                "mixed_diagonal,{theta_f64},{eps},{eps},{},{},{}",
                to_f64(&dd),
                to_f64(&max_od),
                dd <= eps_fb
            );

            clear_caches();
            let sin_alpha = eps / 4.0;
            match synth_fallback(theta_f64, eps, q.clone(), sin_alpha, seed, false) {
                Some(result) => {
                    let theta = to_fbig(theta_f64);
                    let eps_fb = to_fbig(eps);
                    let v_mat = matrix_from_gates(&result.projective_gates);
                    let z = v_mat[0][0].clone();
                    let v = v_mat[1][0].clone();
                    let (bound, off_diag) = fallback_upper_bound(&result, &z, &v, &theta);
                    println!(
                        "fallback,{theta_f64},{eps},{eps},{},{},{}",
                        to_f64(&bound),
                        to_f64(&off_diag),
                        bound <= eps_fb
                    );
                }
                None => println!("fallback,{theta_f64},{eps},{eps},ERROR,NotFound,false"),
            }

            clear_caches();
            match synth_mixed_fallback(theta_f64, eps, q.clone(), seed, false) {
                Some(MixedFallbackResult::Exact { gates }) => {
                    let theta = to_fbig(theta_f64);
                    let eps_fb = to_fbig(eps);
                    let m = matrix_from_gates(&gates);
                    let (dd, off_diag) =
                        pauli_diamond_distance_from_branches(&[(fone(), m)], &theta);
                    println!(
                        "mixed_fallback,{theta_f64},{eps},{eps},{},{},{}",
                        to_f64(&dd),
                        to_f64(&off_diag),
                        dd <= eps_fb
                    );
                }
                Some(MixedFallbackResult::Mixed { lo, hi, p, .. }) => {
                    let theta = to_fbig(theta_f64);
                    let eps_fb = to_fbig(eps);
                    let lo_v = matrix_from_gates(&lo.projective_gates);
                    let hi_v = matrix_from_gates(&hi.projective_gates);
                    let lo_input = MixedFallbackSideInput {
                        side: &lo,
                        z: lo_v[0][0].clone(),
                        v: lo_v[1][0].clone(),
                    };
                    let hi_input = MixedFallbackSideInput {
                        side: &hi,
                        z: hi_v[0][0].clone(),
                        v: hi_v[1][0].clone(),
                    };
                    let (bound, off_diag) =
                        mixed_fallback_total_bound(&lo_input, &hi_input, &p, &theta);
                    println!(
                        "mixed_fallback,{theta_f64},{eps},{eps},{},{},{}",
                        to_f64(&bound),
                        to_f64(&off_diag),
                        bound <= eps_fb
                    );
                }
                None => println!("mixed_fallback,{theta_f64},{eps},{eps},ERROR,NotFound,false"),
            }
        }
    }
}

#[cfg(test)]
mod debug_tests {
    use super::*;
    use rsgridsynth::common::reset_prec_bits;
    use serial_test::serial;

    // `PREC_BITS` (see `rsgridsynth::common`) is a single process-global atomic that every
    // `cos_fbig`/`sin_fbig`/`fb_with_prec` call in this file reads. `cargo test` runs tests in
    // this module concurrently by default, and `plain_diagonal_epsilon_convention_calibration`
    // below mutates `PREC_BITS` (via `config_from_theta_epsilon`) across a range of epsilons --
    // so without `#[serial]` + an explicit `reset_prec_bits()`, a test's own trig calls can
    // observe a precision far lower than intended mid-computation, corrupting a cancellation-
    // heavy result (this is exactly what flagged an apparent bug in
    // `rotation_mixture_matches_quadratic_cancellation_residual` during development: the
    // formula was fine, but a concurrently-running test had transiently dropped `PREC_BITS` to
    // ~36 bits). Every test below follows this crate's own established convention (see
    // `tests/integration_test.rs`, `src/protocol/mixing.rs`) of `#[serial]` plus
    // `reset_prec_bits()` for exactly this reason.

    #[test]
    #[serial]
    fn self_consistency_target_vs_itself() {
        reset_prec_bits();
        let theta = to_fbig(0.7_f64);
        let m = target_matrix(&theta);
        let (dd, max_od) = pauli_diamond_distance_from_branches(&[(fone(), m)], &theta);
        eprintln!("dd={}, max_od={}", to_f64(&dd), to_f64(&max_od));
        assert!(
            to_f64(&dd) < 1e-9,
            "self-comparison should give ~0, got {dd}"
        );
    }

    /// Locks in the derivation `delta_q_from_ptm_diagonal_difference` is built from: a pure
    /// identity channel's PTM diagonal is `[1,1,1,1]` (I,X,Y,Z all fixed) and decomposes to
    /// `q_I=1`, others 0; a pure-Z channel's PTM diagonal is `[1,-1,-1,1]` (Z fixes I and Z,
    /// negates X and Y) and decomposes to `q_Z=1`, others 0.
    #[test]
    #[serial]
    fn q_from_ptm_diagonal_matches_known_pauli_channels() {
        fn q_from_ptm_diagonal(diag: [f64; 4]) -> [f64; 4] {
            let (xx, yy, zz) = (diag[1], diag[2], diag[3]);
            [
                (1.0 + xx + yy + zz) / 4.0,
                (1.0 + xx - yy - zz) / 4.0,
                (1.0 - xx + yy - zz) / 4.0,
                (1.0 - xx - yy + zz) / 4.0,
            ]
        }
        let identity = q_from_ptm_diagonal([1.0, 1.0, 1.0, 1.0]);
        assert!((identity[0] - 1.0).abs() < 1e-12);
        assert!(
            identity[1].abs() < 1e-12 && identity[2].abs() < 1e-12 && identity[3].abs() < 1e-12
        );

        let pure_z = q_from_ptm_diagonal([1.0, -1.0, -1.0, 1.0]);
        assert!((pure_z[3] - 1.0).abs() < 1e-12);
        assert!(pure_z[0].abs() < 1e-12 && pure_z[1].abs() < 1e-12 && pure_z[2].abs() < 1e-12);
    }

    #[test]
    #[serial]
    fn identity_channel_vs_target_theta_zero() {
        reset_prec_bits();
        let theta = to_fbig(0.0_f64);
        let id = [[cx_real(fone()), cx_zero()], [cx_zero(), cx_real(fone())]];
        let (dd, max_od) = pauli_diamond_distance_from_branches(&[(fone(), id)], &theta);
        eprintln!("dd={}, max_od={}", to_f64(&dd), to_f64(&max_od));
        assert!(
            to_f64(&dd) < 1e-9,
            "identity vs target(theta=0) should give ~0, got {dd}"
        );
    }

    /// `single_rotation_diamond_distance` cross-checked against the Pauli-diagonal formula
    /// for the special case where they SHOULD agree.
    /// `rotation_mixture_diamond_distance`'s single-weight special case must agree exactly
    /// with `single_rotation_diamond_distance` (both independently derived, one via the
    /// Pauli-transfer/witness argument, one via the direct `2*|Im(w)|` formula).
    #[test]
    #[serial]
    fn rotation_mixture_matches_single_rotation_special_case() {
        reset_prec_bits();
        let theta = to_fbig(0.5_f64);
        let z_theta = to_fbig(0.52_f64); // theta + a small angular error, unit modulus target
        let z = target_matrix(&z_theta)[0][0].clone();
        let phase_sq = phase_sq_from_z(&z);
        let via_mixture = rotation_mixture_diamond_distance(&[(fone(), phase_sq)], &theta);
        let via_single = single_rotation_diamond_distance(&z, &theta);
        eprintln!(
            "via_mixture={}, via_single={}",
            to_f64(&via_mixture),
            to_f64(&via_single)
        );
        assert!((to_f64(&via_mixture) - to_f64(&via_single)).abs() < 1e-9);
    }

    /// A symmetric straddling pair (equal and opposite angular error `delta`, equal weight
    /// `p=0.5`) exactly cancels the FIRST-order (linear-in-delta) term -- the whole point of
    /// the mixing theorem -- leaving only the theorem's own predicted SECOND-order residual
    /// `2*(p*sin^2(delta_lo) + (1-p)*sin^2(delta_hi))` (with q_lo=q_hi=1, unit-modulus
    /// idealized rotations), NOT exactly zero: mixing reduces the error from linear to
    /// quadratic order in delta, it does not eliminate it. An earlier version of this test
    /// wrongly asserted exactly-zero and "failed" on a numerically-correct result.
    #[test]
    #[serial]
    fn rotation_mixture_matches_quadratic_cancellation_residual() {
        reset_prec_bits();
        let theta = to_fbig(0.5_f64);
        let delta = 0.001_f64; // paper's delta_k = Arg(z_k) - theta/2
        let z_lo = target_matrix(&to_fbig(0.5 - 2.0 * delta))[0][0].clone(); // psi_lo = theta/2 - delta
        let z_hi = target_matrix(&to_fbig(0.5 + 2.0 * delta))[0][0].clone(); // psi_hi = theta/2 + delta
        let phase_sq_lo = phase_sq_from_z(&z_lo);
        let phase_sq_hi = phase_sq_from_z(&z_hi);
        let half = to_fbig(0.5);
        let dd = rotation_mixture_diamond_distance(
            &[(half.clone(), phase_sq_lo), (half, phase_sq_hi)],
            &theta,
        );
        let expected = 2.0 * (0.5 * delta.sin().powi(2) + 0.5 * delta.sin().powi(2));
        eprintln!(
            "symmetric-mixing dd={}, expected quadratic residual={expected}",
            to_f64(&dd)
        );
        assert!(
            (to_f64(&dd) - expected).abs() < 1e-12,
            "got {}, expected {expected}",
            to_f64(&dd)
        );
    }

    #[test]
    #[serial]
    fn single_rotation_distance_matches_known_small_case() {
        reset_prec_bits();
        // R_z(theta) vs R_z(theta + delta): known exact diamond distance 2*|sin(delta/2)|.
        let theta = to_fbig(0.5_f64);
        let delta = 0.02_f64;
        let z = target_matrix(&to_fbig(0.5 + delta))[0][0].clone();
        let got = single_rotation_diamond_distance(&z, &theta);
        let expected = 2.0 * (delta / 2.0).sin().abs();
        eprintln!("got={}, expected={expected}", to_f64(&got));
        assert!(
            (to_f64(&got) - expected).abs() < 1e-9,
            "got {}, expected {expected}",
            to_f64(&got)
        );
    }

    /// Calibration check: what does the EXISTING, unmodified plain-diagonal protocol's
    /// `epsilon` parameter actually correspond to in diamond-norm terms?
    #[test]
    #[serial]
    fn plain_diagonal_epsilon_convention_calibration() {
        use rsgridsynth::config::config_from_theta_epsilon;
        use rsgridsynth::gridsynth::gridsynth_gates;

        for &theta_f64 in &[0.1_f64, 0.7, 1.3, 2.5] {
            for &eps in &[1e-3_f64, 1e-4, 1e-5, 1e-6] {
                clear_caches();
                let mut config = config_from_theta_epsilon(theta_f64, eps, 42, false, false);
                let result = gridsynth_gates(&mut config);
                let theta = to_fbig(theta_f64);
                let m = matrix_from_gates(&result.gates);
                let (dd, _max_od) = pauli_diamond_distance_from_branches(&[(fone(), m)], &theta);
                let dd_f64 = to_f64(&dd);
                eprintln!(
                    "theta={theta_f64}, epsilon_param={eps}, achieved_diamond_distance={dd_f64}, \
                     ratio={:.4}",
                    dd_f64 / eps
                );
            }
        }
    }
}
