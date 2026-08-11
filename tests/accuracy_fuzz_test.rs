//! Fuzz-style accuracy tests.
//!
//! These generate many random target angles across a range of epsilons -- down to 1e-15 -- and
//! check that the synthesized gate string is actually within the requested tolerance of the
//! ideal rotation. Accuracy is computed on demand via
//! `GridSynthResult::achieved_diamond_error` (`AchievedDiamondError`), not cached eagerly during
//! synthesis, and cross-checked against a genuinely different derivation
//! (`independent_operator_error` below): it rebuilds the exact unitary represented by the
//! returned gate string (via `DOmegaUnitary::from_gates`) and the ideal target rotation (via
//! `cos_fbig`/`sin_fbig` at the same working precision) and computes the *operator*-norm
//! distance between them from the full matrix eigenvalue formula -- a different code path from
//! `achieved_diamond_error`'s `WFrame`-based shortcut, related by the well-known
//! `diamond = 2 * operator_norm` identity for this special (SU(2)-with-phase) matrix form. That
//! way a bug in either derivation would show up as a disagreement, not just as both being wrong
//! in the same way.

use dashu_base::{Approximation, SquareRoot};
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use num::Complex;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rsgridsynth::accuracy::AchievedDiamondError;
use rsgridsynth::clear_caches;
use rsgridsynth::common::{cos_fbig, fb_with_prec, get_prec_bits, sin_fbig};
use rsgridsynth::config::config_from_theta_epsilon;
use rsgridsynth::gate::Gate;
use rsgridsynth::gridsynth::gridsynth_gates;
use rsgridsynth::unitary::DOmegaUnitary;
use serial_test::serial;

fn to_fbig(x: f64) -> FBig<HalfEven> {
    FBig::<HalfEven>::try_from(x)
        .unwrap()
        .with_precision(get_prec_bits())
        .value()
}

fn fbig_to_f64(x: &FBig<HalfEven>) -> f64 {
    match x.to_f64() {
        Approximation::Exact(v) => v,
        Approximation::Inexact(v, _) => v,
    }
}

/// Recomputes the *operator*-norm distance between the ideal z-rotation by `theta` and the
/// unitary represented by `gates`, entirely from public API, via the full matrix eigenvalue
/// formula -- a different derivation from `achieved_diamond_error`'s `WFrame`-based shortcut,
/// not a copy of it. `shifted` selects whether the synthesized unitary should be compared up to
/// the extra global phase `e^{i pi/8}` (this crate's `PhaseMode`), matching
/// `GridSynthResult::global_phase`.
fn independent_operator_error(gates: &[Gate], theta: &FBig<HalfEven>, shifted: bool) -> f64 {
    let two = fb_with_prec(FBig::try_from(2.0).unwrap());
    let neg_theta_half = -fb_with_prec(theta / &two);
    let z_x = fb_with_prec(cos_fbig(&neg_theta_half));
    let z_y = fb_with_prec(sin_fbig(&neg_theta_half));

    let synthesized = DOmegaUnitary::from_gates(gates).to_complex_matrix();
    let mut u = synthesized[(0, 0)].clone();
    if shifted {
        let p = to_fbig(std::f64::consts::PI / 8.);
        let phase = Complex::new(cos_fbig(&p), sin_fbig(&p));
        u = &u * &phase;
    }

    // Squared operator norm of (expected - synthesized), via the shared eigenvalue formula:
    // ||A^* A|| for A = expected - synthesized, both being 2x2 unitaries with the same
    // (0,0)-entry phase relationship.
    let eig: FBig<HalfEven> = 2 - 2 * (&z_x * &u.re + &z_y * &u.im);
    let eig = eig.max(FBig::from(0));
    let norm = fb_with_prec(eig.sqrt());
    match norm.to_f64() {
        Approximation::Inexact(v, _) => v,
        Approximation::Exact(v) => v,
    }
}

/// Runs the fuzzer for a given `up_to_phase` setting across a spread of epsilons -- from coarse
/// (1e-2) down to 1e-15 -- and many random target angles per epsilon, checking that:
///  - the on-demand `achieved_diamond_error` is within the requested (diamond-norm, i.e.
///    `2*epsilon`) budget,
///  - an independently derived operator-norm error is *also* within budget (`epsilon`),
///  - the two error computations -- different derivations, related by `diamond = 2*operator` --
///    agree with each other.
fn run_accuracy_fuzz(up_to_phase: bool, thetas_per_epsilon: usize, seeds: &[u64]) {
    // From coarse tolerances down to 1e-15, spanning the precision regimes the algorithm has to
    // handle differently (see `config_from_theta_epsilon`'s `calculated_prec_bits`).
    let epsilons = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-15];

    let mut rng = StdRng::seed_from_u64(0xACC0_FA22);

    for &epsilon in &epsilons {
        for _ in 0..thetas_per_epsilon {
            let theta = rng.random_range(0.0..std::f64::consts::TAU);
            for &seed in seeds {
                clear_caches();
                let mut config =
                    config_from_theta_epsilon(theta, epsilon, seed, false, up_to_phase);
                let res = gridsynth_gates(&mut config);

                let diamond_error = fbig_to_f64(&res.achieved_diamond_error(&config.theta));
                assert!(
                    diamond_error <= 2.0 * epsilon,
                    "achieved diamond error {diamond_error:e} exceeds requested budget \
                     2*epsilon={:e} for theta={theta}, epsilon={epsilon:e}, seed={seed}, \
                     up_to_phase={up_to_phase}, gates={}",
                    2.0 * epsilon,
                    res.gates
                );

                let independent_error =
                    independent_operator_error(&res.gates, &config.theta, res.global_phase);
                assert!(
                    independent_error <= epsilon,
                    "independently computed operator error {independent_error:e} exceeds \
                     requested epsilon {epsilon:e} for theta={theta}, seed={seed}, \
                     up_to_phase={up_to_phase}, gates={}",
                    res.gates
                );

                // Two different derivations (WFrame-shortcut diamond error vs. full-matrix
                // eigenvalue operator norm), related by `diamond = 2*operator_norm` for this
                // special matrix form -- should agree up to the last couple of bits of f64
                // rounding.
                let independent_diamond_error = 2.0 * independent_error;
                let diff = (diamond_error - independent_diamond_error).abs();
                let scale = diamond_error.max(independent_diamond_error).max(1e-300);
                assert!(
                    diff <= scale * 1e-6,
                    "achieved diamond error {diamond_error:e} and 2x independently computed \
                     operator error {independent_diamond_error:e} disagree for theta={theta}, \
                     epsilon={epsilon:e}, seed={seed}, up_to_phase={up_to_phase}"
                );
            }
        }
    }
}

#[test]
#[serial]
fn fuzz_accuracy_exact_phase() {
    run_accuracy_fuzz(false, 6, &[0, 1234, 987654321]);
}

#[test]
#[serial]
fn fuzz_accuracy_up_to_phase() {
    run_accuracy_fuzz(true, 6, &[0, 1234, 987654321]);
}

/// Dedicated, larger sweep specifically at the 1e-15 tolerance boundary this crate's README
/// flags as not fully supported for the CLI's f64-based epsilon parsing -- exercising many more
/// random angles at exactly that precision to build confidence the library entry point
/// (`config_from_theta_epsilon`/`gridsynth_gates`) still produces accurate results there.
#[test]
#[serial]
fn fuzz_accuracy_at_1e_minus_15() {
    let epsilon = 1e-15;
    let mut rng = StdRng::seed_from_u64(0x1E_15FA22);

    for _ in 0..40 {
        let theta = rng.random_range(0.0..std::f64::consts::TAU);
        clear_caches();
        let mut config = config_from_theta_epsilon(theta, epsilon, 42, false, false);
        let res = gridsynth_gates(&mut config);

        let diamond_error = fbig_to_f64(&res.achieved_diamond_error(&config.theta));
        assert!(
            diamond_error <= 2.0 * epsilon,
            "achieved diamond error {diamond_error:e} exceeds requested budget 2*epsilon={:e} \
             for theta={theta}",
            2.0 * epsilon
        );

        let independent_error =
            independent_operator_error(&res.gates, &config.theta, res.global_phase);
        assert!(
            independent_error <= epsilon,
            "independently computed operator error {independent_error:e} exceeds requested \
             epsilon {epsilon:e} for theta={theta}, gates={}",
            res.gates
        );
    }
}
