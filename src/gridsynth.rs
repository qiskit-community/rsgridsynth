// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::common::{cos_fbig, fb_with_prec, get_prec_bits, ib_to_bf_prec, sin_fbig};
use crate::config::{GridSynthConfig, GridSynthResult};
use crate::diophantine::diophantine_dyadic;
use crate::grid_op::GridOp;
use crate::math::solve_quadratic;
use crate::math::sqrt_fbig;
use crate::region::{Ellipse, Rectangle};
use crate::ring::{DOmega, DRootTwo, ZOmega, ZRootTwo};
use crate::synthesis_of_clifford_t::decompose_domega_unitary;
use crate::tdgp::solve_tdgp;
use crate::tdgp::Region;
use crate::to_upright::to_upright_set_pair;
use crate::unitary::DOmegaUnitary;
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;

//use log::{debug, info};
use log::debug;

use nalgebra::{Matrix2, Vector2};
use std::cmp::Ordering;
use std::time::{Duration, Instant};

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

// PhaseMode::Exact synthesize gate including exact phase
// PhaseMode::Shifted synthesize gate with a fixed phase factor of `exp(i pi/8)`

// If we don't care about phase, then it is enough to check both `U` and `exp(i pi/8) U`.
//
// To synthesize up to a phase, we run both `PhaseMode::Exact` and
// `PhaseMode::Shifted` and keep the one with lower T count. We first compute the best
// exact solution and then the best with the phase factor. An optimization would be to  interleave candidates
// from each to avoid doing more work than necessary.
//
// The following comments assume we are checking `exp(i pi/8) U`.
// The pair 2 ± √2 enter in some places as scale factors.
//
// omega = exp(-i pi/4)
// delta = 1 + omega
// |delta|^2 = 2 + 2cos(pi/4) = 2 + √2
// From Lemma 9.6 in R + S, we must scale the epsilon region by
// |delta| = √(2 + √2)
//
// We scale the UnitDisk by the root-2 conjugate of |delta|:
// |delta^●| = √(2 - √2)
// See Algorithm 9.8, page 20 of R+S.
#[derive(Debug, Clone, Copy)]
pub enum PhaseMode {
    Exact,   // no scaling
    Shifted, // do scaling
}

#[derive(Debug)]
pub struct EpsilonRegion {
    _theta: FBig<HalfEven>,
    _epsilon: FBig<HalfEven>,
    scale: ZRootTwo,
    d: FBig<HalfEven>,
    z_x: FBig<HalfEven>,
    z_y: FBig<HalfEven>,
    ellipse: Ellipse,
}

impl EpsilonRegion {
    pub fn new(theta: FBig<HalfEven>, epsilon: FBig<HalfEven>, scale: ZRootTwo) -> Self {
        let two = fb_with_prec(FBig::try_from(2.0).unwrap());
        let theta_half = fb_with_prec(&theta / &two);
        let neg_theta_half = -fb_with_prec(theta_half);
        let z_x: FBig<HalfEven> = fb_with_prec(cos_fbig(&neg_theta_half));
        let z_y: FBig<HalfEven> = fb_with_prec(sin_fbig(&neg_theta_half));
        Self::from_target_direction_impl(z_x, z_y, epsilon, scale, theta)
    }

    /// Builds the same region as [`EpsilonRegion::new`], but from the target direction's
    /// half-angle cosine/sine directly (`z_x = cos(-phi/2)`, `z_y = sin(-phi/2)` for whatever
    /// angle `phi` the caller wants as the target), instead of from a raw angle `theta`.
    ///
    /// This exists so a caller that already has `phi` expressed as an exact `(cos(phi/2),
    /// sin(phi/2))` pair -- e.g. derived algebraically via angle-addition/half-angle identities
    /// from another region's synthesized unitary, as the fallback protocol's residual-angle
    /// correction step needs -- can build the region without a lossy angle round-trip through
    /// an `atan2`-style inverse-trig call (which this crate does not otherwise implement) and a
    /// second `cos_fbig`/`sin_fbig` evaluation. `(z_x, z_y)` must satisfy `z_x^2 + z_y^2 == 1`
    /// (a unit vector); this is not checked.
    pub fn from_target_direction(
        z_x: FBig<HalfEven>,
        z_y: FBig<HalfEven>,
        epsilon: FBig<HalfEven>,
        scale: ZRootTwo,
    ) -> Self {
        // No real `theta` is available here; the `_theta`/`_epsilon` fields are stored only for
        // `Debug` output and are never read by any `Region` method, so a placeholder is fine.
        let placeholder_theta: FBig<HalfEven> = ib_to_bf_prec(IBig::ZERO);
        Self::from_target_direction_impl(z_x, z_y, epsilon, scale, placeholder_theta)
    }

    fn from_target_direction_impl(
        z_x: FBig<HalfEven>,
        z_y: FBig<HalfEven>,
        epsilon: FBig<HalfEven>,
        scale: ZRootTwo,
        theta_for_debug: FBig<HalfEven>,
    ) -> Self {
        let one = fb_with_prec(FBig::try_from(1.0).unwrap());
        let four = fb_with_prec(FBig::try_from(4.0).unwrap());
        let epsilon_squared = fb_with_prec(&epsilon * &epsilon);
        let half_eps_sq = fb_with_prec(&epsilon_squared / &four);
        // `epsilon` >= 2 (or a derived epsilon that overshoots it, e.g. the fallback
        // protocol's rescaled correction-step epsilon on a near-degenerate candidate) is
        // already past the point where any point of the disk fails to qualify: the sane
        // mathematical limit of the formula below is `d = 0` (no angular restriction at
        // all), not a negative radicand. Clamp rather than let a tiny-precision rounding
        // artifact (or a legitimately oversized derived epsilon) panic in `sqrt_fbig`.
        let one_minus_half_eps_sq = (one - half_eps_sq).max(ib_to_bf_prec(IBig::ZERO));
        let scale_to_real = scale.to_real();
        let d = fb_with_prec(sqrt_fbig(&one_minus_half_eps_sq) * sqrt_fbig(&scale_to_real));

        let neg_z_y: FBig<HalfEven> = -fb_with_prec(z_y.clone());
        let zero: FBig<HalfEven> = ib_to_bf_prec(IBig::ZERO);
        let epsilon_neg4: FBig<HalfEven> = fb_with_prec(epsilon.clone().powi(IBig::from(-4)));
        let epsilon_neg2: FBig<HalfEven> = fb_with_prec(epsilon.clone().powi(IBig::from(-2)));
        let d1: Matrix2<FBig<HalfEven>> =
            Matrix2::new(z_x.clone(), neg_z_y.clone(), z_y.clone(), z_x.clone());
        let d2: Matrix2<FBig<HalfEven>> = Matrix2::new(
            fb_with_prec(64 * epsilon_neg4 / &scale_to_real),
            zero.clone(),
            zero.clone(),
            fb_with_prec(4 * epsilon_neg2 / &scale_to_real),
        );
        let d3: Matrix2<FBig<HalfEven>> =
            Matrix2::new(z_x.clone(), z_y.clone(), neg_z_y, z_x.clone());
        let px = fb_with_prec(&d * &z_x);
        let py = fb_with_prec(&d * &z_y);
        let p = Vector2::new(px, py);
        let m1: Matrix2<FBig<HalfEven>> = matrix_multiply_2x2(&d1, &d2);
        let m: Matrix2<FBig<HalfEven>> = matrix_multiply_2x2(&m1, &d3);
        let ellipse = Ellipse::new(m, p);
        Self {
            _theta: theta_for_debug,
            _epsilon: epsilon,
            scale,
            d,
            z_x,
            z_y,
            ellipse,
        }
    }
}

impl Region for EpsilonRegion {
    fn ellipse(&self) -> Ellipse {
        self.ellipse.clone()
    }

    // Return true if `u` is inside shaded region in figure in eq (14) in R + S
    // The radius is 1 in the figure.
    // For "up to phase" it is scaled by |δ|^2 = 2 + √2.
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

#[derive(Debug)]
pub struct UnitDisk {
    scale: ZRootTwo,
    ellipse: Ellipse,
}

impl UnitDisk {
    pub fn new(scale: ZRootTwo) -> Self {
        let s_inv: FBig<HalfEven> = 1 / scale.to_real();
        let ellipse = Ellipse::from(
            s_inv.clone(),
            ib_to_bf_prec(IBig::ZERO),
            ib_to_bf_prec(IBig::ZERO),
            s_inv.clone(),
            ib_to_bf_prec(IBig::ZERO),
            ib_to_bf_prec(IBig::ZERO),
        );
        Self { scale, ellipse }
    }

    pub fn ellipse(&self) -> &Ellipse {
        &self.ellipse
    }
}

impl Region for UnitDisk {
    fn ellipse(&self) -> Ellipse {
        self.ellipse.clone()
    }
    fn inside(&self, u: &DOmega) -> bool {
        DRootTwo::from_domega(u.conj() * u) <= DRootTwo::from_zroottwo(self.scale.clone())
    }

    fn intersect(&self, u0: &DOmega, v: &DOmega) -> Option<(FBig<HalfEven>, FBig<HalfEven>)> {
        let a = v.conj() * v;
        let b = 2 * (v.conj() * u0);
        let c = u0.conj() * u0 - DOmega::from_zroottwo(&self.scale);
        solve_quadratic(a.real(), b.real(), c.real())
    }
}

/// Lightweight counters over the (lazy) candidate stream examined during the
/// Diophantine search inside [`search_for_solution`]. Intentionally does not require
/// materializing the candidate iterator (e.g. via `.len()`/`.count()`) -- each field is
/// incremented one candidate at a time, inside the loop that is already consuming the
/// iterator lazily.
///
/// Not yet wired into any public entry point; this exists so a future stage can surface
/// search diagnostics (e.g. for tuning mixed/fallback region predicates) without having
/// to re-plumb the search loop.
#[derive(Debug, Default)]
pub struct CandidateStats {
    /// Number of candidates pulled from the `solve_tdgp` iterator.
    pub examined: usize,
    /// Number of times `diophantine_dyadic` was invoked on a candidate.
    pub diophantine_attempts: usize,
    /// Number of candidates for which `diophantine_dyadic` found a solution.
    pub solved: usize,
}

/// The output of [`to_upright_set_pair`](crate::to_upright::to_upright_set_pair): the grid
/// operator that maps the original region pair to an "upright" pair, the transformed
/// ellipses, and their axis-aligned bounding boxes.
pub struct UprightTransform {
    pub op_g: GridOp,
    pub ellipse_a: Ellipse,
    pub ellipse_b: Ellipse,
    pub bbox_a: Rectangle,
    pub bbox_b: Rectangle,
}

pub(crate) fn process_solution_candidate(
    mut z: DOmega,
    mut w: DOmega,
    phase: PhaseMode,
) -> DOmegaUnitary {
    z = z.reduce_denomexp();
    w = w.reduce_denomexp();

    match z.k.cmp(&w.k) {
        Ordering::Greater => {
            w = w.renew_denomexp(z.k);
        }
        Ordering::Less => {
            z = z.renew_denomexp(w.k);
        }
        Ordering::Equal => {}
    }

    match phase {
        // Question: this is a bit different from pygridsynth
        PhaseMode::Exact => {
            if (z.clone() + w.clone()).reduce_denomexp().k < z.k {
                DOmegaUnitary::new(z, w, 0, None)
            } else {
                DOmegaUnitary::new(z, w.mul_by_omega(), 0, None)
            }
        }
        PhaseMode::Shifted => {
            // todo: remove clones
            let k1 = (z.clone() + w.clone()).reduce_denomexp().k;
            let k2 = (z.clone() + w.mul_by_omega()).reduce_denomexp().k;
            let k3 = (z.clone() + w.mul_by_omega_inv()).reduce_denomexp().k;

            if k1 <= k2.min(k3) {
                DOmegaUnitary::new(z, w, 7, None)
            } else {
                DOmegaUnitary::new(z, w.mul_by_omega_inv(), 7, None)
            }
        }
    }
}

pub(crate) fn process_solutions<I>(
    config: &mut GridSynthConfig,
    solutions: I,
    time_of_diophantine_dyadic: &mut Duration,
    phase: PhaseMode,
    mut stats: Option<&mut CandidateStats>,
) -> Option<DOmegaUnitary>
where
    I: Iterator<Item = DOmega>,
{
    let start_diophantine = if config.measure_time {
        Some(Instant::now())
    } else {
        None
    };

    for z in solutions {
        if let Some(s) = stats.as_deref_mut() {
            s.examined += 1;
        }

        if (&z * z.conj()).residue() == 0 {
            continue;
        }

        let z_with_phase = match phase {
            PhaseMode::Exact => z.clone(),
            // todo: make constant
            PhaseMode::Shifted => {
                &z * &DOmega::new(
                    ZOmega::new(IBig::from(0), IBig::from(-1), IBig::from(1), IBig::from(0)),
                    1,
                )
            }
        };

        let xi = DRootTwo::from_int(IBig::ONE)
            - DRootTwo::from_domega(z_with_phase.conj() * &z_with_phase);
        if let Some(s) = stats.as_deref_mut() {
            s.diophantine_attempts += 1;
        }
        if let Some(w_val) = diophantine_dyadic(xi, &mut config.diophantine_data) {
            if let Some(start) = start_diophantine {
                *time_of_diophantine_dyadic += start.elapsed();
                if config.measure_time {
                    debug!(
                        "time of diophantine_dyadic: {:.3} ms",
                        time_of_diophantine_dyadic.as_secs_f64() * 1000.0
                    );
                }
            }
            if config.verbose {
                debug!("------------------");
            }
            if let Some(s) = stats.as_deref_mut() {
                s.solved += 1;
            }
            return Some(process_solution_candidate(z_with_phase, w_val, phase));
        }
    }

    if let Some(start) = start_diophantine {
        *time_of_diophantine_dyadic += start.elapsed();
    }
    None
}

/// Constructs the concrete `EpsilonRegion`/`UnitDisk` pair for the (theta, epsilon, phase)
/// triple. This is the only place that needs to know the concrete region types used by
/// `gridsynth`; `setup_regions_and_transform` and `search_for_solution` below are generic
/// over the first region so a future sibling `Region` implementation can reuse them.
pub(crate) fn epsilon_region_and_unit_disk(
    theta: FBig<HalfEven>,
    epsilon: FBig<HalfEven>,
    phase: PhaseMode,
) -> (EpsilonRegion, UnitDisk) {
    let epsilon_region_scale = match phase {
        PhaseMode::Exact => ZRootTwo {
            a: IBig::from(1),
            b: IBig::from(0),
        },
        PhaseMode::Shifted => ZRootTwo {
            a: IBig::from(2),
            b: IBig::from(1),
        },
    };

    let unit_disk_scale = match phase {
        PhaseMode::Exact => ZRootTwo {
            a: IBig::from(1),
            b: IBig::from(0),
        },
        PhaseMode::Shifted => ZRootTwo {
            a: IBig::from(2),
            b: IBig::from(-1),
        },
    };

    let epsilon_region = EpsilonRegion::new(theta, epsilon, epsilon_region_scale);
    let unit_disk = UnitDisk::new(unit_disk_scale);
    (epsilon_region, unit_disk)
}

pub(crate) fn setup_regions_and_transform<A: Region + std::fmt::Debug>(
    set_a: &A,
    set_b: &UnitDisk,
    verbose: bool,
    measure_time: bool,
) -> UprightTransform {
    let start_upright = if measure_time {
        Some(Instant::now())
    } else {
        None
    };
    let (op_g, ellipse_a, ellipse_b, bbox_a, bbox_b) = to_upright_set_pair(set_a, set_b, verbose);
    if let Some(start) = start_upright {
        if measure_time {
            debug!(
                "to_upright_set_pair: {:.3} s",
                start.elapsed().as_secs_f64()
            );
        }
    }

    if verbose {
        debug!("------------------");
    }

    UprightTransform {
        op_g,
        ellipse_a,
        ellipse_b,
        bbox_a,
        bbox_b,
    }
}

pub(crate) fn search_for_solution<A: Region + std::fmt::Debug>(
    epsilon_region: &A,
    unit_disk: &UnitDisk,
    transformed: &UprightTransform,
    config: &mut GridSynthConfig,
    phase: PhaseMode,
    mut stats: Option<&mut CandidateStats>,
) -> Option<DOmegaUnitary> {
    // A generous upper bound on the number of grid-refinement steps `k` the search will
    // try before giving up. Working precision is already set as a multiple of
    // log2(1/epsilon), and legitimate solutions are found at `k` roughly proportional to
    // log2(1/epsilon) too, so this bound is effectively unreachable for any well-formed
    // (theta, epsilon, region) input -- it exists only to fail loudly, in finite time,
    // if a `Region` predicate is ever incorrect (in which case no solution may exist and
    // the loop below would otherwise spin forever).
    let max_k = 4 * get_prec_bits() as i64;

    let mut k = 0;
    let mut time_of_solve_tdgp = Duration::ZERO;
    let mut time_of_diophantine_dyadic = Duration::ZERO;

    while k <= max_k {
        let start_tdgp = if config.measure_time {
            Some(Instant::now())
        } else {
            None
        };
        let solutions = solve_tdgp(
            epsilon_region,
            unit_disk,
            &transformed.op_g,
            &transformed.bbox_a,
            &transformed.bbox_b,
            k,
            config.verbose,
        );
        // TODO: Reenable
        // if config.verbose {
        //     // Warning! Printing the length will materialize a potentially large iterator.
        //     let lensol = match &solutions {
        //         None => 0,
        //         Some(sols) => sols.len(),
        //     };
        //     info!("k = {}, found {} candidates", k, lensol);
        // }
        if let Some(start) = start_tdgp {
            time_of_solve_tdgp += start.elapsed();
        }
        if let Some(solutions) = solutions {
            if let Some(result) = process_solutions(
                config,
                solutions,
                &mut time_of_diophantine_dyadic,
                phase,
                stats.as_deref_mut(),
            ) {
                if config.measure_time {
                    debug!(
                        "time of solve_TDGP: {:.3} ms",
                        time_of_solve_tdgp.as_secs_f64() * 1000.0
                    );
                }
                return Some(result);
            }
        }
        k += 1;
    }
    None
}

/// Core gridsynth algorithm that finds an optimal Clifford+T approximation.
///
/// # Arguments
/// * `theta` - The rotation angle to approximate
/// * `epsilon` - The approximation tolerance
/// * `diophantine_timeout` - Timeout for diophantine equation solving (ms)
/// * `factoring_timeout` - Timeout for integer factoring (ms)
/// * `verbose` - Enable verbose output
/// * `measure_time` - Enable timing measurements
///
/// # Returns
/// `Some(DOmegaUnitary)` representing the optimal Clifford+T approximation, or `None` if
/// the search exceeded its (very generous) internal bound on `k` without finding a
/// solution -- see [`search_for_solution`] for why that can only happen if a `Region`
/// predicate is broken.
pub(crate) fn gridsynth(config: &mut GridSynthConfig, phase: PhaseMode) -> Option<DOmegaUnitary> {
    let (epsilon_region, unit_disk) =
        epsilon_region_and_unit_disk(config.theta.clone(), config.epsilon.clone(), phase);
    let transformed = setup_regions_and_transform(
        &epsilon_region,
        &unit_disk,
        config.verbose,
        config.measure_time,
    );

    search_for_solution(
        &epsilon_region,
        &unit_disk,
        &transformed,
        config,
        phase,
        None,
    )
}

/// Public wrapper around the core gridsynth search that returns the synthesized
/// `DOmegaUnitary` directly, rather than a decomposed Clifford+T gate string. Intended for
/// callers that need the raw unitary (e.g. future mixed-diagonal/fallback synthesis
/// protocols composing multiple `gridsynth` calls).
///
/// # Panics
/// Panics if the internal search exceeds its generous bound on `k` without finding a
/// solution. This is not expected to trigger for any well-formed `(theta, epsilon)` input;
/// see [`search_for_solution`] for details. A panic here is preferable to the unbounded
/// hang this replaces.
pub fn gridsynth_unitary(config: &mut GridSynthConfig, phase: PhaseMode) -> DOmegaUnitary {
    gridsynth(config, phase).expect(
        "gridsynth: exceeded max_k without finding a solution — region predicate is likely incorrect",
    )
}

pub fn gridsynth_gates(config: &mut GridSynthConfig) -> GridSynthResult {
    // let start_total = if config.measure_time {
    //     Some(Instant::now())
    // } else {
    //     None
    // };

    // let start_decompose = if config.measure_time {
    //     Some(Instant::now())
    // } else {
    //     None
    // };

    if !config.up_to_phase {
        // exact synthesis only
        let u_approx = gridsynth_unitary(config, PhaseMode::Exact);
        let gates = decompose_domega_unitary(u_approx);

        GridSynthResult {
            gates,
            global_phase: false,
        }
    } else {
        // exact synthesis
        let u_approx = gridsynth_unitary(config, PhaseMode::Exact);
        let gates_exact = decompose_domega_unitary(u_approx);
        let t_count_exact = gates_exact.chars().filter(|&c| c == 'T').count();

        // also shifted synthesis
        let u_approx = gridsynth_unitary(config, PhaseMode::Shifted);
        let gates_shifted = decompose_domega_unitary(u_approx);
        let t_count_shifted = gates_shifted.chars().filter(|&c| c == 'T').count();

        if t_count_exact <= t_count_shifted {
            GridSynthResult {
                gates: gates_exact,
                global_phase: false,
            }
        } else {
            GridSynthResult {
                gates: gates_shifted,
                global_phase: true,
            }
        }
    }
}

impl crate::accuracy::AchievedDiamondError for GridSynthResult {
    /// Diamond-norm distance between this result's synthesized channel and the ideal
    /// Z-rotation by `theta`, decoded on demand from `self.gates` (accounting for the extra
    /// `e^{i pi/8}` phase `self.global_phase` records having been used, per `PhaseMode`).
    fn achieved_diamond_error(&self, theta: &FBig<HalfEven>) -> FBig<HalfEven> {
        crate::accuracy::gate_string_diamond_error(theta, &self.gates, self.global_phase)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::config_from_theta_epsilon;
    use serial_test::serial;

    /// Exercises the `Option<&mut CandidateStats>` path threaded through
    /// `search_for_solution`/`process_solutions` (not yet reachable from any public API).
    /// Confirms the counters are populated sensibly for a real search, without ever
    /// materializing the underlying lazy candidate iterator.
    #[test]
    #[serial]
    fn candidate_stats_are_populated_when_provided() {
        crate::clear_caches();
        let theta = std::f64::consts::PI / 8.0;
        let epsilon = 1e-10;
        let mut config = config_from_theta_epsilon(theta, epsilon, 1234, false, false);

        let (epsilon_region, unit_disk) = epsilon_region_and_unit_disk(
            config.theta.clone(),
            config.epsilon.clone(),
            PhaseMode::Exact,
        );
        let transformed = setup_regions_and_transform(
            &epsilon_region,
            &unit_disk,
            config.verbose,
            config.measure_time,
        );

        let mut stats = CandidateStats::default();
        let result = search_for_solution(
            &epsilon_region,
            &unit_disk,
            &transformed,
            &mut config,
            PhaseMode::Exact,
            Some(&mut stats),
        );

        assert!(result.is_some(), "expected a solution to be found");
        assert!(
            stats.examined > 0,
            "expected at least one candidate to be examined"
        );
        assert!(
            stats.diophantine_attempts > 0,
            "expected at least one diophantine attempt"
        );
        assert_eq!(
            stats.solved, 1,
            "expected exactly one solved candidate on success"
        );
    }
}
