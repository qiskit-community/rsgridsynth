// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::common::{cos_fbig, fb_with_prec, get_prec_bits, ib_to_bf_prec, sin_fbig};
use crate::config::{GridSynthConfig, GridSynthResult};
use crate::diophantine::diophantine_dyadic;
use crate::math::solve_quadratic;
use crate::math::sqrt_fbig;
use crate::region::Ellipse;
use crate::ring::{DOmega, DRootTwo, ZOmega, ZRootTwo};
use crate::synthesis_of_clifford_t::decompose_domega_unitary;
use crate::tdgp::solve_tdgp;
use crate::tdgp::Region;
use crate::to_upright::to_upright_set_pair;
use crate::unitary::DOmegaUnitary;
use dashu_base::SquareRoot;
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;

//use log::{debug, info};
use log::debug;

use nalgebra::{Matrix2, Vector2};
use num::Complex;
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
// `PhaseMode::Shifted` and keep the one with lower T count. (We don't compute the best
// exact solution and then the best with the phase factor. Rather we interleave candidates
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

fn rotation_mat(theta: &FBig<HalfEven>) -> Matrix2<Complex<FBig<HalfEven>>> {
    let two = fb_with_prec(FBig::try_from(2.0).unwrap());
    let theta_half = fb_with_prec(theta / &two);
    let neg_theta_half = -fb_with_prec(theta_half);
    let z_x: FBig<HalfEven> = fb_with_prec(cos_fbig(&neg_theta_half));
    let z_y: FBig<HalfEven> = fb_with_prec(sin_fbig(&neg_theta_half));
    let neg_z_y: FBig<HalfEven> = -fb_with_prec(z_y.clone());
    let zero: FBig<HalfEven> = ib_to_bf_prec(IBig::ZERO);

    Matrix2::new(
        Complex::new(z_x.clone(), z_y.clone()),
        Complex::new(zero.clone(), zero.clone()),
        Complex::new(zero.clone(), zero.clone()),
        Complex::new(z_x.clone(), neg_z_y.clone()),
    )
}

fn mult_complex_nums(
    u: &Complex<FBig<HalfEven>>,
    v: &Complex<FBig<HalfEven>>,
) -> Complex<FBig<HalfEven>> {
    let re = &u.re * &v.re - &u.im * &v.im;
    let im = &u.re * &v.im + &u.im * &v.re;
    Complex::new(re, im)
}

fn to_fbig(x: f64) -> FBig<HalfEven> {
    FBig::<HalfEven>::try_from(x)
        .unwrap()
        .with_precision(get_prec_bits())
        .value()
}

/// Checks correctness of the synthesized circuit.
fn check_solution(
    gates: &str,
    theta: &FBig<HalfEven>,
    epsilon: &FBig<HalfEven>,
    phase: PhaseMode,
) -> bool {
    let expected = rotation_mat(theta);
    let synthesized = DOmegaUnitary::from_gates(gates).to_complex_matrix();

    // x = e^{-i theta / 2}
    let x = expected[(0, 0)].clone();
    let u = match phase {
        PhaseMode::Exact => synthesized[(0, 0)].clone(),
        PhaseMode::Shifted => {
            let p = to_fbig(std::f64::consts::PI / 8.);
            let phase = Complex::new(cos_fbig(&p), sin_fbig(&p));
            mult_complex_nums(&synthesized[(0, 0)].clone(), &phase)
        }
    };

    // This computes the eigenvalues of A^* A, where A = expected - synthesized.
    // The operator norm is the square root of this.
    let eig: FBig<HalfEven> = 2 - 2 * (&x.re * &u.re + &x.im * &u.im);

    // Due to small numerical imprecisions when computing cos(\theta / 2), it may happen that we
    // get a slightly negative value.
    let eig = eig.max(FBig::from(0));

    // Compute the norm.
    let norm = fb_with_prec(eig.sqrt());

    let res = norm < *epsilon;

    println!("solution is correct: {:?}", res);

    res
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
        let one = fb_with_prec(FBig::try_from(1.0).unwrap());
        let two = fb_with_prec(FBig::try_from(2.0).unwrap());
        let four = fb_with_prec(FBig::try_from(4.0).unwrap());
        let epsilon_squared = fb_with_prec(&epsilon * &epsilon);
        let half_eps_sq = fb_with_prec(&epsilon_squared / &four);
        let one_minus_half_eps_sq = one - half_eps_sq;
        let scale_to_real = scale.to_real();
        let d = fb_with_prec(sqrt_fbig(&one_minus_half_eps_sq) * sqrt_fbig(&scale_to_real));

        let theta_half = fb_with_prec(&theta / &two);
        let neg_theta_half = -fb_with_prec(theta_half);
        let z_x: FBig<HalfEven> = fb_with_prec(cos_fbig(&neg_theta_half));
        let z_y: FBig<HalfEven> = fb_with_prec(sin_fbig(&neg_theta_half));
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
            _theta: theta,
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

fn process_solution_candidate(mut z: DOmega, mut w: DOmega, phase: PhaseMode) -> DOmegaUnitary {
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

fn process_solutions<I>(
    config: &mut GridSynthConfig,
    solutions: I,
    time_of_diophantine_dyadic: &mut Duration,
    phase: PhaseMode,
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
            return Some(process_solution_candidate(z_with_phase, w_val, phase));
        }
    }

    if let Some(start) = start_diophantine {
        *time_of_diophantine_dyadic += start.elapsed();
    }
    None
}

fn setup_regions_and_transform(
    theta: FBig<HalfEven>,
    epsilon: FBig<HalfEven>,
    verbose: bool,
    measure_time: bool,
    phase: PhaseMode,
) -> (
    EpsilonRegion,
    UnitDisk,
    (
        crate::grid_op::GridOp,
        crate::region::Ellipse,
        crate::region::Ellipse,
        crate::region::Rectangle,
        crate::region::Rectangle,
    ),
) {
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

    let start_upright = if measure_time {
        Some(Instant::now())
    } else {
        None
    };
    let transformed = to_upright_set_pair(&epsilon_region, &unit_disk, verbose);
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

    (epsilon_region, unit_disk, transformed)
}

fn search_for_solution(
    epsilon_region: &EpsilonRegion,
    unit_disk: &UnitDisk,
    transformed: &(
        crate::grid_op::GridOp,
        crate::region::Ellipse,
        crate::region::Ellipse,
        crate::region::Rectangle,
        crate::region::Rectangle,
    ),
    config: &mut GridSynthConfig,
    phase: PhaseMode,
) -> DOmegaUnitary {
    let mut k = 0;
    let mut time_of_solve_tdgp = Duration::ZERO;
    let mut time_of_diophantine_dyadic = Duration::ZERO;

    loop {
        let start_tdgp = if config.measure_time {
            Some(Instant::now())
        } else {
            None
        };
        let solutions = solve_tdgp(
            epsilon_region,
            unit_disk,
            &transformed.0,
            &transformed.3,
            &transformed.4,
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
            if let Some(result) =
                process_solutions(config, solutions, &mut time_of_diophantine_dyadic, phase)
            {
                if config.measure_time {
                    debug!(
                        "time of solve_TDGP: {:.3} ms",
                        time_of_solve_tdgp.as_secs_f64() * 1000.0
                    );
                }
                return result;
            }
        }
        k += 1;
    }
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
/// A DOmegaUnitary representing the optimal Clifford+T approximation
fn gridsynth(config: &mut GridSynthConfig, phase: PhaseMode) -> DOmegaUnitary {
    let (epsilon_region, unit_disk, transformed) = setup_regions_and_transform(
        config.theta.clone(),
        config.epsilon.clone(),
        config.verbose,
        config.measure_time,
        phase,
    );

    search_for_solution(&epsilon_region, &unit_disk, &transformed, config, phase)
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
        let u_approx = gridsynth(config, PhaseMode::Exact);
        let gates = decompose_domega_unitary(u_approx);

        // Peform validation check, if required.
        let is_correct = if config.check_solution {
            Some(check_solution(
                &gates,
                &config.theta,
                &config.epsilon,
                PhaseMode::Exact,
            ))
        } else {
            None
        };

        GridSynthResult {
            gates,
            global_phase: false,
            is_correct,
        }
    } else {
        // exact synthesis
        let u_approx = gridsynth(config, PhaseMode::Exact);
        let gates_exact = decompose_domega_unitary(u_approx);

        // also shifted synthesis
        let u_approx = gridsynth(config, PhaseMode::Shifted);
        let gates_shifted = decompose_domega_unitary(u_approx);

        println!(
            "gates_exact = {:?} of len {:?}",
            gates_exact,
            gates_exact.len()
        );
        println!(
            "gates_shifted = {:?} of len {:?}",
            gates_shifted,
            gates_shifted.len()
        );

        // Peform validation check, if required.
        // For now, make sure to check both sequences.
        let is_correct = if config.check_solution {
            Some(
                check_solution(
                    &gates_exact,
                    &config.theta,
                    &config.epsilon,
                    PhaseMode::Exact,
                ) && check_solution(
                    &gates_shifted,
                    &config.theta,
                    &config.epsilon,
                    PhaseMode::Shifted,
                ),
            )
        } else {
            None
        };

        if gates_exact.len() < gates_shifted.len() {
            GridSynthResult {
                gates: gates_exact,
                global_phase: false,
                is_correct,
            }
        } else {
            GridSynthResult {
                gates: gates_shifted,
                global_phase: true,
                is_correct,
            }
        }
    }

    // if let Some(start) = start_decompose {
    //     if config.measure_time {
    //         debug!(
    //             "time of decompose_domega_unitary: {:.3} ms",
    //             start.elapsed().as_secs_f64() * 1000.0
    //         );
    //     }
    // }
    // if let Some(start) = start_total {
    //     if config.measure_time {
    //         debug!(
    //             "total time: {:.3} ms",
    //             start.elapsed().as_secs_f64() * 1000.0
    //         );
    //     }
    // }
}
