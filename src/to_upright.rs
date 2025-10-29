// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use log::debug;
use std::fmt::Debug;

use crate::common::ib_to_bf_prec;
use crate::grid_op::{EllipsePair, GridOp};
use crate::math::{floorsqrt, log};
use crate::region::{Ellipse, Rectangle};
use crate::ring::{ZOmega, ZRootTwo};
use crate::tdgp::Region;
use dashu_int::IBig;
use num_traits::Pow;
const LAMBDA: ZRootTwo = ZRootTwo {
    a: IBig::ONE,
    b: IBig::ONE,
};

fn reduction(
    ellipse_pair: EllipsePair,
    op_g_l: &GridOp,
    op_g_r: &GridOp,
    new_op_g: &GridOp,
) -> (EllipsePair, GridOp, GridOp, bool) {
    (
        new_op_g.clone() * ellipse_pair,
        op_g_l.clone(),
        new_op_g.clone() * op_g_r.clone(),
        false,
    )
}

fn shift_ellipse_pair(mut ep: EllipsePair, n: i32) -> EllipsePair {
    let lambda_n = LAMBDA.pow(&IBig::from(n));
    let lambda_inv_n = LAMBDA.pow(&IBig::from(-n));
    let lambda_n_f = lambda_n.to_real();
    let lambda_inv_n_f = lambda_inv_n.to_real();

    ep.a.d[(0, 0)] = &ep.a.d[(0, 0)] * &lambda_inv_n_f;
    ep.a.d[(1, 1)] = &ep.a.d[(1, 1)] * &lambda_n_f;
    ep.b.d[(0, 0)] = &ep.b.d[(0, 0)] * &lambda_n_f;
    ep.b.d[(1, 1)] = &ep.b.d[(1, 1)] * &lambda_inv_n_f;
    if n & 1 == 1 {
        let b_val = ep.b.d[(0, 1)].clone();
        ep.b.d[(0, 1)] = -b_val;
    }
    ep
}

pub fn step_lemma(
    mut ellipse_pair: EllipsePair,
    op_g_l: &GridOp,
    op_g_r: &GridOp,
    verbose: bool,
) -> (EllipsePair, GridOp, GridOp, bool) {
    let a = &ellipse_pair.a;
    let b = &ellipse_pair.b;
    let check0 = ib_to_bf_prec(IBig::from(33971)) / ib_to_bf_prec(IBig::from(1000)); // 33.971
    let check1 = ib_to_bf_prec(IBig::from(29437)) / ib_to_bf_prec(IBig::from(1000000)); // 0.029437
    let check2 = ib_to_bf_prec(IBig::from(58285)) / ib_to_bf_prec(IBig::from(10000)); // 5.8285
    let check3 = ib_to_bf_prec(IBig::from(17157)) / ib_to_bf_prec(IBig::from(100000)); // 0.17157
    let check4 = ib_to_bf_prec(IBig::from(2441)) / ib_to_bf_prec(IBig::from(10000)); // 0.2441
    let check5 = ib_to_bf_prec(IBig::from(40968)) / ib_to_bf_prec(IBig::from(10000)); // 4.0968
    let check6 = ib_to_bf_prec(IBig::from(16969)) / ib_to_bf_prec(IBig::from(10000)); // 1.6969
                                                                                      // println!("ellipse.bias(): {}", ellipse_pair.bias());
    if b.b() < &ib_to_bf_prec(IBig::ZERO) {
        if verbose {
            debug!("Z");
        }
        let op_z = GridOp::new(
            ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
            ZOmega::new(IBig::ZERO, IBig::NEG_ONE, IBig::ZERO, IBig::ZERO),
        );
        reduction(ellipse_pair, op_g_l, op_g_r, &op_z)
    } else if a.bias() * b.bias() < ib_to_bf_prec(IBig::ONE) {
        if verbose {
            debug!("X");
        }
        let op_x = GridOp::new(
            ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
            ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
        );
        reduction(ellipse_pair, op_g_l, op_g_r, &op_x)
    } else if ellipse_pair.bias() > check0 || ellipse_pair.bias() < check1 {
        let n: IBig = ((log(ellipse_pair.bias()) / log(LAMBDA.to_real()))
            / ib_to_bf_prec(IBig::from(8)))
        .round()
        .try_into()
        .unwrap();
        let op_s = GridOp::new(
            ZOmega::new(IBig::NEG_ONE, IBig::ZERO, IBig::ONE, IBig::ONE),
            ZOmega::new(IBig::ONE, IBig::NEG_ONE, IBig::ONE, IBig::ZERO),
        )
        .pow(n.clone());
        if verbose {
            debug!("S (n={})", n);
        }
        reduction(ellipse_pair, op_g_l, op_g_r, &op_s)
    } else if ellipse_pair.skew() <= ib_to_bf_prec(IBig::from(15)) {
        (ellipse_pair, op_g_l.clone(), op_g_r.clone(), true)
    } else if ellipse_pair.bias() > check2 || ellipse_pair.bias() < check3 {
        let n: i32 = ((log(ellipse_pair.bias()) / log(LAMBDA.to_real()))
            / ib_to_bf_prec(IBig::from(4)))
        .round()
        .try_into()
        .unwrap();
        ellipse_pair = shift_ellipse_pair(ellipse_pair, n);
        if verbose {
            debug!("sigma (n={})", n);
        }
        let op_sigma_l = if n >= 0 {
            GridOp::new(
                ZOmega::new(IBig::NEG_ONE, IBig::ZERO, IBig::ONE, IBig::ONE),
                ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
            )
            .pow(n)
        } else {
            GridOp::new(
                ZOmega::new(IBig::NEG_ONE, IBig::ZERO, IBig::ONE, IBig::NEG_ONE),
                ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
            )
            .pow(-n)
        };
        let op_sigma_r = if n >= 0 {
            GridOp::new(
                ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
                ZOmega::new(IBig::ONE, IBig::NEG_ONE, IBig::ONE, IBig::ZERO),
            )
            .pow(n)
        } else {
            GridOp::new(
                ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
                ZOmega::new(IBig::ONE, IBig::ONE, IBig::ONE, IBig::ZERO),
            )
            .pow(-n)
        };
        (
            ellipse_pair,
            op_g_l.clone() * op_sigma_l,
            op_sigma_r * op_g_r.clone(),
            false,
        )
    } else if (check4.clone()..=check5.clone()).contains(&a.bias())
        && (check4..=check5).contains(&b.bias())
    {
        if verbose {
            debug!("R");
        }
        let op_r = GridOp::new(
            ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ONE, IBig::ZERO),
            ZOmega::new(IBig::ONE, IBig::ZERO, IBig::ZERO, IBig::ZERO),
        );
        reduction(ellipse_pair, op_g_l, op_g_r, &op_r)
    } else if a.b() >= &ib_to_bf_prec(IBig::ZERO) && a.bias() <= check6.clone() {
        if verbose {
            debug!("K");
        }
        let op_k = GridOp::new(
            ZOmega::new(IBig::NEG_ONE, IBig::NEG_ONE, IBig::ZERO, IBig::ZERO),
            ZOmega::new(IBig::ZERO, IBig::NEG_ONE, IBig::ONE, IBig::ZERO),
        );
        reduction(ellipse_pair, op_g_l, op_g_r, &op_k)
    } else if a.b() >= &ib_to_bf_prec(IBig::ZERO) && b.bias() <= check6 {
        if verbose {
            debug!("K_conj_sq2");
        }
        let op_kc = GridOp::new(
            ZOmega::new(IBig::ONE, IBig::NEG_ONE, IBig::ZERO, IBig::ZERO),
            ZOmega::new(IBig::ZERO, IBig::NEG_ONE, IBig::NEG_ONE, IBig::ZERO),
        );
        reduction(ellipse_pair, op_g_l, op_g_r, &op_kc)
    } else if a.b() >= &ib_to_bf_prec(IBig::ZERO) {
        let n: IBig =
            floorsqrt((a.bias().min(b.bias())) / ib_to_bf_prec(IBig::from(4))).max(IBig::ONE);
        if verbose {
            debug!("A (n={})", n);
        }
        let op_a = GridOp::new(
            ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
            ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, 2 * n),
        );
        reduction(ellipse_pair, op_g_l, op_g_r, &op_a)
    } else {
        let n: IBig =
            floorsqrt((a.bias().min(b.bias())) / ib_to_bf_prec(IBig::from(2))).max(IBig::ONE);
        if verbose {
            debug!("B (n={})", n);
        }
        let op_b = GridOp::new(
            ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
            ZOmega::new(n.clone(), IBig::ONE, -n, IBig::ZERO),
        );
        reduction(ellipse_pair, op_g_l, op_g_r, &op_b)
    }
}

fn to_upright_ellipse_pair(a: &Ellipse, b: &Ellipse, verbose: bool) -> GridOp {
    if verbose {
        debug!(
            "Converting ellipse pair to upright: A = {:?}, B = {:?}",
            a, b
        );
    }
    let mut ep = EllipsePair::new(a.normalize(), b.normalize());
    let op_i = GridOp::new(
        ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
        ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
    );
    let mut op_g_l = op_i.clone();
    let mut op_g_r = op_i;

    loop {
        let (next_ep, next_l, next_r, end) = step_lemma(ep, &op_g_l, &op_g_r, verbose);
        ep = next_ep;
        op_g_l = next_l;
        op_g_r = next_r;
        if end {
            break;
        }
    }
    op_g_l * op_g_r
}

pub fn to_upright_set_pair(
    set_a: &(impl Region + Debug),
    set_b: &(impl Region + Debug),
    verbose: bool,
) -> (GridOp, Ellipse, Ellipse, Rectangle, Rectangle) {
    let op_g = to_upright_ellipse_pair(&set_a.ellipse(), &set_b.ellipse(), verbose);
    let ellipse_pair = op_g.clone() * EllipsePair::new(set_a.ellipse(), set_b.ellipse());
    let ellipse_a_upright = ellipse_pair.a;
    let ellipse_b_upright = ellipse_pair.b;
    let bbox_a = ellipse_a_upright.bbox();
    let bbox_b = ellipse_b_upright.bbox();

    (op_g, ellipse_a_upright, ellipse_b_upright, bbox_a, bbox_b)
}
