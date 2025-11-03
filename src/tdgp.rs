// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;

use crate::common::{fb_with_prec, ib_to_bf_prec};
use crate::grid_op::GridOp;
use crate::odgp::{solve_scaled_odgp, solve_scaled_odgp_with_parity_k_ne_0};
use crate::region::{Ellipse, Interval, Rectangle};
use crate::ring::{DOmega, DRootTwo};

pub trait Region {
    fn ellipse(&self) -> Ellipse;
    fn inside(&self, u: &DOmega) -> bool;
    fn intersect(&self, u: &DOmega, v: &DOmega) -> Option<(FBig<HalfEven>, FBig<HalfEven>)>;
}

pub fn solve_tdgp(
    set_a: &impl Region,
    set_b: &impl Region,
    op_g: &GridOp,
    bbox_a: &Rectangle,
    bbox_b: &Rectangle,
    k: i64,
    _verbose: bool,
) -> Vec<DOmega> {
    let mut sol_sufficient = Vec::with_capacity(100); // Pre-allocate reasonable capacity

    let mut sol_x = solve_scaled_odgp(bbox_a.x.clone(), bbox_b.x.clone(), k + 1);

    let sol_y = solve_scaled_odgp(
        bbox_a
            .y
            .fatten(&(bbox_a.y.width() / ib_to_bf_prec(IBig::from(10000)))),
        bbox_b
            .y
            .fatten(&(bbox_b.y.width() / ib_to_bf_prec(IBig::from(10000)))),
        k + 1,
    );

    let alpha0 = match sol_x.next() {
        Some(val) => val,
        None => return vec![],
    };

    let droot_zero = DRootTwo::from_int(IBig::ZERO);
    let _k_ibig = IBig::from(k);
    let dx = DRootTwo::power_of_inv_sqrt2(k);
    let op_g_inv_result = op_g.inv();
    let op_g_inv = op_g_inv_result.as_ref().unwrap();
    let zero_droottwo = DRootTwo::from_int(IBig::ZERO);

    let v = op_g_inv.clone() * DOmega::from_droottwo_vector(&dx, &zero_droottwo, k);
    let v_conj_sq2 = v.conj_sq2();
    for beta in sol_y {
        let dx = DRootTwo::power_of_inv_sqrt2(k);
        let z0 = op_g.inv().as_ref().unwrap().clone()
            * DOmega::from_droottwo_vector(&alpha0, &beta, k + 1);
        let v = op_g.inv().as_ref().unwrap().clone()
            * DOmega::from_droottwo_vector(&dx, &droot_zero, k);
        let t_a = set_a.intersect(&z0, &v);
        let t_b = set_b.intersect(z0.conj_sq2(), v_conj_sq2);
        if t_a.is_none() || t_b.is_none() {
            continue;
        }
        let (t_a, t_b) = (t_a.unwrap(), t_b.unwrap());

        let parity = (&beta - &alpha0).mul_by_sqrt2_power_renewing_denomexp(k);
        let (mut int_a, mut int_b) = (Interval::new(t_a.0, t_a.1), Interval::new(t_b.0, t_b.1));
        let dt_a = {
            let ten = ib_to_bf_prec(IBig::from(10));
            let max_val = {
                let shift_k = IBig::ONE << (k as usize);
                let width_product = shift_k * int_b.width();
                if ten > width_product {
                    ten.clone()
                } else {
                    width_product
                }
            };
            fb_with_prec(&ten / &max_val)
        };
        let dt_b = {
            let ten = ib_to_bf_prec(IBig::from(10));
            let max_val = {
                let shift_k = IBig::from(1) << (k as usize);
                let width_product = shift_k * int_a.width();
                if ten > width_product {
                    ten.clone()
                } else {
                    width_product
                }
            };
            fb_with_prec(&ten / &max_val)
        };

        int_a = int_a.fatten(&dt_a);
        int_b = int_b.fatten(&dt_b);

        let sol_t = solve_scaled_odgp_with_parity_k_ne_0(int_a, int_b, 1, &parity);
        let sol_x = sol_t.map(|alpha| alpha * dx.clone() + alpha0.clone());
        for alpha in sol_x {
            sol_sufficient.push(DOmega::from_droottwo_vector(&alpha, &beta, k));
        }
    }

    let op_g_inv = op_g.inv().unwrap();

    sol_sufficient
        .into_iter()
        .map(|z| op_g_inv.clone() * z)
        .filter(|z| set_a.inside(z) && set_b.inside(z.conj_sq2()))
        .collect()
}
