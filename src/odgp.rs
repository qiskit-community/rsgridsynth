// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::common::ib_to_bf_prec;
use crate::math::{floorlog, pow_sqrt2, sqrt2};
use crate::region::Interval;
use crate::ring::z_root_two::LAMBDA;
use crate::ring::{DRootTwo, ZRootTwo};
use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use dashu_int::IBig;
use std::iter;

pub fn solve_odgp(i: Interval, j: Interval) -> impl Iterator<Item = ZRootTwo> {
    // Can't return two different iterator types. So we can't do this check.
    // I checked with dbg! to confirm that omitting this check is ok:
    // If the condition is true, then an empty iterator is returned at the end.
    //    if i.width() < ib_to_bf_prec(IBig::ZERO) || j.width() < ib_to_bf_prec(IBig::ZERO) {
    //        return vec![].into_iter();
    //    }

    let sum = &i.l + &j.l;
    let div_result1: FBig<HalfEven> = sum / 2;
    let a = div_result1.floor().try_into().unwrap();

    let diff = &i.l - &j.l;
    let div_result2: FBig<HalfEven> = (sqrt2() * diff) / 4;
    let b = div_result2.floor().try_into().unwrap();
    let alpha = ZRootTwo::new(a, b);
    let sub_i = i.sub_ref(&alpha.to_real());
    let sub_j = j.sub_ref(&alpha.conj_sq2().to_real());
    let sol = solve_odgp_internal(sub_i, sub_j);
    sol.into_iter()
        .map(move |beta| &beta + &alpha)
        .filter(move |beta| {
            let real = beta.to_real();
            let real_conj = beta.conj_sq2().to_real();
            i.within(&real) && j.within(&real_conj)
        })
}

fn solve_odgp_internal(i: Interval, j: Interval) -> Box<dyn Iterator<Item = ZRootTwo>> {
    let bfzero = ib_to_bf_prec(IBig::ZERO);
    if i.width() < bfzero || j.width() < bfzero {
        return Box::new(vec![].into_iter());
    } else if i.width() > bfzero && j.width() <= bfzero {
        return Box::new(solve_odgp_internal(j, i).map(|beta| beta.conj_sq2()));
    }

    // Compute the scaling factor LAMBDA^n
    let n = if j.width() <= bfzero {
        IBig::ZERO
    } else {
        floorlog(j.width(), LAMBDA.to_real()).0
    };

    let lambda_n = LAMBDA.pow(&n);
    let lambda_inv_n = LAMBDA.pow(&(-&n));
    let lambda_n_f = lambda_n.to_real();

    let lambda_conj_sq2_n = LAMBDA.conj_sq2().pow(&n);
    let lambda_conj_sq2_n_f = lambda_conj_sq2_n.to_real();

    // Here we replace the intervals (i, j) by their scaled versions
    // (avoiding an extra recursion).
    let i = i.scale(&lambda_n_f);
    let j = j.scale(&lambda_conj_sq2_n_f);

    let sum_min = &i.l + &j.l;
    let div_min: FBig<HalfEven> = sum_min / 2;
    let a_min: IBig = div_min.ceil().try_into().unwrap();

    let sum_max = &i.r + &j.r;
    let div_max: FBig<HalfEven> = sum_max / 2;
    let a_max = div_max.floor().try_into().unwrap();

    let sol_iter = iter::successors(Some(a_min.clone()), move |a| {
        let mut a_next = a.clone();
        a_next += 1;
        if a_next <= a_max {
            Some(a_next)
        } else {
            None
        }
    })
    .flat_map(move |a| {
        let a_real = ib_to_bf_prec(a.clone()); // 明示的に clone して消費
        let tmp1: FBig<HalfEven> = sqrt2() * (&a_real - &j.r) / 2;
        let b_min: IBig = tmp1.ceil().try_into().unwrap();

        let tmp2: FBig<HalfEven> = sqrt2() * (&a_real - &j.l) / 2;
        let b_max: IBig = tmp2.floor().try_into().unwrap();

        iter::successors(Some(b_min.clone()), move |b| {
            let mut b_next = b.clone();
            b_next += 1;
            if b_next <= b_max {
                Some(b_next)
            } else {
                None
            }
        })
        .map({
            let value = lambda_inv_n.clone();
            move |b| ZRootTwo::new(a.clone(), b.clone()) * &value
        })
    });

    Box::new(sol_iter)
}

pub fn solve_odgp_with_parity(
    i: Interval,
    j: Interval,
    beta: &DRootTwo,
) -> impl Iterator<Item = ZRootTwo> {
    let p = beta.parity();
    let scale_factor1 = sqrt2() / 2;
    let scale_factor2 = -sqrt2() / 2;
    let scaled_i = (i - p.clone()).scale(&scale_factor1);
    let scaled_j = (j - p.clone()).scale(&scale_factor2);
    let sol = solve_odgp(scaled_i, scaled_j);

    sol.into_iter()
        .map(move |alpha| (alpha * ZRootTwo::new(IBig::ZERO, IBig::ONE)) + &p)
}

pub fn first_solve_scaled_odgp(i: &Interval, j: &Interval, k: i64) -> Option<DRootTwo> {
    solve_scaled_odgp(i, j, k).next()
}

pub fn solve_scaled_odgp(i: &Interval, j: &Interval, k: i64) -> impl Iterator<Item = DRootTwo> {
    let scale = pow_sqrt2(k);
    let neg_scale = -scale.clone();
    let scaled_j = if k & 1 == 0 {
        j.scale(&scale)
    } else {
        j.scale(&neg_scale)
    };
    solve_odgp(i.scale(&scale), scaled_j).map(move |alpha| DRootTwo::new(alpha, k))
}

pub fn solve_scaled_odgp_with_parity_k_ne_0(
    i: Interval,
    j: Interval,
    k: i64,
    beta: &DRootTwo,
) -> impl Iterator<Item = DRootTwo> {
    // Can't do this because the iterators are not of the same type.
    // But this function is only called with k == 1. So we don't need the k == 0 branch.
    // if k == 0 {
    //     let base = beta.renew_denomexp(0);
    //     return solve_odgp_with_parity(i, j, &base)
    //          .map(DRootTwo::from_zroottwo);
    //  }

    let p = beta.renew_denomexp(k).parity();
    let offset = if p == IBig::ZERO {
        DRootTwo::from_int(IBig::ZERO)
    } else {
        DRootTwo::power_of_inv_sqrt2(k)
    };

    let sub_i = i - offset.to_real();
    let sub_j = j - offset.conj_sq2().to_real();
    let sol = solve_scaled_odgp(&sub_i, &sub_j, k - 1);
    sol.map(move |a| a + offset.clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_empty_interval() -> (Interval, Interval) {
        let inti = Interval::new(FBig::<HalfEven>::from(4), FBig::<HalfEven>::from(2));
        let intj = Interval::new(FBig::<HalfEven>::from(2), FBig::<HalfEven>::from(4));
        (inti, intj)
    }

    #[test]
    fn test_empty_interval() {
        let (inti, intj) = create_empty_interval();
        let mut result = solve_odgp(inti, intj);
        assert!(result.next().is_none());
    }

    #[test]
    fn test_use_empty_interval() {
        let (inti, intj) = create_empty_interval();
        let mut result = solve_scaled_odgp(&inti, &intj, 2);
        assert!(result.next().is_none());
    }
}
