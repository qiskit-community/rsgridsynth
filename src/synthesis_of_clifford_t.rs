// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use dashu_int::IBig;

use crate::normal_form::NormalForm;
use crate::ring::ZOmega;
use crate::unitary::DOmegaUnitary;

const BIT_SHIFT: [i32; 16] = [0, 0, 1, 0, 2, 0, 1, 3, 3, 3, 0, 2, 2, 1, 0, 0];
const BIT_COUNT: [i32; 16] = [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4];

fn reduce_denomexp(mut unitary: DOmegaUnitary) -> (String, DOmegaUnitary) {
    let t_power_and_h = ["H", "TH", "SH", "TSH"];
    let residue_z = unitary.z.residue();
    let residue_w = unitary.w.residue();
    let residue_squared_z = (unitary.z.u.clone() * unitary.z.conj().u.clone()).residue();

    let mut m = BIT_SHIFT[residue_w as usize] - BIT_SHIFT[residue_z as usize];
    if m < 0 {
        m += 4;
    }
    if residue_squared_z == 0b0000 {
        unitary = unitary
            .mul_by_h_and_t_power_from_left(0)
            .renew_denomexp(unitary.k() - 1);
        (t_power_and_h[0].to_string(), unitary)
    } else if residue_squared_z == 0b1010 {
        unitary = unitary
            .mul_by_h_and_t_power_from_left(-m)
            .renew_denomexp(unitary.k() - 1);
        (t_power_and_h[m as usize].to_string(), unitary)
    } else if residue_squared_z == 0b0001 {
        if BIT_COUNT[residue_z as usize] == BIT_COUNT[residue_w as usize] {
            unitary = unitary
                .mul_by_h_and_t_power_from_left(-m)
                .renew_denomexp(unitary.k() - 1);
            (t_power_and_h[m as usize].to_string(), unitary)
        } else {
            unitary = unitary.mul_by_h_and_t_power_from_left(-m);
            (t_power_and_h[m as usize].to_string(), unitary)
        }
    } else {
        panic!("Invalid residue");
    }
}

pub fn decompose_domega_unitary(mut unitary: DOmegaUnitary) -> String {
    let omega_power: [ZOmega; 8] = [
        ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::ONE),
        ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ONE, IBig::ZERO),
        ZOmega::new(IBig::ZERO, IBig::ONE, IBig::ZERO, IBig::ZERO),
        ZOmega::new(IBig::ONE, IBig::ZERO, IBig::ZERO, IBig::ZERO),
        ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::ZERO, IBig::NEG_ONE),
        ZOmega::new(IBig::ZERO, IBig::ZERO, IBig::NEG_ONE, IBig::ZERO),
        ZOmega::new(IBig::ZERO, IBig::NEG_ONE, IBig::ZERO, IBig::ZERO),
        ZOmega::new(IBig::NEG_ONE, IBig::ZERO, IBig::ZERO, IBig::ZERO),
    ];
    let mut gates = String::new();
    while unitary.k() > 0 {
        let (g, next_unitary) = reduce_denomexp(unitary);
        gates += &g;
        unitary = next_unitary;
    }
    if unitary.n & 1 != 0 {
        gates += "T";
        unitary = unitary.mul_by_t_inv_from_left();
    }

    if unitary.z == IBig::ZERO {
        gates += "X";
        unitary = unitary.mul_by_x_from_left();
    }

    let mut m_w = 0;
    for (m, omega) in omega_power.iter().enumerate() {
        if unitary.z == *omega {
            m_w = m;
            unitary = unitary.mul_by_w_power_from_left(-(m as i32));
            break;
        }
    }

    let m_s = (unitary.n >> 1) as usize;
    gates += &"S".repeat(m_s);
    unitary = unitary.mul_by_s_power_from_left(-(m_s as i32));
    gates += &"W".repeat(m_w);
    assert_eq!(
        unitary,
        DOmegaUnitary::identity(),
        "decomposition failed..."
    );
    NormalForm::from_gates(&gates).to_gates()
}
