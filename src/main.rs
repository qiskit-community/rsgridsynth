// Copyright (c) 2025 IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use clap::{Arg, Command};
use dashu_int::UBig;
use log::info;
use rand::{rngs::StdRng, SeedableRng};

pub mod common;
pub mod config;
pub mod diophantine;
pub mod grid_op;
pub mod gridsynth;
pub mod math;
pub mod normal_form;
pub mod odgp;
pub mod region;
pub mod ring;
pub mod synthesis_of_clifford_t;
pub mod tdgp;
pub mod to_upright;
pub mod unitary;
use std::{f32::consts::LOG2_10, time::Instant};

use crate::common::{ib_to_bf_prec, set_prec_bits};
use crate::config::parse_decimal_with_exponent;
use crate::config::{DiophantineData, GridSynthConfig};
use gridsynth::gridsynth_gates;

fn main() {
    let matches = build_command().get_matches();

    let verbose = matches.get_flag("verbose");
    let time = matches.get_flag("time");

    if verbose || time {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    } else {
        env_logger::init();
    }

    let mut args = parse_arguments(&matches);

    let start = if args.measure_time {
        Some(Instant::now())
    } else {
        None
    };

    let res = gridsynth_gates(&mut args);

    if let Some(start_time) = start {
        let elapsed = start_time.elapsed();
        info!("Elapsed time: {:.3?}", elapsed);
    }

    println!("{}", res.gates);
}

fn build_command() -> Command {
    Command::new("rsgridsynth")
        .arg(Arg::new("theta").required(true))
        .arg(Arg::new("epsilon").required(true))
        .arg(Arg::new("dps").long("dps").default_value(None))
        .arg(Arg::new("seed").long("seed").short('s').default_value("1"))
        .arg(
            Arg::new("dtimeout")
                .long("dtimeout")
                .short('d')
                .default_value("200"),
        )
        .arg(
            Arg::new("ftimeout")
                .long("ftimeout")
                .short('f')
                .default_value("50"),
        )
        .arg(
            Arg::new("verbose")
                .long("verbose")
                .short('v')
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("time")
                .long("time")
                .short('t')
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("showgraph")
                .long("showgraph")
                .short('g')
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("check")
                .long("check")
                .action(clap::ArgAction::SetTrue),
        )
}

fn parse_arguments(matches: &clap::ArgMatches) -> GridSynthConfig {
    let theta_str = matches.get_one::<String>("theta").unwrap();
    let (theta_num, theta_den) = parse_decimal_with_exponent(theta_str).unwrap();
    let theta = ib_to_bf_prec(theta_num) / ib_to_bf_prec(theta_den);
    let epsilon_str = matches.get_one::<String>("epsilon").unwrap();
    let (epsilon_num, epsilon_den) = parse_decimal_with_exponent(epsilon_str).unwrap();

    let dps: Option<u32> = matches
        .get_one::<String>("dps")
        .and_then(|s| s.parse().ok());
    // The magic number 12 safely overapproximates the bits of precision.
    let calculated_prec_bits =
        12 * (epsilon_den.ilog(&UBig::from(10u8)) - epsilon_num.ilog(&UBig::from(10u8)));
    let prec_bits: usize = if let Some(dps_val) = dps {
        (dps_val as f64 * LOG2_10 as f64) as usize
    } else {
        calculated_prec_bits
    };
    set_prec_bits(prec_bits);
    let epsilon = ib_to_bf_prec(epsilon_num) / ib_to_bf_prec(epsilon_den);
    let diophantine_timeout = matches
        .get_one::<String>("dtimeout")
        .unwrap()
        .parse()
        .unwrap();
    let factoring_timeout = matches
        .get_one::<String>("ftimeout")
        .unwrap()
        .parse()
        .unwrap();
    let verbose = matches.get_flag("verbose");
    let measure_time = matches.get_flag("time");
    let check_solution = matches.get_flag("check");

    let seed = matches.get_one::<String>("seed").unwrap().parse().unwrap();
    let rng: StdRng = SeedableRng::seed_from_u64(seed);
    let diophantine_data = DiophantineData {
        diophantine_timeout,
        factoring_timeout,
        rng,
    };

    GridSynthConfig {
        theta,
        epsilon,
        verbose,
        measure_time,
        diophantine_data,
        check_solution,
    }
}
