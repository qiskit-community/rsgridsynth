// Copyright (c) 2025 IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use clap::{Arg, Command};
use dashu_float::{round::mode::HalfEven, FBig};
use dashu_int::{IBig, UBig};
use log::info;

pub mod common;
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
use std::{f32::consts::LOG2_10, str::FromStr, time::Instant};

use gridsynth::gridsynth_gates;

use crate::common::{ib_to_bf_prec, set_prec_bits};

pub fn parse_decimal_with_exponent(input: &str) -> Option<(IBig, IBig)> {
    let input = input.trim();
    let (sign, body) = if let Some(s) = input.strip_prefix('-') {
        (-1, s)
    } else if let Some(s) = input.strip_prefix('+') {
        (1, s)
    } else {
        (1, input)
    };

    let (base_str, exp_str) = match body.split_once(['e', 'E']) {
        Some((b, e)) => (b, e),
        None => (body, "0"),
    };

    let parts: Vec<&str> = base_str.split('.').collect();
    if parts.len() > 2 {
        return None;
    }
    let int_part = parts[0];
    let frac_part = if parts.len() == 2 { parts[1] } else { "" };
    let digits = format!("{}{}", int_part, frac_part);
    let decimal_digits = frac_part.len() as i32;

    let exponent: i32 = exp_str.parse().ok()?;
    let scale = exponent - decimal_digits;

    let mut numerator = IBig::from_str(&digits).ok()? * sign;
    let mut denominator = IBig::from(1);

    match scale.cmp(&0) {
        std::cmp::Ordering::Greater => {
            numerator *= IBig::from(10u8).pow(scale as usize);
        }
        std::cmp::Ordering::Less => {
            denominator = IBig::from(10u8).pow((-scale) as usize);
        }
        std::cmp::Ordering::Equal => {}
    }

    Some((numerator, denominator))
}

fn main() {
    let matches = build_command().get_matches();

    let verbose = matches.get_flag("verbose");
    let time = matches.get_flag("time");

    if verbose || time {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    } else {
        env_logger::init();
    }

    let args = parse_arguments(&matches);

    let start = if args.time {
        Some(Instant::now())
    } else {
        None
    };

    let gates = gridsynth_gates(
        args.theta,
        args.epsilon,
        args.dtimeout,
        args.factoring_timeout,
        args.verbose,
        args.time,
    );

    if let Some(start_time) = start {
        let elapsed = start_time.elapsed();
        info!("Elapsed time: {:.3?}", elapsed);
    }

    println!("{}", gates);
}

fn build_command() -> Command {
    Command::new("rsgridsynth")
        .arg(Arg::new("theta").required(true))
        .arg(Arg::new("epsilon").required(true))
        .arg(Arg::new("dps").long("dps").default_value(None))
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
}

struct Args {
    theta: FBig<HalfEven>,
    epsilon: FBig<HalfEven>,
    dtimeout: u128,
    factoring_timeout: u128,
    verbose: bool,
    time: bool,
}

fn parse_arguments(matches: &clap::ArgMatches) -> Args {
    let theta_str = matches.get_one::<String>("theta").unwrap();
    let (theta_num, theta_den) = parse_decimal_with_exponent(theta_str).unwrap();
    let theta = ib_to_bf_prec(theta_num) / ib_to_bf_prec(theta_den);
    let epsilon_str = matches.get_one::<String>("epsilon").unwrap();
    let (epsilon_num, epsilon_den) = parse_decimal_with_exponent(epsilon_str).unwrap();

    let dps: Option<u32> = matches
        .get_one::<String>("dps")
        .and_then(|s| s.parse().ok());
    let calculated_prec_bits =
        12 * (epsilon_den.ilog(&UBig::from(10u8)) - epsilon_num.ilog(&UBig::from(10u8)));
    let prec_bits: usize = if let Some(dps_val) = dps {
        (dps_val as f64 * LOG2_10 as f64) as usize
    } else {
        calculated_prec_bits
    };
    set_prec_bits(prec_bits);
    let epsilon = ib_to_bf_prec(epsilon_num) / ib_to_bf_prec(epsilon_den);
    let dtimeout = matches
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
    let time = matches.get_flag("time");

    Args {
        theta,
        epsilon,
        dtimeout,
        factoring_timeout,
        verbose,
        time,
    }
}
