use dashu_float::round::mode::HalfEven;
use dashu_float::FBig;
use rsgridsynth::accuracy::AchievedDiamondError;
use rsgridsynth::clear_caches;
use rsgridsynth::config::config_from_theta_epsilon;
use rsgridsynth::gridsynth::gridsynth_gates;
use serial_test::serial;

fn fbig_to_f64(x: &FBig<HalfEven>) -> f64 {
    match x.to_f64() {
        dashu_base::Approximation::Exact(v) => v,
        dashu_base::Approximation::Inexact(v, _) => v,
    }
}

#[test]
#[serial]
fn simple_test() {
    let pi = std::f64::consts::PI;
    let theta = pi / 8.0; // ≈ 0.39269908169872414
    let epsilon = 1e-10;

    let gates_1234 = "HTHTSHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTSHTHTSHTHTSHTXSSWWW";

    let gates_101 = "HTSHTSHTSHTHTHTSHTHTHTSHTSHTHTHTHTHTHTSHTHTHTHTSHTSHTHTSHTHTHTHTHTHTHTHTSHTHTHTHTSHTHTSHTSHTSHTSHTHTSHTSHTHTSHTSHTHTSHTHTHTSHTSHTHTHTHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTHTSHTHTSHTSHTHTHTHTSHTHTHTSHTHTHTHTSHTHTSHTSHTHTSHTHTSHTHTHTHTHTHTHTHTSHTHTHTSHTSSSWW";

    let gates_1 = "HTSHTHTHTSHTHTHTHTHTHTHTHTSHTHTSHTHTSHTSHTHTSHTHTHTHTSHTHTHTSHTHTHTHTSHTSHTHTSHTHTHTHTHTSHTSHTHTHTSHTHTSHTHTSHTHTHTHTSHTSHTHTHTSHTHTSHTSHTHTSHTSHTHTSHTSHTSHTSHTHTSHTHTHTHTSHTHTHTHTHTHTHTHTSHTHTSHTSHTHTHTHTSHTHTHTHTHTHTSHTSHTHTHTSHTHTHTSHTSHTSHTSSSWW";

    let test_inputs = vec![(1234, gates_1234), (101, gates_101), (1, gates_1)];

    let verbose = false;
    let up_to_phase = false;
    for (seed, expected_gates) in test_inputs {
        clear_caches();
        let mut gridsynth_config =
            config_from_theta_epsilon(theta, epsilon, seed, verbose, up_to_phase);
        let gates = gridsynth_gates(&mut gridsynth_config).gates;
        assert_eq!(
            gates.to_string(),
            expected_gates,
            "Test failed for seed: {}",
            seed
        );
    }
}

#[test]
#[serial]
fn pi_over_two_test() {
    let pi = std::f64::consts::PI;
    let theta = pi / 2.0;

    let epsilons = vec![1e-2, 1e-3, 1e-10];

    let verbose = false;
    let up_to_phase = false;
    for epsilon in epsilons {
        let seeds = 10..50;
        for seed in seeds {
            clear_caches();
            let mut gridsynth_config =
                config_from_theta_epsilon(theta, epsilon, seed, verbose, up_to_phase);
            let gates = gridsynth_gates(&mut gridsynth_config).gates;
            let expected_gates = "SWWWWWWW";
            assert_eq!(gates.to_string(), expected_gates);
        }
    }
}

#[test]
#[serial]
fn pi_over_4_exact_test() {
    let pi = std::f64::consts::PI;
    let theta = pi / 4.0;
    let epsilon = 1e-10;
    let seed = 1234;
    let up_to_phase = false;
    let verbose = false;
    clear_caches();
    let mut gridsynth_config =
        config_from_theta_epsilon(theta, epsilon, seed, verbose, up_to_phase);
    let gates = gridsynth_gates(&mut gridsynth_config).gates;
    let expected_gates = "SHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTHTSHTSHTSHTHTHTSHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTSHTHTHTHTHTHTSHTSHTHTHTSHTSHTHTHTHTSHTSHTHTSHTSHTSHTHTSHTHTHTSHTSHTHTSHTSHTSHTHTSHTHTHTSHTHTHTSHTSHTSHTSHTHTHTHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTHTHTHTSHTSHTHTSHTSHTSHTHTSHTSHTSHTHTSHTHTHTSHTHTHTSHTSHTHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTHTHTSHTSHTHTSHSWW";

    assert_eq!(gates.to_string(), expected_gates);
}

#[test]
#[serial]
fn pi_over_4_with_phase_test() {
    let pi = std::f64::consts::PI;
    let theta = pi / 4.0;
    let epsilon = 1e-10;
    let seed = 1234;
    let up_to_phase = true;
    let verbose = false;
    clear_caches();
    let mut gridsynth_config =
        config_from_theta_epsilon(theta, epsilon, seed, verbose, up_to_phase);
    let gates = gridsynth_gates(&mut gridsynth_config).gates;
    let expected_gates = "TWWWWWWW";
    assert_eq!(gates.to_string(), expected_gates);
}

/// Pins the user-visible contract that a synthesized result's printed (`Display`) form is
/// always re-readable via `GateSeq::from_str` -- the CLI's whole output contract rests on this.
#[test]
#[serial]
fn public_api_gates_round_trip_through_display_and_parse() {
    let theta = std::f64::consts::PI / 8.0;
    let epsilon = 1e-10;
    clear_caches();
    let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, 1234, false, false);
    let res = gridsynth_gates(&mut gridsynth_config);

    let round_tripped: rsgridsynth::gate::GateSeq = res.gates.to_string().parse().unwrap();
    assert_eq!(round_tripped, res.gates);
}

#[test]
#[serial]
fn test_correct_decomposition_exact() {
    let epsilon = 1e-10;

    let verbose = false;
    let seed = 0u64;
    let up_to_phase = false;

    let thetas = (0..64).map(|k| k as f64 * std::f64::consts::PI / 32.0);

    for theta in thetas {
        clear_caches();
        let mut gridsynth_config =
            config_from_theta_epsilon(theta, epsilon, seed, verbose, up_to_phase);

        let res = gridsynth_gates(&mut gridsynth_config);
        let error = res.achieved_diamond_error(&gridsynth_config.theta);
        let is_correct = fbig_to_f64(&error) < fbig_to_f64(&gridsynth_config.epsilon) * 2.0;

        // not printed, unless cargo test is run with -- -no-capture
        println!(
            "theta = {theta}, gates = {}, error = {:.6e}, correct = {is_correct:?}",
            res.gates,
            fbig_to_f64(&error),
        );

        // Check that the diamond-norm error is within the requested (doubled, per the
        // diamond-vs-operator-norm convention) budget.
        assert!(is_correct);
    }
}

#[test]
#[serial]
fn test_correct_decomposition_up_to_phase() {
    let epsilon = 1e-10;

    let verbose = false;
    let seed = 0u64;
    let up_to_phase = true;

    let thetas = (0..64).map(|k| k as f64 * std::f64::consts::PI / 32.0);

    for theta in thetas {
        clear_caches();
        let mut gridsynth_config =
            config_from_theta_epsilon(theta, epsilon, seed, verbose, up_to_phase);

        let res = gridsynth_gates(&mut gridsynth_config);
        let error = res.achieved_diamond_error(&gridsynth_config.theta);
        let is_correct = fbig_to_f64(&error) < fbig_to_f64(&gridsynth_config.epsilon) * 2.0;

        // not printed, unless cargo test is run with -- -no-capture
        println!(
            "theta = {theta}, gates = {}, error = {:.6e}, correct = {is_correct:?}",
            res.gates,
            fbig_to_f64(&error),
        );
        // Check that the diamond-norm error is within the requested (doubled, per the
        // diamond-vs-operator-norm convention) budget.
        assert!(is_correct);
    }
}

#[test]
#[serial]
fn test_low_precision_bug() {
    let pi = std::f64::consts::PI;
    let theta = pi / 2.0;
    let epsilon = 1e-1;
    let verbose = false;
    let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, 1234, verbose, false);
    gridsynth_gates(&mut gridsynth_config);
}

#[test]
#[serial]
fn test_shared_cache_across_denomexp_no_panic() {
    let epsilon = 1e-3;
    let verbose = false;
    let seed = 0u64;
    let up_to_phase = true;
    let mut panic_thetas = Vec::new();
    let mut incorrect_thetas = Vec::new();
    let mut max_error = 0.0_f64;

    clear_caches();
    for i in 1..=320 {
        let theta = i as f64 * 0.01;
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let mut gridsynth_config =
                config_from_theta_epsilon(theta, epsilon, seed, verbose, up_to_phase);
            let res = gridsynth_gates(&mut gridsynth_config);
            let error = fbig_to_f64(&res.achieved_diamond_error(&gridsynth_config.theta));
            let is_correct = error < fbig_to_f64(&gridsynth_config.epsilon) * 2.0;
            (error, is_correct)
        }));

        match result {
            Ok((error, is_correct)) => {
                max_error = max_error.max(error);
                if !is_correct {
                    incorrect_thetas.push(theta);
                }
            }
            Err(_) => panic_thetas.push(theta),
        }
    }

    println!(
        "shared-cache scan: n=320 epsilon={epsilon} panics={} incorrect={} max_error={max_error:.6e}",
        panic_thetas.len(),
        incorrect_thetas.len()
    );
    assert!(panic_thetas.is_empty(), "panic_thetas={panic_thetas:?}");
    assert!(
        incorrect_thetas.is_empty(),
        "incorrect_thetas={incorrect_thetas:?}"
    );
}

#[test]
#[serial]
fn test_timeouts_preserved_after_synthesis() {
    let mut gridsynth_config =
        config_from_theta_epsilon(std::f64::consts::PI / 8.0, 1e-10, 1234, false, true);
    gridsynth_config.diophantine_data.diophantine_timeout = 237;
    gridsynth_config.diophantine_data.factoring_timeout = 61;

    clear_caches();
    let result = gridsynth_gates(&mut gridsynth_config);

    assert!(!result.gates.is_empty());
    assert_eq!(gridsynth_config.diophantine_data.diophantine_timeout, 237);
    assert_eq!(gridsynth_config.diophantine_data.factoring_timeout, 61);
}
