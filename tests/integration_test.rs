use ntest::timeout;
use rsgridsynth::clear_caches;
use rsgridsynth::config::config_from_theta_epsilon;
use rsgridsynth::gridsynth::gridsynth_gates;
use serial_test::serial;

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
        assert_eq!(gates, expected_gates, "Test failed for seed: {}", seed);
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
            assert_eq!(gates, expected_gates);
        }
    }
}

#[timeout(1000)]
#[test]
#[serial]
fn slowness_from_allocation_test() {
    let pi = std::f64::consts::PI;
    let theta = pi / 4.0;
    let epsilon = 1e-10;
    let seed = 1234;
    let up_to_phase = false;
    let verbose = false;
    let mut gridsynth_config =
        config_from_theta_epsilon(theta, epsilon, seed, verbose, up_to_phase);
    let gates = gridsynth_gates(&mut gridsynth_config).gates;
    let expected_gates = "SHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTHTSHTSHTSHTHTHTSHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTSHTHTHTHTHTHTSHTSHTHTHTSHTSHTHTHTHTSHTSHTHTSHTSHTSHTHTSHTHTHTSHTSHTHTSHTSHTSHTHTSHTHTHTSHTHTHTSHTSHTSHTSHTHTHTHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTHTHTHTSHTSHTHTSHTSHTSHTHTSHTSHTSHTHTSHTHTHTSHTHTHTSHTSHTHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTHTHTSHTSHTHTSHSWW";

    assert_eq!(gates, expected_gates);
}

#[test]
#[serial]
fn test_correct_decomposition() {
    let epsilon = 1e-8;

    let verbose = false;
    let seed = 0u64;
    let up_to_phase = false;

    let thetas = (0..32).map(|k| k as f64 * std::f64::consts::PI / 16.0);

    for theta in thetas {
        clear_caches();
        let mut gridsynth_config =
            config_from_theta_epsilon(theta, epsilon, seed, verbose, up_to_phase)
                .with_check_solution(true);

        let res = gridsynth_gates(&mut gridsynth_config);

        // not printed, unless cargo test is run with -- -no-capture
        println!(
            "theta = {theta}, gates = {}, correct = {:?}",
            res.gates, res.is_correct
        );

        // Check that the check result exits and is valid.
        assert!(res.is_correct.is_some_and(|v| v));
    }
}
