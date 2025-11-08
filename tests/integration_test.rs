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
    for (seed, expected_gates) in test_inputs {
        clear_caches();
        let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, seed, verbose);
        let gates = gridsynth_gates(&mut gridsynth_config);
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
    for epsilon in epsilons {
        let seeds = 10..50;
        for seed in seeds {
            clear_caches();
            let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, seed, verbose);
            let gates = gridsynth_gates(&mut gridsynth_config);
            let expected_gates = "SWWWWWWW";
            assert_eq!(gates, expected_gates);
        }
    }
}

#[timeout(1000)]
#[test]
#[should_panic(expected = "timeout")]
#[serial]
fn slowness_from_allocation_test() {
    let pi = std::f64::consts::PI;
    let theta = pi / 4.0;
    let epsilon = 1e-10;
    let seed = 1234;
    let verbose = false;
    let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, seed, verbose);
    let gates = gridsynth_gates(&mut gridsynth_config);
    let expected_gates = "SHTHTSHTSHTHTSHTHTSHTSHTSHTSHTSHTHTHTHTSHTHTSHTHTHTHTSHTSHTHTSHTSHTHTSHTSHTSHTSHTSHTSHTSHTSHTHTSHTHTSHTSHTHTHTHTHTHTHTSHTSHTSHTSHTHTHTHTHTHTHTHTSHTSHTSHTHTSHTHTSHTHTHTHTSHTHTHTSHTHTHTHTHTSHTSHTHTSHTSHTSHTSHTSHTHTSHTHTSHTSHTHTSHTHTSHTHTSHTSHTSHTSHTHTSHTHTSHTSHTHTHTSHTHTSHTHTSHTHTHTSHTHTHTHTHTSHTHTHTHTSHTSHTSHTSHTHTSHTSHTHTHTHTSHTHTHTHTHTHTHXSWWW";
    assert_eq!(gates, expected_gates);
}
