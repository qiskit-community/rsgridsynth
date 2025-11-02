use rsgridsynth::config::config_from_theta_epsilon;
use rsgridsynth::gridsynth::gridsynth_gates;

#[test]
#[ignore]
fn simple_test() {
    let pi = std::f64::consts::PI;
    let theta = pi / 8.0; // ≈ 0.39269908169872414
    let epsilon = 1e-10;
    let seed = 1234;
    let verbose = false;
    let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, seed, verbose);
    let gates = gridsynth_gates(&mut gridsynth_config);
    let expected_gates = "HTHTSHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTSHTHTSHTHTSHTXSSWWW";
    assert_eq!(gates, expected_gates);
}

#[test]
#[ignore]
fn simple_test2() {
    let pi = std::f64::consts::PI;
    let theta = pi / 8.0; // ≈ 0.39269908169872414
    let epsilon = 1e-10;

    let gates1 = "HTHTSHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTSHTHTSHTHTSHTXSSWWW";

    let test_inputs = vec![(1234, gates1), (101, gates1), (1, gates1)];

    let verbose = false;
    for (seed, expected_gates) in test_inputs {
        let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, seed, verbose);
        let gates = gridsynth_gates(&mut gridsynth_config);
        assert_eq!(gates, expected_gates, "Test failed for seed: {}", seed);
    }
}

//#[ignore]
#[test]
fn pi_over_two_test() {
    let pi = std::f64::consts::PI;
    let theta = pi / 2.0;

    let epsilons = vec![1e-2, 1e-3, 1e-10];

    let verbose = false;
    for epsilon in epsilons {
        let seeds = 100..300;
        for seed in seeds {
            let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, seed, verbose);
            let gates = gridsynth_gates(&mut gridsynth_config);
            let expected_gates = "SWWWWWWW";
            assert_eq!(gates, expected_gates);
        }
    }
}
