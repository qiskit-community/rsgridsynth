use rsgridsynth::config::config_from_theta_epsilon;
use rsgridsynth::gridsynth::gridsynth_gates;
//use rsgridsynth::common::{PI, get_prec_bits};

#[test]
fn simple_test() {
    // let theta = PI.clone() / to_fbig(16.0);
    // let epsilon = to_fbig(1e-10);
    let pi = std::f64::consts::PI;
    let theta = pi / 8.0;
    let epsilon = 1e-10;
    let seed = 1234;
    let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, seed);
    let gates = gridsynth_gates(&mut gridsynth_config);
    let expected_gates = "HTHTSHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTHTHTHTHTHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTSHTSHTSHTSHTSHTHTHTHTHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTSHTHTSHTHTSHTXSSWWW";
    assert_eq!(gates, expected_gates);
}
