//use dashu_float::{round::mode::HalfEven, FBig};

use rsgridsynth::config::config_from_theta_epsilon;
use rsgridsynth::gridsynth::gridsynth_gates;
// use rsgridsynth::common::{PI, get_prec_bits};
// fn to_fbig(x: f64) -> FBig<HalfEven> {
//     FBig::<HalfEven>::try_from(x)
//         .unwrap()
//         .with_precision(get_prec_bits())
//         .value()
// }

fn main() {
    // let theta = PI.clone() / to_fbig(16.0);
    // let epsilon = to_fbig(1e-10);
    let pi = std::f64::consts::PI;
    let theta = pi / 8.0;
    let epsilon = 1e-10;
    let seed = 1234;
    // This flag is currently ignored
    let up_to_phase = false;
    let mut gridsynth_config = config_from_theta_epsilon(theta, epsilon, seed, up_to_phase);
    let gates = gridsynth_gates(&mut gridsynth_config);
    println!("{}", gates);
}
