// Regression test for the SHT-syllable merge formula in NormalForm::append_gate
// (src/normal_form.rs). The merge branch for `Syllable::SHT` must right-multiply
// the running Clifford part by S*H*S (matching pygridsynth's
// `CLIFFORD_SHS = CLIFFORD_S * CLIFFORD_H * CLIFFORD_S`), not H*S*H -- these are
// different elements of the single-qubit Clifford group (H and S do not commute).
//
// The wrong formula (H*S*H) leaves NormalForm::from_gates("SHTT").to_gates()
// unitarily different from the input, while reporting the SAME T-count as the
// correct formula would -- a T-count-only check cannot catch this.

use rsgridsynth::normal_form::NormalForm;

type C64 = (f64, f64);
type Mat2 = [[C64; 2]; 2];

fn c_mul(a: C64, b: C64) -> C64 {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}
fn c_add(a: C64, b: C64) -> C64 {
    (a.0 + b.0, a.1 + b.1)
}

fn mat_mul(a: &Mat2, b: &Mat2) -> Mat2 {
    let mut out = [[(0.0, 0.0); 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            out[i][j] = c_add(c_mul(a[i][0], b[0][j]), c_mul(a[i][1], b[1][j]));
        }
    }
    out
}

fn gate_matrix(g: char) -> Mat2 {
    const R2: f64 = std::f64::consts::FRAC_1_SQRT_2;
    const W_RE: f64 = std::f64::consts::FRAC_1_SQRT_2; // cos(pi/4)
    const W_IM: f64 = std::f64::consts::FRAC_1_SQRT_2; // sin(pi/4)
    match g {
        'I' => [[(1.0, 0.0), (0.0, 0.0)], [(0.0, 0.0), (1.0, 0.0)]],
        'H' => [[(R2, 0.0), (R2, 0.0)], [(R2, 0.0), (-R2, 0.0)]],
        'S' => [[(1.0, 0.0), (0.0, 0.0)], [(0.0, 0.0), (0.0, 1.0)]],
        'T' => [[(1.0, 0.0), (0.0, 0.0)], [(0.0, 0.0), (R2, R2)]],
        'X' => [[(0.0, 0.0), (1.0, 0.0)], [(1.0, 0.0), (0.0, 0.0)]],
        'W' => [[(W_RE, W_IM), (0.0, 0.0)], [(0.0, 0.0), (W_RE, W_IM)]],
        other => panic!("unexpected gate char '{}'", other),
    }
}

/// Left-to-right gate string -> unitary, applying gates in the SAME temporal
/// order convention NormalForm itself uses internally (append_gate right-
/// multiplies the running Clifford: `self.c = self.c * new_gate`).
fn unitary_of(s: &str) -> Mat2 {
    let mut u: Mat2 = [[(1.0, 0.0), (0.0, 0.0)], [(0.0, 0.0), (1.0, 0.0)]];
    for g in s.chars() {
        u = mat_mul(&u, &gate_matrix(g));
    }
    u
}

fn max_abs_diff(a: &Mat2, b: &Mat2) -> f64 {
    let mut m: f64 = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            let dre = a[i][j].0 - b[i][j].0;
            let dim = a[i][j].1 - b[i][j].1;
            let d = (dre * dre + dim * dim).sqrt();
            if d > m {
                m = d;
            }
        }
    }
    m
}

/// Minimal case that exercises the SHT-merge branch directly: "S","H" opens an
/// SHT syllable on the first T, and the SECOND T must merge into it.
#[test]
fn sht_merge_preserves_unitary() {
    let input = "SHTT";
    let nf = NormalForm::from_gates(input);
    let output = nf.to_gates();

    let u_in = unitary_of(input);
    let u_out = unitary_of(&output);

    let diff = max_abs_diff(&u_in, &u_out);
    assert!(
        diff < 1e-9,
        "NormalForm::from_gates(\"{input}\").to_gates() = \"{output}\" is not the \
         same unitary as the input (max diff = {diff}). This reproduces the \
         SHT-merge formula bug: the merge branch for Syllable::SHT must use \
         S*H*S, not H*S*H."
    );
}

/// A second, slightly longer case that also exercises the SHT-merge branch,
/// with a trailing T that must in turn re-open a fresh syllable.
#[test]
fn sht_merge_preserves_unitary_longer() {
    let input = "SHTTHT";
    let nf = NormalForm::from_gates(input);
    let output = nf.to_gates();

    let u_in = unitary_of(input);
    let u_out = unitary_of(&output);

    let diff = max_abs_diff(&u_in, &u_out);
    assert!(
        diff < 1e-9,
        "NormalForm::from_gates(\"{input}\").to_gates() = \"{output}\" is not the \
         same unitary as the input (max diff = {diff})."
    );
}

/// Control case: does NOT exercise the SHT-merge branch (no SHT syllable is
/// ever formed), so it must pass regardless of the H*S*H/S*H*S formula. Kept
/// here to make clear the bug is isolated to that one branch.
#[test]
fn control_case_without_sht_merge_is_unaffected() {
    let input = "THT";
    let nf = NormalForm::from_gates(input);
    let output = nf.to_gates();

    let u_in = unitary_of(input);
    let u_out = unitary_of(&output);

    let diff = max_abs_diff(&u_in, &u_out);
    assert!(
        diff < 1e-9,
        "NormalForm::from_gates(\"{input}\").to_gates() = \"{output}\" is not the \
         same unitary as the input (max diff = {diff})."
    );
}
