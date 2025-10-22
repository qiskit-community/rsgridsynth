# gridsynth

`gridsynth` is a high-precision Rust-based synthesizer that decomposes single-qubit Z-axis rotations into Clifford+T gate sequences, using number-theoretic algorithms and geometry-of-numbers methods.


## Build

Make sure you have [Rust](https://www.rust-lang.org/tools/install) installed.

```bash
cargo build --release
```

The binary will be at:

```
target/release/gridsynth
```

## Usage

```bash
./target/release/gridsynth <theta> <epsilon> [OPTIONS]
```

* `<theta>`: The rotation angle in radians (e.g. `0.6`)
* `<epsilon>`: The target approximation error(still do not support precision higher than 1e-8.) (e.g. `1e-8`)

### Options

| Option                | Description                                             | Default |
| --------------------- | ------------------------------------------------------- | ------- |
| `--dps <dps>`         | Decimal precision used in floating-point computations   | `425`   |
| `-d, --dtimeout <ms>` | Timeout for the Diophantine step (in milliseconds)      | `200`   |
| `-f, --ftimeout <ms>` | Timeout for integer factoring (in milliseconds)         | `50`    |
| `-v, --verbose`       | Enable verbose output (intermediate steps, diagnostics) | `false` |
| `-t, --time`          | Show elapsed wall-clock time after computation          | `false` |

---

## Example

```bash
./target/release/gridsynth 0.6 1e-8 -v -t
```

Sample output:

```
to_upright_set_pair: 0.001 s
time of diophantine_dyadic: 0.036 ms
time of decompose_domega_unitary: 0.089 ms
total time: 1.294 ms
Elapsed time: 1.308ms
HTSHTHTSHTSHTSHTHTSHTHTHTSHTHTSHTHTSHTHTSHTSHTSHTSHTSHTHTSHTHTSHTSHTHTHTSHTHTHTHTSHTSHTSHTSHTSHTSHTHTHTSHTSHTSHTSHTSHTHTSHTHTSHTSHTHTHTSHTSHTSHTHTHTSHTHTSHTSHTSHTSHTHTSHTHTSHTHTHTHTSHTHTSHTSHTSHTHTSHTSHTSHTHTHTHTSHSSSWWWWW
```


## Acknowledgement

This package is a reimplementation of the [pygridsynth](https://github.com/quantum-programming/pygridsynth) python package by Shuntaro Yamato and Noboyuki Yoshioka.
This, in turn, is based on "Optimal ancilla-free Clifford+T approximation of z-rotations" by Neil J. Ross and Peter Selinger ([arXiv:1403.2975](https://arxiv.org/abs/1403.2975)) and its implementation, [newsynth](https://www.mathstat.dal.ca/~selinger/newsynth/).
Please consider citing these previous works.
