// C5: the order of operations, (AB)x against A(Bx).
// Run in release mode, or the timings mean nothing:
//   cargo run -q --release --example convergence_cost_c5
// Matrix products with nalgebra (its own optimised kernel, one core, no BLAS).
// Times depend on your machine; the random numbers differ from R and Python.
use nalgebra::{DMatrix, DVector};
use rand::{rngs::StdRng, RngExt, SeedableRng};
use std::time::Instant;

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let n = 2000;
    let a = DMatrix::from_fn(n, n, |_, _| rng.random::<f64>());
    let b = DMatrix::from_fn(n, n, |_, _| rng.random::<f64>());
    let x = DVector::from_fn(n, |_, _| rng.random::<f64>());
    let start = Instant::now();
    let y1 = (&a * &b) * &x; // `&`: borrow the matrices, don't move them
    println!("{:.3}", start.elapsed().as_secs_f64());
    let start = Instant::now();
    let y2 = &a * (&b * &x);
    println!("{:.3}", start.elapsed().as_secs_f64());
    println!("{:e}", (&y1 - &y2).amax() / y1.amax()); // amax: largest absolute value
}
