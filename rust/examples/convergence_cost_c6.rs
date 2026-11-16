// C6: dot product with a loop, with an iterator, and with a library call.
// Run in release mode, or the timings mean nothing:
//   cargo run -q --release --example convergence_cost_c6
// Times depend on your machine; the random numbers differ from R and Python.
use nalgebra::DVector;
use rand::{rngs::StdRng, RngExt, SeedableRng};
use std::hint::black_box;
use std::time::Instant;

fn dot_loop(x: &[f64], y: &[f64]) -> f64 {
    let mut s = 0.0;
    for i in 0..x.len() {
        s += x[i] * y[i];
    }
    s
}

// seconds for one call of f(); black_box keeps the optimiser from skipping it
fn elapsed(f: impl Fn() -> f64) -> f64 {
    let start = Instant::now();
    black_box(f());
    start.elapsed().as_secs_f64()
}

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let x: Vec<f64> = (0..10_000_000).map(|_| rng.random::<f64>()).collect();
    let y: Vec<f64> = (0..10_000_000).map(|_| rng.random::<f64>()).collect();
    let (xv, yv) = (DVector::from_vec(x.clone()), DVector::from_vec(y.clone()));
    println!("{:.3}", elapsed(|| dot_loop(black_box(&x), black_box(&y))));
    println!("{:.3}", elapsed(|| black_box(&x).iter().zip(&y).map(|(a, b)| a * b).sum()));
    println!("{:.3}", elapsed(|| black_box(&xv).dot(&yv))); // nalgebra
}
