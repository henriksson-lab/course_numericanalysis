// Timing and measuring memory in Rust.
// Run in release mode for honest timings:
//   cargo run -q --release --example convergence_cost_background
// Times depend on your machine; the random numbers differ from R and Python.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use std::hint::black_box;
use std::time::Instant;

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let mut x: Vec<f64> = (0..1_000_000).map(|_| rng.random::<f64>()).collect();
    let bytes = x.len() * std::mem::size_of::<f64>(); // bytes of the data
    let start = Instant::now();
    x.sort_by(f64::total_cmp); // f64 has no plain sort(): NaN has no order
    black_box(&x); // stop the optimiser from skipping the sort
    println!("{}", start.elapsed().as_secs_f64()); // seconds, wall clock
    println!("{bytes}"); // bytes
}
