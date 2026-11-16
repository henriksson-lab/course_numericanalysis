// Background: one step of Gillespie's algorithm, with a seeded generator.
use rand::{rngs::StdRng, RngExt, SeedableRng};

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let a = [10.0, 10.0]; // propensities at m = 10: make, degrade (per hour)
    let a0: f64 = a.iter().sum();
    let tau = -(1.0 - rng.random::<f64>()).ln() / a0; // waiting time, Exp(a0) by inverse transform
    let r = rng.random::<f64>() * a0; // reaction j with probability a_j / a0
    let j = if r < a[0] { 0 } else { 1 };
    println!("tau = {}  reaction = {}", tau, ["make", "degrade"][j]);
}
