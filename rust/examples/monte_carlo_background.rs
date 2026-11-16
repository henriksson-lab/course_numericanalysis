// Background: a seeded generator, uniform numbers, and a Monte Carlo mean with its SE.
use rand::{rngs::StdRng, RngExt, SeedableRng};

fn main() {
    let mut rng = StdRng::seed_from_u64(42); // fix the stream
    let u: Vec<f64> = (0..5).map(|_| rng.random::<f64>()).collect(); // uniform in [0, 1)
    let rounded: Vec<f64> = u.iter().map(|v| (v * 1e4).round() / 1e4).collect();
    println!("{:?}", rounded);

    let n = 10000;
    let fx: Vec<f64> = (0..n).map(|_| rng.random::<f64>().powi(2)).collect(); // f(U), f(x) = x^2
    let mean = fx.iter().sum::<f64>() / n as f64;
    let var = fx.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
    println!("{} {}", mean, (var / n as f64).sqrt()); // exact: 1/3
}
