// C3: rejection sampling from Beta(2, 4) with a uniform proposal.
use rand::{rngs::StdRng, RngExt, SeedableRng};

fn f(x: f64) -> f64 {
    20.0 * x * (1.0 - x).powi(3)
}

// returns (acceptance rate, accepted draws)
fn rejection(n: usize, m: f64, rng: &mut StdRng) -> (f64, Vec<f64>) {
    let mut acc = Vec::new();
    for _ in 0..n {
        let x = rng.random::<f64>();
        if rng.random::<f64>() < f(x) / m {
            acc.push(x);
        }
    }
    (acc.len() as f64 / n as f64, acc)
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    for m in [f(0.25), 10.0] {
        let (rate, x) = rejection(100000, m, &mut rng);
        println!("M = {:6.3}  acceptance {:.4}  (1/M = {:.4})  mean {:.4}", m, rate, 1.0 / m, mean(&x));
    }
    let (_, x) = rejection(100000, 1.5, &mut rng); // M too small
    println!("{}", mean(&x));
}
