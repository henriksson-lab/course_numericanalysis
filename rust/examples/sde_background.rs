// Background: Euler-Maruyama for a bead in an optical trap (Ornstein-Uhlenbeck), 2000 beads.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::StandardNormal;

fn main() {
    let (theta, d) = (1.0f64, 430.0f64); // trap: per ms, nm^2/ms
    let sigma = (2.0 * d).sqrt(); // nm/sqrt(ms)
    let mut rng = StdRng::seed_from_u64(1);
    let h = 0.01;
    let mut x = vec![60.0f64; 2000]; // 2000 beads, all at 60 nm
    for _ in 0..100 {
        // 100 EM steps to t = 1 ms
        for xi in x.iter_mut() {
            let z: f64 = rng.sample(StandardNormal);
            *xi += -theta * *xi * h + sigma * h.sqrt() * z;
        }
    }
    let n = x.len() as f64;
    let mean = x.iter().sum::<f64>() / n;
    let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
    println!("mean {:.4}  var {:.4}", mean, var);
    println!("exact mean {:.4}  var {:.4}", 60.0 * (-theta).exp(), d / theta * (1.0 - (-2.0 * theta).exp()));
}
