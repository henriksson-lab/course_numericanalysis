// C4: geometric Brownian motion: the mean grows, the median (typical population) shrinks.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::StandardNormal;
use statrs::distribution::{ContinuousCDF, Normal};

fn median(v: &[f64]) -> f64 {
    let mut w = v.to_vec();
    w.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = w.len();
    if n % 2 == 0 { (w[n / 2 - 1] + w[n / 2]) / 2.0 } else { w[n / 2] }
}

fn main() {
    let (r, n0, h, m) = (0.1f64, 1000.0f64, 0.05f64, 10000usize);
    let phi = Normal::new(0.0, 1.0).unwrap();
    let mut rng = StdRng::seed_from_u64(1);
    for s in [0.6f64, 0.4] {
        println!("sigma = {:.1}", s);
        let mut n = vec![n0; m];
        let mut at10 = Vec::new();
        for i in 1..=(50.0 / h).round() as usize {
            for x in n.iter_mut() {
                let z: f64 = rng.sample(StandardNormal);
                *x += r * *x * h + s * *x * h.sqrt() * z;
            }
            if i == (10.0 / h).round() as usize {
                at10 = n.clone();
            }
        }
        for (yr, v) in [(10.0f64, &at10), (50.0, &n)] {
            let mean = v.iter().sum::<f64>() / m as f64;
            let low = v.iter().filter(|&&x| x < 10.0).count() as f64 / m as f64;
            let p_ex = phi.cdf(((10.0f64 / 1000.0).ln() - (r - s * s / 2.0) * yr) / (s * yr.sqrt()));
            println!("  t = {:2}  mean {:9.0} (exact {:9.0})  median {:6.1} (exact {:6.1})  P(N < 10) {:.4} (exact {:.4})",
                     yr, mean, n0 * (r * yr).exp(), median(v), n0 * ((r - s * s / 2.0) * yr).exp(), low, p_ex);
        }
    }
}
