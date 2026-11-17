// C7: mean-square stability: Euler-Maruyama vs drift-implicit EM for noisy drug clearance.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::StandardNormal;

fn main() {
    let (k, mu, c0, p) = (0.8f64, 1.0f64, 10.0f64, 10000usize); // per hour; 10 000 patients, 48 h
    let mut rng = StdRng::seed_from_u64(1);
    for h in [0.5f64, 1.0, 2.0] {
        let (mut ce, mut ci) = (vec![c0; p], vec![c0; p]);
        for _ in 0..(48.0 / h).round() as usize {
            for j in 0..p {
                let dw = h.sqrt() * rng.sample::<f64, _>(StandardNormal);
                ce[j] += -k * ce[j] * h + mu * ce[j] * dw; // Euler-Maruyama
                ci[j] = (ci[j] + mu * ci[j] * dw) / (1.0 + k * h); // drift-implicit EM
            }
        }
        let max = |v: &Vec<f64>| v.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
        let m2 = ce.iter().map(|c| c * c).sum::<f64>() / p as f64;
        println!("h = {:<3} EM: factor {:.3} max {:.1e}  mean C^2 {:.1e} | implicit: factor {:.3} max {:.1e}",
                 h, (1.0 - k * h).powi(2) + mu * mu * h, max(&ce), m2,
                 (1.0 + mu * mu * h) / (1.0 + k * h).powi(2), max(&ci));
    }
    println!("exact mean C^2 {:.5e}", c0 * c0 * ((-2.0 * k + mu * mu) * 48.0).exp());
}
