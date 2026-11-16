// C4: Monte Carlo AUC of the logistic OD curve, its SE, and the 1/sqrt(N) law.
use rand::{rngs::StdRng, RngExt, SeedableRng};

const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

fn mean_sd(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    (m, (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0)).sqrt())
}

fn main() {
    let exact = K / R * (((R * 10.0).exp() + K / N0 - 1.0) / (K / N0)).ln();
    let mut rng = StdRng::seed_from_u64(1);
    let fx: Vec<f64> = (0..1000).map(|_| 10.0 * od(rng.random_range(0.0..10.0))).collect();
    let (est, s) = mean_sd(&fx);
    let se = s / 1000f64.sqrt();
    println!("{} {} {} {} {}", est, se, est - 2.0 * se, est + 2.0 * se, exact);

    let ns = [1e2, 1e3, 1e4, 1e5];
    let mut rms = Vec::new();
    for &n in &ns {
        let mut sq = 0.0;
        for _ in 0..20 {
            let est: f64 = (0..n as usize).map(|_| 10.0 * od(rng.random_range(0.0..10.0))).sum::<f64>() / n;
            sq += (est - exact).powi(2);
        }
        rms.push((sq / 20.0).sqrt());
    }
    println!("{:.4?}", rms);
    // least-squares slope of log(rms) against log(N)
    let lx: Vec<f64> = ns.iter().map(|n| n.ln()).collect();
    let ly: Vec<f64> = rms.iter().map(|e| e.ln()).collect();
    let (mx, my) = (lx.iter().sum::<f64>() / 4.0, ly.iter().sum::<f64>() / 4.0);
    let slope = lx.iter().zip(&ly).map(|(x, y)| (x - mx) * (y - my)).sum::<f64>()
        / lx.iter().map(|x| (x - mx).powi(2)).sum::<f64>();
    println!("{}", slope);
    println!("{}", (2.0 * s / 0.01).powi(2).ceil()); // N for 2 SE = 0.01
}
