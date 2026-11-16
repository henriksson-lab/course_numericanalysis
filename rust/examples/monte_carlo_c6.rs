// C6: the curse of dimensionality. Marginal likelihood Z in d dimensions,
// grid (n points per axis at prior quantiles) vs Monte Carlo from the prior.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::StandardNormal;
use statrs::distribution::{ContinuousCDF, Normal};

fn log_dnorm(x: f64, mu: f64, sd: f64) -> f64 {
    -0.5 * ((x - mu) / sd).powi(2) - sd.ln() - 0.5 * (2.0 * std::f64::consts::PI).ln()
}

// likelihood of y_i = 1 (noise sd 1) at one parameter vector th
fn lik(th: &[f64]) -> f64 {
    th.iter().map(|&t| log_dnorm(1.0, t, 1.0)).sum::<f64>().exp()
}

fn main() {
    let std_normal = Normal::new(0.0, 1.0).unwrap();
    let mut rng = StdRng::seed_from_u64(1);
    println!("   d       n       N  grid_err   mc_err mc_relse");
    for d in [1usize, 2, 3, 4, 6, 8, 10] {
        let z = log_dnorm(1.0, 0.0, 2f64.sqrt()).exp().powi(d as i32);
        let n = (1e5f64.powf(1.0 / d as f64) + 1e-9).floor() as usize;
        let nodes: Vec<f64> = (1..=n).map(|i| std_normal.inverse_cdf((i as f64 - 0.5) / n as f64)).collect();
        let big_n = n.pow(d as u32);
        // grid: decode point index k into d digits in base n
        let mut th = vec![0.0; d];
        let mut grid_sum = 0.0;
        for k in 0..big_n {
            let mut q = k;
            for j in 0..d {
                th[j] = nodes[q % n];
                q /= n;
            }
            grid_sum += lik(&th);
        }
        // Monte Carlo with the same number of points
        let l: Vec<f64> = (0..big_n).map(|_| {
            let th: Vec<f64> = (0..d).map(|_| rng.sample::<f64, _>(StandardNormal)).collect();
            lik(&th)
        }).collect();
        let m = l.iter().sum::<f64>() / big_n as f64;
        let s = (l.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (big_n - 1) as f64).sqrt();
        println!("{:4} {:7} {:7} {:9.4} {:8.4} {:8.4}", d, n, big_n,
                 (grid_sum / big_n as f64 / z - 1.0).abs(), (m / z - 1.0).abs(),
                 s / (big_n as f64).sqrt() / m);
    }
}
