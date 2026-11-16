// C8: importance sampling for a rare event, and a proposal that is too narrow.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::StandardNormal;
use statrs::distribution::{ContinuousCDF, Normal};

fn dnorm(x: f64, mu: f64, sd: f64) -> f64 {
    (-0.5 * ((x - mu) / sd).powi(2)).exp() / (sd * (2.0 * std::f64::consts::PI).sqrt())
}

fn mean_sd(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    (m, (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0)).sqrt())
}

fn main() {
    let n = 100000;
    let sqn = (n as f64).sqrt();
    let mut rng = StdRng::seed_from_u64(1);
    let mut z = || -> f64 { rng.sample(StandardNormal) };
    let hits = (0..n).filter(|_| z() > 5.0).count();
    println!("{} {:e}", hits as f64 / n as f64, Normal::new(0.0, 1.0).unwrap().sf(5.0));
    for mu in [3.0, 5.0, 7.0] {
        let hw: Vec<f64> = (0..n).map(|_| {
            let x = mu + z();
            if x > 5.0 { dnorm(x, 0.0, 1.0) / dnorm(x, mu, 1.0) } else { 0.0 }
        }).collect();
        let (m, s) = mean_sd(&hw);
        println!("proposal N({}, 1): estimate {:.4e}  relative se {:.4}", mu, m, s / sqn / m);
    }
    for sd in [0.5, 1.0, 2.0, 5.0] {
        let w: Vec<f64> = (0..n).map(|_| { let x = sd * z(); dnorm(x, 0.0, 1.0) / dnorm(x, 0.0, sd) }).collect();
        let (m, s) = mean_sd(&w);
        let sum: f64 = w.iter().sum();
        let ess = sum * sum / w.iter().map(|x| x * x).sum::<f64>();
        let maxw = w.iter().cloned().fold(0.0, f64::max);
        println!("proposal sd {:.1}: estimate {:.4}  se {:.4}  ESS {:6.0}  max share {:.5}",
                 sd, m, s / sqn, ess, maxw / sum);
    }
}
