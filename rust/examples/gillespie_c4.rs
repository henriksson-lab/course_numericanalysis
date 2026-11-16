// C4: tau-leaping for mRNA birth-death; variance and negative counts versus tau.
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, Poisson};

fn rpois(lambda: f64, rng: &mut StdRng) -> i64 {
    if lambda <= 0.0 { 0 } else { Poisson::new(lambda).unwrap().sample(rng) as i64 }
}

fn tau_leap(tau: f64, rng: &mut StdRng) -> (i64, bool) {
    let (t_end, k, delta) = (10.0, 10.0, 1.0);
    let (mut m, mut went_negative) = (0i64, false);
    for _ in 0..(t_end / tau as f64).round() as usize {
        m += rpois(k * tau, rng) - rpois(delta * m as f64 * tau, rng); // many reactions per leap
        if m < 0 { went_negative = true; m = 0; } // impossible: clip
    }
    (m, went_negative)
}

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    for tau in [0.01, 0.1, 0.5, 1.0, 1.5] {
        let r: Vec<(i64, bool)> = (0..2000).map(|_| tau_leap(tau, &mut rng)).collect();
        let n = r.len() as f64;
        let mean = r.iter().map(|x| x.0 as f64).sum::<f64>() / n;
        let var = r.iter().map(|x| (x.0 as f64 - mean).powi(2)).sum::<f64>() / (n - 1.0);
        let neg = r.iter().filter(|x| x.1).count() as f64 / n;
        println!("tau {:4.2}  mean {:6.3}  var {:6.3}  theory var {:6.3}  negative {:.4}",
                 tau, mean, var, 20.0 / (2.0 - tau), neg);
    }
}
