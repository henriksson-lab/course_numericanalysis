// C3: averaging over events instead of over time biases the mean.
use rand::{rngs::StdRng, RngExt, SeedableRng};

fn main() {
    let (k, delta) = (10.0, 1.0);
    let mut rng = StdRng::seed_from_u64(1);
    let (mut t, mut m) = (0.0f64, 0i64);
    let (mut n_ev, mut sum_ev, mut area) = (0u64, 0i64, 0.0f64); // after a 10 h burn-in
    while t < 5000.0 {
        let (a1, a2) = (k, delta * m as f64);
        let a0 = a1 + a2;
        let tau = -(1.0 - rng.random::<f64>()).ln() / a0;
        if t > 10.0 { area += m as f64 * tau; } // m holds for a time tau
        t += tau;
        m += if rng.random::<f64>() * a0 < a1 { 1 } else { -1 };
        if t > 10.0 { n_ev += 1; sum_ev += m; } // one record per event
    }
    println!("per event {:.4}   time-weighted {:.4}   ({} events)",
             sum_ev as f64 / n_ev as f64, area / (t - 10.0), n_ev);
}
