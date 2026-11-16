// C2: inverse-transform sampling, continuous (exponential) and discrete (bases).
use rand::{rngs::StdRng, RngExt, SeedableRng};

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let lambda = 0.5;
    let n = 100000;
    // random::<f64>() is in [0, 1): use 1 - U so the log never sees 0
    let x: Vec<f64> = (0..n).map(|_| -(1.0 - rng.random::<f64>()).ln() / lambda).collect();
    let mean = x.iter().sum::<f64>() / n as f64;
    let sd = (x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64).sqrt();
    let p6 = x.iter().filter(|&&v| v > 6.0).count() as f64 / n as f64;
    println!("{} {} {} {}", mean, sd, p6, (-lambda * 6.0).exp());

    let bases = ['A', 'C', 'G', 'T'];
    let cum = [0.3, 0.5, 0.7, 1.0]; // cumulative probabilities
    let mut counts = [0usize; 4];
    for _ in 0..n {
        let u = rng.random::<f64>();
        let k = cum.partition_point(|&c| c <= u); // first k with cum[k] > u
        counts[k] += 1;
    }
    for k in 0..4 {
        println!("{} {:.3}", bases[k], counts[k] as f64 / n as f64);
    }
}
