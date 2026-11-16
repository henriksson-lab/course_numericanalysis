// C5: finite vs infinite variance; the SE is only trustworthy for the first.
use rand::{rngs::StdRng, RngExt, SeedableRng};

fn mean_se(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0);
    (m, (var / n).sqrt())
}

fn main() {
    for seed in 1..=5 {
        let mut rng = StdRng::seed_from_u64(seed);
        // 1 - U is in (0, 1], so x^-0.75 never divides by 0
        let x: Vec<f64> = (0..100000).map(|_| 1.0 - rng.random::<f64>()).collect();
        let fa: Vec<f64> = x.iter().map(|v| v.powf(-0.75)).collect(); // exact 4
        let fb: Vec<f64> = x.iter().map(|v| v.powf(-0.25)).collect(); // exact 4/3
        let (ma, sa) = mean_se(&fa);
        let (mb, sb) = mean_se(&fb);
        println!("seed {}   x^-0.75: {:.4} (se {:.4})   x^-0.25: {:.5} (se {:.5})", seed, ma, sa, mb, sb);
    }
}
