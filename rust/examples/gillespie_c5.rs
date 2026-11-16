// C5: Euler-Maruyama for the chemical Langevin equation of mRNA birth-death.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::StandardNormal;

fn em(h: f64, k: f64, rng: &mut StdRng) -> f64 {
    let (t_end, delta) = (10.0, 1.0);
    let mut m = 0.0f64;
    for _ in 0..(t_end / h as f64).round() as usize {
        let z: f64 = rng.sample(StandardNormal);
        m += (k - delta * m) * h + ((k + delta * m) * h).sqrt() * z; // sqrt(<0) is NaN
    }
    m
}

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let m: Vec<f64> = (0..2000).map(|_| em(0.01, 10.0, &mut rng)).collect();
    let n = m.len() as f64;
    let mean = m.iter().sum::<f64>() / n;
    let var = m.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
    println!("mean {:.4}  var {:.4}  fano {:.4}", mean, var, var / mean);
    let m: Vec<f64> = (0..2000).map(|_| em(0.01, 0.5, &mut rng)).collect();
    let nan = m.iter().filter(|x| x.is_nan()).count();
    let ok: Vec<&f64> = m.iter().filter(|x| !x.is_nan()).collect();
    let neg = ok.iter().filter(|&&&x| x < 0.0).count() as f64 / ok.len() as f64;
    println!("nan {:.4}  negative {:.4}", nan as f64 / n, neg);
}
