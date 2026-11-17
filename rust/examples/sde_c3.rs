// C3: the fresh-noise bug: each step size gets its own new Brownian path.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::StandardNormal;

fn main() {
    let (r, s, n0, t) = (0.8f64, 0.5f64, 1.0f64, 1.0f64);
    let (m, nf) = (10000usize, 1024usize);
    let hf = t / nf as f64;
    let rs: Vec<usize> = (0..9).map(|k| 1usize << k).collect();
    let mut rng = StdRng::seed_from_u64(1);
    let (mut strong, mut diff) = (vec![0.0; rs.len()], vec![0.0; rs.len()]);
    for _ in 0..m {
        let w: f64 = (0..nf).map(|_| hf.sqrt() * rng.sample::<f64, _>(StandardNormal)).sum();
        let nt = n0 * ((r - s * s / 2.0) * t + s * w).exp(); // exact N(T) on the fine path
        for (j, &rr) in rs.iter().enumerate() {
            let h = rr as f64 * hf;
            let mut n = n0;
            for _ in 0..nf / rr {
                let dw = h.sqrt() * rng.sample::<f64, _>(StandardNormal); // BUG: fresh increment
                n += r * n * h + s * n * dw;
            }
            strong[j] += (n - nt).abs() / m as f64;
            diff[j] += (n - nt) / m as f64;
        }
    }
    for j in 0..rs.len() {
        println!("h = 1/{:<4}  strong {:.4}  weak {:.4}", nf / rs[j], strong[j], diff[j].abs());
    }
}
