// C2: strong and weak error of Euler-Maruyama on geometric Brownian motion, same Brownian path.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::StandardNormal;

// least-squares slope of log(err) against log(h)
fn slope(h: &[f64], e: &[f64]) -> f64 {
    let (x, y): (Vec<f64>, Vec<f64>) = h.iter().zip(e).map(|(a, b)| (a.ln(), b.ln())).unzip();
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let sxy: f64 = x.iter().zip(&y).map(|(a, b)| (a - mx) * (b - my)).sum();
    sxy / x.iter().map(|a| (a - mx).powi(2)).sum::<f64>()
}

fn main() {
    let (r, s, n0, t) = (0.8f64, 0.5f64, 1.0f64, 1.0f64);
    let (m, nf) = (10000usize, 1024usize);
    let hf = t / nf as f64;
    let rs: Vec<usize> = (0..9).map(|k| 1usize << k).collect(); // coarse step h = R hf
    let mut rng = StdRng::seed_from_u64(1);
    let (mut strong, mut diff) = (vec![0.0; rs.len()], vec![0.0; rs.len()]);
    let mut dwf = vec![0.0f64; nf];
    for _ in 0..m {
        for w in dwf.iter_mut() {
            *w = hf.sqrt() * rng.sample::<f64, _>(StandardNormal); // fine increments
        }
        let nt = n0 * ((r - s * s / 2.0) * t + s * dwf.iter().sum::<f64>()).exp(); // exact N(T)
        for (j, &rr) in rs.iter().enumerate() {
            let h = rr as f64 * hf;
            let mut n = n0;
            for chunk in dwf.chunks(rr) {
                let dw: f64 = chunk.iter().sum(); // sum R fine increments
                n += r * n * h + s * n * dw;
            }
            strong[j] += (n - nt).abs() / m as f64;
            diff[j] += (n - nt) / m as f64;
        }
    }
    let h: Vec<f64> = rs.iter().map(|&rr| rr as f64 * hf).collect();
    let weak: Vec<f64> = diff.iter().map(|d| d.abs()).collect();
    for j in 0..rs.len() {
        println!("h = {:.7}  strong {:.5}  weak {:.6}", h[j], strong[j], weak[j]);
    }
    println!("slopes: strong {:.4}  weak {:.4}", slope(&h, &strong), slope(&h, &weak));
}
