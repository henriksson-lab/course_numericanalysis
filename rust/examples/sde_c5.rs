// C5: Milstein against Euler-Maruyama on geometric Brownian motion: strong error, same path.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::StandardNormal;

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
    let rs: Vec<usize> = (0..9).map(|k| 1usize << k).collect();
    let mut rng = StdRng::seed_from_u64(1);
    let (mut e_em, mut e_mil) = (vec![0.0; rs.len()], vec![0.0; rs.len()]);
    let mut dwf = vec![0.0f64; nf];
    for _ in 0..m {
        for w in dwf.iter_mut() {
            *w = hf.sqrt() * rng.sample::<f64, _>(StandardNormal);
        }
        let nt = n0 * ((r - s * s / 2.0) * t + s * dwf.iter().sum::<f64>()).exp();
        for (j, &rr) in rs.iter().enumerate() {
            let h = rr as f64 * hf;
            let (mut ne, mut nm) = (n0, n0);
            for chunk in dwf.chunks(rr) {
                let dw: f64 = chunk.iter().sum();
                ne += r * ne * h + s * ne * dw;
                nm += r * nm * h + s * nm * dw + 0.5 * s * s * nm * (dw * dw - h); // + (1/2) g g' (dW^2 - h)
            }
            e_em[j] += (ne - nt).abs() / m as f64;
            e_mil[j] += (nm - nt).abs() / m as f64;
        }
    }
    let h: Vec<f64> = rs.iter().map(|&rr| rr as f64 * hf).collect();
    for j in 0..rs.len() {
        println!("h = {:.7}  EM {:.5}  Milstein {:.6}", h[j], e_em[j], e_mil[j]);
    }
    println!("slopes: EM {:.4}  Milstein {:.4}", slope(&h, &e_em), slope(&h, &e_mil));
}
