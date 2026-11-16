// C6: smoothing a noisy bioluminescence trace by low-pass filtering, raw and detrended.
// The random noise differs from R and Python, so the digits differ a little.
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, Normal};
use rustfft::{num_complex::Complex, FftPlanner};
use std::f64::consts::PI;

fn lowpass(y: &[f64], kmax: usize) -> Vec<f64> {
    let n = y.len();
    let mut planner = FftPlanner::<f64>::new();
    let mut buf: Vec<Complex<f64>> = y.iter().map(|&v| Complex::new(v, 0.0)).collect();
    planner.plan_fft_forward(n).process(&mut buf);
    for k in 0..n {
        let signed = if k <= n / 2 { k } else { n - k };   // |frequency| of entry k
        if signed > kmax { buf[k] = Complex::new(0.0, 0.0); }
    }
    planner.plan_fft_inverse(n).process(&mut buf);         // unnormalised
    buf.iter().map(|z| z.re / n as f64).collect()
}

fn rms(e: &[f64]) -> f64 {
    (e.iter().map(|v| v * v).sum::<f64>() / e.len() as f64).sqrt()
}

fn main() {
    let n = 240;                                          // every 30 min for 5 days
    let t: Vec<f64> = (0..n).map(|j| 0.5 * j as f64).collect();
    let truth: Vec<f64> = t.iter().map(|&s| 20.0 + 60.0 * (-s / 40.0).exp()
        + 30.0 * (-s / 80.0).exp() * (2.0 * PI * (s - 4.0) / 24.0).cos()).collect();
    let mut rng = StdRng::seed_from_u64(3);
    let noise = Normal::new(0.0, 6.0).unwrap();
    let y: Vec<f64> = truth.iter().map(|&v| v + noise.sample(&mut rng)).collect();
    // straight-line fit by the textbook formulas for slope and intercept
    let (tm, ym) = (t.iter().sum::<f64>() / n as f64, y.iter().sum::<f64>() / n as f64);
    let slope = (0..n).map(|j| (t[j] - tm) * (y[j] - ym)).sum::<f64>()
        / t.iter().map(|&s| (s - tm).powi(2)).sum::<f64>();
    let trend: Vec<f64> = t.iter().map(|&s| ym + slope * (s - tm)).collect();
    let ends: Vec<usize> = (0..10).chain(n - 10..n).collect();   // first and last 5 h
    for kmax in [5, 10, 15, 20, 30, 60] {
        let s = lowpass(&y, kmax);
        let resid: Vec<f64> = (0..n).map(|j| y[j] - trend[j]).collect();
        let d: Vec<f64> = lowpass(&resid, kmax).iter().zip(&trend).map(|(a, b)| a + b).collect();
        let es: Vec<f64> = (0..n).map(|j| s[j] - truth[j]).collect();
        let ed: Vec<f64> = (0..n).map(|j| d[j] - truth[j]).collect();
        let pick = |e: &[f64]| ends.iter().map(|&j| e[j]).collect::<Vec<f64>>();
        println!("kmax {:2} (periods >= {:5.1} h): raw {:.2} (ends {:.2})  detrended {:.2} (ends {:.2})",
                 kmax, 120.0 / kmax as f64, rms(&es), rms(&pick(&es)), rms(&ed), rms(&pick(&ed)));
    }
    let e: Vec<f64> = (0..n).map(|j| y[j] - truth[j]).collect();
    println!("noise: {:.2}", rms(&e));
}
