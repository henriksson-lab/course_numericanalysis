// C5: leakage of a 25.5-h rhythm, the Hann window, and zero-padding.
use rustfft::{num_complex::Complex, FftPlanner};
use std::f64::consts::PI;

fn fft(x: &[f64]) -> Vec<Complex<f64>> {
    let mut buf: Vec<Complex<f64>> = x.iter().map(|&v| Complex::new(v, 0.0)).collect();
    FftPlanner::<f64>::new().plan_fft_forward(x.len()).process(&mut buf);
    buf
}

fn main() {
    let n = 240;                                       // 5 days, every 30 min
    let t: Vec<f64> = (0..n).map(|j| 0.5 * j as f64).collect();
    let y: Vec<f64> = t.iter().map(|&s| (2.0 * PI * s / 25.5).cos()).collect();
    let w: Vec<f64> = (0..n).map(|j| 0.5 - 0.5 * (2.0 * PI * j as f64 / n as f64).cos()).collect();
    let sw: f64 = w.iter().sum();
    let wy: Vec<f64> = (0..n).map(|j| w[j] * y[j]).collect();
    let amp: Vec<f64> = fft(&y).iter().map(|z| 2.0 * z.norm() / n as f64).collect();
    let ampw: Vec<f64> = fft(&wy).iter().map(|z| 2.0 * z.norm() / sw).collect();
    for k in [1, 2, 3, 4, 5, 6, 7, 8, 10, 20] {
        println!("k = {:2}  no window {:.3}  Hann {:.3}", k, amp[k], ampw[k]);
    }
    // (b) zero-padding to 10x the length: a finer grid of frequencies
    let big_n = 10 * n;
    let mut padded = wy.clone();
    padded.resize(big_n, 0.0);
    let yp = fft(&padded);
    let kp = (1..big_n / 2).max_by(|&a, &b| yp[a].norm().partial_cmp(&yp[b].norm()).unwrap()).unwrap();
    println!("{}", 1200.0 / kp as f64);                // period estimate (h)
    // (c) add a weak 12-h rhythm (amplitude 0.05)
    let y2: Vec<f64> = (0..n).map(|j| y[j] + 0.05 * (2.0 * PI * t[j] / 12.0).cos()).collect();
    let wy2: Vec<f64> = (0..n).map(|j| w[j] * y2[j]).collect();
    println!("{:.4} {:.4}", 2.0 * fft(&y2)[10].norm() / n as f64, 2.0 * fft(&wy2)[10].norm() / sw);
    println!("{:.4} {:.4}", amp[10], ampw[10]);        // leakage alone at k = 10
}
