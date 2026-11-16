// C4: aliasing. The same gene sampled every 2 h and every 4 h; a rhythm at Nyquist.
use rustfft::{num_complex::Complex, FftPlanner};
use std::f64::consts::PI;

fn fft(x: &[f64]) -> Vec<Complex<f64>> {
    let mut buf: Vec<Complex<f64>> = x.iter().map(|&v| Complex::new(v, 0.0)).collect();
    FftPlanner::<f64>::new().plan_fft_forward(x.len()).process(&mut buf);
    buf
}

fn gene(t: f64) -> f64 {
    8.0 + 1.5 * (2.0 * PI * (t - 6.0) / 24.0).cos() + (2.0 * PI * t / 4.8).cos()
}

fn main() {
    for dt in [2usize, 4] {
        let x: Vec<f64> = (0..48 / dt).map(|j| gene((j * dt) as f64)).collect();
        let n = x.len();
        let xf = fft(&x);
        let amp: Vec<f64> = xf.iter().map(|z| 2.0 * z.norm() / n as f64).collect();
        let kmax = (1..n / 2).max_by(|&a, &b| amp[a].partial_cmp(&amp[b]).unwrap()).unwrap();
        let k24 = n * dt / 24;                          // the 24-h bin
        let peak = (-xf[k24].arg() / (2.0 * PI) * 24.0).rem_euclid(24.0);
        println!("every {} h: strongest k = {} (period {:.1} h); 24-h amplitude {:.3}, peak {:.2} h",
                 dt, kmax, 48.0 / kmax as f64, amp[k24], peak);
    }
    // (c) an 8-h rhythm at the Nyquist frequency of the 4-h design
    for p in [0.0, 1.0, 2.0] {
        let x: Vec<f64> = (0..12).map(|j| (2.0 * PI * (4.0 * j as f64 - p) / 8.0).cos()).collect();
        println!("peak at {} h: amplitude seen {:.3}", p, fft(&x)[6].norm() / 12.0);
    }
}
