// C3: amplitude spectrum, cosinor least squares (nalgebra) and the trig interpolant.
use nalgebra::{DMatrix, DVector};
use rustfft::{num_complex::Complex, FftPlanner};
use std::f64::consts::PI;

fn fft(x: &[f64]) -> Vec<Complex<f64>> {
    let mut buf: Vec<Complex<f64>> = x.iter().map(|&v| Complex::new(v, 0.0)).collect();
    FftPlanner::<f64>::new().plan_fft_forward(x.len()).process(&mut buf);
    buf
}

fn trig_interp(x: &[f64], big_t: f64, tt: f64) -> f64 {
    let n = x.len();
    let xf: Vec<Complex<f64>> = fft(x).iter().map(|z| z / n as f64).collect();
    let mut p = xf[0].re + xf[n / 2].re * (PI * n as f64 * tt / big_t).cos();
    for k in 1..n / 2 {
        let a = 2.0 * PI * k as f64 * tt / big_t;
        p += 2.0 * xf[k].re * a.cos() - 2.0 * xf[k].im * a.sin();
    }
    p
}

fn main() {
    let x = [7.81, 8.81, 9.05, 9.98, 9.40, 8.50, 8.15, 7.47, 6.87, 6.41, 7.15, 7.37,
             7.81, 8.09, 9.64, 9.49, 9.29, 9.03, 8.25, 7.43, 6.98, 6.73, 6.72, 6.65];
    let n = x.len();
    let t: Vec<f64> = (0..n).map(|j| 2.0 * j as f64).collect();
    let big_t = 48.0;
    let xf = fft(&x);
    for k in 0..8 {                                   // (a) amplitude spectrum
        let a = if k == 0 { xf[0].norm() / n as f64 } else { 2.0 * xf[k].norm() / n as f64 };
        println!("k = {}  period {:5.1} h  amplitude {:.3}", k, big_t / k as f64, a);
    }
    // (b) cosinor: least squares with mean + 24-h cosine + 24-h sine
    let a = DMatrix::from_fn(n, 3, |i, j| match j {
        0 => 1.0,
        1 => (2.0 * PI * t[i] / 24.0).cos(),
        _ => (2.0 * PI * t[i] / 24.0).sin(),
    });
    let coef = a.svd(true, true).solve(&DVector::from_column_slice(&x), 1e-12).unwrap();
    println!("{:.4} {:.4} {:.4}", coef[0], coef[1], coef[2]);
    println!("{:.4} {:.4} {:.4}", xf[0].re / n as f64, 2.0 * xf[2].re / n as f64,
             -2.0 * xf[2].im / n as f64);
    // (c) the interpolant, on a fine grid
    let err = (0..n).map(|j| (trig_interp(&x, big_t, t[j]) - x[j]).abs()).fold(0.0, f64::max);
    println!("{:e}", err);                            // passes through the data
    let (mut pmax, mut tmax) = (f64::MIN, 0.0);
    for i in 0..=4800 {
        let tt = i as f64 * 0.01;
        let p = trig_interp(&x, big_t, tt);
        if p > pmax { pmax = p; tmax = tt; }
    }
    println!("{:.4} {:.4}", pmax, tmax);              // highest point of the curve
}
